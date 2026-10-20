"""
Shared results database: a species-level cache of computed energies plus a
bond-type index of computed BDEs (target workflow Aim 3).

Why SQLite, on Lustre: one file, no server process to start/stop inside the
Slurm block, microsecond primary-key lookups, real SQL for the retrieval
queries the planner uses (storage.rag), and stdlib-only. Writes are rare (one
per computed species, minutes apart), so lock contention is negligible.
Orion is mounted with `flock`, so SQLite's POSIX locks are honoured across
compute nodes. The journal stays in rollback (DELETE) mode on purpose: WAL
needs a shared-memory index that does not work across nodes on a network
filesystem.

Enabled by pointing COMPCHEM_DB at a file path; when it is unset get_db()
returns None and every caller falls back to computing as before.

Cross-pipeline deduplication: a worker that misses inserts a 'running' claim
row (atomic - the primary key makes exactly one inserter win); other workers
that need the same species wait for it instead of computing it again. While
computing, the claimant refreshes the claim every COMPCHEM_DB_CLAIM_TTL/4
seconds; a claim not refreshed for COMPCHEM_DB_CLAIM_TTL seconds (default
600) belongs to a dead worker (SIGKILL, OOM, walltime) and is taken over.
Each claim carries a unique token, and only its holder can complete it - a
worker whose claim was taken over cannot overwrite the new claimant's row.
"""
from __future__ import annotations

import hashlib
import json
import os
import socket
import sqlite3
import threading
import time
import uuid
from contextlib import contextmanager
from typing import Any, Callable, Iterator

SCHEMA_VERSION = 1

_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT);
CREATE TABLE IF NOT EXISTS species (
    key              TEXT PRIMARY KEY,
    canonical_smiles TEXT NOT NULL,
    charge           INTEGER NOT NULL,
    multiplicity     INTEGER NOT NULL,
    method           TEXT NOT NULL,
    status           TEXT NOT NULL,          -- running | done | failed
    energy_hartree   REAL,
    converged        INTEGER,
    result_json      TEXT,                   -- full compute_species_* result
    error            TEXT,
    wall_s           REAL,                   -- compute time of the run that produced it
    owner            TEXT,                   -- host:pid[:run] of the claimant
    run_id           TEXT,
    claimed_at       REAL,
    finished_at      REAL,
    hits             INTEGER NOT NULL DEFAULT 0
);
CREATE INDEX IF NOT EXISTS species_smiles ON species (canonical_smiles, charge, multiplicity);
CREATE TABLE IF NOT EXISTS bde (
    id                 INTEGER PRIMARY KEY AUTOINCREMENT,
    parent_key         TEXT NOT NULL,
    fragment_keys      TEXT NOT NULL,        -- JSON list, sorted
    parent_smiles      TEXT NOT NULL,
    fragment_smiles    TEXT NOT NULL,        -- JSON list, sorted canonical
    bond_type          TEXT,                 -- e.g. 'O-H', NULL if undetermined
    method             TEXT NOT NULL,
    bde_hartree        REAL,
    bde_kcal_mol       REAL,
    qc_flag            TEXT,
    source_question    TEXT,
    run_id             TEXT,
    created_at         REAL NOT NULL,
    UNIQUE (parent_key, fragment_keys)
);
CREATE INDEX IF NOT EXISTS bde_bond ON bde (bond_type);
CREATE INDEX IF NOT EXISTS bde_parent ON bde (parent_smiles);
"""


# --------------------------------------------------------------------------
# Keys
# --------------------------------------------------------------------------

def canonical_smiles(smiles: str) -> str:
    """RDKit canonical SMILES, so 'O', '[H]O[H]' and 'OH2'-style spellings of
    one species share a cache entry. Falls back to the raw string (prefixed,
    so it can never collide with a real canonical form) if RDKit can't parse
    it - the compute will most likely fail too, and that failure is cached
    under the same raw key."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return f"raw:{smiles}"
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
        except Exception:
            return f"raw:{smiles}"
    return Chem.MolToSmiles(mol)


def normalize_basis(basis: str) -> str:
    """One spelling per basis set: planners write 6-31G*, 6-31g(d), 6-31G(d)
    and def2SVP for the same thing (all seen on Frontier), which would
    otherwise be separate cache entries."""
    b = basis.lower().replace(" ", "")
    b = b.replace("(d,p)", "**").replace("(d)", "*")
    if b.startswith("def2") and not b.startswith("def2-"):
        b = "def2-" + b[4:]
    return b


def method_string(engine: str, xc: str, basis: str, protocol: str) -> str:
    """Everything besides the species itself that changes the energy. The
    protocol distinguishes a single-embed xtb start (BDE species) from a
    conformer-search start (single_point), which can land in different
    minima."""
    return f"{engine.lower()}:{xc.lower()}/{normalize_basis(basis)}|xtb-gfn2|{protocol}"


def species_key(smiles: str, charge: int, multiplicity: int, method: str) -> tuple[str, str]:
    """(key, canonical_smiles). The key is a short stable hash so it can also
    name files/tasks."""
    can = canonical_smiles(smiles)
    raw = f"{can}|{int(charge)}|{int(multiplicity)}|{method}"
    return hashlib.sha1(raw.encode()).hexdigest()[:20], can


# --------------------------------------------------------------------------
# Bond typing
# --------------------------------------------------------------------------

_BOND_SYMBOL = {1.0: "-", 1.5: ":", 2.0: "=", 3.0: "#"}


def _fragment_signature(mol, ignore_electrons: bool) -> str:
    """Canonical SMILES of a fragment with explicit hydrogens - keeping them
    is what tells a C-H cut from an O-H cut. ignore_electrons drops radicals
    and charges (connectivity-only comparison for heterolysis)."""
    from rdkit import Chem

    m = Chem.Mol(mol)
    if ignore_electrons:
        for a in m.GetAtoms():
            a.SetNumRadicalElectrons(0)
            a.SetFormalCharge(0)
    return Chem.MolToSmiles(m, allHsExplicit=True)


def _signature_of_smiles(smiles: str, ignore_electrons: bool) -> str | None:
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    for a in mol.GetAtoms():
        a.SetNoImplicit(True)
    return _fragment_signature(mol, ignore_electrons)


def detect_broken_bond(parent_smiles: str, fragment_smiles: list[str]) -> str | None:
    """Which bond of the parent the fragments come from, as 'O-H', 'C-C',
    'C=O', ... (heavier element first). Tries every bond of the parent
    (hydrogens explicit): cut it, give each end the electron it keeps
    (homolysis), and compare the pieces with the planned fragments. If no
    homolytic match exists (e.g. heterolysis, or the planner put charges on
    the fragments), compares connectivity only. None if nothing matches or
    the case isn't a single-bond, two-fragment break."""
    from collections import Counter
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    if len(fragment_smiles) != 2:
        return None
    parent = Chem.MolFromSmiles(parent_smiles)
    if parent is None:
        return None
    parent = Chem.AddHs(parent)
    for a in parent.GetAtoms():
        a.SetNoImplicit(True)

    for ignore in (False, True):
        targets = [_signature_of_smiles(s, ignore) for s in fragment_smiles]
        if any(t is None for t in targets):
            return None
        want = Counter(targets)
        matches: set[str] = set()
        for bond in parent.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            rw = Chem.RWMol(parent)
            rw.RemoveBond(i, j)
            for idx in (i, j):
                atom = rw.GetAtomWithIdx(idx)
                atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() + 1)
            pieces = Chem.GetMolFrags(rw.GetMol(), asMols=True, sanitizeFrags=False)
            if len(pieces) != 2:
                continue  # ring bond: cutting it doesn't split the molecule
            got = Counter(_fragment_signature(p, ignore) for p in pieces)
            if got == want:
                a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
                if a1.GetAtomicNum() < a2.GetAtomicNum():
                    a1, a2 = a2, a1
                sym = _BOND_SYMBOL.get(bond.GetBondTypeAsDouble(), "~")
                matches.add(f"{a1.GetSymbol()}{sym}{a2.GetSymbol()}")
        if matches:
            # equivalent bonds (water's two O-H) agree; anything else is ambiguous
            return matches.pop() if len(matches) == 1 else None
    return None


# --------------------------------------------------------------------------
# Database
# --------------------------------------------------------------------------

def _new_claim_token() -> str:
    """host:pid:uuid - unique per claim, since every thread of one process
    shares host:pid."""
    return f"{socket.gethostname()}:{os.getpid()}:{uuid.uuid4().hex[:12]}"


class ResultsDB:
    def __init__(self, path: str, claim_ttl_s: float | None = None):
        self.path = os.path.abspath(path)
        self.claim_ttl_s = float(claim_ttl_s or os.environ.get("COMPCHEM_DB_CLAIM_TTL", 600))
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        with self._conn() as c:
            c.executescript(_SCHEMA)
            c.execute("INSERT OR IGNORE INTO meta VALUES ('schema_version', ?)", (str(SCHEMA_VERSION),))

    # One short-lived connection per operation: safe across the threads
    # Academy runs agent actions on, and across processes/nodes.
    @contextmanager
    def _conn(self) -> Iterator[sqlite3.Connection]:
        c = sqlite3.connect(self.path, timeout=60, isolation_level=None)
        try:
            c.row_factory = sqlite3.Row
            c.execute("PRAGMA journal_mode=DELETE")
            c.execute("PRAGMA busy_timeout=60000")
            yield c
        finally:
            c.close()

    @contextmanager
    def _write_txn(self) -> Iterator[sqlite3.Connection]:
        with self._conn() as c:
            c.execute("BEGIN IMMEDIATE")
            try:
                yield c
                c.execute("COMMIT")
            except BaseException:
                c.execute("ROLLBACK")
                raise

    # ---------------- species ----------------

    def get_species(self, key: str) -> dict | None:
        with self._conn() as c:
            row = c.execute("SELECT * FROM species WHERE key = ?", (key,)).fetchone()
        return dict(row) if row else None

    def lookup(self, key: str) -> dict | None:
        """A finished, successful entry or None. Does not count as a hit."""
        row = self.get_species(key)
        return row if row and row["status"] == "done" else None

    def claim(self, key: str, canonical: str, charge: int, multiplicity: int, method: str) -> tuple[str, Any]:
        """Atomically decide who computes `key`. Returns one of
        ('done', row)      - already computed: reuse it
        ('running', row)   - someone else is computing it right now: wait
        ('claimed', token) - the caller must compute it, keep the claim alive
                             with heartbeat(), and finish with
                             store_result()/store_failure() passing `token`
        A previously failed entry is re-claimed (failures aren't reused: the
        next attempt may have different inputs upstream, e.g. a re-plan)."""
        now = time.time()
        token = _new_claim_token()
        with self._write_txn() as c:
            row = c.execute("SELECT * FROM species WHERE key = ?", (key,)).fetchone()
            if row is not None:
                row = dict(row)
                if row["status"] == "done":
                    c.execute("UPDATE species SET hits = hits + 1 WHERE key = ?", (key,))
                    return "done", row
                if row["status"] == "running" and now - (row["claimed_at"] or 0) < self.claim_ttl_s:
                    return "running", row
            c.execute(
                """INSERT INTO species (key, canonical_smiles, charge, multiplicity, method,
                                        status, owner, run_id, claimed_at)
                   VALUES (?, ?, ?, ?, ?, 'running', ?, ?, ?)
                   ON CONFLICT(key) DO UPDATE SET status='running', owner=excluded.owner,
                       run_id=excluded.run_id, claimed_at=excluded.claimed_at, error=NULL""",
                (key, canonical, int(charge), int(multiplicity), method,
                 token, os.environ.get("COMPCHEM_RUN_ID"), now),
            )
        return "claimed", token

    def heartbeat(self, key: str, token: str) -> bool:
        """Refresh a live claim; False if it was taken over."""
        with self._write_txn() as c:
            cur = c.execute(
                "UPDATE species SET claimed_at=? WHERE key=? AND owner=? AND status='running'",
                (time.time(), key, token),
            )
            return cur.rowcount == 1

    def store_result(self, key: str, token: str, result: dict, wall_s: float) -> bool:
        """Complete our claim; False (nothing written) if it was taken over."""
        refinement = result.get("refinement") or {}
        with self._write_txn() as c:
            cur = c.execute(
                """UPDATE species SET status='done', energy_hartree=?, converged=?, result_json=?,
                       error=NULL, wall_s=?, finished_at=? WHERE key=? AND owner=? AND status='running'""",
                (refinement.get("energy_hartree"), int(bool(refinement.get("converged"))),
                 json.dumps(result, default=str), wall_s, time.time(), key, token),
            )
            return cur.rowcount == 1

    def store_failure(self, key: str, token: str, error: str, wall_s: float) -> bool:
        with self._write_txn() as c:
            cur = c.execute(
                """UPDATE species SET status='failed', error=?, wall_s=?, finished_at=?
                   WHERE key=? AND owner=? AND status='running'""",
                (error, wall_s, time.time(), key, token),
            )
            return cur.rowcount == 1

    def release_run_claims(self, run_id: str, reason: str) -> int:
        """Mark claims a run left 'running' (its workers died) as failed, so
        nothing waits on them for the full TTL. Called by the staged driver
        after its compute stage."""
        with self._write_txn() as c:
            cur = c.execute(
                "UPDATE species SET status='failed', error=?, finished_at=? WHERE run_id=? AND status='running'",
                (reason, time.time(), run_id),
            )
            return cur.rowcount

    def wait(self, key: str, poll_s: float = 5.0) -> dict | None:
        """Block until another worker's claim on `key` finishes. Returns the
        row once it is 'done' or 'failed'; None if the claim went stale (no
        heartbeat for a TTL - the caller should then claim it itself). A live
        claimant keeps heartbeating, so a long compute is waited out."""
        while True:
            row = self.get_species(key)
            if row is None:
                return None
            if row["status"] in ("done", "failed"):
                if row["status"] == "done":
                    with self._write_txn() as c:
                        c.execute("UPDATE species SET hits = hits + 1 WHERE key = ?", (key,))
                return row
            if time.time() - (row["claimed_at"] or 0) >= self.claim_ttl_s:
                return None
            time.sleep(poll_s)

    # ---------------- bde ----------------

    def record_bde(self, *, parent_key: str, fragment_keys: list[str], parent_smiles: str,
                   fragment_smiles: list[str], bond_type: str | None, method: str,
                   bde_hartree: float | None, bde_kcal_mol: float | None, qc_flag: str | None,
                   source_question: str | None) -> int:
        """Insert or refresh the BDE for this exact (parent, fragments,
        method) combination - recomputing an identical BDE from cached
        species yields the same number, so one row per combination."""
        keys = sorted(fragment_keys)
        smis = sorted(fragment_smiles)
        with self._write_txn() as c:
            c.execute(
                """INSERT INTO bde (parent_key, fragment_keys, parent_smiles, fragment_smiles,
                                    bond_type, method, bde_hartree, bde_kcal_mol, qc_flag,
                                    source_question, run_id, created_at)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                   ON CONFLICT(parent_key, fragment_keys) DO UPDATE SET
                       bde_hartree=excluded.bde_hartree, bde_kcal_mol=excluded.bde_kcal_mol,
                       qc_flag=excluded.qc_flag, bond_type=excluded.bond_type""",
                (parent_key, json.dumps(keys), parent_smiles, json.dumps(smis), bond_type, method,
                 bde_hartree, bde_kcal_mol, qc_flag, source_question,
                 os.environ.get("COMPCHEM_RUN_ID"), time.time()),
            )
            row = c.execute(
                "SELECT id FROM bde WHERE parent_key = ? AND fragment_keys = ?",
                (parent_key, json.dumps(keys)),
            ).fetchone()
        return int(row["id"])

    def query_bdes(self, *, bond_types: list[str] | None = None, parent_smiles: list[str] | None = None,
                   qc_flag: str | None = "PASS", limit: int = 50) -> list[dict]:
        # bond-type and parent filters are alternatives (OR); qc_flag always applies
        filters, args = [], []
        if bond_types:
            filters.append(f"bond_type IN ({','.join('?' * len(bond_types))})")
            args += bond_types
        if parent_smiles:
            filters.append(f"parent_smiles IN ({','.join('?' * len(parent_smiles))})")
            args += parent_smiles
        clauses = ["(" + " OR ".join(filters) + ")"] if filters else []
        if qc_flag:
            clauses.append("qc_flag = ?")
            args.append(qc_flag)
        sql = "SELECT * FROM bde" + (" WHERE " + " AND ".join(clauses) if clauses else "")
        sql += " ORDER BY created_at DESC LIMIT ?"
        with self._conn() as c:
            rows = c.execute(sql, (*args, limit)).fetchall()
        return [dict(r) for r in rows]

    def species_by_smiles(self, canonical: list[str]) -> list[dict]:
        if not canonical:
            return []
        with self._conn() as c:
            rows = c.execute(
                f"""SELECT key, canonical_smiles, charge, multiplicity, method, energy_hartree, hits
                    FROM species WHERE status='done' AND canonical_smiles IN ({','.join('?' * len(canonical))})""",
                canonical,
            ).fetchall()
        return [dict(r) for r in rows]

    def stats(self) -> dict:
        with self._conn() as c:
            s = c.execute(
                """SELECT COUNT(*) n, SUM(status='done') done, SUM(status='failed') failed,
                          SUM(status='running') running, COALESCE(SUM(hits),0) hits,
                          COALESCE(SUM(hits * wall_s),0) saved_s
                   FROM species"""
            ).fetchone()
            b = c.execute("SELECT COUNT(*) n FROM bde").fetchone()
        return {**dict(s), "bde_rows": b["n"]}


_DB_CACHE: dict[str, ResultsDB] = {}


def get_db() -> ResultsDB | None:
    """The process-wide ResultsDB named by $COMPCHEM_DB, or None (caching off)."""
    path = os.environ.get("COMPCHEM_DB")
    if not path:
        return None
    if path not in _DB_CACHE:
        _DB_CACHE[path] = ResultsDB(path)
    return _DB_CACHE[path]


def get_or_compute_species(
    *, smiles: str, charge: int, multiplicity: int, engine: str, xc: str, basis: str,
    protocol: str, compute: Callable[[], dict], db: ResultsDB | None = None,
) -> dict:
    """Return a compute_species_* style result for this species, reusing the
    database when possible. `compute()` must follow the never-raises
    contract of skills._common (status 'ok' / 'failed').

    The returned dict carries a 'cache' entry: {'key', 'outcome', 'wait_s',
    'saved_wall_s'} where outcome is 'hit' (already in the DB), 'waited'
    (another worker computed it while we waited), 'computed' (we did), or
    'disabled' (no database configured)."""
    db = db if db is not None else get_db()
    if db is None:
        result = compute()
        result["cache"] = {"key": None, "outcome": "disabled", "wait_s": 0.0, "saved_wall_s": 0.0}
        return result

    method = method_string(engine, xc, basis, protocol)
    key, can = species_key(smiles, charge, multiplicity, method)

    def _from_row(row: dict, outcome: str, wait_s: float) -> dict:
        if row["status"] == "done":
            result = json.loads(row["result_json"])
            saved = row["wall_s"] or 0.0
        else:  # failed while we waited: report it rather than recompute it
            result = {"smiles": smiles, "charge": charge, "multiplicity": multiplicity,
                      "status": "failed", "error": row["error"], "refinement_engine": engine,
                      "xtb": None, "refinement": None}
            saved = row["wall_s"] or 0.0
        result["smiles"] = smiles  # keep the planner's spelling in reports
        result["cache"] = {"key": key, "outcome": outcome, "wait_s": round(wait_s, 2),
                           "saved_wall_s": round(saved, 2)}
        return result

    t0 = time.time()
    for _ in range(3):  # at most: stale claim -> re-claim -> compute
        state, got = db.claim(key, can, charge, multiplicity, method)
        if state == "done":
            return _from_row(got, "hit", 0.0)
        if state == "running":
            row = db.wait(key)
            if row is not None:
                return _from_row(row, "waited", time.time() - t0)
            continue
        # claimed: we compute, heartbeating so waiters know we're alive
        token, c0 = got, time.time()
        stop = threading.Event()

        def _beat():
            while not stop.wait(db.claim_ttl_s / 4):
                try:
                    if not db.heartbeat(key, token):
                        return
                except Exception:
                    pass  # transient DB error: the next beat retries

        beater = threading.Thread(target=_beat, daemon=True)
        beater.start()
        try:
            result = compute()
        except BaseException as exc:  # compute() shouldn't raise, but never leave a claim dangling
            stop.set()
            db.store_failure(key, token, f"{type(exc).__name__}: {exc}", time.time() - c0)
            raise
        finally:
            stop.set()
        wall = time.time() - c0
        if result.get("status") == "ok":
            db.store_result(key, token, result, wall)
        else:
            db.store_failure(key, token, result.get("error") or "unknown error", wall)
        result["cache"] = {"key": key, "outcome": "computed", "wait_s": round(c0 - t0, 2),
                           "saved_wall_s": 0.0, "wall_s": round(wall, 2)}
        return result

    result = compute()  # repeated stale claims: give up on coordination, just compute
    result["cache"] = {"key": key, "outcome": "computed", "wait_s": round(time.time() - t0, 2),
                       "saved_wall_s": 0.0}
    return result

"""
Structure rendering via PyMOL open-source (Phase 7 - architecture plan
Sec.17), run as an isolated subprocess (revised during Phase 9 debugging -
architecture plan Sec.18).

Renders one species' 3D structure - the same symbols/positions shape every
skill's xtb result already carries (src/skills/_common.py) - to a PNG, for
embedding in ReporterAgent's markdown reports.

Originally this called PyMOL in-process (first via the pymol2 multi-
instance API, which failed outright on this repo's conda-forge build - see
Sec.17 - then via PyMOL's classic `pymol.finish_launching()` singleton
pattern, which DID render correctly but then caused the whole Python
process to silently exit (code 0, no traceback, no further output)
sometime after `cmd.png()` returned - discovered while testing Phase 9's
provenance capture, whose real symptom (an always-empty Flowcept buffer)
turned out to be a consequence of ReporterAgent never actually finishing,
not a Flowcept bug at all. PyMOL's classic embedding pattern was designed
assuming PyMOL owns the whole process (a standalone script that calls
`cmd.quit()` when done); something in its process-level teardown apparently
doesn't coexist safely with being called partway through a larger async,
multi-agent application's own lifecycle - the exact mechanism was not
pinned down further, since isolating it in a subprocess sidesteps the
problem regardless of the precise cause.

So: this module now writes the geometry to a temporary .xyz file and a
temporary PyMOL script (.pml), then runs the real `pymol` command-line
program (installed by conda-forge's `pymol-open-source` package) as a
throwaway subprocess via `subprocess.run(...)`. Whatever PyMOL does to ITS
process on exit no longer has any way to affect this one. The PNG this
process writes to `out_path` is picked up as this function's actual result.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile

from storage.flowcept_store import flowcept_task


@flowcept_task
def render_structure_png(
    symbols: list[str],
    positions: list[list[float]],
    out_path: str,
    label: str | None = None,
    width: int = 640,
    height: int = 480,
    style: str = "sticks",
    timeout_s: int = 60,
) -> dict:
    """Render one species' geometry to a PNG via a standalone `pymol -cq`
    subprocess.

    Returns a JSON-safe dict: {"status": "ok"/"failed", "path": out_path or
    None, "error": str or None} - the same graceful-failure contract used
    throughout this repo (architecture plan Sec.7): a missing PyMOL install,
    a malformed geometry, or any other rendering hiccup comes back as an
    explicit failure rather than raising, since a missing structure image
    should never take down the whole report.
    """
    if shutil.which("pymol") is None:
        return {
            "status": "failed", "path": None,
            "error": "pymol executable not found on PATH (conda install -c conda-forge pymol-open-source)",
        }

    tmp_xyz_path = None
    tmp_script_path = None
    try:
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

        xyz_lines = [str(len(symbols)), label or ""]
        for sym, coord in zip(symbols, positions):
            x, y, z = coord
            xyz_lines.append(f"{sym} {float(x):.8f} {float(y):.8f} {float(z):.8f}")
        xyz_text = "\n".join(xyz_lines) + "\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as tmp_xyz:
            tmp_xyz.write(xyz_text)
            tmp_xyz_path = tmp_xyz.name

        abs_out_path = os.path.abspath(out_path)
        script_lines = [
            f"load {tmp_xyz_path}, species",
            "hide everything",
            f"show_as {style}, species",
            "util.cbaw species",
            "bg_color white",
            "set ray_opaque_background, 1",
            "orient species",
            "zoom species, buffer=2",
            f"png {abs_out_path}, width={width}, height={height}, dpi=150, ray=1",
            "quit",
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pml", delete=False) as tmp_script:
            tmp_script.write("\n".join(script_lines) + "\n")
            tmp_script_path = tmp_script.name

        result = subprocess.run(
            ["pymol", "-cq", tmp_script_path],
            capture_output=True, text=True, timeout=timeout_s,
        )

        if not os.path.exists(abs_out_path):
            stderr_tail = (result.stderr or "")[-500:]
            return {
                "status": "failed", "path": None,
                "error": (
                    f"pymol subprocess did not produce an output file "
                    f"(returncode={result.returncode}): {stderr_tail}"
                ),
            }
        return {"status": "ok", "path": out_path, "error": None}
    except subprocess.TimeoutExpired:
        return {"status": "failed", "path": None, "error": f"pymol subprocess timed out after {timeout_s}s"}
    except Exception as exc:
        return {"status": "failed", "path": None, "error": f"{type(exc).__name__}: {exc}"}
    finally:
        for p in (tmp_xyz_path, tmp_script_path):
            if p is not None and os.path.exists(p):
                os.remove(p)

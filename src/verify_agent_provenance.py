"""
Phase 9 verification script (architecture plan Sec.18): read back the local
Flowcept JSONL buffer after running a campaign and show the PROV-AGENT-
tagged agent-decision records (subtype 'agent_tool' or
'ai_model_invocation') captured alongside all the plain computational-task
records that were already being captured before this phase.

Run AFTER a campaign, from src/:

    python -m agents.run_agentic_campaign
    python verify_agent_provenance.py

or the LLM-backed --question variant, to also see one or more
'ai_model_invocation' records for the planner's actual Claude API call(s)
(two of them if the bounded self-correction retry - architecture plan
Sec.17/Sec.18 - kicked in):

    python -m agents.run_agentic_campaign --question "..." --model claude-haiku-4-5-20251001
    python verify_agent_provenance.py
"""
from __future__ import annotations

from storage.flowcept_store import read_captured_tasks


def main() -> None:
    df = read_captured_tasks()
    if df is None or df.empty:
        print("No captured tasks found - is flowcept_buffer.jsonl present in this directory?")
        return

    print(f"Total captured task records: {len(df)}")

    if "subtype" not in df.columns:
        print(
            "No 'subtype' column in the captured records - this flowcept "
            "install may predate agent-decision provenance (needs >= "
            "1.0.2). All captured records' columns:"
        )
        print(sorted(df.columns))
        return

    agent_df = df[df["subtype"].isin(["agent_tool", "ai_model_invocation"])]
    if agent_df.empty:
        print(
            "No PROV-AGENT-tagged records found (subtype agent_tool/"
            "ai_model_invocation) - agent-decision provenance may not have "
            "been captured this run. All subtypes seen: "
            f"{sorted(df['subtype'].dropna().unique().tolist())}"
        )
        return

    cols = [c for c in ("agent_id", "subtype", "activity_id", "started_at", "status") if c in agent_df.columns]
    print(f"\nFound {len(agent_df)} PROV-AGENT-tagged agent-decision records:\n")
    print(agent_df[cols].to_string(index=False))


if __name__ == "__main__":
    main()

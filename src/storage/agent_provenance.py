"""
PROV-AGENT-style agent-decision provenance (Phase 9 - architecture plan
Sec.18).

Flowcept's plain @flowcept_task (storage/flowcept_store.py) captures
COMPUTATIONAL task provenance - which function ran, on what inputs,
producing what outputs. PROV-AGENT (arXiv:2508.02866, ORNL/ANL) is a
purpose-built extension of that same W3C-PROV-derived model specifically
for AI-AGENT DECISIONS: it types Flowcept records with two additional
subtypes - AgentTool (one tool/action an AI agent invoked) and
AIModelInvocation (one LLM prompt->response cycle) - so a provenance query
can distinguish "the ExecutorAgent ran this skill" from "the PlannerAgent
asked Claude this question and got this answer", not just "some function
ran".

Flowcept ships this as flowcept.instrumentation.task_capture.FlowceptTask
(a per-call provenance record you drive explicitly with begin/end, not a
decorator) plus the PROV_AGENT enum (flowcept.commons.vocabulary.PROV_AGENT).
This module wraps both in two small context managers used directly inside
each Academy agent's @action methods (agents/*.py) and inside the LLM
planner's actual API call (llm/planner_llm.py) - the natural attribution
points for "which agent did this" and "what did the model actually see and
say".

Deliberately NOT using Flowcept's own `@agent_flowcept_task` decorator
(researched before choosing this design - see architecture plan Sec.18):
that decorator (a) requires langchain_core importable even though its own
logic never touches LangChain, a dependency this repo doesn't otherwise
need; (b) is implemented sync-only (`result = func(*args, **kwargs)`, no
await/asyncio.iscoroutinefunction check) while every Academy @action method
here is `async def`; and, decisively, (c) reads `agent_id` from a single
process-global class attribute (`BaseAgentContextManager.agent_id`), so
every task tagged through it would collapse onto the SAME agent identity -
useless for a genuinely multi-agent pipeline where Planner/Executor/Critic/
Reporter each need their own distinct identity in the provenance graph.
FlowceptTask takes `agent_id` as a plain per-call keyword argument instead,
which is exactly what a 4-agent pipeline needs, has zero LangChain
dependency, and works identically from sync or async code since it's a
plain object with explicit begin/end, not a decorator.

Same never-crash-the-pipeline philosophy as storage/flowcept_store.py: if
flowcept isn't installed, or the installed version predates these agent-
provenance primitives, both context managers below no-op cleanly rather
than raising - a missing/older Flowcept install should never take down a
real calculation.

Requires: flowcept >= 1.0.2 (confirmed to already include FlowceptTask and
PROV_AGENT.{AI_MODEL_INVOCATION,AGENT_TOOL} - this is not a new pip
dependency, just a version floor on the one already in environment.yml).
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Any, Iterator

try:
    from flowcept import FlowceptTask as _FlowceptTask
    from flowcept.commons.vocabulary import PROV_AGENT as _PROV_AGENT
    _AGENT_PROVENANCE_AVAILABLE = True
except ImportError:
    _FlowceptTask = None
    _PROV_AGENT = None
    _AGENT_PROVENANCE_AVAILABLE = False


class _NullTask:
    """Stand-in for FlowceptTask when agent-provenance capture isn't
    available, so callers can always write `t.end(generated=...)`
    unconditionally without an availability check at every call site."""

    def end(self, *args: Any, **kwargs: Any) -> None:
        pass


@contextmanager
def agent_tool_task(
    agent_id: str, activity_id: str, used: dict[str, Any] | None = None
) -> Iterator[Any]:
    """Record one Academy agent action as a PROV-AGENT 'AgentTool' activity
    - one tool execution invoked by an AI agent (see module docstring).

    Usage:
        with agent_tool_task("ExecutorAgent", "run_analysis", used={...}) as t:
            result = await self.agent_run_sync(...)
            t.end(generated={...})
        return result

    On an uncaught exception inside the `with` block, records the task's
    stderr with the error before re-raising - the same graceful-failure-
    visibility principle used throughout this repo (architecture plan
    Sec.7), applied to the provenance record rather than the calculation
    itself (the calculation's own status/error fields are untouched).
    """
    if not _AGENT_PROVENANCE_AVAILABLE:
        yield _NullTask()
        return
    task = _FlowceptTask(
        agent_id=agent_id, subtype=_PROV_AGENT.AGENT_TOOL,
        activity_id=activity_id, used=used or {},
    )
    try:
        yield task
    except Exception as exc:
        task.end(stderr=f"{type(exc).__name__}: {exc}")
        raise


@contextmanager
def agent_llm_invocation_task(
    agent_id: str, activity_id: str, used: dict[str, Any] | None = None
) -> Iterator[Any]:
    """Record one LLM prompt->response cycle as a PROV-AGENT
    'AIModelInvocation' activity (see module docstring). Same usage and
    failure-visibility pattern as agent_tool_task."""
    if not _AGENT_PROVENANCE_AVAILABLE:
        yield _NullTask()
        return
    task = _FlowceptTask(
        agent_id=agent_id, subtype=_PROV_AGENT.AI_MODEL_INVOCATION,
        activity_id=activity_id, used=used or {},
    )
    try:
        yield task
    except Exception as exc:
        task.end(stderr=f"{type(exc).__name__}: {exc}")
        raise

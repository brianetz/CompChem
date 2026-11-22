"""Learning from the plan board's mistakes.

task_gen   BDE questions with RDKit-computed answers (parent, fragments, spins)
oracle     grades a plan against those answers, independently of the board's own checks
episodes   exports Flowcept planning records into attempt-level episodes
memory     lessons from past mistakes, retrieved for the Planner's prompt
"""

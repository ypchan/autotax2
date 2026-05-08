"""Consensus placeholders."""

from __future__ import annotations


def choose_near_best_hits(scores: dict[str, float], delta: float) -> list[str]:
    """Return hit IDs whose scores are within ``delta`` of the best score."""
    if not scores:
        return []
    best = max(scores.values())
    return sorted(hit_id for hit_id, score in scores.items() if best - score <= delta)

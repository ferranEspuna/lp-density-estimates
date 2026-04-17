#!/usr/bin/env python3
"""
Reproduce two-parameter LP bounds documented in latex/5_experiments_and_bounds.tex.

Run: python3 experiments_two_parameter.py
"""

from __future__ import annotations

import itertools
import time

from two_parameter_lp import evaluate_two_parameter_bound, make_direction_family, pentad_family


def plus_minus_four(m):
    """u, v, (u+v)/m, (u-v)/m with shift 0 on each."""
    return (
        make_direction_family((1, 0), (0,), prefix="u")
        + make_direction_family((0, 1), (0,), prefix="v")
        + make_direction_family((1, 1), (0,), dilation=m, prefix="uv")
        + make_direction_family((1, -1), (0,), dilation=m, prefix="umv")
    )


def search_best_pentad_from_pool(m):
    """
    Fix u, v and choose three extra directions from a fixed integer pool.
    Returns (best_bound, best_index_tuple).
    """
    dirs = [
        ((1, 0), 1),
        ((0, 1), 1),
        ((1, 1), m),
        ((1, -1), m),
        ((2, 1), m),
        ((1, 2), m),
        ((2, -1), m),
        ((1, -2), m),
        ((3, 1), m),
        ((1, 3), m),
        ((3, -1), m),
        ((1, -3), m),
        ((2, 3), m),
        ((3, 2), m),
        ((4, 1), m),
        ((1, 4), m),
    ]

    def family_from_indices(indices):
        forms = []
        for i in indices:
            (a, b), dil = dirs[i]
            forms += make_direction_family((a, b), (0,), dilation=dil, prefix=str(i))
        return forms

    rest = list(range(2, len(dirs)))
    best = 1.0
    best_idx = None
    t0 = time.time()
    for comb in itertools.combinations(rest, 3):
        idxs = (0, 1) + comb
        val = evaluate_two_parameter_bound(m, family_from_indices(idxs))
        if val is None:
            continue
        if val < best - 1e-12:
            best = val
            best_idx = idxs
    return best, best_idx, time.time() - t0


def main():
    print("Two-parameter LP — reproducible benchmarks\n")

    m = 3
    b_quad = evaluate_two_parameter_bound(m, plus_minus_four(m))
    b_pent = evaluate_two_parameter_bound(m, pentad_family(m))
    print(f"m={m}  four directions (u,v,(u+v)/m,(u-v)/m):  {b_quad}")
    print(f"m={m}  pentad_family(m):                         {b_pent}")
    if b_pent is not None:
        print(f"       (8/15 = {8/15})")

    best, idxs, elapsed = search_best_pentad_from_pool(m)
    print(f"\nExhaustive 5-form search (u,v + 3 dirs from pool), m={m}:")
    print(f"  best bound: {best}")
    print(f"  index tuple: {idxs}")
    print(f"  time: {elapsed:.2f}s")
    print("  (indices 6,7 correspond to (2u-v)/m and (u-2v)/m when m=3.)")

    print("\nPentad family for other m (same code path):")
    for m2 in (4, 5, 7):
        val = evaluate_two_parameter_bound(m2, pentad_family(m2))
        print(f"  m={m2}: {val}")


if __name__ == "__main__":
    main()

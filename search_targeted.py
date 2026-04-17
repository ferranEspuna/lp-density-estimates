"""
Targeted search for the two-parameter LP. Includes:

1. Symmetric enlargement of ``pentad_family``.
2. Exhaustive small-pool size-5/6/7 searches with progress printed live.
3. Search over nonzero shifts for the pentad.

Run: python3 search_targeted.py
"""

from __future__ import annotations

import itertools
import time
from fractions import Fraction

from two_parameter_lp import (
    evaluate_two_parameter_bound,
    linear_form_2d,
    make_direction_family,
    pentad_family,
)


def make_family(m: int, directions: list[tuple[tuple[int, int], int]]) -> list[linear_form_2d]:
    forms: list[linear_form_2d] = []
    for idx, ((a, b), dil) in enumerate(directions):
        forms += make_direction_family((a, b), (0,), dilation=dil, prefix=f"d{idx}")
    return forms


def pentad_plus(m: int, extra: list[tuple[tuple[int, int], int]]) -> list[linear_form_2d]:
    base = [
        ((1, 0), 1),
        ((0, 1), 1),
        ((1, 1), m),
        ((2, -1), m),
        ((1, -2), m),
    ]
    return make_family(m, base + extra)


def experiment_symmetric_pentad(m: int):
    """Augment pentad with symmetric sign-flipped partners."""
    print(f"\n--- Symmetric augmentations of pentad_family({m}) ---")
    candidates = [
        ("(u+2v)/m", [((1, 2), m)]),
        ("(2u+v)/m", [((2, 1), m)]),
        ("(u+2v)/m and (2u+v)/m", [((1, 2), m), ((2, 1), m)]),
        ("(u-v)/m", [((1, -1), m)]),
        ("(u+v), (u-v) /m", [((1, 1), m), ((1, -1), m)]),  # already got uv
        ("(3u-v)/m, (u-3v)/m", [((3, -1), m), ((1, -3), m)]),
        ("(u+2v)/m,(2u+v)/m,(u-v)/m", [((1, 2), m), ((2, 1), m), ((1, -1), m)]),
        ("full 8-form: pentad + (u-v)/m + (u+2v)/m + (2u+v)/m",
         [((1, -1), m), ((1, 2), m), ((2, 1), m)]),
    ]
    for label, extra in candidates:
        forms = pentad_plus(m, extra)
        if len(forms) > 10:
            print(f"  [skip] {label}: too many forms ({len(forms)})")
            continue
        t0 = time.time()
        bound = evaluate_two_parameter_bound(m, forms)
        elapsed = time.time() - t0
        print(f"  n={len(forms)} {label:60s}  bound={bound}  t={elapsed:.2f}s")


def experiment_small_pools(m: int):
    """Exhaustive search over compact direction pools, sizes 5..7."""
    pool_integers = [
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
        ((3, 2), m),
        ((2, 3), m),
        ((3, -2), m),
        ((2, -3), m),
    ]
    base = [((1, 0), 1), ((0, 1), 1)]
    for size in (3, 4, 5):
        print(f"\n--- m={m}, size = 2+{size} over pool of {len(pool_integers)} extras ---")
        best = 1.0
        best_idx = None
        start = time.time()
        count = 0
        total = sum(1 for _ in itertools.combinations(range(len(pool_integers)), size))
        for combo in itertools.combinations(range(len(pool_integers)), size):
            chosen = base + [pool_integers[i] for i in combo]
            forms = make_family(m, chosen)
            bound = evaluate_two_parameter_bound(m, forms)
            count += 1
            if bound is None:
                continue
            if bound < best - 1e-12:
                best = bound
                best_idx = combo
                elapsed = time.time() - start
                dirs_labels = [f"({a},{b})/{d}" for (a, b), d in chosen]
                print(f"  [{count}/{total} t={elapsed:.1f}s] new best {best:.6f}  idx={combo}  dirs={dirs_labels}", flush=True)
        print(f"  done: best {best:.6f}  idx={best_idx}  tried={count}/{total}  t={time.time() - start:.1f}s")


def experiment_shifts_on_pentad(m: int):
    """Try shift pools on the pentad."""
    print(f"\n--- Shift pools on pentad_family({m}) ---")

    shift_options = [(0,), (0, 1), (0, 1, 2)]
    directions = [
        ((1, 0), 1),
        ((0, 1), 1),
        ((1, 1), m),
        ((2, -1), m),
        ((1, -2), m),
    ]

    best = 1.0
    best_case = None
    start = time.time()
    tried = 0
    for shift_combo in itertools.product(shift_options, repeat=len(directions)):
        forms = []
        for (direction, dil), shifts in zip(directions, shift_combo):
            forms += make_direction_family(direction, shifts, dilation=dil, prefix=f"{direction}")
        if len(forms) > 9:
            continue
        bound = evaluate_two_parameter_bound(m, forms)
        tried += 1
        if bound is None:
            continue
        if bound < best - 1e-12:
            best = bound
            best_case = shift_combo
            elapsed = time.time() - start
            print(f"  [t={elapsed:.1f}s tried={tried}] new best {best:.6f} shifts={shift_combo}", flush=True)
    print(f"  done: best {best:.6f} case={best_case} tried={tried} t={time.time() - start:.1f}s")


def main():
    experiment_symmetric_pentad(3)
    experiment_small_pools(3)
    experiment_shifts_on_pentad(3)


if __name__ == "__main__":
    main()

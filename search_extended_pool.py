"""
Extended exhaustive searches over direction pools of arbitrary size, for the
two-parameter LP.

Main functions:

- ``pool_directions(max_abs=1, include_scale=True)``: generate a pool of 2D
  integer directions and the canonical dilation factor ``m`` for each.
- ``search_size(m, extra_pool, size, base=('u','v'), verbose=False)``: fix the
  base forms ``u`` and ``v`` and choose ``size`` additional forms from
  ``extra_pool``; run every combination and report the best LP upper bound.
- ``search_with_shifts(m, directions, shift_pool, verbose=False)``: given a list
  of directions (integer pairs with a canonical dilation ``m``), pick shifts
  for each direction from ``shift_pool`` and solve the resulting LP.

All routines short-circuit whenever ``evaluate_two_parameter_bound`` returns
``None`` (solver failure), so they never silently record fake bounds.
"""

from __future__ import annotations

import itertools
import time
from fractions import Fraction
from typing import Iterable

from two_parameter_lp import (
    evaluate_two_parameter_bound,
    linear_form_2d,
    make_direction_family,
)


Direction = tuple[tuple[int, int], int]


def pool_directions(m: int, max_abs: int = 2, include_plain: bool = True) -> list[Direction]:
    """Return a de-duplicated list of ``((a,b), dilation)`` with ``|a|,|b|<=max_abs``.

    Directions with ``gcd(a,b,dil)`` nontrivial are kept once in their reduced
    form. ``(1,0)`` and ``(0,1)`` (with ``dil=1``) are always present.
    """

    seen = set()
    result: list[Direction] = []

    def add(a: int, b: int, dil: int) -> None:
        g = Fraction(a, dil), Fraction(b, dil)
        key = (g[0], g[1])
        if key in seen:
            return
        seen.add(key)
        result.append(((a, b), dil))

    if include_plain:
        add(1, 0, 1)
        add(0, 1, 1)

    for a in range(-max_abs, max_abs + 1):
        for b in range(-max_abs, max_abs + 1):
            if a == 0 and b == 0:
                continue
            add(a, b, m)

    return result


def forms_from_indices(directions: list[Direction], indices: Iterable[int], shifts_per: int = 1) -> list[linear_form_2d]:
    """Build the witness family by taking ``shifts_per`` copies with shift 0 (default)."""

    shifts = tuple(range(shifts_per))
    forms: list[linear_form_2d] = []
    for idx in indices:
        (a, b), dil = directions[idx]
        forms += make_direction_family((a, b), shifts, dilation=dil, prefix=f"dir{idx}")
    return forms


def search_size(
    m: int,
    extra_pool: list[Direction],
    size: int,
    fixed_prefix: tuple[Direction, ...] = (((1, 0), 1), ((0, 1), 1)),
    verbose: bool = False,
):
    """Return ``(best_bound, best_indices, elapsed_seconds)``.

    The witness family is ``fixed_prefix`` followed by ``size`` directions drawn
    from ``extra_pool`` (each with shift 0).
    """

    best = 1.0
    best_indices: tuple[int, ...] | None = None

    start = time.time()
    all_dirs = list(fixed_prefix) + list(extra_pool)
    fixed_len = len(fixed_prefix)
    extras = list(range(fixed_len, len(all_dirs)))
    tried = 0

    for combo in itertools.combinations(extras, size):
        idxs = tuple(range(fixed_len)) + combo
        forms = forms_from_indices(all_dirs, idxs)
        bound = evaluate_two_parameter_bound(m, forms)
        tried += 1
        if bound is None:
            continue
        if bound < best - 1e-12:
            best = bound
            best_indices = idxs
            if verbose:
                elapsed = time.time() - start
                print(f"[{elapsed:6.1f}s tried={tried}] m={m} size={len(idxs)} bound={bound:.6f} idx={idxs}", flush=True)

    return best, best_indices, time.time() - start


def search_with_shifts(
    m: int,
    directions: list[Direction],
    shift_options,
    verbose: bool = False,
):
    """Try every assignment of shifts from ``shift_options`` to ``directions``.

    ``shift_options`` is a list of tuples of integer shifts; one entry per
    direction.
    """

    best = 1.0
    best_assignment = None

    start = time.time()
    tried = 0
    for combo in itertools.product(shift_options, repeat=len(directions)):
        forms: list[linear_form_2d] = []
        for (direction, dil), shifts in zip(directions, combo):
            forms += make_direction_family(direction, shifts, dilation=dil, prefix=f"dir_{direction[0]}_{direction[1]}")
        if not forms:
            continue
        bound = evaluate_two_parameter_bound(m, forms)
        tried += 1
        if bound is None:
            continue
        if bound < best - 1e-12:
            best = bound
            best_assignment = combo
            if verbose:
                elapsed = time.time() - start
                print(f"[{elapsed:6.1f}s tried={tried}] m={m} bound={bound:.6f} shifts={combo}", flush=True)

    return best, best_assignment, time.time() - start


if __name__ == "__main__":
    print("Extended pool search, m=3\n")

    pool = pool_directions(3, max_abs=2)
    print("Directions in pool:")
    for idx, ((a, b), dil) in enumerate(pool):
        print(f"  {idx:2d}: ({a:2d},{b:2d}) / {dil}")

    print("\nSweeping size = 5 .. 8 over the full pool (u, v fixed)...\n")
    for size in (3, 4, 5, 6):  # 3..6 extras in addition to u, v
        bound, idxs, elapsed = search_size(
            3,
            [d for d in pool if d != ((1, 0), 1) and d != ((0, 1), 1)],
            size,
            verbose=True,
        )
        print(f"  size={size + 2}: best={bound:.6f}  idx={idxs}  time={elapsed:.1f}s\n", flush=True)

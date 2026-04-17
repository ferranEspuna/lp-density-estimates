"""Reproduce the three-parameter LP benchmarks.

Run: ``python3 experiments_three_parameter.py``

Each call prints the LP relaxation upper bound on ``delta(A)`` for a set
``A subset F_p`` avoiding ``x+y=mz``, under the assumption that ``m >= 3``.

These bounds are derived with ``three_parameter_lp.py``; they are valid LP
relaxation bounds (i.e.\ upper bounds on the LP optimum, which is an upper
bound on ``delta(A)`` itself).
"""

from __future__ import annotations

import time
from fractions import Fraction

from three_parameter_lp import (
    evaluate_three_parameter_bound,
    linear_form_3d,
)


def F(a, b, c, d=0, label=""):
    return linear_form_3d(a, b, c, d, label)


def basic_pair_family(m: int):
    """``u, v, w, (u+v)/m, (u+w)/m, (v+w)/m``."""
    return [
        F(1, 0, 0, 0, "u"),
        F(0, 1, 0, 0, "v"),
        F(0, 0, 1, 0, "w"),
        F(Fraction(1, m), Fraction(1, m), 0, 0, "u+v/m"),
        F(Fraction(1, m), 0, Fraction(1, m), 0, "u+w/m"),
        F(0, Fraction(1, m), Fraction(1, m), 0, "v+w/m"),
    ]


def signed_pair_family(m: int):
    """``basic_pair_family`` augmented with the three signed differences."""
    return basic_pair_family(m) + [
        F(Fraction(1, m), Fraction(-1, m), 0, 0, "u-v/m"),
        F(Fraction(1, m), 0, Fraction(-1, m), 0, "u-w/m"),
        F(0, Fraction(1, m), Fraction(-1, m), 0, "v-w/m"),
    ]


def symmetric_12_family(m: int):
    """``signed_pair_family`` augmented with the three ``(2a-b-c)/m`` forms."""
    return signed_pair_family(m) + [
        F(Fraction(2, m), Fraction(-1, m), Fraction(-1, m), 0, "2u-v-w/m"),
        F(Fraction(-1, m), Fraction(2, m), Fraction(-1, m), 0, "-u+2v-w/m"),
        F(Fraction(-1, m), Fraction(-1, m), Fraction(2, m), 0, "-u-v+2w/m"),
    ]


FAMILIES = {
    "6-form basic_pair_family (u,v,w + pairwise /m)": basic_pair_family,
    "9-form signed_pair_family (+ signed differences)": signed_pair_family,
    "12-form symmetric_12_family (+ (2a-b-c)/m)": symmetric_12_family,
}


def main():
    print("Three-parameter LP benchmarks\n")
    print(f"{'family':55s}{'m':>4}{'bound':>20}")
    for label, family in FAMILIES.items():
        for m in (3, 4, 5, 7):
            t0 = time.time()
            bound = evaluate_three_parameter_bound(m, family(m))
            t = time.time() - t0
            print(f"{label:55s}{m:>4}{bound:>20}   t={t:.1f}s")
        print()


if __name__ == "__main__":
    main()

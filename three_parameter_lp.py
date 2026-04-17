"""
Three-parameter linear-form LP for ``x + y = m z``.

This generalizes ``two_parameter_lp.py`` by letting each witness form depend
on three (rather than two) random parameters:

    L_r(u, v, w) = a_r * u + b_r * v + c_r * w + d_r.

Using more parameters enlarges the affine-symmetry group from
``GL_2(F_p) ltimes F_p^2`` to ``GL_3(F_p) ltimes F_p^3`` and makes more
identities ``L_i + L_j - m L_k = 0`` possible. The LP size is still ``2^n`` in
the number of forms (so practical only for ``n`` up to about 10).

All rationals are kept as ``fractions.Fraction`` for exact symmetry detection.
Canonical subset shapes are computed by trying every choice of basis (of
rank 1, 2, or 3) and minimizing over normalized representatives.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from fractions import Fraction

from atomic_density_lp import (
    lp_equation,
    lp_objective,
    lp_problem,
    lp_variable,
    optimize,
    relevant_set,
)
from two_parameter_lp import (
    linear_form_2d,
    canonical_subset_shape as canonical_subset_shape_2d,
)
from affine_form_lp import (
    affine_form,
    canonical_subset_shape as canonical_subset_shape_1d,
)


class linear_form_3d:
    def __init__(self, a, b, c, d=0, label: str = ""):
        self.a = Fraction(a)
        self.b = Fraction(b)
        self.c = Fraction(c)
        self.d = Fraction(d)
        self.label = label

    def coeffs(self):
        return (self.a, self.b, self.c)

    def is_constant(self) -> bool:
        return self.a == 0 and self.b == 0 and self.c == 0

    def __repr__(self) -> str:
        return f"linear_form_3d({self.a}, {self.b}, {self.c}, {self.d}, {self.label!r})"


def _mat_vec(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def _transpose(matrix):
    return tuple(zip(*matrix))


def _invert_3x3(matrix):
    a = [list(row) for row in matrix]
    n = 3
    inv = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for i in range(n):
        pivot = Fraction(0)
        pivot_row = -1
        for r in range(i, n):
            if a[r][i] != 0:
                pivot = a[r][i]
                pivot_row = r
                break
        if pivot_row < 0:
            return None
        if pivot_row != i:
            a[i], a[pivot_row] = a[pivot_row], a[i]
            inv[i], inv[pivot_row] = inv[pivot_row], inv[i]
        for j in range(n):
            a[i][j] = a[i][j] / pivot
            inv[i][j] = inv[i][j] / pivot
        for r in range(n):
            if r == i:
                continue
            factor = a[r][i]
            if factor == 0:
                continue
            for j in range(n):
                a[r][j] -= factor * a[i][j]
                inv[r][j] -= factor * inv[i][j]
    return tuple(tuple(row) for row in inv)


def _rank_of_vectors(vectors):
    if not vectors:
        return 0
    nonzero = [v for v in vectors if any(x != 0 for x in v)]
    if not nonzero:
        return 0
    basis = [nonzero[0]]
    for v in nonzero[1:]:
        reduced = list(v)
        for b in basis:
            pivot_idx = next((i for i, x in enumerate(b) if x != 0), None)
            if pivot_idx is not None and reduced[pivot_idx] != 0:
                scale = reduced[pivot_idx] / b[pivot_idx]
                reduced = [reduced[k] - scale * b[k] for k in range(len(b))]
        if any(x != 0 for x in reduced):
            basis.append(tuple(reduced))
            if len(basis) == 3:
                return 3
    return len(basis)


def _normalize_rank_three(forms, idx_triple):
    """Normalize forms by mapping forms[idx_triple] to (u,v,w) respectively."""

    i, j, k = idx_triple
    rows = [forms[i].coeffs(), forms[j].coeffs(), forms[k].coeffs()]
    matrix = tuple(tuple(row) for row in rows)
    inverse = _invert_3x3(matrix)
    if inverse is None:
        return None
    translation = _mat_vec(inverse, (-forms[i].d, -forms[j].d, -forms[k].d))

    normalized = []
    for form in forms:
        new_coeffs = _mat_vec(_transpose(inverse), form.coeffs())
        new_constant = (
            form.a * translation[0]
            + form.b * translation[1]
            + form.c * translation[2]
            + form.d
        )
        normalized.append(new_coeffs + (new_constant,))
    return tuple(sorted(normalized))


def _normalize_rank_two(forms, idx_pair):
    """Normalize forms assuming two independent directions; use them as u,v."""

    i, j = idx_pair
    vi = forms[i].coeffs()
    vj = forms[j].coeffs()

    coord_idx = None
    for kk in range(3):
        sub = [[vi[a] for a in range(3) if a != kk], [vj[a] for a in range(3) if a != kk]]
        det = sub[0][0] * sub[1][1] - sub[0][1] * sub[1][0]
        if det != 0:
            coord_idx = kk
            inv2 = [[sub[1][1] / det, -sub[0][1] / det], [-sub[1][0] / det, sub[0][0] / det]]
            break
    if coord_idx is None:
        return None

    other_coords = [a for a in range(3) if a != coord_idx]
    a_idx, b_idx = other_coords
    inverse_full = [
        [Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(0)],
    ]
    inverse_full[a_idx][0] = inv2[0][0]
    inverse_full[a_idx][1] = inv2[0][1]
    inverse_full[b_idx][0] = inv2[1][0]
    inverse_full[b_idx][1] = inv2[1][1]
    inverse_full[coord_idx][2] = Fraction(1)

    t_solve = _mat_vec(inv2, (-forms[i].d, -forms[j].d))
    translation = [Fraction(0), Fraction(0), Fraction(0)]
    translation[a_idx] = t_solve[0]
    translation[b_idx] = t_solve[1]

    normalized = []
    for form in forms:
        coeffs = form.coeffs()
        new_a = sum(coeffs[c] * inverse_full[c][0] for c in range(3))
        new_b = sum(coeffs[c] * inverse_full[c][1] for c in range(3))
        new_c = sum(coeffs[c] * inverse_full[c][2] for c in range(3))
        new_constant = form.d + sum(form.coeffs()[c] * translation[c] for c in range(3))
        normalized.append((new_a, new_b, new_c, new_constant))
    return tuple(sorted(normalized))


def _normalize_rank_one(forms, idx):
    ref = forms[idx]
    ref_vec = ref.coeffs()
    pivot = next((c for c in ref_vec if c != 0), None)
    if pivot is None:
        return None

    normalized = []
    for form in forms:
        if form.is_constant():
            normalized.append((Fraction(0), Fraction(0), Fraction(0), form.d))
            continue
        lam = next(
            (fc / rc for fc, rc in zip(form.coeffs(), ref_vec) if rc != 0),
            None,
        )
        if lam is None:
            return None
        ok = all(
            fc == lam * rc for fc, rc in zip(form.coeffs(), ref_vec)
        )
        if not ok:
            return None
        normalized.append(
            (lam, Fraction(0), Fraction(0), form.d - lam * ref.d)
        )
    return tuple(sorted(normalized))


def _project_to_2d(forms, idx_pair):
    """Given forms whose coefficient vectors are rank 2, express each form in the
    basis of ``forms[i]`` and ``forms[j]`` (the two coefficient vectors) and
    return a list of ``linear_form_2d`` representatives with the residual
    constant term.

    Specifically, for each form L_r we solve
        (a_r, b_r, c_r) = alpha_r * (a_i, b_i, c_i) + beta_r * (a_j, b_j, c_j)
    and return a 2D form alpha_r * u + beta_r * v + (d_r - alpha_r*d_i - beta_r*d_j).
    Returns ``None`` if the decomposition fails (should not happen for rank 2).
    """

    i, j = idx_pair
    vi = forms[i].coeffs()
    vj = forms[j].coeffs()
    di = forms[i].d
    dj = forms[j].d

    # Find a 2x2 minor with nonzero determinant to invert.
    for idx_a in range(3):
        for idx_b in range(idx_a + 1, 3):
            det = vi[idx_a] * vj[idx_b] - vi[idx_b] * vj[idx_a]
            if det == 0:
                continue
            # Columns at idx_a, idx_b form a 2x2 invertible block.
            twod_forms = []
            for form in forms:
                c = form.coeffs()
                # Solve alpha * vi_sub + beta * vj_sub = c_sub
                a_val = (c[idx_a] * vj[idx_b] - c[idx_b] * vj[idx_a]) / det
                b_val = (vi[idx_a] * c[idx_b] - vi[idx_b] * c[idx_a]) / det
                # Verify third coordinate matches
                remaining = 3 - idx_a - idx_b
                expected_remaining = a_val * vi[remaining] + b_val * vj[remaining]
                if expected_remaining != c[remaining]:
                    return None  # rank actually > 2, should not happen
                twod_forms.append(
                    linear_form_2d(a_val, b_val, form.d - a_val * di - b_val * dj)
                )
            return twod_forms
    return None


def _project_to_1d(forms, idx):
    """Given rank-1 forms, express each in multiples of ``forms[idx]``.

    Returns affine_form objects with slope = lambda_r and intercept = residual
    constant term ``d_r - lambda_r * d_ref``.
    """

    ref = forms[idx]
    ref_vec = ref.coeffs()
    pivot_pos = next((k for k, v in enumerate(ref_vec) if v != 0), None)
    if pivot_pos is None:
        return None

    result = []
    for form in forms:
        if form.is_constant():
            result.append(affine_form(slope=Fraction(0), intercept=form.d, label=""))
            continue
        lam = form.coeffs()[pivot_pos] / ref_vec[pivot_pos]
        # Verify proportionality
        ok = all(fc == lam * rc for fc, rc in zip(form.coeffs(), ref_vec))
        if not ok:
            return None
        intercept = form.d - lam * ref.d
        # ensure slope is nonzero (we use intercept-only forms for constant inputs)
        if lam == 0:
            result.append(affine_form(slope=Fraction(0), intercept=intercept, label=""))
        else:
            result.append(affine_form(slope=lam, intercept=intercept, label=""))
    return result


def canonical_subset_shape(forms):
    """Return a rank-tagged canonical shape of the subset under the 3D affine
    group ``GL_3(F_p) ltimes F_p^3``.

    The rank of the coefficient matrix classifies the subset into affine-type
    cases. For rank ``r`` we project into an ``r``-parameter affine model and
    reuse the canonical-shape logic already in ``affine_form_lp`` (rank 1) and
    ``two_parameter_lp`` (rank 2). Rank 3 uses the dedicated rank-3 normalization.
    """

    if not forms:
        return ("r0",)

    rank = _rank_of_vectors([form.coeffs() for form in forms])
    if rank == 0:
        return ("r0",) + tuple(sorted(f.d for f in forms))

    n = len(forms)
    options = []
    if rank == 1:
        for idx in range(n):
            if forms[idx].is_constant():
                continue
            projected = _project_to_1d(forms, idx)
            if projected is None:
                continue
            options.append(canonical_subset_shape_1d(projected))
        if not options:
            return ("r1_fallback",) + tuple(sorted((f.a, f.b, f.c, f.d) for f in forms))
        return ("r1",) + min(options)

    if rank == 2:
        for pair in itertools.permutations(range(n), 2):
            if forms[pair[0]].is_constant() or forms[pair[1]].is_constant():
                continue
            vi = forms[pair[0]].coeffs()
            vj = forms[pair[1]].coeffs()
            det_nonzero = False
            for ia in range(3):
                for ib in range(ia + 1, 3):
                    if vi[ia] * vj[ib] - vi[ib] * vj[ia] != 0:
                        det_nonzero = True
                        break
                if det_nonzero:
                    break
            if not det_nonzero:
                continue
            projected = _project_to_2d(forms, pair)
            if projected is None:
                continue
            options.append(canonical_subset_shape_2d(projected))
        if not options:
            return ("r2_fallback",) + tuple(sorted((f.a, f.b, f.c, f.d) for f in forms))
        return ("r2",) + min(options)

    # rank == 3
    for triple in itertools.permutations(range(n), 3):
        vecs = [forms[t].coeffs() for t in triple]
        if _rank_of_vectors(vecs) < 3:
            continue
        opt = _normalize_rank_three(forms, triple)
        if opt is not None:
            options.append(opt)
    if not options:
        return ("r3_fallback",) + tuple(sorted((f.a, f.b, f.c, f.d) for f in forms))
    return ("r3",) + min(options)


def form_relation_holds(form_i, form_j, form_k, m: int) -> bool:
    return (
        form_i.a + form_j.a - m * form_k.a == 0
        and form_i.b + form_j.b - m * form_k.b == 0
        and form_i.c + form_j.c - m * form_k.c == 0
        and form_i.d + form_j.d - m * form_k.d == 0
    )


def _mask_contains(mask: int, submask: int) -> bool:
    return (mask & submask) == submask


def build_three_parameter_lp(m: int, forms, verbose: bool = False):
    if not forms:
        raise ValueError("Need at least one form.")
    if all(f.is_constant() for f in forms):
        raise ValueError("Need at least one nonconstant form.")

    base_sets = [
        relevant_set(dilation=1, shift=i + 1, label=f.label or f"form_{i}")
        for i, f in enumerate(forms)
    ]

    invalid_triples = []
    for i in range(len(forms)):
        for j in range(len(forms)):
            for k in range(len(forms)):
                if form_relation_holds(forms[i], forms[j], forms[k], m):
                    invalid_triples.append((i, j, k))

    valid_vars = []
    for bits in itertools.product([0, 1], repeat=len(forms)):
        valid = True
        for i, j, k in invalid_triples:
            if bits[i] and bits[j] and bits[k]:
                valid = False
                break
        if not valid:
            continue
        relevant_sets = []
        for idx, bit in enumerate(bits):
            base = base_sets[idx]
            relevant_sets.append(base if bit else base.complement())
        valid_vars.append(lp_variable(relevant_sets=tuple(relevant_sets)))

    equations = [
        lp_equation(
            variables=valid_vars,
            coefficients=[1] * len(valid_vars),
            independent_term=1,
            equals=True,
        )
    ]

    masks_of_valid_vars = []
    for var in valid_vars:
        mask = 0
        for idx, rs in enumerate(var.relevant_sets):
            if not rs.complementary:
                mask |= 1 << idx
        masks_of_valid_vars.append((mask, var))

    subsets_by_shape = defaultdict(list)
    for mask in range(1, 1 << len(forms)):
        subset_forms = [forms[i] for i in range(len(forms)) if mask & (1 << i)]
        subsets_by_shape[canonical_subset_shape(subset_forms)].append(mask)

    num_shape_eqs = 0
    for masks in subsets_by_shape.values():
        if len(masks) < 2:
            continue
        rep_mask = masks[0]
        rep_vars = [var for mask, var in masks_of_valid_vars if _mask_contains(mask, rep_mask)]
        for current_mask in masks[1:]:
            current_vars = [var for mask, var in masks_of_valid_vars if _mask_contains(mask, current_mask)]
            coeff_map = defaultdict(int)
            for var in current_vars:
                coeff_map[var] += 1
            for var in rep_vars:
                coeff_map[var] -= 1
            eq_vars = []
            eq_coeffs = []
            for var, coeff in coeff_map.items():
                if coeff:
                    eq_vars.append(var)
                    eq_coeffs.append(coeff)
            if eq_vars:
                equations.append(
                    lp_equation(
                        variables=eq_vars,
                        coefficients=eq_coeffs,
                        independent_term=0,
                        equals=True,
                    )
                )
                num_shape_eqs += 1

    objective_idx = next(i for i, f in enumerate(forms) if not f.is_constant())
    density_mask = 1 << objective_idx
    density_vars = [var for mask, var in masks_of_valid_vars if _mask_contains(mask, density_mask)]
    objective = lp_objective(variables=density_vars, coefficients=[-1] * len(density_vars))
    problem = lp_problem(sets=base_sets, equations=equations, objective=objective)
    result = optimize(problem)

    if verbose:
        print(f"[3D LP] forms={len(forms)} forbidden_triples={len(invalid_triples)} atoms={len(valid_vars)} shape_eqs={num_shape_eqs}")

    return {
        "problem": problem,
        "result": result,
        "forms": forms,
        "invalid_triples": invalid_triples,
        "shape_equations": num_shape_eqs,
    }


def evaluate_three_parameter_bound(m: int, forms, verbose: bool = False):
    built = build_three_parameter_lp(m=m, forms=forms, verbose=verbose)
    result = built["result"]
    if result.success:
        return -result.fun
    if verbose:
        print(f"LP solver did not succeed: {getattr(result, 'message', '')}")
    return None

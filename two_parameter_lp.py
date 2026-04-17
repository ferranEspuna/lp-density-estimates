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


class linear_form_2d:
    def __init__(self, a, b, c=0, label=""):
        self.a = Fraction(a)
        self.b = Fraction(b)
        self.c = Fraction(c)
        self.label = label

    def coeffs(self):
        return (self.a, self.b)

    def is_constant(self):
        return self.a == 0 and self.b == 0

    def __repr__(self):
        return f"linear_form_2d({self.a}, {self.b}, {self.c}, {self.label!r})"


def determinant(col_1, col_2):
    return col_1[0] * col_2[1] - col_1[1] * col_2[0]


def invert_2x2_columns(col_1, col_2):
    det = determinant(col_1, col_2)
    if det == 0:
        raise ValueError("Matrix is singular.")
    return (
        (col_2[1] / det, -col_2[0] / det),
        (-col_1[1] / det, col_1[0] / det),
    )


def mat_vec_mul(matrix, vector):
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def rank_of_forms(forms):
    vectors = [form.coeffs() for form in forms if not form.is_constant()]
    if not vectors:
        return 0
    ref = vectors[0]
    for vector in vectors[1:]:
        if determinant(ref, vector) != 0:
            return 2
    return 1


def scalar_multiple(vector, ref):
    if ref[0] != 0:
        return vector[0] / ref[0]
    if ref[1] != 0:
        return vector[1] / ref[1]
    raise ValueError("Reference vector must be nonzero.")


def normalize_rank_one(forms, ref_idx):
    ref = forms[ref_idx]
    ref_vector = ref.coeffs()
    normalized = []
    for form in forms:
        if form.is_constant():
            normalized.append((Fraction(0), Fraction(0), form.c))
            continue
        lam = scalar_multiple(form.coeffs(), ref_vector)
        normalized.append((lam, Fraction(0), form.c - lam * ref.c))
    return tuple(sorted(normalized))


def normalize_rank_two(forms, idx_1, idx_2):
    ref_1 = forms[idx_1]
    ref_2 = forms[idx_2]
    basis = (ref_1.coeffs(), ref_2.coeffs())
    inverse = invert_2x2_columns(*basis)

    translation = mat_vec_mul(
        (
            (inverse[0][0], inverse[1][0]),
            (inverse[0][1], inverse[1][1]),
        ),
        (-ref_1.c, -ref_2.c),
    )

    normalized = []
    for form in forms:
        coeffs = mat_vec_mul(inverse, form.coeffs())
        constant = form.a * translation[0] + form.b * translation[1] + form.c
        normalized.append((coeffs[0], coeffs[1], constant))
    return tuple(sorted(normalized))


def canonical_subset_shape(forms):
    if not forms:
        return tuple()

    subset_rank = rank_of_forms(forms)
    if subset_rank == 0:
        return tuple(sorted((Fraction(0), Fraction(0), form.c) for form in forms))

    normalized_options = []
    if subset_rank == 1:
        for idx, form in enumerate(forms):
            if form.is_constant():
                continue
            normalized_options.append(normalize_rank_one(forms, idx))
        return min(normalized_options)

    for idx_1 in range(len(forms)):
        if forms[idx_1].is_constant():
            continue
        for idx_2 in range(len(forms)):
            if idx_1 == idx_2 or forms[idx_2].is_constant():
                continue
            if determinant(forms[idx_1].coeffs(), forms[idx_2].coeffs()) == 0:
                continue
            normalized_options.append(normalize_rank_two(forms, idx_1, idx_2))
    return min(normalized_options)


def form_relation_holds(form_i, form_j, form_k, m):
    return (
        form_i.a + form_j.a - m * form_k.a == 0
        and form_i.b + form_j.b - m * form_k.b == 0
        and form_i.c + form_j.c - m * form_k.c == 0
    )


def make_direction_family(direction, shifts, dilation=1, prefix="F"):
    a, b = direction
    return [
        linear_form_2d(
            Fraction(a, dilation),
            Fraction(b, dilation),
            Fraction(shift, dilation),
            label=f"{prefix}_{idx}",
        )
        for idx, shift in enumerate(shifts)
    ]


def standard_family(m, shifts_u, shifts_v, shifts_uv):
    return (
        make_direction_family((1, 0), shifts_u, prefix="u")
        + make_direction_family((0, 1), shifts_v, prefix="v")
        + make_direction_family((1, 1), shifts_uv, dilation=m, prefix="uv")
    )


def pentad_family(m):
    """
    Five witness directions: u, v, (u+v)/m, (2u-v)/m, (u-2v)/m, each with a single
    shift at 0. For m=3 the two-parameter LP gives 8/15; for sample choices
    m in {4, 5, 7} the same family currently returns 4/7 in floating-point solves
    (see experiments script and latex notes).
    """
    return (
        make_direction_family((1, 0), (0,), prefix="u")
        + make_direction_family((0, 1), (0,), prefix="v")
        + make_direction_family((1, 1), (0,), dilation=m, prefix="u_plus_v")
        + make_direction_family((2, -1), (0,), dilation=m, prefix="two_u_minus_v")
        + make_direction_family((1, -2), (0,), dilation=m, prefix="u_minus_two_v")
    )


def _mask_contains(mask, submask):
    return (mask & submask) == submask


def build_two_parameter_lp(m, forms, verbose=False):
    if not forms:
        raise ValueError("Need at least one form.")
    if all(form.is_constant() for form in forms):
        raise ValueError("Need at least one nonconstant form.")

    base_sets = [
        relevant_set(dilation=1, shift=idx + 1, label=form.label or f"form_{idx}")
        for idx, form in enumerate(forms)
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
        for idx, rel_set in enumerate(var.relevant_sets):
            if not rel_set.complementary:
                mask |= 1 << idx
        masks_of_valid_vars.append((mask, var))

    subsets_by_shape = defaultdict(list)
    for mask in range(1, 1 << len(forms)):
        subset_forms = [forms[idx] for idx in range(len(forms)) if mask & (1 << idx)]
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

    objective_idx = next(idx for idx, form in enumerate(forms) if not form.is_constant())
    density_mask = 1 << objective_idx
    density_vars = [var for mask, var in masks_of_valid_vars if _mask_contains(mask, density_mask)]
    objective = lp_objective(variables=density_vars, coefficients=[-1] * len(density_vars))
    problem = lp_problem(sets=base_sets, equations=equations, objective=objective)
    result = optimize(problem)

    if verbose:
        print(f"Forbidden triples: {len(invalid_triples)}")
        print(f"Valid atoms: {len(valid_vars)}")
        print(f"Affine-shape equations: {num_shape_eqs}")

    return {
        "problem": problem,
        "result": result,
        "forms": forms,
        "invalid_triples": invalid_triples,
        "shape_equations": num_shape_eqs,
    }


def evaluate_two_parameter_bound(m, forms, verbose=False):
    built = build_two_parameter_lp(m=m, forms=forms, verbose=verbose)
    result = built["result"]
    if result.success:
        return -result.fun
    if verbose:
        print(f"LP solver did not succeed: {getattr(result, 'message', '')}")
    return None

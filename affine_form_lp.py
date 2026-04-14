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


class affine_form:
    def __init__(self, slope: Fraction, intercept: Fraction, label: str):
        self.slope = Fraction(slope)
        self.intercept = Fraction(intercept)
        self.label = label

    @classmethod
    def from_relevant_set(cls, dilation: int, shift: int, label: str):
        return cls(
            slope=Fraction(1, dilation),
            intercept=Fraction(shift, dilation),
            label=label,
        )

    def normalized_against(self, ref):
        scale = Fraction(1, ref.slope)
        shift = -ref.intercept / ref.slope
        return (
            self.slope * scale,
            self.intercept + self.slope * shift,
        )

    def __repr__(self):
        return f"affine_form({self.slope}, {self.intercept}, {self.label!r})"


def canonical_subset_shape(forms):
    if not forms:
        return tuple()

    normalized_options = []
    for ref in forms:
        normalized = [f.normalized_against(ref) for f in forms]
        normalized_options.append(tuple(sorted(normalized)))
    return min(normalized_options)


def form_relation_holds(form_i, form_j, form_k, m):
    return (
        form_i.slope + form_j.slope - m * form_k.slope == 0
        and form_i.intercept + form_j.intercept - m * form_k.intercept == 0
    )


def _mask_contains(mask, submask):
    return (mask & submask) == submask


def build_affine_lp(m, layers, verbose=False):
    base_sets = []
    forms = []
    for dilation, shifts in layers:
        for shift in shifts:
            label = f"{dilation}A-{shift}"
            base_sets.append(relevant_set(dilation=dilation, shift=shift))
            forms.append(affine_form.from_relevant_set(dilation, shift, label))

    if not base_sets:
        raise ValueError("Need at least one base set.")

    invalid_triples = []
    for i in range(len(forms)):
        for j in range(len(forms)):
            for k in range(len(forms)):
                if form_relation_holds(forms[i], forms[j], forms[k], m):
                    invalid_triples.append((i, j, k))

    valid_vars = []
    for bits in itertools.product([0, 1], repeat=len(base_sets)):
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

    density_vars = [
        var
        for mask, var in masks_of_valid_vars
        if _mask_contains(mask, 1)
    ]
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


def evaluate_affine_bound(m, layers, verbose=False):
    built = build_affine_lp(m=m, layers=layers, verbose=verbose)
    result = built["result"]
    if result.success:
        return -result.fun
    return 0.5

import itertools
import random
from collections import defaultdict
from fractions import Fraction
import numpy as np
from atomic_density_lp import relevant_set, lp_variable, lp_equation, lp_objective, lp_problem, optimize
import sys

def evaluate_bound(X, Y, verbose=False):
    base_sets_X = [relevant_set(dilation=1, shift=x) for x in X]
    base_sets_Y = [relevant_set(dilation=2, shift=y) for y in Y]

    invalid_triples = []
    for i in range(len(X)):
        for j in range(len(X)):
            for k in range(len(Y)):
                if X[i] + X[j] - 2 * Y[k] == 0:
                    invalid_triples.append((i, j, k))

    valid_vars = []
    for bits in itertools.product([0, 1], repeat=len(X) + len(Y)):
        valid = True
        for (i, j, k) in invalid_triples:
            if bits[i] == 1 and bits[j] == 1 and bits[len(X) + k] == 1:
                valid = False
                break
        if valid:
            rs_x = [base_sets_X[m] if bits[m] else base_sets_X[m].complement() for m in range(len(X))]
            rs_y = [base_sets_Y[m] if bits[len(X) + m] else base_sets_Y[m].complement() for m in range(len(Y))]
            valid_vars.append(lp_variable(relevant_sets=rs_x + rs_y))

    if not valid_vars:
        return 0

    equations = []
    equations.append(lp_equation(
        variables=valid_vars,
        coefficients=[1] * len(valid_vars),
        independent_term=1,
        equals=True
    ))

    masks_of_valid_vars = []
    for v in valid_vars:
        mask = 0
        for i in range(len(X)):
            if not v.relevant_sets[i].complementary:
                mask |= (1 << i)
        for i in range(len(Y)):
            if not v.relevant_sets[len(X) + i].complementary:
                mask |= (1 << (len(X) + i))
        masks_of_valid_vars.append((mask, v))

    subsets_by_shape = defaultdict(list)
    for mask in range(1, 1 << (len(X) + len(Y))):
        S_X = tuple(X[i] for i in range(len(X)) if (mask & (1 << i)))
        S_Y = tuple(Y[i] for i in range(len(Y)) if (mask & (1 << (len(X) + i))))
        m = min(S_X + S_Y)
        shape = (tuple(x - m for x in S_X), tuple(y - m for y in S_Y))
        subsets_by_shape[shape].append(mask)

    for shape, S_list in subsets_by_shape.items():
        if len(S_list) > 1:
            rep_mask_val = S_list[0]
            rep_vars = [v for m, v in masks_of_valid_vars if (m & rep_mask_val) == rep_mask_val]
            for S_mask_val in S_list[1:]:
                S_vars = [v for m, v in masks_of_valid_vars if (m & S_mask_val) == S_mask_val]
                coeff_map = defaultdict(int)
                for v in S_vars: coeff_map[v] += 1
                for v in rep_vars: coeff_map[v] -= 1
                eq_vars = []
                eq_coeffs = []
                for v, c in coeff_map.items():
                    if c != 0:
                        eq_vars.append(v)
                        eq_coeffs.append(c)
                if eq_vars:
                    equations.append(lp_equation(variables=eq_vars, coefficients=eq_coeffs, independent_term=0, equals=True))

    mask_A = (1 << 0)
    vars_A = [v for m, v in masks_of_valid_vars if (m & mask_A) == mask_A]
    mask_2A = (1 << len(X))
    vars_2A = [v for m, v in masks_of_valid_vars if (m & mask_2A) == mask_2A]
    coeff_map = defaultdict(int)
    for v in vars_2A: coeff_map[v] += 1
    for v in vars_A: coeff_map[v] -= 1
    eq_vars = []
    eq_coeffs = []
    for v, c in coeff_map.items():
        if c != 0:
            eq_vars.append(v)
            eq_coeffs.append(c)
    if eq_vars:
        equations.append(lp_equation(variables=eq_vars, coefficients=eq_coeffs, independent_term=0, equals=True))

    # Dilation invariance mapping cleanly across sets A to 2A natively
    for shape, masks in subsets_by_shape.items():
        S_X = shape[0]
        S_Y = shape[1]
        if S_Y or not S_X: continue
        dilated_S_Y = tuple(sorted(2 * x for x in S_X))
        m = min(dilated_S_Y)
        dilated_shape = (tuple(), tuple(y - m for y in dilated_S_Y))
        if dilated_shape in subsets_by_shape:
            rep_mask_X = masks[0]
            rep_mask_Y = subsets_by_shape[dilated_shape][0]
            vars_X = [v for m, v in masks_of_valid_vars if (m & rep_mask_X) == rep_mask_X]
            vars_Y = [v for m, v in masks_of_valid_vars if (m & rep_mask_Y) == rep_mask_Y]
            coeff_map = defaultdict(int)
            for v in vars_X: coeff_map[v] += 1
            for v in vars_Y: coeff_map[v] -= 1
            eq_vars = []
            eq_coeffs = []
            for v, c in coeff_map.items():
                if c != 0:
                    eq_vars.append(v)
                    eq_coeffs.append(c)
            if eq_vars:
                equations.append(lp_equation(variables=eq_vars, coefficients=eq_coeffs, independent_term=0, equals=True))

    # 4. Scaled Homogeneous Shift Dilation Invariance
    # Since the optimal bounding set A* avoiding an equation can be symmetrized across all dilations c \in F_p*, 
    # its atomic densities natively satisfy invariance if we dilate the identical structural shape's shifts by c. 
    # We group all derived shapes structurally related by a scalar multiplier.
    dilation_groups = defaultdict(list)
    for shape in subsets_by_shape.keys():
        S_X, S_Y = shape
        all_vals = S_X + S_Y
        nonzero_vals = [v for v in all_vals if v != 0]
        if not nonzero_vals:
            # Pure anchors natively group alone
            dilation_groups[shape].append(shape)
            continue
        
        # Build canonical rational representation via dividing by the first non-zero element
        f = nonzero_vals[0]
        canon_X = tuple(Fraction(x, f) for x in S_X)
        canon_Y = tuple(Fraction(y, f) for y in S_Y)
        canon = (canon_X, canon_Y)
        dilation_groups[canon].append(shape)

    num_hom_eqs = 0
    for canon, shapes in dilation_groups.items():
        if len(shapes) > 1:
            rep_shape = shapes[0]
            rep_mask = subsets_by_shape[rep_shape][0]
            rep_vars = [v for m, v in masks_of_valid_vars if (m & rep_mask) == rep_mask]
            
            for S_shape in shapes[1:]:
                S_mask = subsets_by_shape[S_shape][0]
                S_vars = [v for m, v in masks_of_valid_vars if (m & S_mask) == S_mask]
                
                coeff_map = defaultdict(int)
                for v in S_vars: coeff_map[v] += 1
                for v in rep_vars: coeff_map[v] -= 1
                    
                eq_vars = []
                eq_coeffs = []
                for v, c in coeff_map.items():
                    if c != 0:
                        eq_vars.append(v)
                        eq_coeffs.append(c)
                        
                if eq_vars:
                    equations.append(lp_equation(
                        variables=eq_vars,
                        coefficients=eq_coeffs,
                        independent_term=0,
                        equals=True
                    ))
                    num_hom_eqs += 1
    
    print(f"Added {num_hom_eqs} scalar shift dilation equations.")

    x0_vars = [v for m, v in masks_of_valid_vars if (m & mask_A) == mask_A]
    if not x0_vars: return 0
    objective = lp_objective(variables=x0_vars, coefficients=[-1]*len(x0_vars))
    prob = lp_problem(sets=base_sets_X+base_sets_Y, equations=equations, objective=objective)
    res = optimize(prob)
    
    if res.success and verbose:
        print("\nNonzero atoms in optimal solution:")
        for idx, val in enumerate(res.x):
            if val > 1e-8:
                v = res.lp_variables[idx]
                bits = []
                for i in range(len(X) + len(Y)):
                    bits.append('0' if v.relevant_sets[i].complementary else '1')
                print(f"Atom X:{''.join(bits[:len(X)])} Y:{''.join(bits[len(X):])} : {val:.5f}")

    return -res.fun if res.success else 0.5

def random_search():
    best_bound = 0.5
    size_X = 6
    size_Y = 6
    max_val = 15
    print("Starting randomized search with full dilation equations...")
    
    attempts = 0
    while True:
        rem_X = random.sample(range(1, max_val + 1), size_X - 1)
        X = [0] + sorted(rem_X)
        Y = sorted(random.sample(range(0, max_val + 1), size_Y))
        
        has_overlap = False
        for i in range(size_X):
            for j in range(size_X):
                for k in range(size_Y):
                    if X[i] + X[j] == 2 * Y[k]:
                        has_overlap = True
                        break
        if not has_overlap:
            continue
            
        bound = evaluate_bound(X, Y)
        attempts += 1
        
        if bound < best_bound - 1e-5:
            best_bound = bound
            print(f"NEW BEST! {best_bound:.5f} with X={X}, Y={Y}", flush=True)

        if attempts % 20 == 0:
            print(f"Attempt {attempts}, best so far: {best_bound:.5f}", flush=True)

if __name__ == "__main__":
    random_search()

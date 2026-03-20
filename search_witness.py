import itertools
import random
from collections import defaultdict
import numpy as np
from atomic_density_lp import relevant_set, lp_variable, lp_equation, lp_objective, lp_problem, optimize
import sys

def evaluate_bound(X, Y):
    base_sets_X = [relevant_set(dilation=1, shift=x) for x in X]
    base_sets_Y = [relevant_set(dilation=2, shift=y) for y in Y]

    invalid_triples = []
    for i in range(len(X)):
        for j in range(len(X)):
            for k in range(len(Y)):
                if X[i] + X[j] - 2 * Y[k] == 0:
                    invalid_triples.append((i, j, k))

    valid_vars = []
    # optimize generation
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
                    equations.append(lp_equation(
                        variables=eq_vars,
                        coefficients=eq_coeffs,
                        independent_term=0,
                        equals=True
                    ))

    # delta(A) = delta(2A)
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
        equations.append(lp_equation(
            variables=eq_vars,
            coefficients=eq_coeffs,
            independent_term=0,
            equals=True
        ))

    x0_vars = [v for m, v in masks_of_valid_vars if (m & mask_A) == mask_A]
    if not x0_vars:
        return 0

    objective = lp_objective(
        variables=x0_vars,
        coefficients=[-1] * len(x0_vars)
    )

    prob = lp_problem(sets=base_sets_X + base_sets_Y, equations=equations, objective=objective)
    
    # We can pass method="highs" natively or suppress output
    from scipy.optimize import linprog
    # The atomic_density_lp script runs it. 
    res = optimize(prob)
    return -res.fun if res.success else 0.5


def main():
    best_bound = 0.5
    size_X = 5
    size_Y = 5
    max_val = 15
    print("Starting randomized search for witness sets...")
    
    attempts = 0
    while True:
        # keep X[0] = 0 as translational anchor
        # sample rest
        rem_X = random.sample(range(1, max_val + 1), size_X - 1)
        X = [0] + sorted(rem_X)
        Y = sorted(random.sample(range(0, max_val + 1), size_Y))
        
        # We need X[i] + X[j] - 2*Y[k] == 0 to trigger invalid triples, otherwise bound will trivially be 0.5
        # So we should ensure at least SOME overlaps
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
            print(f"[{attempts}] New best bound: {best_bound:.5f} with X={X}, Y={Y}", flush=True)
            if best_bound < 0.49:
                # We found a non-trivial improvement!
                pass
        
        if attempts % 50 == 0:
            print(f"Attempt {attempts}, best so far: {best_bound:.5f}", flush=True)

if __name__ == "__main__":
    main()

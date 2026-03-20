import itertools
from collections import defaultdict
from atomic_density_lp import relevant_set, lp_variable, lp_equation, lp_objective, lp_problem, optimize

X = [0, 2, 3, 4, 7, 8, 9, 10]
Y = [1, 2, 3, 4, 5, 6, 7, 8]

base_sets_X = [relevant_set(dilation=1, shift=x) for x in X]
base_sets_Y = [relevant_set(dilation=2, shift=y) for y in Y]

# find all (i, j, k) where x_i + x_j - 2*x_k == d
invalid_triples = []
for i in range(len(X)):
    for j in range(len(X)):
        for k in range(len(Y)):
            if X[i] + X[j] - 2*Y[k] == 0:
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

print(f"Number of valid atoms: {len(valid_vars)} out of {2**(len(X) + len(Y))}")

equations = []

# 1. sum of all variables is 1
equations.append(lp_equation(
    variables=valid_vars,
    coefficients=[1]*len(valid_vars),
    independent_term=1,
    equals=True
))

# 2. translation invariance
# OPTIMIZATION: Instead of re-evaluating shapes element by element iteratively, we process intersections 
# through simple deterministic bitmasks mapping valid subsets directly inside valid_vars.
# $m$ represents a generic binary shift overlay across all available permutations simultaneously spanning X and Y.
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

# We aggregate structurally identical translation configurations spanning combinations
# Canonical formats anchor combinations against their smallest minimal common element `m`.
subsets_by_shape = defaultdict(list)
for mask in range(1, 1 << (len(X) + len(Y))):
    S_X = tuple(X[i] for i in range(len(X)) if (mask & (1 << i)))
    S_Y = tuple(Y[i] for i in range(len(Y)) if (mask & (1 << (len(X) + i))))
    m = min(S_X + S_Y) # Shifter constant uniquely defining subset clusters
    shape = (tuple(x - m for x in S_X), tuple(y - m for y in S_Y))
    subsets_by_shape[shape].append(mask)

num_invariance_eqs = 0
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
                num_invariance_eqs += 1

print(f"Added {num_invariance_eqs} translate-invariance equations.")

# Also equate delta(A) and delta(2A) mapping shapes directly
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
    print("Added delta(A) = delta(2A) equation.")


# 3. Dilation invariance equations
# NON-TRIVIAL LOGIC: For any structural configuration `S_X` evaluated purely under dilation array 1 (A - x), 
# enforcing a strict scalar multiplication directly onto the shape theoretically projects its target intersections onto 2A intersections.
# Because \delta(A) functionally maps equitably to \delta(2A) natively without overlap modifications due to cardinality bijections inside F_p, 
# identifying instances identically corresponding locally mapped subsets `2*S_X == S_Y` enables strict equality constraint formulations natively!

num_dilation_eqs = 0
for shape, masks in subsets_by_shape.items():
    S_X = shape[0]
    S_Y = shape[1]
    
    # We strictly isolate target subsets evaluating purely structurally inside X 
    if S_Y: continue
    if not S_X: continue
    
    # Scale X subset parameters directly scaling 2*A offsets
    dilated_S_Y = tuple(sorted(2 * x for x in S_X))
    m = min(dilated_S_Y)
    dilated_shape = (tuple(), tuple(y - m for y in dilated_S_Y))
    
    if dilated_shape in subsets_by_shape:
        # We can structurally equate this configuration internally matching densities across the independent witness vectors
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
            equations.append(lp_equation(
                variables=eq_vars,
                coefficients=eq_coeffs,
                independent_term=0,
                equals=True
            ))
            num_dilation_eqs += 1

print(f"Added {num_dilation_eqs} dilation-invariance equations.")

x0_vars = [v for m, v in masks_of_valid_vars if (m & mask_A) == mask_A]
objective = lp_objective(
    variables=x0_vars,
    coefficients=[-1]*len(x0_vars)
)

prob = lp_problem(sets=base_sets_X + base_sets_Y, equations=equations, objective=objective)
print("Executing optimization...")
res = optimize(prob)
print("Optimization finished.")

print(f"Max delta(A) bounded by: {-res.fun:.6f}")
print("Success:", res.success, res.message)
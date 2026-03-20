import itertools
from collections import defaultdict
from atomic_density_lp import relevant_set, lp_variable, lp_equation, lp_objective, lp_problem, optimize

X = [0, 2, 3, 4, 7, 8, 9, 10]
d = 6

base_sets = [relevant_set(dilation=1, shift=x) for x in X]

# find all (i, j, k) where x_i + x_j - 2*x_k == d
invalid_triples = []
for i in range(len(X)):
    for j in range(len(X)):
        for k in range(len(X)):
            if X[i] + X[j] - 2*X[k] == d:
                invalid_triples.append((i, j, k))

valid_vars = []
for bits in itertools.product([0, 1], repeat=len(X)):
    valid = True
    for (i, j, k) in invalid_triples:
        if bits[i] == 1 and bits[j] == 1 and bits[k] == 1:
            valid = False
            break
    if valid:
        rs = [base_sets[m] if bits[m] else base_sets[m].complement() for m in range(len(X))]
        valid_vars.append(lp_variable(relevant_sets=rs))

print(f"Number of valid atoms: {len(valid_vars)} out of {2**len(X)}")

equations = []

# 1. sum of all variables is 1
equations.append(lp_equation(
    variables=valid_vars,
    coefficients=[1]*len(valid_vars),
    independent_term=1,
    equals=True
))

# 2. translation invariance
subsets_by_shape = defaultdict(list)
for r in range(1, len(X)+1):
    for S_tuple in itertools.combinations(X, r):
        m = min(S_tuple)
        shape = frozenset(x - m for x in S_tuple)
        subsets_by_shape[shape].append(set(S_tuple))

num_invariance_eqs = 0
for shape, S_list in subsets_by_shape.items():
    if len(S_list) > 1:
        S_rep = S_list[0]
        rep_vars = []
        for v in valid_vars:
            match = True
            for x_val in S_rep:
                idx = X.index(x_val)
                if v.relevant_sets[idx].complementary:
                    match = False
                    break
            if match:
                rep_vars.append(v)
                
        for S in S_list[1:]:
            S_vars = []
            for v in valid_vars:
                match = True
                for x_val in S:
                    idx = X.index(x_val)
                    if v.relevant_sets[idx].complementary:
                        match = False
                        break
                if match:
                    S_vars.append(v)
                    
            coeff_map = defaultdict(int)
            for v in S_vars:
                coeff_map[v] += 1
            for v in rep_vars:
                coeff_map[v] -= 1
                
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

x0_vars = []
for v in valid_vars:
    if not v.relevant_sets[0].complementary:
        x0_vars.append(v)

objective = lp_objective(
    variables=x0_vars,
    coefficients=[-1]*len(x0_vars)
)

prob = lp_problem(sets=base_sets, equations=equations, objective=objective)
print("Executing optimization...")
res = optimize(prob)
print("Optimization finished.")

print(f"Max delta(A) bounded by: {-res.fun:.6f}")
print("Expected around 2/7 =", 2/7)
print("Success:", res.success, res.message)

from scipy.optimize import linprog
import numpy as np
from pydantic import BaseModel, model_validator

# sets of either $dA-x$ or $\mathbb{F}_p \setminus (dA-x)$
class relevant_set(BaseModel):
    model_config = {"frozen": True}
    dilation: int = 1
    shift: int = 1
    complementary: bool = False

    @model_validator(mode="after")
    def check_nonzero_dilation(self):
        if self.dilation == 0:
            raise ValueError("Dilation cannot be zero.")
        return self

    def complement(self):
        return relevant_set(
            dilation=self.dilation,
            shift=self.shift,
            complementary=not self.complementary
        )

# this is the density of the intersection of a family of relevant sets
class lp_variable(BaseModel):
    model_config = {"frozen": True}
    relevant_sets: tuple[relevant_set, ...]

    @model_validator(mode="after")
    def check_unique_sets(self):
        if len(set(self.relevant_sets)) != len(self.relevant_sets):
            raise ValueError("The list of relevant sets must not contain duplicates.")
        return self

    @model_validator(mode="after")
    def check_no_complement_pairs(self):
        if set(self.relevant_sets) & set(r.complement() for r in self.relevant_sets) != set():
            raise ValueError("The list of relevant sets must not contain complement pairs.")
        return self
    
# relationships between variables
class lp_equation(BaseModel):
    variables: list[lp_variable]
    coefficients: list[int]
    independent_term: int = 0
    equals: bool = False # if true, the equation is of the form \sum c_i x_i = independent_term, otherwise it is \sum c_i x_i \ge independent_term

    @model_validator(mode="after")
    def check_length(self):
        if len(self.variables) != len(self.coefficients):
            raise ValueError("The number of variables must equal the number of coefficients.")
        return self

class lp_objective(BaseModel):
    variables: list[lp_variable]
    coefficients: list[int]

    @model_validator(mode="after")
    def check_length(self):
        if len(self.variables) != len(self.coefficients):
            raise ValueError("The number of variables must equal the number of coefficients.")
        return self

    

class lp_problem(BaseModel):
    sets: list[relevant_set]
    equations: list[lp_equation]
    objective: lp_objective

    def variables(self):
        return list(set(v for eq in self.equations for v in eq.variables))

    # we assume the sets are of the form $dA - x$ and not $\mathbb{F}_p \setminus (dA - x)$
    @model_validator(mode="after")
    def check_no_complementary_sets(self):
        for r in self.sets:
            if r.complementary:
                raise ValueError("Base sets must not be complementary")
        return self

    # make sure the variables in the objective are in the list of variables
    @model_validator(mode="after")
    def check_objective_variables_in_sets(self):
        problem_vars = self.variables()
        for v in self.objective.variables:
            if v not in problem_vars:
                raise ValueError("Objective functions must only contain variables used in the equations")
        return self

    # make sure the variables are composed of the sets
    @model_validator(mode="after")
    def check_variables_in_sets(self):
        allowed_sets = set(self.sets) | set(r.complement() for r in self.sets)
        for v in self.variables():
            if not all(r in allowed_sets for r in v.relevant_sets):
                raise ValueError(f"Variable contains relevant sets not in the problem base sets.")
        return self

    # for all the base sets, each variable should either include it or its complement
    @model_validator(mode="after")
    def check_variables_are_long(self):
        expected_len = len(self.sets)
        for v in self.variables():
            if len(v.relevant_sets) != expected_len:
                raise ValueError(f"Variable does not have the correct length of relevant sets.")
        return self


def optimize(lp: lp_problem):
    all_vars = set()
    for eq in lp.equations:
        all_vars.update(eq.variables)
    all_vars.update(lp.objective.variables)
    variables = list(all_vars)
    
    if not variables:
        raise ValueError("LP problem has no variables")
        
    var_idx = {v: i for i, v in enumerate(variables)}
    n_vars = len(variables)
    
    c = np.zeros(n_vars)
    for v, coef in zip(lp.objective.variables, lp.objective.coefficients):
        c[var_idx[v]] += coef
        
    A_eq = []
    b_eq = []
    A_ub = []
    b_ub = []
    
    for eq in lp.equations:
        row = np.zeros(n_vars)
        for v, coef in zip(eq.variables, eq.coefficients):
            row[var_idx[v]] += coef
            
        if eq.equals:
            A_eq.append(row)
            b_eq.append(eq.independent_term)
        else:
            A_ub.append(-row)
            b_ub.append(-eq.independent_term)
            
    bounds = [(0, 1) for _ in range(n_vars)]
    
    A_eq_arr = np.array(A_eq) if A_eq else None
    b_eq_arr = np.array(b_eq) if b_eq else None
    A_ub_arr = np.array(A_ub) if A_ub else None
    b_ub_arr = np.array(b_ub) if b_ub else None
    
    res = linprog(c, A_ub=A_ub_arr, b_ub=b_ub_arr, A_eq=A_eq_arr, b_eq=b_eq_arr, bounds=bounds)
    return res


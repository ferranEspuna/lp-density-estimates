# LP Investigation Notes

These notes are based on [current_status.tex](/home/fespuna/lp-density-estimates/current_status.tex) and the current scripts
[search_witness.py](/home/fespuna/lp-density-estimates/search_witness.py),
[search_witness_dil.py](/home/fespuna/lp-density-estimates/search_witness_dil.py),
and [run_proof.py](/home/fespuna/lp-density-estimates/run_proof.py).

## 1. What the current obstruction really is

The issue described in `current_status.tex` is not just that the witness sets use too few shifts. The deeper problem is that the present LP is built from a **one-parameter family of events**

\[
t \in dA-x \iff \frac{t+x}{d} \in A.
\]

So each base set should be viewed as an affine form in one random parameter:

\[
L_{d,x}(t) := \frac{t+x}{d}.
\]

Then an atom is recording the joint membership pattern of the values `1_A(L_{d,x}(t))`.

In this language, a forbidden triple appears exactly when

\[
L_i(t) + L_j(t) - m L_k(t) \equiv 0
\]

as an affine identity. Equivalently,

\[
\frac{1}{d_i} + \frac{1}{d_j} = \frac{m}{d_k},
\qquad
\frac{x_i}{d_i} + \frac{x_j}{d_j} = \frac{m x_k}{d_k}.
\]

This is the clean generalization of the `m=4` condition `x_i + x_j = 2 y_k`.

## 2. What changes in the new prototype

The new file [affine_form_lp.py](/home/fespuna/lp-density-estimates/affine_form_lp.py) implements this formulation directly.

It does two things differently from `search_witness_dil.py`:

1. It treats `dA-x` as the affine form `L_{d,x}(t) = (t+x)/d`, so the forbidden triples are derived from the exact identity `L_i + L_j - m L_k = 0`.
2. It groups subset constraints by **full affine reparameterization** of the parameter `t`, rather than by ad hoc translation plus a few extra dilation equalities.

So this is a mathematically cleaner LP, and it works for any `m`.

## 3. What small experiments say

I tested small witnesses with the new builder.

- `m=3`, layers `[(1, [0,1]), (2, [0,1])]`: bound `0.5`
- `m=3`, layers `[(1, [0,1,2]), (2, [0,1,2])]`: bound `0.5`
- `m=4`, layers `[(1, [0,1]), (2, [0,1])]`: bound `0.5`
- `m=4`, layers `[(1, [0,1,2]), (2, [0,1,2])]`: bound `0.5`
- small brute-force searches over dilations like `(1,2,3)` and `(1,2,4)` with 2-3 shifts per layer also stayed at `0.5`

So the conclusion is:

The present failure is not caused by the particular translation/dilation equalities in `search_witness_dil.py`. Even after replacing them by the exact one-parameter affine symmetry, the LP still has too much room.

## 4. Promising LP directions from here

### A. Two-parameter linear-form LP

This is the most promising extension.

Instead of sampling one random parameter `t`, sample two parameters `(u,v)` and look at linear forms

\[
L_r(u,v) = a_r u + b_r v + c_r.
\]

Then a forbidden triple is encoded whenever

\[
L_i + L_j - m L_k \equiv 0.
\]

This removes the main obstruction from `current_status.tex`: cancellation no longer has to come from the single common variable `t`; it can happen in the coefficient vectors `(a_r,b_r)`.

Example for general `m`:

\[
L_1(u,v)=u,\qquad
L_2(u,v)=v,\qquad
L_3(u,v)=\frac{u+v}{m}.
\]

Then `L_1 + L_2 - m L_3 = 0` identically.

Why this matters:

- it gives a genuine LP generalization for every `m`, not just special even cases;
- it uses direct evaluations of `1_A` rather than switching between `A`, `2A`, `3A`, ... as separate layers;
- it is closer in spirit to flag algebras / local-distribution hierarchies.

### B. Extension-hierarchy LP

Keep the one-parameter affine-form model, but do not optimize on a single witness family in isolation.

Introduce variables for a larger witness family `F'` and require that its atom distribution marginalizes to the atom distribution on `F`. This gives a hierarchy of local consistency constraints.

This is still linear, but it is stronger than solving one witness LP at a time.

### C. LP plus additive-combinatorial cuts

If the pure atom LP keeps admitting fake `1/2` pseudodistributions, add inequalities coming from external theorems rather than symmetry alone.

The natural candidates are:

- Pollard-type popular sum inequalities
- Cauchy-Davenport lower bounds on sumsets of dense slices
- robust `3k-4` style inequalities on structured pieces

These are not automatic from the atomic formalism, but they are valid linear cuts once written in terms of the chosen marginals.

## 5. Recommended next step

If the goal is to find a formulation that has a realistic chance to beat `1/2`, the next thing to implement is:

1. a **two-parameter linear-form LP** for `1_A(L_r(u,v))`;
2. optionally with a small extension hierarchy on top of that.

If the goal is only to continue searching within the current family, use `affine_form_lp.py`, not `search_witness_dil.py`, because it is the cleaner formulation and already supports arbitrary `m`.

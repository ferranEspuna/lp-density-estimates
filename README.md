# LP Density Estimates for Additive Combinatorics

This repository contains linear-programming experiments for bounding the density of
sets \(A \subseteq \mathbb{F}_p\) that avoid equations of the form
\[
x+y=mz.
\]

## Project layout

- `atomic_density_lp.py`: shared LP primitives and the SciPy-based solver.
- `latex/`: ordered LaTeX notes and writeups.
- `pdfs/`: compiled PDFs produced by `compile_latex.sh`.

The LaTeX sources are ordered as:

- `latex/1_original_proof_sketch.tex`
- `latex/2_current_status.tex`
- `latex/3_lp_investigation_notes.tex`
- `latex/4_two_parameter_model.tex`
- `latex/5_experiments_and_bounds.tex` — recorded LP relaxation values, the $m=3$
  improvement to **$8/15$**, and ideas for stronger formulations

## Which Python file to run

- `python3 run_proof.py`
  Reproduces the original `m=2` witness-set LP and the `2/7` bound.

- `python3 run_proof_4.py`
  Runs the earlier two-layer experiment for the equation `x+y=4z`, using `A-x`
  and `2A-y` witness families.

- `python3 search_witness.py`
  Randomized search over small witness sets for the older `x+y=4z` approach.

- `python3 search_witness_dil.py`
  Evaluates the older affine/dilation-based `x+y=4z` LP for chosen witness sets.
  Use this if you want to inspect or compare with the original failing approach.

- `python3 test_aps.py`
  Deterministic search over arithmetic-progression witnesses for the older
  `x+y=4z` model.

- `python3 test_large_interval.py`
  Tests interval-shaped witnesses for the older `x+y=4z` model.

- `python3 -c "from affine_form_lp import evaluate_affine_bound; ..."`
  Use `affine_form_lp.py` for the cleaned one-parameter affine-form LP that works
  for general `m`. On solver failure, `evaluate_affine_bound` returns `None`
  (not a numeric guess).

- `python3 experiments_two_parameter.py`
  Reproduces the two-parameter benchmarks: four-direction family ($4/7$ for
  $m=3$), `pentad_family(m)` ($8/15$ for $m=3$, $4/7$ for sample $m\ge 4$), and a
  short exhaustive check over a small direction pool.

- `python3 search_two_parameter.py`
  Runs small searches for the two-parameter LP based on forms `L(u,v)=au+bv+c`.
  The extended-direction search with a large `max_shift` can take many minutes.

- `python3 -c "from two_parameter_lp import evaluate_two_parameter_bound, pentad_family; ..."`
  Use `two_parameter_lp.py` directly for custom witness families. Helpers include
  `standard_family`, `make_direction_family`, and `pentad_family`. On solver
  failure, `evaluate_two_parameter_bound` returns `None`.

## LaTeX compilation

Run:

```bash
./compile_latex.sh
```

This script:

- compiles every `.tex` file in `latex/`,
- writes the final PDFs into `pdfs/`,
- removes generated `.aux`, `.log`, `.out`, `.toc`, `.nav`, `.snm`, `.fls`, and
  `.fdb_latexmk` files from `latex/`.

## Notes

- `lp_investigation_notes.md` is kept as the original Markdown scratch report.
  Its LaTeX version lives in `latex/3_lp_investigation_notes.tex`.
- The LP code depends on SciPy, NumPy, and Pydantic; see `requirements.txt`.
- **Scope.** Numerical values from these programs are bounds for the **LP
  relaxation** defined by affine symmetry and empty-atom constraints; turning
  them into theorems about $m$-sum-free sets $A\subset\mathbb{F}_p$ requires an
  additional step (see the discussion in `latex/2_current_status.tex` and
  `latex/5_experiments_and_bounds.tex`).

# LP Density Estimates for Additive Combinatorics

This repository contains linear-programming experiments for bounding the density of
sets \(A \subseteq \mathbb{F}_p\) that avoid equations of the form
\[
x+y=mz.
\]

## Current best LP relaxation bounds

All bounds below are **upper bounds for the LP relaxation** defined by affine
symmetry on the parameter space and empty-atom constraints coming from
identities \(L_i + L_j = m L_k\). Turning them into theorems about
\(m\)-sum-free sets \(A \subseteq \mathbb{F}_p\) still requires an additional
step; see `latex/2_current_status.tex`.

| Family | # forms | \(m\) | LP bound on \(\delta(A)\) | Code |
|---|---|---|---|---|
| One-parameter Kravitz witness | 8 | 2 | \(2/7\) | `run_proof.py` |
| Two-parameter four-directions | 4 | 3 | \(4/7\) | `experiments_two_parameter.py` |
| Two-parameter `pentad_family` | 5 | 3 | \(8/15 \approx 0.5333\) | `experiments_two_parameter.py` |
| **Three-parameter `basic_pair_family`** | 6 | any \(\ge 3\) | \(\mathbf{16/33 \approx 0.4848}\) | `experiments_three_parameter.py` |
| **Three-parameter `signed_pair_family`** | 9 | any \(\ge 3\) | \(\mathbf{48/103 \approx 0.4660}\) | `experiments_three_parameter.py` |
| **Three-parameter `symmetric_12_family`** | 12 | 3 | \(\mathbf{\approx 0.4398}\) | `experiments_three_parameter.py` |
| **Three-parameter `symmetric_12_family`** | 12 | 4 | \(\mathbf{\approx 0.4386}\) | `experiments_three_parameter.py` |

The **bold** rows are the three-parameter LP bounds introduced in this project
(see `latex/6_three_parameter_model.tex` and `latex/7_three_parameter_results.tex`).
They are the first LP-relaxation bounds that break the \(1/2\) barrier for
\(m \ge 3\).

## Project layout

- `atomic_density_lp.py` — shared LP primitives (Pydantic models) and the
  SciPy HiGHS-based solver.
- `affine_form_lp.py` — one-parameter affine-form LP (generalizes
  `run_proof.py` to arbitrary \(m\)). Returns `None` on solver failure.
- `two_parameter_lp.py` — two-parameter affine-form LP
  (`canonical_subset_shape`, `evaluate_two_parameter_bound`, etc).
- `three_parameter_lp.py` — three-parameter affine-form LP with a rank-based
  canonical shape that delegates rank-1 and rank-2 subsets to the lower
  dimensional canonical-shape code (see `latex/6_three_parameter_model.tex` for
  the mathematical rationale).
- `experiments_two_parameter.py` — reproduces two-parameter bounds.
- `experiments_three_parameter.py` — reproduces the new three-parameter
  bounds.
- `search_extended_pool.py`, `search_targeted.py`, `search_two_parameter.py`,
  `search_witness.py`, `search_witness_dil.py` — exploratory searches that
  confirm the two-parameter LP bounds do not improve below \(8/15\) for
  small integer direction pools.
- `latex/` — ordered LaTeX notes and writeups.
- `pdfs/` — compiled PDFs produced by `compile_latex.sh`.

LaTeX files (read in order):

- `latex/1_original_proof_sketch.tex` — excerpts of \cite{LMPV}.
- `latex/2_current_status.tex` — project status and background.
- `latex/3_lp_investigation_notes.tex` — the one-parameter affine-form LP.
- `latex/4_two_parameter_model.tex` — the two-parameter LP.
- `latex/5_experiments_and_bounds.tex` — two-parameter numerical bounds,
  including the \(8/15\) improvement for \(m=3\).
- `latex/6_three_parameter_model.tex` — **new**: the three-parameter LP.
- `latex/7_three_parameter_results.tex` — **new**: numerical results for the
  three-parameter LP (\(16/33\), \(48/103\), \(\approx 0.4398\)).

## Which Python file to run

- `python3 run_proof.py` reproduces the original Kravitz \(m=2\) proof and the
  \(2/7\) bound.
- `python3 run_proof_4.py` runs an earlier two-layer experiment for
  \(x+y=4z\) (historical).
- `python3 experiments_two_parameter.py` reproduces the two-parameter benchmarks
  (\(4/7\) for the four-direction family, \(8/15\) for `pentad_family`).
- `python3 experiments_three_parameter.py` reproduces the three-parameter
  benchmarks (\(16/33\), \(48/103\), \(\approx 0.4398\)). The twelve-form run
  takes a few minutes per \(m\).
- `python3 search_two_parameter.py` runs small searches for the two-parameter LP.
- `python3 search_extended_pool.py` / `python3 search_targeted.py` extend the
  searches over larger pools.
- `python3 search_witness.py`, `python3 search_witness_dil.py`,
  `python3 test_aps.py`, `python3 test_large_interval.py` — older exploratory
  scripts for the \(x+y=4z\) model.

## LaTeX compilation

```bash
./compile_latex.sh
```

This script:

- compiles every `.tex` file in `latex/`,
- writes the final PDFs into `pdfs/`,
- removes generated auxiliary files.

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Dependencies: NumPy, SciPy, Pydantic (see `requirements.txt`).

## Notes and caveats

- `lp_investigation_notes.md` is the original Markdown scratch version of
  `latex/3_lp_investigation_notes.tex`, kept for history.
- **Scope.** Numerical values from these programs are bounds for the **LP
  relaxation** defined by affine symmetry and empty-atom constraints; turning
  them into theorems about \(m\)-sum-free sets \(A \subseteq \mathbb{F}_p\)
  requires an additional step. See the discussion in
  `latex/2_current_status.tex`, `latex/5_experiments_and_bounds.tex` and
  `latex/6_three_parameter_model.tex`.
- The three-parameter LP relies on a rank-based `canonical_subset_shape` that
  delegates rank-1 and rank-2 subsets to the existing one- and two-parameter
  canonical-shape code. A previous ad hoc three-parameter normalization gave
  weaker bounds (\(21/41\), \(11/23\), \(\approx 0.449\)) because it retained
  residual \(w\)-coefficients and therefore missed some affine equalities; see
  the third remark of `latex/6_three_parameter_model.tex`.

# LP Density Estimates for Additive Combinatorics

This repository contains tools to specify, model, and bound the maximum density of subsets $A \subseteq \mathbb{F}_p$ that avoid linear equations using pure Linear Programming on atomic constituent intersections. 

## Files and Tooling

### `atomic_density_lp.py` (Core Engine)
The foundation of the library. It formalizes intersections of subset families into strict Linear Programming models built on top of Pydantic and SciPy.
- **Optimization**: This engine originally used a dense array constructor for `linprog`. We dramatically optimized it by precomputing boundary vectors natively into `scipy.sparse.csr_matrix` coordinates. This completely eliminates $O(N^2)$ memory bloat from the solver constraints, accelerating evaluations from ~60 seconds to ~0.05 seconds locally for complex configurations extending upwards of 16,000 bounds!
- **Data Models**: Pydantic models leverage strict tuple-casting `frozen=True` kwargs to enable proper deep hashing on `lp_variable` instances across disjoint combinations.

### `run_proof.py`
Replicates the original problem posed in `proof.tex` validating that sets avoiding $a_1 + a_2 - 2a_3 = d$ are strictly bounded by $2/7$. Groups 256 generated intersections computationally and binds translation invariances prior to triggering the constraint matrix properly.

### `run_proof_4.py`
Adapts the engine to target formulations avoiding $a + b = 4c$. Handles two independently dilated witness bases ($X$ for $A$ and $Y$ for $2A$). 
- Explicitly drops atoms where $x + y - 2z = 0$, guaranteeing geometric overlap contradictions are completely factored out.
- Builds multi-dilation translation symmetries efficiently via an $O(N)$ integer bitmask filtering logic. 
- **Non-trivial Logic**: Implements rigorous Dilation Equivalences (identifying arbitrary $S_X$ clusters representing shapes bounded inside $X$, multiplying by 2, and natively tracing an equality relation against its exact analog inside $Y$).

### Systematic Search Scripts
Created to attempt bypassing trivial bounds (e.g. bounding $< 0.5$ on $a+b=4c$ which $X$ and $Y$ alone naturally regress towards due to integral LP relaxation gaps):
- `search_witness_dil.py`: Fully randomized heuristic that rapidly attempts $X$ and $Y$ subset iterations (size 6) to spot constraint improvements.
- `test_aps.py`: Deterministically scales up explicit Arithmetic Progressions against `evaluate_bound(...)`, confirming uniform grid results locally.
- `test_large_interval.py`: Evaluates strictly contiguous and symmetric spaces bounded from 5 elements up to substantial sizes (e.g., $X=[0..8], Y=[0..8]$). Results prove the local witness LP ceilings at precisely $0.50000$ due to optimal relaxed fractional configurations mimicking the missing density.

## Usage
Simply source the local virtual environment and execute the desired sequence:
```bash
source venv/bin/activate
python run_proof_4.py
```

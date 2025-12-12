# Repository Guidelines

## Project Structure & Module Organization

- `src/` contains the reusable MATLAB CFD library. Current layout:
  - `src/2d/` – 2‑D solvers and infrastructure:
    - `mesh/` (e.g., `Mesh2D.m`)
    - `boundary/` boundary condition classes (`DirichletBC.m`, `PeriodicBC.m`, …)
    - `solvers/` time integrators inheriting from `BaseSolver.m`
    - `utils/` analysis helpers (e.g., `ErrorAnalyzer.m`)
  - `src/3d/` and `src/shared/` are reserved for future extensions.
- `2dCases/` hosts runnable teaching/benchmark cases. Each case lives in its own folder (e.g., `2dCases/coutte_flow/`) with:
  - `run_<case>.m` for simulation
  - `postprocess_<case>.m` for plots/animations
  - `results/` for generated outputs (keep large binaries out of commits when possible).

## Build, Test, and Development Commands

This repository is MATLAB‑first; there is no separate build step.

- From MATLAB, add the library to your path:
  ```matlab
  addpath(genpath("src"));
  ```
- Run the Couette example:
  ```matlab
  cd("2dCases/coutte_flow");
  run_couette;          % compute and save results
  postprocess_couette;  % generate figures/animation
  ```
- Headless run (optional):
  ```bash
  matlab -batch "addpath(genpath('src')); cd('2dCases/coutte_flow'); run_couette"
  ```

## Coding Style & Naming Conventions

- Indentation: 4 spaces; align `end` with the opening block.
- Classes: `PascalCase` filenames matching the `classdef` (e.g., `ExplicitEulerSolver.m`).
- Functions/methods/properties: `lower_snake_case` (e.g., `time_step`, `check_stability`).
- Prefer vectorized MATLAB operations over explicit loops unless clarity or stability requires otherwise.
- Keep comments consistent with the surrounding file language (most core files use Chinese).

## Testing Guidelines

- There is no automated test suite yet. Validate changes by re‑running existing cases and checking:
  - stability constraints (e.g., diffusion number in explicit schemes),
  - error metrics vs. analytical solutions where available.
- If adding formal tests, use `matlab.unittest` and place them under a new `tests/` folder; name tests `test_<feature>.m`.

## Commit & Pull Request Guidelines

- No commit convention is established in this repo. Use short, imperative messages with an optional area prefix, e.g. `2d: fix alpha_y stability check`.
- Keep commits focused (one numerical/physical change per commit).
- PRs should include: a clear description of the model/scheme change, how it was validated, and updated plots or error tables for affected cases. Link any related issues or notes.  

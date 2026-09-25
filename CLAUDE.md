# CLAUDE.md

## Overview

Multi-language scientific computing library focused on quantum information, random matrix theory (RMT), spin chains, and tensor networks. Published at [github.com/carlospgmat03/libs](https://github.com/carlospgmat03/libs).

## C++ Toolbox (`cpp/`)

### Dependencies

Requires `libitpp-dev` and `libtclap-dev` (Ubuntu packages).

### Build Commands

```bash
cd cpp/

make                    # Clean + build test_spins (default)
make all                # Build all targets + documentation
make test_spins         # Standard spin chain simulation
make test_spins_fast    # Optimized build (-O3 -march=native)
make test_spins2d       # 2D spin chain variant
make test_rmt           # Random matrix theory tests
make testing            # Runs testing.cpp with cfp_math + itpp_ext_math
make documentation      # Generate Doxygen docs into ./doc/
make clean              # Remove all compiled binaries and docs
```

Compiler flags come from `itpp-config --cflags` and `itpp-config --static --libs`. The Makefile also includes `~/makefile` for system-wide configuration.

### Key Source Files

- `spinchain.cpp` / `spinchain.h` — Spin chain Hamiltonian construction and diagonalization
- `RMT.cpp` / `RMT.h` — Random matrix ensembles (GUE, GOE, CUE) and statistics
- `purity_RMT.cpp` / `purity_RMT.h` — Purity and entanglement calculations
- `itpp_ext_math.cpp` / `itpp_ext_math.h` — Extensions to the ITPP library
- `cfp_math.cpp` / `cfp_math.h` — Mathematical utilities
- `dev_random.cpp` — High-quality random number generation
- `test_spins.cpp` — Main entry point for spin chain simulations (TCLAP CLI)

## Mathematica Packages

### Setup

Add to `~/.Wolfram/Kernel/init.m`:
```mathematica
AppendTo[$Path, FileNameJoin[{"/home/" <> ToString[$Username], "libs"}]];
```

Then load with `<<Quantum\`` or `<<SpinChain\`` in any notebook.

### Key Packages

- `Quantum.m` — Core quantum information primitives: partial trace, entanglement entropy, random states/matrices, quantum gates, bitwise operations, RMT ensembles
- `Quantum2.m` — Extended quantum operations
- `SpinChain.m` — Spin chain specific constructions
- `Carlos.m` — General mathematical/plotting utilities
- `TheoreticalRMT.m` — RMT visualization and analytical results
- `test.nb` — Test notebook for verifying the Mathematica libraries

`Quantum.m` is being rewritten into a modular paclet under `Quantum/` (`Kernel/`, `Tests/`, `PacletInfo.wl`) with qudit generalization and `VerificationTest` coverage; `Quantum_legacy.m` is the frozen pre-rewrite snapshot. See the rewrite section below.

## Julia (`julia/`)

- `quantum.jl` — Quantum operations
- `rmt.jl` — Random matrix implementations
- `customs.jl` — Shared utilities

## Architecture Notes

- The same algorithms are implemented in multiple languages: C++ for performance-critical simulations, Mathematica for symbolic/interactive analysis, MATLAB and Julia for alternative implementations
- C++ programs are standalone binaries with TCLAP command-line interfaces — check `--help` flags for available options
- The Fortran 90 modules in `modules/` are legacy code
- Maxima files (in `Maxima/`) were provided by Thomas Gorin for computing integrals of polynomials of unitary matrices with the Haar measure

---

# Quantum.wl Rewrite Project

This file gives Claude Code the context it needs to work on this repository.
Read this before making any changes.

## Project Goal

We are rewriting `Quantum.m`, an existing Mathematica package for quantum
mechanics/quantum information calculations, into a **leaner, modular Wolfram
Language package**. The rewrite has three objectives, in priority order:

1. **Compaction** — remove dead code, redundant helper functions, and
   duplicated logic. Consolidate near-identical functions (e.g. functions
   that differ only by matrix dimension or convention) into single
   parameterized functions wherever this doesn't hurt clarity or performance.
2. **Qudit generalization** — some functions that currently assume
   `d = 2` (qubits) should be generalized to accept an arbitrary local
   dimension `d` (qudits), defaulting to `d = 2` for backward compatibility.
3. **Test coverage** — every public function gets at least one
   `VerificationTest`, using Wolfram's built-in testing framework
   (`TestReport`), covering both the qubit case (regression against old
   behavior) and at least one qudit case (d = 3 or d = 4).

Do not change the public API's semantics for `d = 2` calls unless
explicitly asked — old notebooks that use this package should keep working.

## Repository Structure (target layout)

```
Quantum/
  Kernel/
    init.m                 -- BeginPackage / minimal loader, declares dependencies
    Core.wl                -- core data structures, state/operator representations
    Gates.wl                -- gate definitions (qubit + qudit generalized)
    States.wl                -- state prep (Bell, GHZ, W, generalized versions)
    Measures.wl              -- entropy, entanglement, fidelity, etc.
    Utilities.wl             -- partial trace, tensor/Kronecker helpers, basis conversion
  PacletInfo.wl             -- paclet metadata (if packaged as a paclet)
  Tests/
    CoreTests.wlt
    GatesTests.wlt
    StatesTests.wlt
    MeasuresTests.wlt
  README.md
  CLAUDE.md                 -- this file
```

If the current `Quantum.m` is a single flat file, the first task is to
split it along these lines — propose the split and confirm the function
list before doing a large mechanical move.

## Wolfram Language Conventions

- One `BeginPackage["Quantum`"]` / `EndPackage[]` per package, with each
  file inside using `Begin["`Private`"]` / `End[]`.
- Every exported symbol needs a `::usage` message written *before* the
  function definition, in the standard Mathematica style
  (`Foo::usage = "Foo[x] does thing."`).
- Use `Options[]` + `OptionsPattern[]` for functions with configurable
  behavior (e.g. `Dimension -> 2`) rather than overloading arity.
- Prefer `Module` over `Block` for local variables unless dynamic scoping
  is specifically needed.
- Guard public functions with argument patterns
  (`Foo[state_?VectorQ, d_Integer:2]`) rather than manual `If[!...,
  Message[...]; Return[$Failed]]` where a pattern will do — but keep
  explicit `Message[Foo::baddim, d]`-style errors for user-facing dimension
  mismatches, since these are the most common qudit-generalization bug.
- Numeric vs. symbolic: functions should work symbolically where feasible
  (e.g. symbolic `d`), but it's acceptable to require numeric `d` for
  functions where symbolic generalization isn't practical (say so in the
  usage message if so).
- Avoid global mutable state; if the old package uses global variables for
  configuration, convert to `Options`/optional arguments during the rewrite.

## Qudit Generalization Notes

- Generalized Pauli operators: use the clock (`Z_d`) and shift (`X_d`)
  operator convention (`Z_d |k> = ω^k |k>`, `X_d |k> = |k+1 mod d>`,
  `ω = Exp[2 Pi I / d]`). Name them something like `GeneralizedPauliX[d]`,
  `GeneralizedPauliZ[d]`, and keep `PauliX`, `PauliZ` etc. as `d=2`
  convenience wrappers if the old API used those names.
- Generalized Bell states: `Sum[Exp[2 Pi I j k/d], {k,0,d-1}] |k>|k>` family;
  confirm normalization matches whatever convention the old package used
  for `d=2` before generalizing.
- CNOT generalizes to the qudit SUM gate: `SUM|j,k> = |j, j+k mod d>`.
- When in doubt about which generalization convention to use, ask rather
  than guessing — there are several inequivalent standard conventions in
  the literature (Weyl-Heisenberg vs. Galois-field-based for prime-power
  d, etc.), and picking the wrong one silently breaks downstream physics.

## Testing

- Use `VerificationTest[expr, expected, opts]` in `.wlt` files under `Tests/`.
- Every `.wlt` file should be runnable standalone via:
  ```
  wolframscript -file Tests/RunAll.wl
  ```
  where `RunAll.wl` loads the package and calls `TestReport["Tests/"]`,
  printing a pass/fail summary and exiting non-zero on failure (so it can
  be used in a pre-commit hook or CI).
- For every generalized function, include:
  - a regression test at `d = 2` checked against the known old-package
    output (paste in the literal expected matrix/state where practical),
  - at least one test at `d = 3` and one at `d = 4`,
  - one test for the "bad dimension" error path (e.g. `d` non-integer or
    `d < 2`).
- Keep test data (matrices, states) as explicit literal expressions in the
  test file rather than computed from the same function under test, to
  avoid tautological tests.

## Useful Commands

```bash
# Run the full test suite
wolframscript -file Tests/RunAll.wl

# Load the package interactively for manual checks
wolframscript -code 'Needs["Quantum`"]; Foo[...]'

# Syntax-check a single file without evaluating top-level code
wolframscript -file Kernel/Gates.wl -check   # (if using a checker script)
```

Adjust these once the actual test runner script exists — the first task
in this project should be writing `Tests/RunAll.wl`.

## Workflow Expectations

- Before large mechanical refactors (splitting files, renaming many
  functions), propose the plan and function list rather than doing it
  silently.
- When generalizing a qubit function to qudits, show the `d=2` case
  reduces to the original behavior explicitly (either in a comment or a
  test) — this is the easiest place to introduce a silent bug.
- Prefer several small, reviewable commits/diffs (per-file or per-feature)
  over one large rewrite dump.
- Flag anywhere the original `Quantum.m` code relies on a Mathematica
  version-specific feature or deprecated function.

## Notes
- `Quantum.m` is currently a single flat file.
- Use a normalization/phase convention that is currently there for Bell states and Pauli operators.
- This should become a proper **paclet**
- I will ask you to generalize to qudits one function at a time. 


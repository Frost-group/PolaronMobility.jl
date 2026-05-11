# PolaronMobility.jl

`PolaronMobility.jl` is a Julia package for variational polaron calculations,
mobility estimates, and frequency-dependent response. The package is organized
around one composable pipeline:

```julia
problem = VariationalProblem(model, trial)
result = solve(problem; temperatures, frequencies)
```

The model describes the electron-phonon Hamiltonian. The trial process is the
solvable stochastic process used in the Feynman-Jensen variational objective.
The solver optimizes trial parameters on a temperature grid and then evaluates
model-specific observables on a frequency grid.

## Implemented Pipelines

- `FrohlichModel + GaussianFeynmanTrial`: the literature-pinned one-mode
  Feynman path-integral variational polaron.
- `FrohlichModel + MultiGaussianTrial`: Martin-Frost-style finite fictitious
  oscillator modes.
- `FrohlichModel + ProfileGaussianTrial`: a positive profile-function Gaussian
  trial in the Adamowski, Gerlach, and Leschke functional-integral spirit.
- `FrohlichModel + NonlocalGaussianTrial`: an experimental nonlocal kernel
  scaffold for exploratory response calculations.
- `HolsteinModel + PoissonTrial`: a lattice Holstein polaron trialled by a
  continuous-time Markov chain (CTMC) hopping process.
- `PeierlsModel + PoissonTrial`: a standalone bond-coupled Peierls lattice
  polaron using the same CTMC trial.
- `CompositePolaronModel + PoissonTrial`: compatible lattice influence
  functionals, currently Holstein plus Peierls on one Poisson path measure.

All core calculations are dimensionless. Unitful conversion is deliberately
kept at the API boundary through `FrohlichMaterial`, `HolsteinMaterial`,
`PeierlsMaterial`, `material_to_problem`, and `material_units`.

## The Shared Shape

Fröhlich, Holstein, Peierls, and compatible composite workflows return results with the same top-level
fields:

```julia
result.problem
result.temperatures
result.frequencies
result.solutions
result.mobilities
result.responses
```

Fröhlich results also expose `result.zero_temperature`, because the
zero-temperature variational solution is the standard literature reference and
is useful even when the requested temperature grid starts at finite
temperature.

The objective minimized by `solve_variational(problem, beta)` is always

```julia
free_energy(trial, parameters, beta) +
entropy_cost(trial, parameters, beta) +
interaction_free_energy(model, trial, parameters, beta)
```

This decomposition is the extension contract: a new model/trial pair only needs
to define its parameters, bounds, objective terms, and observable hooks.

## Stochastic Process View

The package treats variational trials as probability measures over paths.

For Fröhlich polarons, phonons are integrated out to produce a retarded
self-interaction of a continuum particle path. Gaussian trials replace the
exact path measure by a Gaussian process with a prescribed memory kernel. The
original Feynman trial is the one-fictitious-oscillator Gaussian process.

For Holstein and Peierls polarons, the electron lives on a lattice. The `PoissonTrial`
replaces the exact lattice path measure by a symmetric CTMC with hopping rate
`rate`. Holstein local phonons generate a retarded attraction between times at
which the CTMC returns to the same lattice site. Peierls bond phonons generate
a retarded bond-order influence evaluated with the full periodic CTMC bridge.

This is why the two models can share the same Julia pipeline even though their
path spaces are different: Fröhlich integrates Gaussian displacements, while
lattice models integrate CTMC site and bond bridge correlations.

## Literature Map

The Fröhlich side follows the older path-integral variational tradition:

- Fröhlich introduced the continuum polar electron-phonon model.
- Feynman introduced the all-coupling path-integral variational solution.
- Schultz provided early self-energy, mass, and mobility comparisons.
- Osaka extended the variational free energy to finite temperature.
- Feynman, Hellwarth, Iddings, and Platzman (FHIP) built the linear-response
  mobility theory around Feynman's trial action.
- Kadanoff and Devreese-related optical absorption work clarify the
  low-temperature and frequency-response limits.
- Hellwarth and Biaggio introduced practical multimode effective-frequency
  schemes for real polar lattices.
- Frost (2017) turned the finite-temperature Feynman/Hellwarth machinery into
  a first-principles mobility workflow for halide perovskites.
- Martin and Frost (2022) generalized the fictitious oscillator sector to
  multiple phonon modes and frequency-dependent multimode response.
- Adamowski, Gerlach, and Leschke developed broader functional-integral
  Gaussian-memory viewpoints that motivate the profile-function trial.

The Holstein and Peierls side follows the lattice-polaron path-integral tradition:

- Holstein introduced molecular-crystal small-polaron transport.
- De Raedt and Lagendijk developed path-integral lattice-polaron treatments.
- Kornilovitch and collaborators showed how phonons can be integrated out to
  give continuous-time worldline Monte Carlo for lattice polarons.
- The current package implements a lightweight CTMC variational analogue of
  that worldline picture, including full-periodic Holstein/Peierls bridge free
  energies and a general-`d` CTMC first-return mobility/response formalism
  with exact Holstein cloud sidebands and Peierls current-vertex sidebands.

Validated Fröhlich kernels are regression-tested against the historical
literature values. Lattice mobility and response tests check mathematical
consistency and limiting behavior; Peierls v1 is explicitly documented as a
controlled closure awaiting stronger external benchmarks.

## Where To Go Next

- Use [Examples](@ref) for executable scripts covering the public workflows.
- Use [Scientific Discussion](@ref) for equations, path-integral structure,
  mobility kernels, and trial-family caveats.
- Use [API Reference](@ref) for exported symbols and docstrings.

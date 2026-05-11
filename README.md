# PolaronMobility.jl

`PolaronMobility.jl` is a compact Julia package for variational polaron calculations, DC mobility estimates, and frequency-dependent response. The core is dimensionless and model-agnostic: build a `VariationalProblem(model, trial)`, solve the variational objective, then evaluate model-specific observables on temperature and frequency grids.

Implemented pipelines:

- `FrohlichModel + GaussianFeynmanTrial` for continuum Fröhlich polarons.
- `FrohlichModel + MultiGaussianTrial` for finite-mode Martin/Frost-style Gaussian variational trials.
- `FrohlichModel + ProfileGaussianTrial` for general profile-function Gaussian variational trials.
- `FrohlichModel + NonlocalGaussianTrial` for experimental finite-basis nonlocal Gaussian kernels that do not validate energy bounds.
- `HolsteinModel + PoissonTrial` for lattice Holstein polarons with a continuous-time hopping-rate CTMC trial.
- `PeierlsModel + PoissonTrial` for bond-coupled Peierls lattice polarons.
- `CompositePolaronModel + PoissonTrial` for compatible lattice influence functionals, currently Holstein plus Peierls.

For continuum Fröhlich models, the production path follows the Feynman/FHIP/Hellwarth Gaussian-trial literature. For lattice Holstein, Peierls, and Holstein-Peierls models, the production path uses a CTMC variational bridge free energy and a CTMC first-return transport kernel, both dressed with exact phonon sidebands.

## Core Pipeline

All implemented pipelines use the same public shape:

```julia
using PolaronMobility

problem = frohlich_feynman_problem(coupling = 3.0)
result = solve(problem; temperatures = [0.0, 1.0], frequencies = [0.0, 1.0])

problem_h = holstein_poisson_problem(coupling = 1.5)
result_h = solve(problem_h; temperatures = [0.25, 1.0], frequencies = [0.0, 1.0])

problem_p = peierls_poisson_problem(coupling = 0.5)
result_p = solve(problem_p; temperatures = [0.25, 1.0], frequencies = [0.0, 1.0])
```

Every full result exposes:

```julia
result.problem
result.temperatures
result.frequencies
result.solutions
result.mobilities
result.responses
```

Fröhlich's results also expose `result.zero_temperature`, even if zero temperature was not included in the requested grid.

The generic variational objective is decomposed as:

```julia
objective(problem, parameters, beta) =
    free_energy(problem.trial, parameters, beta) +
    entropy_cost(problem.trial, parameters, beta) +
    interaction_free_energy(problem.model, problem.trial, parameters, beta)
```

Use `solve_variational(problem, beta)` for a single inverse-temperature optimisation.

## Fröhlich

For dimensionless inputs:

```julia
problem = frohlich_feynman_problem(
    coupling = 3.0,
    phonon_frequency = 1.0,
    dimension = 3,
)

result = solve(problem; temperatures = [0.0, 0.5, 1.0], frequencies = [0.0, 2.0])

v = feynman_v(result.zero_temperature)
w = feynman_w(result.zero_temperature)
mobility = result.mobilities[2].mobility
conductivity = result.responses[2, 2].conductivity
```

`solve_frohlich(alpha; ...)` is a convenience wrapper around `frohlich_feynman_problem(...)` and `solve(...)`.

For material inputs, `material_to_problem` builds the same generic Fröhlich problem, and `solve` interprets temperatures as Kelvin and frequencies as THz:

```julia
material = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
problem = material_to_problem(material)
result = solve(problem; temperatures = [0.0, 300.0], frequencies = [0.0, 3.0])

unitful = material_units(result)
unitful.mobility[2]
unitful.frequency[2]
```

All Fröhlich Gaussian trial families can be selected from a material-derived problem:

```julia
material_to_problem(material; trial = :feynman)
material_to_problem(material; trial = :multi_gaussian, modes = 2)
material_to_problem(material; trial = :profile_gaussian, matsubara_terms = 256)
material_to_problem(material; trial = :nonlocal_gaussian)
```

Holstein material inputs use a sibling `HolsteinMaterial` type. For Rubrene, the convenience constructor encodes the local Holstein parameters from Ordejón et al., Phys. Rev. B 96, 035202 (2017), Table II:

```julia
rubrene = rubrene_holstein_material(
    direction = :high_mobility,
    lattice_constant_angstrom = 7.2, # optional, needed for unitful mobility
)
problem = material_to_problem(rubrene)
result = solve(problem; temperatures = [300.0], frequencies = [0.0, 10.0])

unitful = material_units(result)
unitful.mobility[1]
lambda_holstein(rubrene)
```

This maps `J = 134.0 meV`, `E_H = 106.8 meV`, and `ω_H = 1208.9 cm^-1` to the reduced Holstein model `hopping = J/(hν)`, `phonon_frequency = 1`, and `coupling = sqrt(E_H/(hν))`.

Fröhlich kernels are literature-pinned by tests, including:

- `frohlich_energy`
- `frohlich_memory_function`
- `frohlich_complex_impedance`
- `frohlich_complex_conductivity`
- `frohlich_mobility`
- `fhip_low_temperature_mobility`
- `kadanoff_low_temperature_mobility`
- `hellwarth_mobility`

The material coupling helper is `frohlich_alpha(...)`.

### Multi-Mode Gaussian Trials

`GaussianFeynmanTrial` remains the default, literature-pinned one-mode Fröhlich trial. For exploratory finite-mode Gaussian variational calculations, use:

```julia
problem = frohlich_multi_gaussian_problem(
    coupling = 3.0,
    phonon_frequency = 1.0,
    modes = 2,
)

result = solve(problem; temperatures = [0.0, 0.5], frequencies = [0.0, 1.0])

multi_gaussian_v(result.zero_temperature)
multi_gaussian_w(result.zero_temperature)
```

Parameters are ordered as `w1, delta1, w2, delta2, ...`, with `v_i = w_i + delta_i`. For `modes = 1`, the finite-mode kernels reduce to the current Feynman formulas and are tested against the one-mode implementation. For more than one fictitious mode, this is a Martin/Frost-style finite-mode extension: useful for variational exploration, but not yet as extensively literature-pinned as the one-mode Fröhlich kernels.

### General Profile Gaussian Trials

For functional Gaussian trial-action work, use `ProfileGaussianTrial`:

```julia
problem = frohlich_profile_gaussian_problem(
    coupling = 1.0,
    basis_frequencies = [0.5, 1.0, 2.0, 4.0],
    matsubara_terms = 512,
)

result = solve(problem; temperatures = [0.0, 0.5], frequencies = [0.0, 1.0])
```

The optimised parameters are amplitudes in a positive profile

```math
\Gamma(\omega) = \sum_i \frac{a_i\nu_i^2}{\omega^2+\nu_i^2}.
```

This is the main Adamowski/Gerlach/Leschke-style general Gaussian functional entry point. The one-mode Feynman trial remains the most thoroughly literature-pinned production path.

The profile family contains Feynman's trial exactly. If the one-mode Feynman solution has parameters `v,w`, then

```math
\Gamma_F(\omega) = \frac{v^2-w^2}{\omega^2+w^2}
```

is represented by `basis_frequencies = [w]` and `a1 = (v^2-w^2)/w^2`. The test suite verifies that this embedding reproduces the Feynman displacement kernel, the entropy correction, the interaction, and the total objective. General profile calculations should improve energies only modestly unless the basis has been carefully validated.

### Experimental Nonlocal Gaussian Kernels

For general functional experiments inspired by nonlocal Gaussian trial-action approaches, the package also exposes a finite-basis kernel scaffold:

```julia
problem = frohlich_nonlocal_gaussian_problem(
    coupling = 1.0,
    basis_frequencies = [1.0, 2.0, 4.0],
    regularization = 1e-3,
)

result = solve(
    problem;
    temperatures = [0.0, 0.5],
    frequencies = [0.0, 1.0],
)
```

The optimised parameters are amplitudes `a1, a2, ...` in a positive nonlocal memory-kernel basis. This path is explicitly experimental: the interaction and response use the generic Gaussian mean-square-displacement kernel, while the entropy contribution is a configurable quadratic regularizer rather than the profile-functional log-determinant action. A lower nonlocal pseudo-objective is not evidence of a better variational upper bound; use `ProfileGaussianTrial` for energy comparisons.

## Holstein, Peierls, And Composite Lattice Models

The lattice implementation follows the same pipeline as the continuum Fröhlich implementation, but the underlying theory differs. The lattice trial is a continuous-time Markov chain (CTMC) with variational nearest-neighbour hopping rate `κ`. The free-energy side uses a Feynman-Jensen/relative-entropy CTMC bridge objective. The transport side uses the optimised `κ` in a first-return current-blip kernel dressed by exact Holstein and/or Peierls phonon sidebands.

```julia
problem = holstein_poisson_problem(
    hopping = 1.0,
    phonon_frequency = 1.0,
    coupling = 1.5,
    dimension = 1,
)

result = solve(problem; temperatures = [0.25, 0.5, 1.0], frequencies = [0.0, 1.0])

rate = result.solutions[1].rate
einstein = result.mobilities[1].mobility_einstein
projected_mobility = result.mobilities[1].mobility
mobility_factor = result.mobilities[1].mobility_factor
complex_mobility = result.responses[2, 1].mobility
```

`solve_holstein(...)` is the convenience wrapper.

### Lattice Examples

For rubrene Holstein parameters from Ordejón Table II:

```julia
rubrene_h = rubrene_holstein_material(lattice_constant_angstrom = 7.2)
problem_h = material_to_problem(rubrene_h)
lambda_holstein(rubrene_h)
result_h = solve(problem_h; temperatures = [300.0], frequencies = [0.0, 10.0])
unitful_h = material_units(result_h)
unitful_h.mobility[1]
```

Peierls support is a standalone bond-coupled lattice model:

```julia
problem = peierls_poisson_problem(
    hopping = 1.0,
    phonon_frequency = 0.8,
    coupling = 0.35,
    dimension = 1,
)

result = solve(problem; temperatures = [0.5, 1.0], frequencies = [0.0, 0.5])
rate = result.solutions[1].rate
```

For rubrene Peierls parameters from Ordejón Table II:

```julia
rubrene_p = rubrene_peierls_material(lattice_constant_angstrom = 7.2)
problem_p = material_to_problem(rubrene_p)
lambda_peierls(rubrene_p)
result_p = solve(problem_p; temperatures = [300.0], frequencies = [0.0, 10.0])
unitful_p = material_units(result_p)
unitful_p.mobility[1]
```

This maps `J = 134.0 meV`, `E_P = 21.9 meV`, and `ω_P = 117.9 cm^-1` to `coupling = sqrt(E_P/(hν_P))`.

Compatible lattice models can be combined on one CTMC path space:

```julia
holstein = holstein_poisson_problem(hopping = 1.0, phonon_frequency = 1.0, coupling = 0.4)
peierls = peierls_poisson_problem(hopping = 1.0, phonon_frequency = 0.8, coupling = 0.35)

problem = VariationalProblem(
    combine_models(holstein.model, peierls.model),
    PoissonTrial(bare_hopping = 1.0),
)
result = solve(problem; temperatures = [1.0], frequencies = [0.0, 0.5])
```

For rubrene, `rubrene_holstein_peierls_problem(...)` builds both materials, converts them to a shared Holstein phonon reference frequency, combines the models, and uses one `PoissonTrial`.

```julia
problem_hp = rubrene_holstein_peierls_problem(lattice_constant_angstrom = 7.2)
result_hp = solve(problem_hp; temperatures = [300.0], frequencies = [0.0, 10.0])
unitful_hp = material_units(result_hp)
unitful_hp.mobility[1]
```

Composing incompatible path spaces, such as Fröhlich Gaussian plus Holstein Poisson, results in an `ArgumentError`.

### Holstein-Peierls Sidebands

For independent site and bond phonons, the transport cloud is the product

```math
C_{HP}(t)=C_H(t)\mathcal C_P(t).
```

The sideband list is therefore the convolution of the Holstein Franck-Condon sidebands with the Peierls vertex sidebands.

For transport-focused studies, the package also provides a dedicated guide-style
sweeps:

```julia
rows = holstein_peierls_transport_sweep(
    problem;
    temperatures = 0.05:0.05:0.5,
    frequencies = [0.0, 0.25, 0.5],
    kappa_source = :zero_temperature,
)
```

These helpers reuse one frozen `T = 0` optimised rate by default, cache repeated first-return evaluations over reused frequency shifts, and return flat
rows containing `mobility`, `mobility_einstein`, `mobility_factor`,
`conductivity_real`, `conductivity_imag`, and sideband normalization
diagnostics.

## Limiting Regimes

The examples and tests include coded checks for simple limits:

- Fröhlich weak coupling has `E ≈ -α`, while strong coupling drives `w -> ω` and `v ∼ 4α²/(9π)`.
- Fröhlich high-temperature Hellwarth mobility decreases as thermal phonon occupation grows; adiabaticity is interpreted through the reduced phonon scale.
- Holstein zero coupling returns `κ = J`; weak coupling gives the bare-rate interaction correction; stronger coupling suppresses `κ` and approaches a localised shift `-g²/ω0` plus exponentially narrowed transport.
- Holstein high temperature broadens the finite-temperature sideband distribution and modifies the CTMC-projected mobility through the \(\beta\kappa_\star\) prefactor and thermal sideband weights.
- Peierls zero coupling recovers the bare Poisson walk or the other composite component; weak Peierls corrections scale as `-g_P²`.
- Peierls antiadiabatic phonons produce short bond memory, while adiabatic phonons produce long-lived bond modulation and dynamic-disorder-like behaviour.
- Holstein-Peierls sidebands are a convolution of the Holstein Franck-Condon ladder and the Peierls current-vertex channels.

## Theory And Literature Guide

The Fröhlich implementation follows the historical path-integral arc: Fröhlich's continuum Hamiltonian, Feynman's all-coupling Gaussian variational solution, Schultz's early self-energy/mass/mobility comparisons, Osaka's finite-temperature free energy, FHIP mobility, Devreese optical absorption and memory-function work, Hellwarth and Biaggio's multimode effective-frequency treatment, Frost's 2017 first-principles halide-perovskite workflow, and Martin and Frost's 2022 multimode extension.

The general Gaussian profile trial follows the broader functional-integral viewpoint associated with Adamowski, Gerlach, and Leschke. In the current package, this is implemented as a finite positive profile basis with an analytic rational-kernel displacement decomposition.

The Holstein and Peierls implementations follow the lattice-polaron worldline tradition of Holstein, De Raedt, Lagendijk, Kornilovitch, and collaborators. In this package, the lattice free energy is approximated by lightweight CTMC bridge variational kernels, while lattice mobility and optical response are approximated by CTMC first-return kernels dressed by exact Holstein clouds and/or Peierls current-vertex sidebands. This is a compact deterministic alternative to continuous-time quantum Monte Carlo, not a replacement for full numerically exact lattice-polaron simulation.

## Optimizer Options

`OptimizerOptions` contains only generic numerical controls:

```julia
options = OptimizerOptions(
    lower = [1e-8, 0.0],
    upper = [60.0, 60.0],
    initial_parameters = [2.8, 0.3],
    multistart = true,
    warm_start = true,
    adaptive_bounds = true,
    quadrature_rtol = 1e-5,
)
```

Trial-specific defaults live on trial constructors such as `GaussianFeynmanTrial` and `PoissonTrial`.

## Continuation Sweeps

Continuation helpers run forward and backwards warm-started branches and select the lower-free-energy row:

```julia
continued_frohlich_coupling_sweep(0.5:0.5:5.0; temperature = 1.0)
continued_holstein_adiabaticity_sweep([0.5, 1.0, 2.0]; coupling = 1.5, temperature = 0.5)
```

Rows are flat `NamedTuple`s intended for validation tables, CSV export, or plotting.

Frequency sweeps use the same convention, but optimise once per temperature and then evaluate the requested response grid:

```julia
rows = frohlich_frequency_sweep(
    0.0:0.25:4.0;
    coupling = 3.0,
    temperatures = [0.5],
)

rows_h = holstein_frequency_sweep(
    0.0:0.25:4.0;
    coupling = 1.5,
    temperatures = [0.5],
)

rows_p = peierls_frequency_sweep(
    0.0:0.25:4.0;
    coupling = 0.5,
    temperatures = [0.5],
)
```

Full solve results can be flattened explicitly:

```julia
solution_rows = solution_table(result)
mobility_rows = mobility_table(result)
response_rows = response_table(result)

write_sweep_csv("response.csv", response_rows)
```

Plot helpers are optional. Install and load `Plots.jl` to activate the package extension:

```julia
using Plots

plot_frequency_sweep(response_rows)
plot_response_components(response_rows)
plot_coupling_sweep(continued_frohlich_coupling_sweep(0.5:0.5:5.0; temperature = 1.0))
```

The core package does not depend on a plotting backend; without `Plots.jl`, the `plot_*` helpers throw a clear error.

## Examples

Executable demos live in `examples/`:

- `frohlich_basic.jl`: one-mode Fröhlich solve, kernels, mobility, response, and tables.
- `frohlich_trials.jl`: Feynman, multi-Gaussian, profile Gaussian, and nonlocal Gaussian trials.
- `material_workflows.jl`: single-mode and multimode materials, units, and all material-derived trial choices.
- `holstein_basic.jl`: Holstein Poisson CTMC solve, return probabilities, mobility, and response.
- `peierls_basic.jl`: Peierls Poisson solve plus Holstein+Peierls composition.
- `full_lattice_free_energy.jl`: full periodic Holstein/Peierls free-energy kernels and bridge correlators.
- `lattice_mobility_response.jl`: CTMC first-return sideband mobility response for Holstein, Peierls, and composites.
- `holstein_peierls_mobility_comparison.jl`: frequency-by-frequency Holstein, Peierls, and composite mobility-factor comparison.
- `model_limits.jl`: coded weak/strong/high-temperature/adiabaticity sanity checks.
- `sweeps_tables_plots.jl`: continuation sweeps, frequency sweeps, CSV output, and optional plots.

## Installation And Tests

From the package directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

The test suite checks literature regressions, API shape, optimiser behaviour, sweep shape, exported docstrings, and README examples.

## Adding A Model

To add a new model/trial pair:

1. Define `struct MyModel <: AbstractPolaronModel`.
2. Define `struct MyTrial <: AbstractTrialProcess`.
3. Implement `parameter_names`, `initial_parameters`, and `parameter_bounds`.
4. Implement `free_energy`, `entropy_cost`, and `interaction_free_energy`.
5. Implement internal hooks `solution_result`, `mobility_result`, and `response_result`.
6. Add a `solve(problem::VariationalProblem{MyModel,MyTrial}; temperatures, frequencies, options)` method that calls the shared grid solver.

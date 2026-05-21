# Examples

The scripts in `examples/*.jl` are executable demos and are also run by the
test suite as smoke tests. They are intentionally compact, but together they
exercise the package's main workflows.

## Fröhlich Basic Workflow

Run:

```julia
include("examples/frohlich_basic.jl")
```

This demonstrates:

- Building `frohlich_feynman_problem`.
- Solving zero and finite reduced temperatures.
- Reading optimized `v,w`, mobility, memory function, impedance, and
  conductivity.
- Calling `solve_variational` directly for a single inverse temperature.
- Flattening typed results with `solution_table`, `mobility_table`, and
  `response_table`.

The core pattern is:

```julia
problem = frohlich_feynman_problem(coupling = 3.0)
result = solve(problem; temperatures = [0.0, 0.5], frequencies = [0.0, 1.0])
```

## Fröhlich Trial Families

Run:

```julia
include("examples/frohlich_trials.jl")
```

This compares:

- `GaussianFeynmanTrial`: the one-fictitious-oscillator production path.
- `MultiGaussianTrial`: several fictitious oscillator modes.
- `ProfileGaussianTrial`: a positive profile-function Gaussian trial.
- `NonlocalGaussianTrial`: an experimental finite-basis kernel scaffold.

The profile trial is useful because it contains the Feynman trial as a special
case. If Feynman's optimized parameters are `v,w`, then

```math
\Gamma_F(\omega) = \frac{v^2-w^2}{\omega^2+w^2}
```

is represented by a one-basis profile with `basis_frequencies = [w]` and
`a1 = (v^2 - w^2)/w^2`. The example prints this embedding as a sanity check.

A more general profile calculation looks like:

```julia
problem = frohlich_profile_gaussian_problem(
    coupling = 2.0,
    basis_frequencies = [0.75, 1.5],
    matsubara_terms = 128,
)
```

`profile_function(trial, parameters, omega)` and
`mean_square_displacement(trial, parameters, tau, beta)` are exported so the
kernel itself can be inspected.

`NonlocalGaussianTrial` is intentionally not shown as a validated lower-energy
variational improvement. It uses a regularized kernel objective rather than the
profile log-determinant entropy, so its pseudo-objective is not directly
comparable to Feynman or profile variational energies.

## Material Workflows

Run:

```julia
include("examples/material_workflows.jl")
```

This demonstrates:

- Single-mode material construction from dielectric constants, effective mass,
  and phonon frequency.
- Multimode material construction from phonon frequencies, infrared
  activities, and cell volume.
- `material_to_problem` for every Fröhlich trial family.
- Rubrene `HolsteinMaterial` construction from local Holstein parameters.
- Rubrene `PeierlsMaterial` construction from bond-coupling parameters.
- Combined Rubrene Holstein + Peierls composition on one Poisson trial.
- Conversion back to common units through `material_units`.

Material-derived problems keep Kelvin and THz at the `solve` boundary:

```julia
material = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
problem = material_to_problem(material; trial = :feynman)
result = solve(problem; temperatures = [0, 300], frequencies = [0, 3])
```

All Fröhlich trial families can be selected from the same material:

```julia
material_to_problem(material; trial = :feynman)
material_to_problem(material; trial = :multi_gaussian, modes = 2)
material_to_problem(material; trial = :profile_gaussian, matsubara_terms = 256)
material_to_problem(material; trial = :nonlocal_gaussian)
```

For profile and nonlocal trials, default basis frequencies are the material's
reduced phonon frequencies. This means multimode material-derived profile
problems solve without manual basis construction.

Rubrene can be mapped into the Holstein-Poisson pipeline from the Ordejón
et al. local-coupling data:

```julia
rubrene = rubrene_holstein_material(lattice_constant_angstrom = 7.2)
problem = material_to_problem(rubrene)
result = solve(problem; temperatures = [300], frequencies = [0, 10])
units = material_units(result)
```

The Holstein mapping uses `J = 134.0 meV`, `E_H = 106.8 meV`, and
`ω_H = 1208.9 cm^-1`. The Peierls mapping uses the same transfer integral
with `E_P = 21.9 meV` and `ω_P = 117.9 cm^-1`:

```julia
rubrene_p = rubrene_peierls_material(lattice_constant_angstrom = 7.2)
problem_p = material_to_problem(rubrene_p)
result_p = solve(problem_p; temperatures = [300], frequencies = [0, 10])
```

To combine local and bond couplings, use the convenience builder:

```julia
problem_hp = rubrene_holstein_peierls_problem(lattice_constant_angstrom = 7.2)
result_hp = solve(problem_hp; temperatures = [300], frequencies = [0, 10])
```

## Peierls And Limits

Run:

```julia
include("examples/peierls_basic.jl")
include("examples/model_limits.jl")
```

These scripts demonstrate standalone Peierls variational solves, Holstein plus
Peierls model composition, and coded limit checks: Fröhlich weak/strong
coupling markers, Holstein zero-coupling and high-temperature behavior,
Peierls `g_P^2` weak-coupling scaling, and adiabatic/antiadiabatic ratios.

## Holstein Basic Workflow

Run:

```julia
include("examples/holstein_basic.jl")
```

This demonstrates:

- Building `holstein_poisson_problem`.
- Optimizing the CTMC hopping rate.
- Reading Einstein and CTMC first-return projected mobilities.
- Evaluating return probabilities.
- Evaluating frequency-dependent complex mobility and mobility factors.

The core pattern is:

```julia
problem = holstein_poisson_problem(coupling = 1.5, hopping = 1.0)
result = solve(problem; temperatures = [0.25, 0.5], frequencies = [0.0, 1.0])
```

At zero coupling, the optimized rate returns the bare hopping. Finite coupling
adds a full-periodic bridge free-energy correction and a CTMC first-return sideband
response through `mobility_factor`, `impedance`, `conductivity`, and
`mobility` fields.

## CTMC Lattice Transport

Run:

```julia
include("examples/full_lattice_free_energy.jl")
include("examples/lattice_mobility_response.jl")
include("examples/holstein_peierls_mobility_comparison.jl")
include("examples/general_d_lattice_transport.jl")
```

These scripts demonstrate the active lattice theory: periodic phonon kernels,
CTMC site/bond bridge correlations, the general-`d` CTMC first-return kernel,
exact Holstein blip sidebands, Peierls current-vertex sidebands, and
Holstein-Peierls sideband convolution.

`general_d_lattice_transport.jl` also demonstrates the guide-style helpers

```julia
holstein_transport_sweep(...)
peierls_transport_sweep(...)
holstein_peierls_transport_sweep(...)
```

which default to `kappa_source = :zero_temperature` so one optimized
zero-temperature rate can be reused across a temperature or frequency study.

## Sweeps, Tables, CSV, And Plotting

Run:

```julia
include("examples/sweeps_tables_plots.jl")
```

This demonstrates:

- Fröhlich coupling sweeps.
- Holstein temperature sweeps.
- Fröhlich and Holstein frequency sweeps.
- CSV export through `write_sweep_csv`.
- Optional plots through the `Plots.jl` extension.

Frequency sweeps optimize once per temperature and then evaluate every
frequency:

```julia
rows = frohlich_frequency_sweep(
    [0.0, 0.5, 1.0];
    coupling = 1.0,
    temperatures = [0.5],
)

write_sweep_csv("frequency.csv", rows)
```

If `Plots.jl` is loaded, the plotting extension activates:

```julia
using Plots

plot_frequency_sweep(rows)
plot_response_components(rows)
```

The plotting API intentionally consumes flat rows rather than result structs.
That keeps the core solver independent from plotting choices and makes CSV,
notebook, and batch-report workflows use the same data representation.

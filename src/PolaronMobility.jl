"""
    PolaronMobility

Variational polaron models, mobility estimates, and frequency-dependent
transport response.

Core calculations are dimensionless and organized around a generic
`VariationalProblem(model, trial)` pipeline. The package currently provides
continuum Fröhlich Gaussian trials, lattice Holstein and Peierls Poisson
trials, and compatible lattice model composition.
"""
module PolaronMobility

using LinearAlgebra: eigvals, eigen, SymTridiagonal
using Optim
using QuadGK
using SpecialFunctions: besselix, besseli, ellipk, gamma, loggamma, coth
using Unitful

export PolaronResult,
       AbstractPolaronModel,
       AbstractContinuumModel,
       AbstractLatticeModel,
       AbstractMaterial,
       AbstractTrialProcess,
       AbstractJumpTrial,
       AbstractGaussianTrial,
       VariationalProblem,
       VariationalResult,
       ObservableResult,
       FrohlichModel,
       GaussianFeynmanTrial,
       MultiGaussianTrial,
       NonlocalGaussianTrial,
       ProfileGaussianTrial,
       HolsteinModel,
       PeierlsModel,
       CompositePolaronModel,
       PoissonTrial,
       FrohlichMaterial,
       HolsteinMaterial,
       PeierlsMaterial,
       FrohlichSolution,
       LatticeSolution,
       FrohlichMobilityResult,
       FrohlichResponseResult,
       LatticeMobilityResult,
       LatticeResponseResult,
       OptimizerOptions,
       EnergyComponents,
       MaterialUnits,
       FrohlichResultUnits,
       LatticeResultUnits,
       solve,
       solve_variational,
       solve_frohlich,
       solve_holstein,
       solve_peierls,
       frohlich_feynman_problem,
       frohlich_multi_gaussian_problem,
       frohlich_nonlocal_gaussian_problem,
       frohlich_profile_gaussian_problem,
       holstein_poisson_problem,
       peierls_poisson_problem,
       combine_models,
       return_probability,
       periodic_phonon_kernel,
       lattice_q0,
       lattice_q1,
       lattice_current_kernel,
       site_return_bridge,
       bond_order_bridge,
       bond_current_bridge,
       holstein_integral_d,
       peierls_integral_d,
       lattice_green_function_d,
       first_return_laplace_d,
       lattice_holstein_phonon_factor,
       lattice_peierls_phonon_factor,
       holstein_transport_sidebands,
       peierls_transport_sidebands,
       holstein_peierls_transport_sidebands,
       lattice_mobility_factor,
       lattice_mobility,
       lattice_conductivity,
       lattice_impedance,
       objective,
       free_energy,
       entropy_cost,
       interaction_free_energy,
       parameter_names,
       feynman_v,
       feynman_w,
       multi_gaussian_v,
       multi_gaussian_w,
       mean_square_displacement,
       profile_function,
       frequency_sweep,
       frohlich_frequency_sweep,
       holstein_frequency_sweep,
       peierls_frequency_sweep,
       lattice_transport_sweep,
       holstein_transport_sweep,
       peierls_transport_sweep,
       holstein_peierls_transport_sweep,
       solution_table,
       mobility_table,
       response_table,
       sweep_table,
       write_sweep_csv,
       plot_coupling_sweep,
       plot_temperature_sweep,
       plot_adiabaticity_sweep,
       plot_frequency_sweep,
       plot_response_components,
       continued_frohlich_temperature_sweep,
       continued_frohlich_coupling_sweep,
       continued_frohlich_adiabaticity_sweep,
       continued_holstein_temperature_sweep,
       continued_holstein_coupling_sweep,
       continued_holstein_adiabaticity_sweep,
       continued_peierls_temperature_sweep,
       continued_peierls_coupling_sweep,
       continued_peierls_adiabaticity_sweep,
       material_to_problem,
       material_units,
       mobility_cm2_per_v_s,
       lattice_mobility_cm2_per_v_s,
       energy_meV,
       frequency_THz,
       wavenumber_THz,
       wavenumber_meV,
       radius_angstrom,
       lambda_holstein,
       lambda_peierls,
       rubrene_holstein_material,
       rubrene_peierls_material,
       rubrene_holstein_peierls_problem,
       frohlich_alpha,
       dielectric_ionic_mode,
       dielectric_ionic_total,
       hellwarth_b_scheme,
       hellwarth_a_scheme,
       frohlich_energy,
       frohlich_memory_function,
       frohlich_complex_impedance,
       frohlich_complex_conductivity,
       frohlich_mobility,
       fhip_low_temperature_mobility,
       kadanoff_low_temperature_mobility,
       hellwarth_mobility

include("Types.jl")
include("Units.jl")
include("HellwarthTheory.jl")
include("ContinuumKernels.jl")
include("VariationalSolver.jl")
include("LatticeKernels.jl")
include("LatticeTransport.jl")
include("Material.jl")
include("FrohlichPolaron.jl")
include("HolsteinPolaron.jl")
include("PeierlsPolaron.jl")
include("Sweeps.jl")

end

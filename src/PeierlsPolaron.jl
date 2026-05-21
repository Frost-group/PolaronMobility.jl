"""
    peierls_poisson_problem(; hopping=1, phonon_frequency=1, coupling, dimension=1)

Build a Peierls variational problem using the same symmetric Poisson hopping
trial used by the Holstein lattice model. The Peierls coupling modulates bonds
instead of site density.
"""
function peierls_poisson_problem(;
    hopping::Real = 1.0,
    phonon_frequency::Real = 1.0,
    coupling::Real,
    dimension::Integer = 1,
    effective_frequency::Real = phonon_frequency,
    material::Union{Nothing,PeierlsMaterial} = nothing,
)
    model = PeierlsModel(
        dimension = dimension,
        hopping = hopping,
        phonon_frequency = phonon_frequency,
        coupling = coupling,
        effective_frequency = effective_frequency,
        material = material,
    )
    trial = PoissonTrial(dimension = dimension, bare_hopping = hopping)
    return VariationalProblem(model, trial)
end

"""
    combine_models(models...)

Combine compatible lattice polaron models so they share one `PoissonTrial`.
The v1 composition path supports Holstein and Peierls components with the same
dimension and reduced hopping. Continuum Fröhlich/Gaussian models deliberately
cannot be combined with lattice Poisson models.
"""
function combine_models(models::AbstractLatticeModel...)
    length(models) >= 2 || throw(ArgumentError("combine_models requires at least two component models."))
    all(model -> model isa Union{HolsteinModel,PeierlsModel}, models) ||
        throw(ArgumentError("v1 model composition supports only HolsteinModel and PeierlsModel components with PoissonTrial."))
    return CompositePolaronModel(Tuple(models))
end

function combine_models(models::AbstractPolaronModel...)
    throw(ArgumentError("combine_models only supports compatible lattice models that share one PoissonTrial path space."))
end

initial_parameters(problem::VariationalProblem{PeierlsModel,PoissonTrial}, beta::Real) = initial_parameters(problem.trial)
initial_parameters(problem::VariationalProblem{<:CompositePolaronModel,PoissonTrial}, beta::Real) = initial_parameters(problem.trial)

"""
    peierls_bond_correlation(model, trial, parameters, tau)

Peierls full-periodic bond-order bridge correlation. At finite `beta` this is
the CTMC loop bridge; at zero temperature it approaches
`2d * (q0(τ) + q1(τ))`.
"""
function peierls_bond_correlation(model::PeierlsModel, trial::PoissonTrial, parameters::AbstractVector, tau::Real)
    rate = poisson_rate(trial, parameters)
    return bond_order_bridge(rate, model.dimension, tau, Inf)
end

"""
    interaction_free_energy(model::PeierlsModel, trial::PoissonTrial, parameters, beta; rtol=1e-4)

Full-periodic retarded Peierls bond influence evaluated against the CTMC
bond-order bridge. At zero temperature this reduces to
`-0.5 g² peierls_integral_d(κ, d, ω)`.
"""
function interaction_free_energy(
    model::PeierlsModel,
    trial::PoissonTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    iszero(model.coupling) && return 0.0
    rate = poisson_rate(trial, parameters)
    if beta == Inf
        return -0.5 * model.coupling^2 * peierls_integral_d(rate, model.dimension, model.phonon_frequency)
    end
    finite_temperature_integrand(tau) = periodic_phonon_kernel(tau, beta, model.phonon_frequency) *
                                        bond_order_bridge(rate, model.dimension, tau, beta)
    return -0.5 * model.coupling^2 * quadgk(finite_temperature_integrand, 0, beta; rtol = rtol)[1]
end

function peierls_interaction_integrand(model::PeierlsModel, trial::PoissonTrial, parameters, tau, beta)
    rate = poisson_rate(trial, parameters)
    beta == Inf && return 0.5 * exp(-model.phonon_frequency * tau) *
                          2 * model.dimension * (lattice_q0(rate, model.dimension, tau) + lattice_q1(rate, model.dimension, tau))
    return 0.5 * periodic_phonon_kernel(tau, beta, model.phonon_frequency) *
           bond_order_bridge(rate, model.dimension, tau, beta)
end

"""
    interaction_free_energy(composite, trial, parameters, beta; rtol=1e-4)

Composite lattice interaction free energy, implemented as the sum of component
influence functionals evaluated on the same Poisson path measure.
"""
function interaction_free_energy(
    model::CompositePolaronModel,
    trial::PoissonTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    return sum(component -> interaction_free_energy(component, trial, parameters, beta; rtol = rtol), model.models)
end

function diagnostics(
    model::PeierlsModel,
    trial::PoissonTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    rate = Float64(poisson_rate(trial, parameters))
    total_jump_rate = 2 * model.dimension * rate
    phonon_time = inv(model.phonon_frequency)
    return (;
        rate = rate,
        peierls_shift = model.coupling^2 / model.phonon_frequency,
        adiabatic_ratio = model.hopping / model.phonon_frequency,
        lambda_peierls = model.coupling^2 / (2 * model.dimension * model.hopping * model.phonon_frequency),
        total_jump_rate = Float64(total_jump_rate),
        mean_waiting_time = inv(max(total_jump_rate, 1e-15)),
        bond_correlation_at_phonon_time = Float64(bond_order_bridge(rate, model.dimension, phonon_time, beta)),
    )
end

function diagnostics(
    model::CompositePolaronModel,
    trial::PoissonTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    rate = Float64(poisson_rate(trial, parameters))
    pairs = Pair{Symbol,Any}[
        :rate => rate,
        :composite_components => length(model.models),
        :total_jump_rate => 2 * model.dimension * rate,
        :mean_waiting_time => inv(max(2 * model.dimension * rate, 1e-15)),
    ]
    for (index, component) in enumerate(model.models)
        component_diagnostics = diagnostics(component, trial, parameters, beta; options = options)
        prefix = Symbol("component", index)
        push!(pairs, Symbol(prefix, :_model) => Symbol(lowercase(string(nameof(typeof(component))))))
        for key in keys(component_diagnostics)
            key in (:rate, :total_jump_rate, :mean_waiting_time) && continue
            push!(pairs, Symbol(prefix, :_, key) => getproperty(component_diagnostics, key))
        end
    end
    return NamedTuple(pairs)
end

"""
    solve(problem::VariationalProblem{PeierlsModel,PoissonTrial}; temperatures, frequencies, options)

Solve a Peierls-Poisson variational problem over temperature and frequency
grids. Direct inputs are dimensionless. Material-derived problems interpret
temperatures as Kelvin and frequencies as THz.
"""
function solve(
    problem::VariationalProblem{PeierlsModel,PoissonTrial};
    temperatures = [1.0],
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
)
    return _solve_lattice_poisson_problem(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

"""
    solve(problem::VariationalProblem{CompositePolaronModel,PoissonTrial}; temperatures, frequencies, options)

Solve a compatible composite lattice model over temperature and frequency
grids using a shared Poisson trial.
"""
function solve(
    problem::VariationalProblem{<:CompositePolaronModel,PoissonTrial};
    temperatures = [1.0],
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
)
    return _solve_lattice_poisson_problem(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

function _solve_lattice_poisson_problem(
    problem::VariationalProblem{M,PoissonTrial};
    temperatures,
    frequencies,
    options::OptimizerOptions,
) where {M<:Union{PeierlsModel,CompositePolaronModel}}
    _assert_poisson_compatible(problem)
    solve_temperatures, solve_frequencies = _lattice_solve_grid(problem.model, temperatures, frequencies)
    temperature_values, frequency_values, zero_solution, solutions, mobilities, responses =
        _solve_grid(problem; temperatures = solve_temperatures, frequencies = solve_frequencies, options = options)
    return PolaronResult(problem, temperature_values, frequency_values, zero_solution, solutions, mobilities, responses)
end

function _assert_poisson_compatible(problem::VariationalProblem{PeierlsModel,PoissonTrial})
    problem.trial.dimension == problem.model.dimension ||
        throw(ArgumentError("PoissonTrial dimension must match PeierlsModel dimension."))
    isapprox(problem.trial.bare_hopping, problem.model.hopping; rtol = 1e-8, atol = 1e-12) ||
        throw(ArgumentError("PoissonTrial bare_hopping must match PeierlsModel hopping."))
    return nothing
end

function _assert_poisson_compatible(problem::VariationalProblem{<:CompositePolaronModel,PoissonTrial})
    problem.trial.dimension == problem.model.dimension ||
        throw(ArgumentError("PoissonTrial dimension must match CompositePolaronModel dimension."))
    isapprox(problem.trial.bare_hopping, problem.model.hopping; rtol = 1e-8, atol = 1e-12) ||
        throw(ArgumentError("PoissonTrial bare_hopping must match CompositePolaronModel hopping."))
    return nothing
end

function _lattice_solve_grid(model::PeierlsModel, temperatures, frequencies)
    model.material === nothing && return temperatures, frequencies
    reduced_T = [iszero(T) ? 0.0 : reduced_temperature(T, model.effective_frequency) for T in _as_vector(temperatures)]
    reduced_Ω = [reduced_frequency(Ω, model.effective_frequency) for Ω in _as_vector(frequencies)]
    return reduced_T, reduced_Ω
end

function _lattice_solve_grid(model::CompositePolaronModel, temperatures, frequencies)
    !_composite_has_material(model) && return temperatures, frequencies
    reduced_T = [iszero(T) ? 0.0 : reduced_temperature(T, model.effective_frequency) for T in _as_vector(temperatures)]
    reduced_Ω = [reduced_frequency(Ω, model.effective_frequency) for Ω in _as_vector(frequencies)]
    return reduced_T, reduced_Ω
end

function _composite_has_material(model::CompositePolaronModel)
    return any(component -> hasproperty(component, :material) && getproperty(component, :material) !== nothing, model.models)
end

"""
    solve_peierls(; hopping=1, phonon_frequency=1, coupling, dimension=1, temperatures=[1], frequencies=[0])

Convenience wrapper around `peierls_poisson_problem(...)` and `solve(...)`.
"""
function solve_peierls(;
    hopping::Real = 1.0,
    phonon_frequency::Real = 1.0,
    coupling::Real,
    dimension::Integer = 1,
    effective_frequency::Real = phonon_frequency,
    temperatures = [1.0],
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
)
    problem = peierls_poisson_problem(
        hopping = hopping,
        phonon_frequency = phonon_frequency,
        coupling = coupling,
        dimension = dimension,
        effective_frequency = effective_frequency,
    )
    return solve(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

"""
    solution_result(problem, temperature, variational)

Convert a generic `VariationalResult` into a Peierls/composite lattice solution
record.
"""
function solution_result(
    problem::VariationalProblem{M,PoissonTrial},
    temperature::Real,
    variational::VariationalResult,
) where {M<:Union{PeierlsModel,CompositePolaronModel}}
    rate = variational.diagnostics.rate
    return LatticeSolution(
        Float64(temperature),
        variational.beta,
        variational,
        Float64(rate),
        variational.free_energy,
        variational.entropy_cost,
        variational.interaction_energy,
        variational.diagnostics,
    )
end

"""
    peierls_mobility(problem, solution)

Compute Einstein and transport-corrected DC mobilities for optimized
Peierls/composite Poisson solutions using the CTMC first-return plus
sideband/vertex transport kernel.
"""
function mobility_result(
    problem::VariationalProblem{M,PoissonTrial},
    solution::LatticeSolution,
    options::OptimizerOptions,
) where {M<:Union{PeierlsModel,CompositePolaronModel}}
    return _lattice_dc_mobility_result(
        problem.model,
        solution.rate,
        solution.beta,
        solution.temperature;
        rtol = options.quadrature_rtol,
    )
end

"""
    peierls_response(problem, solution, frequency)

Frequency-dependent lattice transport response for optimized
Peierls/composite Poisson solutions using the CTMC first-return plus
sideband/vertex transport kernel. At zero temperature and finite frequency,
the reported reduced optical response is the transport kernel itself rather
than `βκ` times that kernel.
"""
function response_result(
    problem::VariationalProblem{M,PoissonTrial},
    solution::LatticeSolution,
    frequency::Real,
    options::OptimizerOptions,
) where {M<:Union{PeierlsModel,CompositePolaronModel}}
    return _lattice_response_result(
        problem.model,
        solution.rate,
        solution.beta,
        solution.temperature,
        frequency;
        rtol = options.quadrature_rtol,
    )
end

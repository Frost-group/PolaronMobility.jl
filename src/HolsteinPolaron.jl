"""
    holstein_poisson_problem(; hopping=1, phonon_frequency=1, coupling, dimension=1)

Build a Holstein variational problem using a Poisson hopping trial process. The
result can be passed to `solve_variational` for a single inverse temperature or
to `solve` for temperature/frequency grids. Material-derived problems pass
physical Kelvin/THz grids through `effective_frequency` and `material`.
"""
function holstein_poisson_problem(;
    hopping::Real = 1.0,
    phonon_frequency::Real = 1.0,
    coupling::Real,
    dimension::Integer = 1,
    effective_frequency::Real = phonon_frequency,
    material::Union{Nothing,HolsteinMaterial} = nothing,
)
    model = HolsteinModel(
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

parameter_names(::PoissonTrial) = [:rate]
initial_parameters(trial::PoissonTrial) = [trial.bare_hopping]
initial_parameters(problem::VariationalProblem{HolsteinModel,PoissonTrial}, beta::Real) = initial_parameters(problem.trial)
parameter_bounds(trial::PoissonTrial, beta::Real) = ([trial.lower_rate], [trial.upper_rate])

"""
    return_probability(trial::PoissonTrial, tau, parameters)

Return probability of the optimized symmetric Poisson walk at imaginary time
`tau`. Uses the exponentially-scaled Bessel function representation for
numerical stability.
"""
function return_probability(trial::PoissonTrial, tau::Real, parameters::AbstractVector)
    rate = poisson_rate(trial, parameters)
    return lattice_q0(rate, trial.dimension, tau)
end

function poisson_rate(trial::PoissonTrial, parameters::AbstractVector)
    length(parameters) == 1 || throw(ArgumentError("PoissonTrial expects one variational parameter: rate."))
    rate = parameters[1]
    rate >= 0 || throw(DomainError(rate, "rate must be non-negative."))
    return rate
end

"""
    free_energy(trial::PoissonTrial, parameters, beta)

Bare continuous-time random-walk energy for the symmetric nearest-neighbor
Poisson trial process.
"""
free_energy(trial::PoissonTrial, parameters::AbstractVector, beta::Real) =
    -2 * trial.dimension * poisson_rate(trial, parameters)

"""
    entropy_cost(trial::PoissonTrial, parameters, beta)

Donsker-Varadhan relative-entropy rate for changing the bare hopping rate to the
variational hopping rate.
"""
function entropy_cost(trial::PoissonTrial, parameters::AbstractVector, beta::Real)
    rate = poisson_rate(trial, parameters)
    iszero(rate) && return zero(rate)
    return 2 * trial.dimension * rate * log(rate / trial.bare_hopping)
end

"""
    interaction_free_energy(model::HolsteinModel, trial::PoissonTrial, parameters, beta; rtol=1e-4)

Full-periodic retarded local Holstein interaction evaluated against the CTMC
site-return bridge. At zero temperature this reduces to
`-g² holstein_integral_d(κ, d, ω)`.
"""
function interaction_free_energy(
    model::HolsteinModel,
    trial::PoissonTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    iszero(model.coupling) && return 0.0
    rate = poisson_rate(trial, parameters)
    if beta == Inf
        return -model.coupling^2 * holstein_integral_d(rate, model.dimension, model.phonon_frequency)
    end
    finite_temperature_integrand(tau) = periodic_phonon_kernel(tau, beta, model.phonon_frequency) *
                                        site_return_bridge(rate, model.dimension, tau, beta)
    return -0.5 * model.coupling^2 * quadgk(finite_temperature_integrand, 0, beta; rtol = rtol)[1]
end

function holstein_interaction_integrand(model::HolsteinModel, trial::PoissonTrial, parameters, tau, beta)
    rate = poisson_rate(trial, parameters)
    beta == Inf && return exp(-model.phonon_frequency * tau) * lattice_q0(rate, model.dimension, tau)
    return 0.5 * periodic_phonon_kernel(tau, beta, model.phonon_frequency) *
           site_return_bridge(rate, model.dimension, tau, beta)
end

function diagnostics(
    model::HolsteinModel,
    trial::PoissonTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    rate = Float64(poisson_rate(trial, parameters))
    total_jump_rate = 2 * model.dimension * rate
    long_time = max(10 / model.phonon_frequency, 10 / max(rate, 1e-12))
    return (;
        rate = rate,
        polaron_shift = model.coupling^2 / model.phonon_frequency,
        adiabatic_ratio = model.hopping / model.phonon_frequency,
        lambda_holstein = model.coupling^2 / (2 * model.dimension * model.hopping * model.phonon_frequency),
        total_jump_rate = Float64(total_jump_rate),
        mean_waiting_time = inv(max(total_jump_rate, 1e-15)),
        long_time_return_probability = Float64(return_probability(trial, long_time, parameters)),
    )
end

"""
    solve(problem::VariationalProblem{HolsteinModel,PoissonTrial}; temperatures, frequencies, options)

Solve a Holstein-Poisson variational problem over temperature and frequency
grids. Direct model inputs are dimensionless. If the model was built by
`material_to_problem(::HolsteinMaterial)`, temperatures are interpreted as
Kelvin and frequencies as THz.
"""
function solve(
    problem::VariationalProblem{HolsteinModel,PoissonTrial};
    temperatures = [1.0],
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
)
    solve_temperatures, solve_frequencies = _holstein_solve_grid(problem.model, temperatures, frequencies)
    temperature_values, frequency_values, zero_solution, solutions, mobilities, responses =
        _solve_grid(problem; temperatures = solve_temperatures, frequencies = solve_frequencies, options = options)
    return PolaronResult(problem, temperature_values, frequency_values, zero_solution, solutions, mobilities, responses)
end

function _holstein_solve_grid(model::HolsteinModel, temperatures, frequencies)
    model.material === nothing && return temperatures, frequencies
    reduced_T = [iszero(T) ? 0.0 : reduced_temperature(T, model.effective_frequency) for T in _as_vector(temperatures)]
    reduced_Ω = [reduced_frequency(Ω, model.effective_frequency) for Ω in _as_vector(frequencies)]
    return reduced_T, reduced_Ω
end

"""
    solve_holstein(; hopping=1, phonon_frequency=1, coupling, dimension=1, temperatures=[1], frequencies=[0])

Convenience wrapper around `holstein_poisson_problem(...)` and `solve(...)`.
"""
function solve_holstein(;
    hopping::Real = 1.0,
    phonon_frequency::Real = 1.0,
    coupling::Real,
    dimension::Integer = 1,
    effective_frequency::Real = phonon_frequency,
    temperatures = [1.0],
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
)
    problem = holstein_poisson_problem(
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

Convert a generic `VariationalResult` into a lattice solution record.
"""
function solution_result(
    problem::VariationalProblem{HolsteinModel,PoissonTrial},
    temperature::Real,
    variational::VariationalResult,
)
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
    holstein_mobility(problem, solution)

Compute Einstein and transport-corrected DC mobilities for an optimized
Holstein-Poisson solution using the CTMC first-return plus exact-cloud
sideband transport kernel.
"""
function mobility_result(problem::VariationalProblem{HolsteinModel,PoissonTrial}, solution::LatticeSolution, options::OptimizerOptions)
    return _lattice_dc_mobility_result(
        problem.model,
        solution.rate,
        solution.beta,
        solution.temperature;
        rtol = options.quadrature_rtol,
    )
end

"""
    holstein_response(problem, solution, frequency)

Frequency-dependent lattice transport response for the optimized Holstein
solution using the CTMC first-return plus exact-cloud sideband transport
kernel. At zero temperature and finite frequency, the reported reduced optical
response is the transport kernel itself rather than `βκ` times that kernel.
"""
function response_result(
    problem::VariationalProblem{HolsteinModel,PoissonTrial},
    solution::LatticeSolution,
    frequency::Real,
    options::OptimizerOptions,
)
    return _lattice_response_result(
        problem.model,
        solution.rate,
        solution.beta,
        solution.temperature,
        frequency;
        rtol = options.quadrature_rtol,
    )
end

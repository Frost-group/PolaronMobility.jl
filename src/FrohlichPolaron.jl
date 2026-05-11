"""
    frohlich_feynman_problem(; coupling, phonon_frequency=1, dimension=3, ...)

Build a `VariationalProblem{FrohlichModel,GaussianFeynmanTrial}`. This is the
preferred Fröhlich entry point when using dimensionless coupling/frequency data
directly.
"""
function frohlich_feynman_problem(;
    coupling,
    phonon_frequency = 1.0,
    dimension::Integer = 3,
    band_mass::Real = 1.0,
    effective_frequency::Real = first(_as_vector(phonon_frequency)),
    material::Union{Nothing,FrohlichMaterial} = nothing,
    initial_v::Union{Nothing,Real} = nothing,
    initial_w::Union{Nothing,Real} = nothing,
)
    model = FrohlichModel(
        coupling,
        phonon_frequency;
        dimension = dimension,
        band_mass = band_mass,
        effective_frequency = effective_frequency,
        material = material,
    )
    guess_v, guess_w = variational_initial_guess(model.alpha)
    trial = GaussianFeynmanTrial(
        dimension = dimension,
        initial_v = initial_v === nothing ? guess_v : initial_v,
        initial_w = initial_w === nothing ? guess_w : initial_w,
    )
    return VariationalProblem(model, trial)
end

"""
    frohlich_multi_gaussian_problem(; coupling, modes, phonon_frequency=1, dimension=3, ...)

Build a Fröhlich variational problem using a finite-mode Gaussian trial. The
trial parameters are ordered as `w1, delta1, w2, delta2, ...`; `modes = 1`
reduces to the Feynman Gaussian trial formulas.
"""
function frohlich_multi_gaussian_problem(;
    coupling,
    modes::Integer,
    phonon_frequency = 1.0,
    dimension::Integer = 3,
    band_mass::Real = 1.0,
    effective_frequency::Real = first(_as_vector(phonon_frequency)),
    material::Union{Nothing,FrohlichMaterial} = nothing,
    initial_v = nothing,
    initial_w = nothing,
)
    model = FrohlichModel(
        coupling,
        phonon_frequency;
        dimension = dimension,
        band_mass = band_mass,
        effective_frequency = effective_frequency,
        material = material,
    )
    guess_v, guess_w = variational_initial_guess(model.alpha)
    if initial_w === nothing
        initial_w = [guess_w * (1 + 0.12 * (index - 1)) for index in 1:modes]
    end
    if initial_v === nothing
        initial_v = [
            index == 1 ? guess_v : initial_w[index] + max(1e-4 / index, (guess_v - guess_w) / (4index))
            for index in 1:modes
        ]
    end
    trial = MultiGaussianTrial(
        modes = modes,
        dimension = dimension,
        initial_v = initial_v,
        initial_w = initial_w,
    )
    return VariationalProblem(model, trial)
end

"""
    frohlich_nonlocal_gaussian_problem(; coupling, basis_frequencies, phonon_frequency=1, dimension=3, ...)

Build an experimental Fröhlich variational problem using a finite nonlocal
Gaussian kernel basis. The optimized parameters are non-negative amplitudes
`a1, a2, ...`; the entropy term is a configurable quadratic regularizer rather
than a literature-pinned closed-form action.
"""
function frohlich_nonlocal_gaussian_problem(;
    coupling,
    basis_frequencies,
    phonon_frequency = 1.0,
    dimension::Integer = 3,
    band_mass::Real = 1.0,
    effective_frequency::Real = first(_as_vector(phonon_frequency)),
    material::Union{Nothing,FrohlichMaterial} = nothing,
    initial_amplitudes = nothing,
    regularization::Real = 1e-3,
)
    model = FrohlichModel(
        coupling,
        phonon_frequency;
        dimension = dimension,
        band_mass = band_mass,
        effective_frequency = effective_frequency,
        material = material,
    )
    trial = NonlocalGaussianTrial(
        basis_frequencies = basis_frequencies,
        dimension = dimension,
        initial_amplitudes = initial_amplitudes,
        regularization = regularization,
    )
    return VariationalProblem(model, trial)
end

"""
    frohlich_profile_gaussian_problem(; coupling, basis_frequencies, phonon_frequency=1, dimension=3, ...)

Build a Fröhlich variational problem using a general Gaussian profile-function
trial. The profile is `Γ(ω) = sum(aᵢνᵢ² / (ω² + νᵢ²))`, with amplitudes
optimized by the generic solver.
"""
function frohlich_profile_gaussian_problem(;
    coupling,
    basis_frequencies,
    phonon_frequency = 1.0,
    dimension::Integer = 3,
    band_mass::Real = 1.0,
    effective_frequency::Real = first(_as_vector(phonon_frequency)),
    material::Union{Nothing,FrohlichMaterial} = nothing,
    initial_amplitudes = nothing,
    matsubara_terms::Integer = 4096,
)
    model = FrohlichModel(
        coupling,
        phonon_frequency;
        dimension = dimension,
        band_mass = band_mass,
        effective_frequency = effective_frequency,
        material = material,
    )
    trial = ProfileGaussianTrial(
        basis_frequencies = basis_frequencies,
        dimension = dimension,
        initial_amplitudes = initial_amplitudes,
        matsubara_terms = matsubara_terms,
    )
    return VariationalProblem(model, trial)
end

"""
    solve(problem::VariationalProblem{FrohlichModel,<:AbstractGaussianTrial}; temperatures=[0], frequencies=[0], options=OptimizerOptions())

Solve a Fröhlich Gaussian variational problem on reduced temperature and
frequency grids, returning a `PolaronResult`.
"""
function solve(
    problem::VariationalProblem{FrohlichModel,T};
    temperatures = [0.0],
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
) where {T<:AbstractGaussianTrial}
    solve_temperatures, solve_frequencies = _frohlich_solve_grid(problem.model, temperatures, frequencies)
    temperature_values, frequency_values, zero_solution, solutions, mobilities, responses =
        _solve_grid(problem; temperatures = solve_temperatures, frequencies = solve_frequencies, options = options)
    return PolaronResult(problem, temperature_values, frequency_values, zero_solution, solutions, mobilities, responses)
end

function _frohlich_solve_grid(model::FrohlichModel, temperatures, frequencies)
    model.material === nothing && return temperatures, frequencies
    reduced_T = [iszero(T) ? 0.0 : reduced_temperature(T, model.effective_frequency) for T in _as_vector(temperatures)]
    reduced_Ω = [reduced_frequency(Ω, model.effective_frequency) for Ω in _as_vector(frequencies)]
    return reduced_T, reduced_Ω
end

"""
    solve_frohlich(alpha, phonon_frequencies=1; temperatures=[0], frequencies=[0], dimension=3, options=OptimizerOptions())

Convenience wrapper around `frohlich_feynman_problem(...)` and `solve(...)`.
Use the builder-first form when composing with generic solver or sweep tools.
"""
function solve_frohlich(
    alpha,
    phonon_frequencies = 1.0;
    temperatures = [0.0],
    frequencies = [0.0],
    dimension::Integer = 3,
    band_mass::Real = 1.0,
    effective_frequency::Real = first(_as_vector(phonon_frequencies)),
    options::OptimizerOptions = OptimizerOptions(),
)
    problem = frohlich_feynman_problem(
        coupling = alpha,
        phonon_frequency = phonon_frequencies,
        dimension = dimension,
        band_mass = band_mass,
        effective_frequency = effective_frequency,
    )
    return solve(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

function mobility_result(
    problem::VariationalProblem{FrohlichModel,T},
    solution::FrohlichSolution,
    options::OptimizerOptions,
) where {T<:AbstractGaussianTrial}
    if solution.beta == Inf
        return FrohlichMobilityResult(solution.temperature, solution.beta, Inf, Inf, Inf, Inf, 0.0, Inf, Inf)
    end

    model = problem.model
    mobility = if problem.trial isa NonlocalGaussianTrial || problem.trial isa ProfileGaussianTrial
        _nonlocal_frohlich_mobility(problem, solution, options)
    else
        variational_v = problem.trial isa MultiGaussianTrial ? multi_gaussian_v(solution) : solution.v
        variational_w = problem.trial isa MultiGaussianTrial ? multi_gaussian_w(solution) : solution.w
        frohlich_mobility(
            variational_v,
            variational_w,
            model.alpha,
            model.phonon_frequencies,
            solution.beta;
            dimension = model.dimension,
            cutoff = options.memory_cutoff,
            rtol = options.quadrature_rtol,
        )
    end
    fhip = problem.trial isa NonlocalGaussianTrial || problem.trial isa ProfileGaussianTrial ? NaN :
        fhip_low_temperature_mobility(solution.v, solution.w, model.alpha, model.phonon_frequencies, solution.beta)
    kadanoff_devreese, kadanoff, relaxation_time = problem.trial isa NonlocalGaussianTrial || problem.trial isa ProfileGaussianTrial ? (NaN, NaN, NaN) :
        kadanoff_low_temperature_mobility(solution.v, solution.w, model.alpha, model.phonon_frequencies, solution.beta)
    hellwarth, hellwarth_b0 = problem.trial isa NonlocalGaussianTrial || problem.trial isa ProfileGaussianTrial ? (NaN, NaN) :
        hellwarth_mobility(solution.v, solution.w, model.alpha, model.phonon_frequencies, solution.beta)

    return FrohlichMobilityResult(
        solution.temperature,
        solution.beta,
        Float64(mobility),
        Float64(fhip),
        Float64(kadanoff_devreese),
        Float64(kadanoff),
        Float64(relaxation_time),
        Float64(hellwarth),
        Float64(hellwarth_b0),
    )
end

function response_result(
    problem::VariationalProblem{FrohlichModel,T},
    solution::FrohlichSolution,
    frequency::Real,
    options::OptimizerOptions,
) where {T<:AbstractGaussianTrial}
    if solution.beta == Inf && iszero(frequency)
        return FrohlichResponseResult(solution.temperature, solution.beta, Float64(frequency), Inf + 0im, 0 + Inf * im, 0 + 0im)
    end

    model = problem.model
    memory = if problem.trial isa NonlocalGaussianTrial || problem.trial isa ProfileGaussianTrial
        _nonlocal_frohlich_memory_function(problem, solution, frequency, options)
    else
        variational_v = problem.trial isa MultiGaussianTrial ? multi_gaussian_v(solution) : solution.v
        variational_w = problem.trial isa MultiGaussianTrial ? multi_gaussian_w(solution) : solution.w
        frohlich_memory_function(
            frequency,
            variational_v,
            variational_w,
            model.alpha,
            model.phonon_frequencies,
            solution.beta;
            dimension = model.dimension,
            cutoff = options.memory_cutoff,
            rtol = options.quadrature_rtol,
        )
    end
    impedance = -im * (frequency + memory)
    conductivity = inv(impedance)
    return FrohlichResponseResult(solution.temperature, solution.beta, Float64(frequency), memory, impedance, conductivity)
end

function _nonlocal_frohlich_structure_factor(problem, solution, t, alpha, omega, profile_decomposition = nothing)
    model = problem.model
    trial = problem.trial
    parameters = solution.variational.parameters
    coupling = frohlich_coupling(1, alpha, omega; dimension = model.dimension) * omega
    displacement = profile_decomposition === nothing ?
        mean_square_displacement(trial, parameters, im * t, solution.beta) :
        _profile_mean_square_displacement(profile_decomposition, im * t, solution.beta)
    propagator = displacement * omega / 2
    integral = ball_surface(model.dimension) / (2π)^model.dimension * sqrt(π) / 4 / propagator^(3 / 2)
    return 2 / model.dimension * coupling * integral * phonon_propagator(im * t, omega, solution.beta)
end

function _nonlocal_frohlich_memory_function(problem, solution, frequency, options)
    model = problem.model
    profile_decomposition = problem.trial isa ProfileGaussianTrial ?
        _profile_decomposition(problem.trial, solution.variational.parameters) :
        nothing
    memory = zero(ComplexF64)
    for (alpha, omega) in zip(model.alpha, model.phonon_frequencies)
        structure_factor(t) = _nonlocal_frohlich_structure_factor(problem, solution, t, alpha, omega, profile_decomposition)
        memory += ComplexF64(memory_integral(frequency, structure_factor; cutoff = options.memory_cutoff, rtol = options.quadrature_rtol))
    end
    return memory
end

function _nonlocal_frohlich_mobility(problem, solution, options)
    inverse_mobility = abs(imag(_nonlocal_frohlich_memory_function(problem, solution, 0.0, options)))
    return iszero(inverse_mobility) ? Inf : inv(inverse_mobility)
end

function solution_result(
    problem::VariationalProblem{FrohlichModel,T},
    temperature::Real,
    variational::VariationalResult,
) where {T<:AbstractGaussianTrial}
    return derived_solution(temperature, variational.beta, problem.model.dimension, variational)
end

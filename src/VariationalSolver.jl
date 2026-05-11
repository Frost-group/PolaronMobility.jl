"""
    parameter_names(trial)

Return the symbolic names of a trial process's variational parameters in solver
order.
"""
parameter_names(::AbstractTrialProcess) = Symbol[]

"""
    initial_parameters(trial)
    initial_parameters(problem, beta)

Return a solver-ready initial parameter vector for a trial process or
temperature-specific variational problem.
"""
initial_parameters(trial::AbstractTrialProcess) = error("initial_parameters is not implemented for $(typeof(trial)).")
initial_parameters(problem::VariationalProblem, beta::Real) = initial_parameters(problem.trial)

"""
    parameter_bounds(trial, beta)

Return `(lower, upper)` parameter bounds for a trial process at inverse
temperature `beta`.
"""
parameter_bounds(trial::AbstractTrialProcess, beta::Real) = error("parameter_bounds is not implemented for $(typeof(trial)).")

"""
    free_energy(trial, parameters, beta)

Reference/free contribution of a trial process to the variational objective.
"""
free_energy(trial::AbstractTrialProcess, parameters::AbstractVector, beta::Real) = error("free_energy is not implemented for $(typeof(trial)).")

"""
    entropy_cost(trial, parameters, beta)

Relative-entropy or trial-correction contribution to the variational objective.
"""
entropy_cost(trial::AbstractTrialProcess, parameters::AbstractVector, beta::Real) = error("entropy_cost is not implemented for $(typeof(trial)).")

"""
    interaction_free_energy(model, trial, parameters, beta; kwargs...)

Model/trial interaction contribution to the variational objective.
"""
interaction_free_energy(model::AbstractPolaronModel, trial::AbstractTrialProcess, parameters::AbstractVector, beta::Real; kwargs...) =
    error("interaction_free_energy is not implemented for $(typeof(model)) with $(typeof(trial)).")
diagnostics(model::AbstractPolaronModel, trial::AbstractTrialProcess, parameters::AbstractVector, beta::Real; kwargs...) = (;)

parameter_names(::GaussianFeynmanTrial) = [:w, :delta]
initial_parameters(trial::GaussianFeynmanTrial) = [trial.initial_w, trial.initial_delta]
parameter_bounds(trial::GaussianFeynmanTrial, beta::Real) =
    ([trial.lower_w, trial.lower_delta], [trial.upper_w, trial.upper_delta])

function parameter_names(trial::MultiGaussianTrial)
    names = Symbol[]
    for index in 1:trial.modes
        push!(names, Symbol("w", index))
        push!(names, Symbol("delta", index))
    end
    return names
end

function initial_parameters(trial::MultiGaussianTrial)
    values = Float64[]
    for index in 1:trial.modes
        push!(values, trial.initial_w[index])
        push!(values, trial.initial_delta[index])
    end
    return values
end

function parameter_bounds(trial::MultiGaussianTrial, beta::Real)
    lower = Float64[]
    upper = Float64[]
    for index in 1:trial.modes
        push!(lower, trial.lower_w[index])
        push!(lower, trial.lower_delta[index])
        push!(upper, trial.upper_w[index])
        push!(upper, trial.upper_delta[index])
    end
    return lower, upper
end

parameter_names(trial::NonlocalGaussianTrial) = [Symbol("a", index) for index in eachindex(trial.basis_frequencies)]
initial_parameters(trial::NonlocalGaussianTrial) = copy(trial.initial_amplitudes)
parameter_bounds(trial::NonlocalGaussianTrial, beta::Real) = (copy(trial.lower_amplitudes), copy(trial.upper_amplitudes))

parameter_names(trial::ProfileGaussianTrial) = [Symbol("a", index) for index in eachindex(trial.basis_frequencies)]
initial_parameters(trial::ProfileGaussianTrial) = copy(trial.initial_amplitudes)
parameter_bounds(trial::ProfileGaussianTrial, beta::Real) = (copy(trial.lower_amplitudes), copy(trial.upper_amplitudes))

function unpack(trial::GaussianFeynmanTrial, parameters::AbstractVector)
    length(parameters) == 2 || throw(ArgumentError("GaussianFeynmanTrial expects parameters [w, delta]."))
    w, delta = parameters
    w > 0 || throw(DomainError(w, "w must be positive."))
    delta >= 0 || throw(DomainError(delta, "delta must be non-negative."))
    return w + delta, w
end

function unpack(trial::MultiGaussianTrial, parameters::AbstractVector)
    expected = 2 * trial.modes
    length(parameters) == expected ||
        throw(ArgumentError("MultiGaussianTrial expects $expected parameters ordered as w1, delta1, w2, delta2, ...."))
    w = [parameters[2 * index - 1] for index in 1:trial.modes]
    delta = [parameters[2 * index] for index in 1:trial.modes]
    all(>(0), w) || throw(DomainError(w, "all w values must be positive."))
    all(>=(0), delta) || throw(DomainError(delta, "all delta values must be non-negative."))
    v = w .+ delta
    _validate_multi_gaussian_separation(trial, v, w)
    return v, w
end

function _validate_multi_gaussian_separation(trial::MultiGaussianTrial, v, w)
    trial.modes == 1 && return nothing
    for left in 1:(trial.modes - 1), right in (left + 1):trial.modes
        abs(w[left] - w[right]) > trial.min_separation ||
            throw(DomainError(w, "MultiGaussianTrial w values must be distinct."))
        abs(v[left] - v[right]) > trial.min_separation ||
            throw(DomainError(v, "MultiGaussianTrial v values must be distinct."))
    end
    return nothing
end

function nonlocal_amplitudes(trial::NonlocalGaussianTrial, parameters::AbstractVector)
    length(parameters) == length(trial.basis_frequencies) ||
        throw(ArgumentError("NonlocalGaussianTrial expects one amplitude per basis frequency."))
    all(>=(0), parameters) || throw(DomainError(parameters, "all nonlocal Gaussian amplitudes must be non-negative."))
    sum(parameters) < 1 || throw(DomainError(parameters, "nonlocal Gaussian amplitudes must sum to less than one for a positive diffusion tail."))
    return parameters
end

function profile_amplitudes(trial::ProfileGaussianTrial, parameters::AbstractVector)
    length(parameters) == length(trial.basis_frequencies) ||
        throw(ArgumentError("ProfileGaussianTrial expects one amplitude per basis frequency."))
    all(>=(0), parameters) || throw(DomainError(parameters, "all Gaussian profile amplitudes must be non-negative."))
    return parameters
end

function variational_initial_guess(alpha::AbstractVector{<:Real})
    alpha_eff = sum(alpha)
    w = 2 + tanh((6 - alpha_eff) / 3)
    v = alpha_eff < 7 ? 3 + alpha_eff / 4 : 4 * alpha_eff^2 / (9π) - 3 / 2 * (2log(2) + 0.5772) - 3 / 4
    return max(v, w + sqrt(eps(Float64))), w
end

function initial_parameters(problem::VariationalProblem{FrohlichModel,GaussianFeynmanTrial}, beta::Real)
    v, w = variational_initial_guess(problem.model.alpha)
    return [w, v - w]
end

function initial_parameters(problem::VariationalProblem{FrohlichModel,MultiGaussianTrial}, beta::Real)
    base_v, base_w = variational_initial_guess(problem.model.alpha)
    trial = problem.trial
    values = Float64[]
    for index in 1:trial.modes
        w = index == 1 ? base_w : max(trial.initial_w[index], base_w * (1 + 0.12 * (index - 1)))
        delta = index == 1 ? base_v - base_w : max(trial.initial_delta[index], 1e-4 / index)
        push!(values, w)
        push!(values, delta)
    end
    return values
end

function free_energy(trial::GaussianFeynmanTrial, parameters::AbstractVector, beta::Real)
    return zero(eltype(parameters))
end

function entropy_cost(trial::GaussianFeynmanTrial, parameters::AbstractVector, beta::Real)
    v, w = unpack(trial, parameters)
    A, C = beta == Inf ? trial_energy(v, w; dimension = trial.dimension) : trial_energy(v, w, beta; dimension = trial.dimension)
    return -(A + C)
end

free_energy(trial::MultiGaussianTrial, parameters::AbstractVector, beta::Real) = zero(eltype(parameters))

function entropy_cost(trial::MultiGaussianTrial, parameters::AbstractVector, beta::Real)
    v, w = unpack(trial, parameters)
    A, C = beta == Inf ? trial_energy(v, w; dimension = trial.dimension) : trial_energy(v, w, beta; dimension = trial.dimension)
    return -(A + C)
end

function interaction_free_energy(
    model::FrohlichModel,
    trial::GaussianFeynmanTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    v, w = unpack(trial, parameters)
    interaction = beta == Inf ?
        frohlich_interaction_energy(v, w, model.alpha, model.phonon_frequencies; dimension = model.dimension, rtol = rtol) :
        frohlich_interaction_energy(v, w, model.alpha, model.phonon_frequencies, beta; dimension = model.dimension, rtol = rtol)
    return -interaction
end

function interaction_free_energy(
    model::FrohlichModel,
    trial::MultiGaussianTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    v, w = unpack(trial, parameters)
    interaction = beta == Inf ?
        frohlich_interaction_energy(v, w, model.alpha, model.phonon_frequencies; dimension = model.dimension, rtol = rtol) :
        frohlich_interaction_energy(v, w, model.alpha, model.phonon_frequencies, beta; dimension = model.dimension, rtol = rtol)
    return -interaction
end

free_energy(trial::NonlocalGaussianTrial, parameters::AbstractVector, beta::Real) = zero(eltype(parameters))

function entropy_cost(trial::NonlocalGaussianTrial, parameters::AbstractVector, beta::Real)
    amplitudes = nonlocal_amplitudes(trial, parameters)
    return trial.regularization * sum(abs2, amplitudes)
end

function interaction_free_energy(
    model::FrohlichModel,
    trial::NonlocalGaussianTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    interaction = sum(_nonlocal_frohlich_interaction_energy(
        trial,
        parameters,
        alpha,
        omega,
        beta;
        dimension = model.dimension,
        rtol = rtol,
    ) for (alpha, omega) in zip(model.alpha, model.phonon_frequencies))
    return -interaction
end

free_energy(trial::ProfileGaussianTrial, parameters::AbstractVector, beta::Real) = zero(eltype(parameters))

function entropy_cost(trial::ProfileGaussianTrial, parameters::AbstractVector, beta::Real)
    profile_amplitudes(trial, parameters)
    if beta == Inf
        integrand(ω) = begin
            Γ = profile_function(trial, parameters, ω)
            log1p(Γ) - Γ / (1 + Γ)
        end
        return trial.dimension / (2π) * _profile_integral_0inf(integrand)
    end
    total = zero(eltype(parameters))
    for n in 1:trial.matsubara_terms
        ωn = 2π * n / beta
        Γ = profile_function(trial, parameters, ωn)
        total += log1p(Γ) - Γ / (1 + Γ)
    end
    return trial.dimension * total / beta
end

function interaction_free_energy(
    model::FrohlichModel,
    trial::ProfileGaussianTrial,
    parameters::AbstractVector,
    beta::Real;
    rtol::Real = 1e-4,
)
    interaction = sum(_profile_frohlich_interaction_energy(
        trial,
        parameters,
        alpha,
        omega,
        beta;
        dimension = model.dimension,
        rtol = rtol,
    ) for (alpha, omega) in zip(model.alpha, model.phonon_frequencies))
    return -interaction
end

function _profile_frohlich_interaction_energy(
    trial::ProfileGaussianTrial,
    parameters,
    alpha,
    omega,
    beta;
    dimension::Integer,
    rtol::Real,
)
    coupling = frohlich_coupling(1, alpha, omega; dimension = dimension)
    upper = beta == Inf ? Inf : beta / 2
    decomposition = _profile_decomposition(trial, parameters)
    integrand(τ) = phonon_propagator(τ, omega, beta) / sqrt(_profile_mean_square_displacement(decomposition, τ, beta) * omega)
    integral = quadgk(integrand, 0, upper; rtol = rtol)[1]
    return coupling * ball_surface(dimension) / (2π)^dimension * sqrt(π / 2) * integral
end

function _nonlocal_frohlich_interaction_energy(
    trial::NonlocalGaussianTrial,
    parameters,
    alpha,
    omega,
    beta;
    dimension::Integer,
    rtol::Real,
)
    coupling = frohlich_coupling(1, alpha, omega; dimension = dimension)
    upper = beta == Inf ? Inf : beta / 2
    integrand(τ) = phonon_propagator(τ, omega, beta) / sqrt(mean_square_displacement(trial, parameters, τ, beta) * omega)
    integral = quadgk(integrand, 0, upper; rtol = rtol)[1]
    return coupling * ball_surface(dimension) / (2π)^dimension * sqrt(π / 2) * integral
end

"""
    objective(problem, parameters, beta; options=OptimizerOptions())

Evaluate the full variational objective for a `VariationalProblem` as
`free_energy + entropy_cost + interaction_free_energy`.
"""
function objective(problem::VariationalProblem, parameters::AbstractVector, beta::Real; options::OptimizerOptions = OptimizerOptions())
    return free_energy(problem.trial, parameters, beta) +
           entropy_cost(problem.trial, parameters, beta) +
           interaction_free_energy(problem.model, problem.trial, parameters, beta; rtol = options.quadrature_rtol)
end

function _objective_or_penalty(problem::VariationalProblem, parameters::AbstractVector, beta::Real; options::OptimizerOptions)
    try
        value = objective(problem, parameters, beta; options = options)
        return isfinite(value) ? value : 1e100
    catch error
        error isa DomainError || error isa ArgumentError || rethrow()
        return 1e100
    end
end

function diagnostics(
    model::FrohlichModel,
    trial::GaussianFeynmanTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    v, w = unpack(trial, parameters)
    energy = frohlich_energy(v, w, model.alpha, model.phonon_frequencies, beta; dimension = model.dimension, rtol = options.quadrature_rtol)
    spring_constant = v^2 - w^2
    reduced_mass = spring_constant / v^2
    return (;
        v = Float64(v),
        w = Float64(w),
        delta = Float64(v - w),
        total_energy = Float64(energy.total),
        trial_free = Float64(energy.trial_free),
        interaction = Float64(energy.interaction),
        trial_correction = Float64(energy.trial_correction),
        spring_constant = Float64(spring_constant),
        fictitious_mass = Float64(spring_constant / w^2),
        asymptotic_mass = Float64(v / w),
        reduced_mass = Float64(reduced_mass),
        radius = Float64(sqrt(3 / (2 * reduced_mass * v))),
    )
end

function diagnostics(
    model::FrohlichModel,
    trial::MultiGaussianTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    v, w = unpack(trial, parameters)
    energy = frohlich_energy(v, w, model.alpha, model.phonon_frequencies, beta; dimension = model.dimension, rtol = options.quadrature_rtol)
    primary_v = first(v)
    primary_w = first(w)
    spring_constant = sum(v .^ 2 .- w .^ 2)
    reduced_mass = spring_constant / sum(v .^ 2)
    mode_pairs = Pair{Symbol,Float64}[]
    for index in eachindex(v)
        push!(mode_pairs, Symbol("v", index) => Float64(v[index]))
        push!(mode_pairs, Symbol("w", index) => Float64(w[index]))
        push!(mode_pairs, Symbol("delta", index) => Float64(v[index] - w[index]))
    end
    return merge(
        (;
            modes = Float64(trial.modes),
            v = Float64(primary_v),
            w = Float64(primary_w),
            delta = Float64(primary_v - primary_w),
            total_energy = Float64(energy.total),
            trial_free = Float64(energy.trial_free),
            interaction = Float64(energy.interaction),
            trial_correction = Float64(energy.trial_correction),
            spring_constant = Float64(spring_constant),
            fictitious_mass = Float64(spring_constant / primary_w^2),
            asymptotic_mass = Float64(primary_v / primary_w),
            reduced_mass = Float64(reduced_mass),
            radius = Float64(sqrt(3 / (2 * max(reduced_mass, eps(Float64)) * primary_v))),
        ),
        NamedTuple(mode_pairs),
    )
end

function diagnostics(
    model::FrohlichModel,
    trial::NonlocalGaussianTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    amplitudes = Float64.(nonlocal_amplitudes(trial, parameters))
    amplitude_pairs = Pair{Symbol,Float64}[]
    for (index, amplitude) in enumerate(amplitudes)
        push!(amplitude_pairs, Symbol("amplitude", index) => Float64(amplitude))
    end
    return merge(
        (;
            kernel_basis_modes = Float64(length(trial.basis_frequencies)),
            amplitude_sum = Float64(sum(amplitudes)),
            regularization = Float64(trial.regularization),
            total_energy = Float64(objective(VariationalProblem(model, trial), parameters, beta; options = options)),
            trial_free = 0.0,
            interaction = Float64(-interaction_free_energy(model, trial, parameters, beta; rtol = options.quadrature_rtol)),
            trial_correction = Float64(entropy_cost(trial, parameters, beta)),
            v = Float64(first(trial.basis_frequencies)),
            w = Float64(first(trial.basis_frequencies)),
            delta = 0.0,
            spring_constant = Float64(sum(amplitudes)),
            fictitious_mass = Float64(sum(amplitudes)),
            asymptotic_mass = 1.0,
            reduced_mass = Float64(max(sum(amplitudes), eps(Float64))),
            radius = Float64(sqrt(3 / (2 * max(sum(amplitudes), eps(Float64)) * first(trial.basis_frequencies)))),
        ),
        NamedTuple(amplitude_pairs),
    )
end

function diagnostics(
    model::FrohlichModel,
    trial::ProfileGaussianTrial,
    parameters::AbstractVector,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
)
    amplitudes = Float64.(profile_amplitudes(trial, parameters))
    amplitude_pairs = Pair{Symbol,Float64}[]
    for (index, amplitude) in enumerate(amplitudes)
        push!(amplitude_pairs, Symbol("amplitude", index) => Float64(amplitude))
    end
    entropy = Float64(entropy_cost(trial, parameters, beta))
    interaction = Float64(-interaction_free_energy(model, trial, parameters, beta; rtol = options.quadrature_rtol))
    total = entropy - interaction
    return merge(
        (;
            kernel_basis_modes = Float64(length(trial.basis_frequencies)),
            amplitude_sum = Float64(sum(amplitudes)),
            total_energy = Float64(total),
            trial_free = 0.0,
            interaction = interaction,
            trial_correction = entropy,
            v = Float64(first(trial.basis_frequencies)),
            w = Float64(first(trial.basis_frequencies)),
            delta = 0.0,
            spring_constant = Float64(sum(amplitudes)),
            fictitious_mass = Float64(sum(amplitudes)),
            asymptotic_mass = 1.0,
            reduced_mass = Float64(max(sum(amplitudes), eps(Float64))),
            radius = Float64(sqrt(3 / (2 * max(sum(amplitudes), eps(Float64)) * first(trial.basis_frequencies)))),
        ),
        NamedTuple(amplitude_pairs),
    )
end

"""
    solve_variational(problem, beta; options=OptimizerOptions(), initial_parameters_override=nothing, use_multistart=options.multistart)

Optimize one `VariationalProblem` at inverse temperature `beta`, returning a
`VariationalResult` with named parameters, objective components, and optimizer
diagnostics.
"""
function solve_variational(
    problem::VariationalProblem,
    beta::Real;
    options::OptimizerOptions = OptimizerOptions(),
    initial_parameters_override = nothing,
    use_multistart::Bool = options.multistart,
)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    lower, upper = _solver_bounds(problem.trial, beta, options)
    initial = initial_parameters_override === nothing ?
        _default_initial_parameters(problem, beta, options, lower, upper) :
        Float64.(collect(initial_parameters_override))
    length(initial) == length(lower) || throw(ArgumentError("initial parameter length does not match bounds."))

    result, final_lower, final_upper, expansions = _solve_with_adaptive_bounds(
        problem,
        beta,
        initial,
        lower,
        upper,
        options;
        use_multistart = use_multistart,
    )

    parameters = Float64.(Optim.minimizer(result))
    reference = free_energy(problem.trial, parameters, beta)
    entropy = entropy_cost(problem.trial, parameters, beta)
    interaction = interaction_free_energy(problem.model, problem.trial, parameters, beta; rtol = options.quadrature_rtol)
    model_diagnostics = diagnostics(problem.model, problem.trial, parameters, beta; options = options)
    optimizer_diagnostics = _optimizer_diagnostics(result, parameters, final_lower, final_upper, expansions, options)
    return VariationalResult(
        parameters,
        parameter_names(problem.trial),
        Float64(beta),
        Float64(reference + entropy + interaction),
        Float64(reference),
        Float64(entropy),
        Float64(interaction),
        merge(model_diagnostics, optimizer_diagnostics),
    )
end

temperature_beta(temperature::Real) = iszero(temperature) ? Inf : inv(Float64(temperature))

function _solve_grid(
    problem::VariationalProblem;
    temperatures,
    frequencies,
    options::OptimizerOptions,
)
    temperature_values = _as_vector(temperatures)
    frequency_values = _as_vector(frequencies)
    !isempty(temperature_values) || throw(ArgumentError("temperatures must not be empty."))
    !isempty(frequency_values) || throw(ArgumentError("frequencies must not be empty."))
    all(>=(0), temperature_values) ||
        throw(DomainError(temperature_values, "temperatures must be non-negative reduced temperatures."))

    zero_variational = solve_variational(problem, Inf; options = options)
    zero_solution = solution_result(problem, 0.0, zero_variational)

    solution_type = typeof(zero_solution)
    mobility_type = typeof(mobility_result(problem, zero_solution, options))
    response_type = typeof(response_result(problem, zero_solution, first(frequency_values), options))

    solutions = Vector{solution_type}(undef, length(temperature_values))
    mobilities = Vector{mobility_type}(undef, length(temperature_values))
    responses = Matrix{response_type}(undef, length(frequency_values), length(temperature_values))

    previous_parameters = zero_variational.parameters
    for (temperature_index, temperature) in pairs(temperature_values)
        beta = temperature_beta(temperature)
        solution = if beta == Inf
            zero_solution
        else
            variational = solve_variational(
                problem,
                beta;
                options = options,
                initial_parameters_override = options.warm_start ? previous_parameters : nothing,
                use_multistart = previous_parameters === nothing && options.multistart,
            )
            previous_parameters = variational.parameters
            solution_result(problem, temperature, variational)
        end
        solutions[temperature_index] = solution
        mobilities[temperature_index] = mobility_result(problem, solution, options)
        for (frequency_index, frequency) in pairs(frequency_values)
            responses[frequency_index, temperature_index] = response_result(problem, solution, frequency, options)
        end
    end

    return temperature_values, frequency_values, zero_solution, solutions, mobilities, responses
end

function _solver_bounds(trial::AbstractTrialProcess, beta::Real, options::OptimizerOptions)
    default_lower, default_upper = parameter_bounds(trial, beta)
    lower = isempty(options.lower) ? default_lower : copy(options.lower)
    upper = isempty(options.upper) ? default_upper : copy(options.upper)
    length(lower) == length(upper) || throw(ArgumentError("lower and upper bounds must have the same length."))
    return Float64.(lower), Float64.(upper)
end

function _default_initial_parameters(problem::VariationalProblem, beta::Real, options::OptimizerOptions, lower, upper)
    guess = options.initial_parameters === nothing ? initial_parameters(problem, beta) : options.initial_parameters
    return _clip_to_bounds(Float64.(guess), lower, upper)
end

function _solve_with_adaptive_bounds(problem, beta, initial, lower, upper, options; use_multistart::Bool)
    current_lower = copy(lower)
    current_upper = copy(upper)
    current_initial = _clip_to_bounds(initial, current_lower, current_upper)
    expansions = 0
    while true
        result = _run_optimizer(problem, beta, current_initial, current_lower, current_upper, options; use_multistart = use_multistart && expansions == 0)
        parameters = Float64.(Optim.minimizer(result))
        expanded_upper = _expanded_upper_bounds(parameters, current_lower, current_upper, options)
        if !options.adaptive_bounds || expanded_upper === nothing || expansions >= options.max_bound_expansions
            return result, current_lower, current_upper, expansions
        end
        current_initial = _clip_to_bounds(parameters, current_lower, expanded_upper)
        current_upper = expanded_upper
        expansions += 1
    end
end

function _run_optimizer(problem, beta, initial, lower, upper, options; use_multistart::Bool)
    seeds = use_multistart ? _seed_candidates(problem, beta, initial, lower, upper, options) : [initial]
    results = map(seeds) do seed
        objective_function(x) = _objective_or_penalty(problem, x, beta; options = options)
        _optimize_seed(problem, objective_function, seed, lower, upper, options)
    end
    return results[argmin([Optim.minimum(result) for result in results])]
end

function _optimize_seed(problem::VariationalProblem, objective_function, seed, lower, upper, options)
    return Optim.optimize(
        Optim.OnceDifferentiable(objective_function, seed; autodiff = :forward),
        lower,
        upper,
        seed,
        Fminbox(BFGS()),
        Optim.Options(g_tol = options.gradient_tolerance),
    )
end

function _optimize_seed(problem::VariationalProblem{FrohlichModel,ProfileGaussianTrial}, objective_function, seed, lower, upper, options)
    return Optim.optimize(
        Optim.OnceDifferentiable(objective_function, seed; autodiff = :finite),
        lower,
        upper,
        seed,
        Fminbox(BFGS(linesearch = Optim.LineSearches.BackTracking())),
        Optim.Options(g_tol = options.gradient_tolerance),
    )
end

function _seed_candidates(problem, beta, initial, lower, upper, options)
    finite_bounds = all(isfinite, lower) && all(isfinite, upper)
    finite_bounds || return [initial]
    candidates = [initial]
    grids = map(eachindex(lower)) do index
        span = upper[index] - lower[index]
        inner_lower = lower[index] + 0.02 * span
        inner_upper = upper[index] - 0.02 * span
        range(inner_lower, inner_upper; length = 5)
    end
    mesh = Iterators.product(grids...)
    scored = vec(map(mesh) do tuple_seed
        seed = Float64[tuple_seed...]
        _objective_or_penalty(problem, seed, beta; options = options), seed
    end)
    sort!(scored; by = first)
    for (_, seed) in first(scored, min(8, length(scored)))
        push!(candidates, _clip_to_bounds(seed, lower, upper))
    end
    unique_candidates = RealVector[]
    seen = Set{Tuple{Vararg{Float64}}}()
    for candidate in candidates
        key = Tuple(round.(candidate; digits = 12))
        key in seen && continue
        push!(seen, key)
        push!(unique_candidates, candidate)
    end
    return unique_candidates
end

function _expanded_upper_bounds(parameters, lower, upper, options)
    expanded = false
    next_upper = copy(upper)
    for index in eachindex(parameters)
        isfinite(upper[index]) || continue
        threshold = _bound_hit_threshold(lower[index], upper[index], options)
        if upper[index] - parameters[index] <= threshold
            span = upper[index] - lower[index]
            next_upper[index] = upper[index] + max(2span, abs(upper[index]), 1.0)
            expanded = true
        end
    end
    return expanded ? next_upper : nothing
end

function _bound_hit_threshold(lower, upper, options)
    span = upper - lower
    return max(1e-8 * max(abs(upper), 1.0), options.bound_hit_fraction * span)
end

function _clip_to_bounds(values, lower, upper)
    return map(eachindex(values)) do index
        lo = lower[index]
        hi = upper[index]
        value = min(max(values[index], lo), hi)
        if isfinite(lo) && isfinite(hi)
            margin = max(eps(Float64) * max(abs(hi), 1.0), 1e-10 * (hi - lo))
            value = min(max(value, lo + margin), hi - margin)
        end
        value
    end
end

function _optimizer_diagnostics(result, parameters, lower, upper, expansions, options)
    base = (;
        optimizer_success = Optim.converged(result),
        optimizer_objective = Float64(Optim.minimum(result)),
        optimizer_bound_expansions = Float64(expansions),
    )
    bound_pairs = Pair{Symbol,Float64}[]
    hit_pairs = Pair{Symbol,Bool}[]
    for index in eachindex(parameters)
        threshold = _bound_hit_threshold(lower[index], upper[index], options)
        push!(bound_pairs, Symbol("optimizer_bound_lower_", index) => Float64(lower[index]))
        push!(bound_pairs, Symbol("optimizer_bound_upper_", index) => Float64(upper[index]))
        push!(bound_pairs, Symbol("optimizer_bound_span_", index) => Float64(upper[index] - lower[index]))
        push!(hit_pairs, Symbol("optimizer_hit_upper_bound_", index) => upper[index] - parameters[index] <= threshold)
    end
    return merge(base, NamedTuple(bound_pairs), NamedTuple(hit_pairs))
end

"""
    feynman_v(result)

Return Feynman's `v` parameter from a Fröhlich `VariationalResult` or
`FrohlichSolution`.
"""
function feynman_v(result::VariationalResult)
    index = findfirst(==(:delta), result.parameter_names)
    w_index = findfirst(==(:w), result.parameter_names)
    (index === nothing || w_index === nothing) && throw(ArgumentError("result does not contain Feynman parameters."))
    return result.parameters[w_index] + result.parameters[index]
end

"""
    feynman_w(result)

Return Feynman's `w` parameter from a Fröhlich `VariationalResult` or
`FrohlichSolution`.
"""
function feynman_w(result::VariationalResult)
    index = findfirst(==(:w), result.parameter_names)
    index === nothing && throw(ArgumentError("result does not contain Feynman w parameter."))
    return result.parameters[index]
end

feynman_v(solution::FrohlichSolution) = feynman_v(solution.variational)
feynman_w(solution::FrohlichSolution) = feynman_w(solution.variational)

"""
    multi_gaussian_v(result)

Return all `v_i = w_i + delta_i` values from a multi-mode Gaussian
`VariationalResult` or `FrohlichSolution`.
"""
function multi_gaussian_v(result::VariationalResult)
    values = Float64[]
    index = 1
    while true
        w_index = findfirst(==(Symbol("w", index)), result.parameter_names)
        delta_index = findfirst(==(Symbol("delta", index)), result.parameter_names)
        (w_index === nothing || delta_index === nothing) && break
        push!(values, result.parameters[w_index] + result.parameters[delta_index])
        index += 1
    end
    isempty(values) && return [feynman_v(result)]
    return values
end

"""
    multi_gaussian_w(result)

Return all `w_i` values from a multi-mode Gaussian `VariationalResult` or
`FrohlichSolution`.
"""
function multi_gaussian_w(result::VariationalResult)
    values = Float64[]
    index = 1
    while true
        w_index = findfirst(==(Symbol("w", index)), result.parameter_names)
        w_index === nothing && break
        push!(values, result.parameters[w_index])
        index += 1
    end
    isempty(values) && return [feynman_w(result)]
    return values
end

multi_gaussian_v(solution::FrohlichSolution) = multi_gaussian_v(solution.variational)
multi_gaussian_w(solution::FrohlichSolution) = multi_gaussian_w(solution.variational)

function energy_components(result::VariationalResult)
    diagnostics = result.diagnostics
    return EnergyComponents(
        diagnostics.total_energy,
        diagnostics.trial_free,
        diagnostics.interaction,
        diagnostics.trial_correction,
    )
end

function derived_solution(temperature::Real, beta::Real, dimension::Integer, result::VariationalResult)
    diagnostics = result.diagnostics
    parameters = (; zip(result.parameter_names, result.parameters)...)
    return FrohlichSolution(
        Float64(temperature),
        Float64(beta),
        Int(dimension),
        result,
        parameters,
        Float64(diagnostics.v),
        Float64(diagnostics.w),
        energy_components(result),
        Float64(diagnostics.spring_constant),
        Float64(diagnostics.fictitious_mass),
        Float64(diagnostics.asymptotic_mass),
        Float64(diagnostics.reduced_mass),
        Float64(diagnostics.radius),
    )
end

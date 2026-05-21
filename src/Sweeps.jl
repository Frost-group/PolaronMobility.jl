function _variational_result_row(
    problem::VariationalProblem{FrohlichModel,T},
    result::VariationalResult;
    temperature::Real,
    sweep_direction::Symbol = :none,
) where {T<:AbstractGaussianTrial}
    model = problem.model
    representative_frequency = length(model.phonon_frequencies) == 1 ? first(model.phonon_frequencies) : model.effective_frequency
    return merge(
        (;
            dimension = model.dimension,
            beta = result.beta,
            temperature = Float64(temperature),
            coupling = sum(model.alpha),
            phonon_frequency = representative_frequency,
            phonon_modes = length(model.phonon_frequencies),
            adiabatic_ratio = inv(representative_frequency),
            free_energy = result.free_energy,
            reference_free_energy = result.reference_free_energy,
            entropy_cost = result.entropy_cost,
            interaction_free_energy = result.interaction_energy,
            sweep_direction = sweep_direction,
        ),
        _parameter_row(result),
        result.diagnostics,
    )
end

function _variational_result_row(
    problem::VariationalProblem{HolsteinModel,PoissonTrial},
    result::VariationalResult;
    temperature::Real,
    sweep_direction::Symbol = :none,
)
    model = problem.model
    return merge(
        (;
            dimension = model.dimension,
            beta = result.beta,
            temperature = Float64(temperature),
            coupling = model.coupling,
            hopping = model.hopping,
            phonon_frequency = model.phonon_frequency,
            adiabatic_ratio = model.hopping / model.phonon_frequency,
            free_energy = result.free_energy,
            reference_free_energy = result.reference_free_energy,
            entropy_cost = result.entropy_cost,
            interaction_free_energy = result.interaction_energy,
            sweep_direction = sweep_direction,
        ),
        _parameter_row(result),
        result.diagnostics,
    )
end

function _variational_result_row(
    problem::VariationalProblem{M,PoissonTrial},
    result::VariationalResult;
    temperature::Real,
    sweep_direction::Symbol = :none,
) where {M<:Union{PeierlsModel,CompositePolaronModel}}
    model = problem.model
    return merge(
        (;
            dimension = model.dimension,
            beta = result.beta,
            temperature = Float64(temperature),
            hopping = model.hopping,
            adiabatic_ratio = _lattice_adiabatic_ratio(model),
            free_energy = result.free_energy,
            reference_free_energy = result.reference_free_energy,
            entropy_cost = result.entropy_cost,
            interaction_free_energy = result.interaction_energy,
            sweep_direction = sweep_direction,
        ),
        _model_row(problem),
        _parameter_row(result),
        result.diagnostics,
    )
end

function _parameter_row(result::VariationalResult)
    pairs = Pair{Symbol,Float64}[]
    for (name, value) in zip(result.parameter_names, result.parameters)
        push!(pairs, name => Float64(value))
    end
    return NamedTuple(pairs)
end

"""
    solution_table(result)

Return flat `NamedTuple` rows for optimized variational solutions in a
`PolaronResult`.
"""
solution_table(result::PolaronResult) = [_solution_row(result.problem, solution) for solution in result.solutions]

"""
    mobility_table(result)

Return flat `NamedTuple` rows for DC mobility observables in a `PolaronResult`.
"""
mobility_table(result::PolaronResult) = [
    merge(_solution_identity_row(result.solutions[index]), _mobility_row(result.problem, result.mobilities[index]))
    for index in eachindex(result.mobilities)
]

"""
    response_table(result)

Return flat `NamedTuple` rows for frequency-dependent response observables in a
`PolaronResult`.
"""
function response_table(result::PolaronResult)
    rows = NamedTuple[]
    for temperature_index in eachindex(result.temperatures), frequency_index in eachindex(result.frequencies)
        push!(rows, merge(
            _solution_row(result.problem, result.solutions[temperature_index]),
            _response_row(result.problem, result.responses[frequency_index, temperature_index]),
        ))
    end
    return rows
end

"""
    sweep_table(result_or_rows)

Normalize a result or existing row vector to flat sweep rows. For solve results
this returns `response_table(result)`, which includes solution parameters at
each temperature/frequency point.
"""
sweep_table(rows::AbstractVector{<:NamedTuple}) = rows
sweep_table(result::PolaronResult) = response_table(result)

"""
    frequency_sweep(problem; temperatures, frequencies, options=OptimizerOptions())

Solve a variational problem once per temperature and evaluate all requested
frequencies without reoptimizing per frequency. Returns flat response rows.
"""
function frequency_sweep(
    problem::VariationalProblem;
    temperatures,
    frequencies,
    options::OptimizerOptions = OptimizerOptions(),
)
    return response_table(solve(problem; temperatures = temperatures, frequencies = frequencies, options = options))
end

"""
    frohlich_frequency_sweep(frequencies; coupling, phonon_frequency=1, temperatures=[0], dimension=3, options=OptimizerOptions())

Convenience wrapper for `frequency_sweep(frohlich_feynman_problem(...))`.
"""
function frohlich_frequency_sweep(
    frequencies;
    coupling,
    phonon_frequency = 1.0,
    temperatures = [0.0],
    dimension::Integer = 3,
    options::OptimizerOptions = OptimizerOptions(),
)
    problem = frohlich_feynman_problem(coupling = coupling, phonon_frequency = phonon_frequency, dimension = dimension)
    return frequency_sweep(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

"""
    holstein_frequency_sweep(frequencies; coupling, hopping=1, phonon_frequency=1, temperatures=[1], dimension=1, options=OptimizerOptions())

Convenience wrapper for `frequency_sweep(holstein_poisson_problem(...))`.
"""
function holstein_frequency_sweep(
    frequencies;
    coupling,
    hopping = 1.0,
    phonon_frequency = 1.0,
    temperatures = [1.0],
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    problem = holstein_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension)
    return frequency_sweep(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

"""
    peierls_frequency_sweep(frequencies; coupling, hopping=1, phonon_frequency=1, temperatures=[1], dimension=1, options=OptimizerOptions())

Convenience wrapper for `frequency_sweep(peierls_poisson_problem(...))`.
"""
function peierls_frequency_sweep(
    frequencies;
    coupling,
    hopping = 1.0,
    phonon_frequency = 1.0,
    temperatures = [1.0],
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    problem = peierls_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension)
    return frequency_sweep(problem; temperatures = temperatures, frequencies = frequencies, options = options)
end

"""
    lattice_transport_sweep(problem; temperatures, frequencies=[0], options=OptimizerOptions(), kappa_source=:zero_temperature, broadening=1e-6, sideband_tolerance=1e-12, laguerre_points=120)

Evaluate guide-style lattice CTMC transport rows for a Holstein, Peierls, or
compatible composite Poisson problem. When `kappa_source = :zero_temperature`,
the transport kernel reuses a single `beta = Inf` optimized rate across the
entire temperature/frequency grid and caches repeated first-return evaluations.
"""
function lattice_transport_sweep(
    problem::VariationalProblem{M,PoissonTrial};
    temperatures,
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
    kappa_source::Symbol = :zero_temperature,
    broadening::Real = default_lattice_broadening,
    sideband_tolerance::Real = default_lattice_sideband_tolerance,
    laguerre_points::Integer = default_lattice_laguerre_points,
) where {M<:Union{HolsteinModel,PeierlsModel,CompositePolaronModel}}
    kappa_source in (:zero_temperature, :per_temperature) ||
        throw(ArgumentError("kappa_source must be :zero_temperature or :per_temperature."))
    solve_temperatures, solve_frequencies = _transport_solve_grid(problem.model, temperatures, frequencies)
    temperature_values = _as_vector(solve_temperatures)
    frequency_values = _as_vector(solve_frequencies)
    rows = NamedTuple[]

    zero_variational = solve_variational(problem, Inf; options = options, use_multistart = options.multistart)
    shared_rate = zero_variational.diagnostics.rate
    previous_parameters = kappa_source === :zero_temperature ? zero_variational.parameters : nothing

    for temperature in temperature_values
        beta = iszero(temperature) ? Inf : inv(temperature)
        if kappa_source === :per_temperature
            variational = solve_variational(
                problem,
                beta;
                options = options,
                initial_parameters_override = previous_parameters,
                use_multistart = previous_parameters === nothing && options.multistart,
            )
            previous_parameters = variational.parameters
            rate = variational.diagnostics.rate
            row_prefix = merge(
                _variational_result_row(problem, variational; temperature = temperature),
                (; kappa_source = kappa_source),
            )
        else
            rate = shared_rate
            row_prefix = merge(_model_row(problem), (; temperature = Float64(temperature), beta = Float64(beta), rate = Float64(rate), kappa_source = kappa_source))
        end

        sidebands = _transport_sidebands(problem.model, beta; response_hopping = default_lattice_response_hopping, tolerance = sideband_tolerance)
        cache = _transport_kernel_cache(sidebands, rate, problem.model.dimension, frequency_values; broadening = broadening, laguerre_points = laguerre_points)
        mobility = _lattice_dc_mobility_result(
            problem.model,
            rate,
            beta,
            temperature;
            broadening = broadening,
            sideband_tolerance = sideband_tolerance,
            laguerre_points = laguerre_points,
            kappa_source = kappa_source,
        )
        for frequency in frequency_values
            response_factor = _cached_transport_sum(cache, sidebands, frequency)
            mobility_einstein = beta == Inf ? Inf : beta * rate
            response = LatticeResponseResult(
                Float64(temperature),
                Float64(beta),
                Float64(frequency),
                _scale_lattice_response(response_factor, mobility_einstein, frequency),
                response_factor,
                lattice_conductivity(_scale_lattice_response(response_factor, mobility_einstein, frequency)),
                lattice_impedance(lattice_conductivity(_scale_lattice_response(response_factor, mobility_einstein, frequency))),
                Float64(_resolve_response_hopping(problem.model, rate, default_lattice_response_hopping)),
                _component_mobilities(
                    problem.model,
                    rate,
                    beta,
                    frequency;
                    broadening = broadening,
                    sideband_tolerance = sideband_tolerance,
                    laguerre_points = laguerre_points,
                ),
                merge(
                    _transport_diagnostics(
                        problem.model,
                        rate,
                        beta,
                        frequency,
                        response_factor;
                        response_hopping = default_lattice_response_hopping,
                        broadening = broadening,
                        sideband_tolerance = sideband_tolerance,
                        laguerre_points = laguerre_points,
                        kappa_source = kappa_source,
                    ),
                    (; cached_frequency_shifts = length(cache)),
                ),
            )
            push!(rows, merge(row_prefix, _mobility_row(problem, mobility), _response_row(problem, response)))
        end
    end
    return rows
end

"""
    holstein_transport_sweep(; coupling, hopping=1, phonon_frequency=1, dimension=1, temperatures, frequencies=[0], ...)

Convenience wrapper for guide-style Holstein lattice transport sweeps.
"""
function holstein_transport_sweep(;
    coupling,
    hopping = 1.0,
    phonon_frequency = 1.0,
    dimension::Integer = 1,
    temperatures,
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
    kwargs...,
)
    problem = holstein_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension)
    return lattice_transport_sweep(problem; temperatures = temperatures, frequencies = frequencies, options = options, kwargs...)
end

"""
    peierls_transport_sweep(; coupling, hopping=1, phonon_frequency=1, dimension=1, temperatures, frequencies=[0], ...)

Convenience wrapper for guide-style Peierls lattice transport sweeps.
"""
function peierls_transport_sweep(;
    coupling,
    hopping = 1.0,
    phonon_frequency = 1.0,
    dimension::Integer = 1,
    temperatures,
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
    kwargs...,
)
    problem = peierls_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension)
    return lattice_transport_sweep(problem; temperatures = temperatures, frequencies = frequencies, options = options, kwargs...)
end

"""
    holstein_peierls_transport_sweep(; holstein_coupling, peierls_coupling, hopping=1, holstein_frequency=1, peierls_frequency=1, dimension=1, temperatures, frequencies=[0], ...)

Guide-style transport sweep for a composite Holstein-Peierls model with a
shared Poisson trial.
"""
function holstein_peierls_transport_sweep(;
    holstein_coupling,
    peierls_coupling,
    hopping = 1.0,
    holstein_frequency = 1.0,
    peierls_frequency = 1.0,
    dimension::Integer = 1,
    temperatures,
    frequencies = [0.0],
    options::OptimizerOptions = OptimizerOptions(),
    kwargs...,
)
    holstein = holstein_poisson_problem(
        hopping = hopping,
        phonon_frequency = holstein_frequency,
        coupling = holstein_coupling,
        dimension = dimension,
    ).model
    peierls = peierls_poisson_problem(
        hopping = hopping,
        phonon_frequency = peierls_frequency,
        coupling = peierls_coupling,
        dimension = dimension,
    ).model
    problem = VariationalProblem(combine_models(holstein, peierls), PoissonTrial(dimension = dimension, bare_hopping = hopping))
    return lattice_transport_sweep(problem; temperatures = temperatures, frequencies = frequencies, options = options, kwargs...)
end

"""
    write_sweep_csv(path, rows)

Write flat sweep/table rows to a comma-separated file without adding a plotting
or table dependency. Returns `path`.
"""
function write_sweep_csv(path::AbstractString, rows::AbstractVector{<:NamedTuple})
    columns = _table_columns(rows)
    open(path, "w") do io
        println(io, join(string.(columns), ","))
        for row in rows
            values = [_csv_escape(hasproperty(row, column) ? getproperty(row, column) : "") for column in columns]
            println(io, join(values, ","))
        end
    end
    return path
end

function _table_columns(rows)
    columns = Symbol[]
    seen = Set{Symbol}()
    for row in rows
        for column in keys(row)
            column in seen && continue
            push!(seen, column)
            push!(columns, column)
        end
    end
    return columns
end

function _csv_escape(value)
    text = string(value)
    if occursin(",", text) || occursin("\"", text) || occursin("\n", text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end
    return text
end

function _model_row(problem::VariationalProblem{FrohlichModel,T}) where {T<:AbstractGaussianTrial}
    model = problem.model
    representative_frequency = length(model.phonon_frequencies) == 1 ? first(model.phonon_frequencies) : model.effective_frequency
    return (;
        model = :frohlich,
        trial = Symbol(nameof(typeof(problem.trial))),
        dimension = model.dimension,
        coupling = sum(model.alpha),
        phonon_frequency = representative_frequency,
        phonon_modes = length(model.phonon_frequencies),
        adiabatic_ratio = inv(representative_frequency),
    )
end

function _model_row(problem::VariationalProblem{HolsteinModel,PoissonTrial})
    model = problem.model
    return (;
        model = :holstein,
        trial = :PoissonTrial,
        dimension = model.dimension,
        coupling = model.coupling,
        hopping = model.hopping,
        phonon_frequency = model.phonon_frequency,
        adiabatic_ratio = model.hopping / model.phonon_frequency,
    )
end

function _model_row(problem::VariationalProblem{PeierlsModel,PoissonTrial})
    model = problem.model
    return (;
        model = :peierls,
        trial = :PoissonTrial,
        dimension = model.dimension,
        coupling = model.coupling,
        hopping = model.hopping,
        phonon_frequency = model.phonon_frequency,
        adiabatic_ratio = model.hopping / model.phonon_frequency,
    )
end

function _model_row(problem::VariationalProblem{<:CompositePolaronModel,PoissonTrial})
    model = problem.model
    return (;
        model = :composite,
        trial = :PoissonTrial,
        dimension = model.dimension,
        components = length(model.models),
        hopping = model.hopping,
        effective_frequency = model.effective_frequency,
        coupling = sum(_component_coupling, model.models),
        phonon_frequency = _representative_lattice_frequency(model),
        adiabatic_ratio = _lattice_adiabatic_ratio(model),
    )
end

_component_coupling(model::HolsteinModel) = model.coupling
_component_coupling(model::PeierlsModel) = model.coupling

_representative_lattice_frequency(model::PeierlsModel) = model.phonon_frequency
_representative_lattice_frequency(model::HolsteinModel) = model.phonon_frequency
_representative_lattice_frequency(model::CompositePolaronModel) =
    sum(_representative_lattice_frequency(component) for component in model.models) / length(model.models)

_lattice_adiabatic_ratio(model::PeierlsModel) = model.hopping / model.phonon_frequency
_lattice_adiabatic_ratio(model::HolsteinModel) = model.hopping / model.phonon_frequency
_lattice_adiabatic_ratio(model::CompositePolaronModel) = model.hopping / _representative_lattice_frequency(model)

function _solution_identity_row(solution)
    return (;
        temperature = solution.temperature,
        beta = solution.beta,
    )
end

function _solution_row(problem::VariationalProblem{FrohlichModel,T}, solution::FrohlichSolution) where {T<:AbstractGaussianTrial}
    result = solution.variational
    return merge(
        _model_row(problem),
        _solution_identity_row(solution),
        (;
            free_energy = result.free_energy,
            reference_free_energy = result.reference_free_energy,
            entropy_cost = result.entropy_cost,
            interaction_free_energy = result.interaction_energy,
            v = solution.v,
            w = solution.w,
            energy_total = solution.energy.total,
            radius = solution.radius,
            asymptotic_mass = solution.asymptotic_mass,
        ),
        _parameter_row(result),
        result.diagnostics,
    )
end

function _solution_row(problem::VariationalProblem{HolsteinModel,PoissonTrial}, solution::LatticeSolution)
    result = solution.variational
    return merge(
        _model_row(problem),
        _solution_identity_row(solution),
        (;
            free_energy = result.free_energy,
            reference_free_energy = result.reference_free_energy,
            entropy_cost = result.entropy_cost,
            interaction_free_energy = result.interaction_energy,
            rate = solution.rate,
        ),
        _parameter_row(result),
        result.diagnostics,
    )
end

function _solution_row(problem::VariationalProblem{M,PoissonTrial}, solution::LatticeSolution) where {M<:Union{PeierlsModel,CompositePolaronModel}}
    result = solution.variational
    return merge(
        _model_row(problem),
        _solution_identity_row(solution),
        (;
            free_energy = result.free_energy,
            reference_free_energy = result.reference_free_energy,
            entropy_cost = result.entropy_cost,
            interaction_free_energy = result.interaction_energy,
            rate = solution.rate,
        ),
        _parameter_row(result),
        result.diagnostics,
    )
end

function _mobility_row(problem::VariationalProblem{FrohlichModel,T}, mobility::FrohlichMobilityResult) where {T<:AbstractGaussianTrial}
    return merge(
        _model_row(problem),
        (;
            mobility = mobility.mobility,
            fhip_low_temperature = mobility.fhip_low_temperature,
            kadanoff_devreese_low_temperature = mobility.kadanoff_devreese_low_temperature,
            kadanoff_low_temperature = mobility.kadanoff_low_temperature,
            relaxation_time = mobility.relaxation_time,
            hellwarth = mobility.hellwarth,
            hellwarth_b0 = mobility.hellwarth_b0,
        ),
    )
end

function _mobility_row(
    problem::VariationalProblem{M,PoissonTrial},
    mobility::LatticeMobilityResult,
) where {M<:Union{HolsteinModel,PeierlsModel,CompositePolaronModel}}
    return merge(
        _model_row(problem),
        (;
            mobility = mobility.mobility,
            mobility_einstein = mobility.mobility_einstein,
            mobility_factor = mobility.mobility_factor,
            diffusion_constant = mobility.diffusion_constant,
            mean_waiting_time = mobility.mean_waiting_time,
            total_jump_rate = mobility.total_jump_rate,
            response_hopping = mobility.response_hopping,
        ),
        _component_mobility_row(mobility.component_mobilities),
        mobility.diagnostics,
    )
end

function _response_row(problem::VariationalProblem{FrohlichModel,T}, response::FrohlichResponseResult) where {T<:AbstractGaussianTrial}
    return merge(
        _model_row(problem),
        (;
            frequency = response.frequency,
            memory_function_real = real(response.memory_function),
            memory_function_imag = imag(response.memory_function),
            memory_function_abs = abs(response.memory_function),
            impedance_real = real(response.impedance),
            impedance_imag = imag(response.impedance),
            impedance_abs = abs(response.impedance),
            conductivity_real = real(response.conductivity),
            conductivity_imag = imag(response.conductivity),
            conductivity_abs = abs(response.conductivity),
        ),
    )
end

function _response_row(
    problem::VariationalProblem{M,PoissonTrial},
    response::LatticeResponseResult,
) where {M<:Union{HolsteinModel,PeierlsModel,CompositePolaronModel}}
    return merge(
        _model_row(problem),
        (;
            frequency = response.frequency,
            mobility_real = real(response.mobility),
            mobility_imag = imag(response.mobility),
            mobility_abs = abs(response.mobility),
            mobility_factor_real = real(response.mobility_factor),
            mobility_factor_imag = imag(response.mobility_factor),
            mobility_factor_abs = abs(response.mobility_factor),
            conductivity_real = real(response.conductivity),
            conductivity_imag = imag(response.conductivity),
            conductivity_abs = abs(response.conductivity),
            impedance_real = real(response.impedance),
            impedance_imag = imag(response.impedance),
            impedance_abs = abs(response.impedance),
            response_hopping = response.response_hopping,
        ),
        _component_mobility_row(response.component_mobilities),
        response.diagnostics,
    )
end

function _component_mobility_row(components::NamedTuple)
    pairs = Pair{Symbol,Float64}[]
    for key in keys(components)
        value = getproperty(components, key)
        push!(pairs, Symbol(:component_mobility_, key, :_real) => Float64(real(value)))
        push!(pairs, Symbol(:component_mobility_, key, :_imag) => Float64(imag(value)))
        push!(pairs, Symbol(:component_mobility_, key, :_abs) => Float64(abs(value)))
    end
    return NamedTuple(pairs)
end

"""
    plot_coupling_sweep(rows; kwargs...)

Plot coupling-sweep rows. Load a plotting backend extension, such as Plots.jl,
to enable this function.
"""
plot_coupling_sweep(args...; kwargs...) = _plotting_extension_error()

"""
    plot_temperature_sweep(rows; kwargs...)

Plot temperature-sweep rows. Load a plotting backend extension, such as
Plots.jl, to enable this function.
"""
plot_temperature_sweep(args...; kwargs...) = _plotting_extension_error()

"""
    plot_adiabaticity_sweep(rows; kwargs...)

Plot adiabaticity-sweep rows. Load a plotting backend extension, such as
Plots.jl, to enable this function.
"""
plot_adiabaticity_sweep(args...; kwargs...) = _plotting_extension_error()

"""
    plot_frequency_sweep(rows; kwargs...)

Plot frequency-sweep rows. Load a plotting backend extension, such as Plots.jl,
to enable this function.
"""
plot_frequency_sweep(args...; kwargs...) = _plotting_extension_error()

"""
    plot_response_components(rows; kwargs...)

Plot frequency-response components such as mobility factor, conductivity,
impedance, or continuum memory-function observables. Load a plotting backend
extension to enable this function.
"""
plot_response_components(args...; kwargs...) = _plotting_extension_error()

function _plotting_extension_error()
    throw(ArgumentError("Plotting requires an optional plotting extension. Install and load Plots.jl before calling plot_* helpers."))
end

"""
    continued_frohlich_temperature_sweep(temperatures; coupling, phonon_frequency=1, dimension=3, options=OptimizerOptions())

Run forward and backward warm-started Fröhlich temperature continuations and
select the lower-free-energy branch at each point.
"""
function continued_frohlich_temperature_sweep(
    temperatures;
    coupling,
    phonon_frequency = 1.0,
    dimension::Integer = 3,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(temperatures, :temperature; positive = true)
    return _continued_model_sweep(
        values;
        build_problem = _ -> frohlich_feynman_problem(coupling = coupling, phonon_frequency = phonon_frequency, dimension = dimension),
        beta_from_value = temperature -> inv(temperature),
        row_temperature = identity,
        sweep_variable = :temperature,
        options = options,
    )
end

"""
    continued_frohlich_coupling_sweep(couplings; phonon_frequency=1, temperature, dimension=3, options=OptimizerOptions())

Run a Fröhlich coupling continuation at fixed reduced temperature.
"""
function continued_frohlich_coupling_sweep(
    couplings;
    phonon_frequency = 1.0,
    temperature::Real,
    dimension::Integer = 3,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(couplings, :coupling; nonnegative = true)
    return _continued_model_sweep(
        values;
        build_problem = coupling -> frohlich_feynman_problem(coupling = coupling, phonon_frequency = phonon_frequency, dimension = dimension),
        beta_from_value = _ -> inv(temperature),
        row_temperature = _ -> temperature,
        sweep_variable = :coupling,
        options = options,
    )
end

"""
    continued_frohlich_adiabaticity_sweep(adiabatic_ratios; coupling, temperature, dimension=3, options=OptimizerOptions())

Run a Fröhlich adiabaticity continuation where
`adiabatic_ratio = 1 / phonon_frequency`.
"""
function continued_frohlich_adiabaticity_sweep(
    adiabatic_ratios;
    coupling,
    temperature::Real,
    dimension::Integer = 3,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(adiabatic_ratios, :adiabatic_ratio; positive = true)
    return _continued_model_sweep(
        values;
        build_problem = ratio -> frohlich_feynman_problem(coupling = coupling, phonon_frequency = inv(ratio), dimension = dimension),
        beta_from_value = _ -> inv(temperature),
        row_temperature = _ -> temperature,
        sweep_variable = :adiabatic_ratio,
        options = options,
    )
end

"""
    continued_holstein_temperature_sweep(temperatures; hopping=1, phonon_frequency=1, coupling, dimension=1, options=OptimizerOptions())

Run forward and backward Holstein temperature continuations and select the
lower-free-energy branch at each point.
"""
function continued_holstein_temperature_sweep(
    temperatures;
    hopping = 1.0,
    phonon_frequency = 1.0,
    coupling,
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(temperatures, :temperature; positive = true)
    return _continued_model_sweep(
        values;
        build_problem = _ -> holstein_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension),
        beta_from_value = temperature -> inv(temperature),
        row_temperature = identity,
        sweep_variable = :temperature,
        options = options,
    )
end

"""
    continued_holstein_coupling_sweep(couplings; hopping=1, phonon_frequency=1, temperature, dimension=1, options=OptimizerOptions())

Run a Holstein coupling continuation at fixed reduced temperature.
"""
function continued_holstein_coupling_sweep(
    couplings;
    hopping = 1.0,
    phonon_frequency = 1.0,
    temperature::Real,
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(couplings, :coupling; nonnegative = true)
    return _continued_model_sweep(
        values;
        build_problem = coupling -> holstein_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension),
        beta_from_value = _ -> inv(temperature),
        row_temperature = _ -> temperature,
        sweep_variable = :coupling,
        options = options,
    )
end

"""
    continued_holstein_adiabaticity_sweep(adiabatic_ratios; hopping=1, coupling, temperature, dimension=1, options=OptimizerOptions())

Run a Holstein adiabaticity continuation where
`adiabatic_ratio = hopping / phonon_frequency`.
"""
function continued_holstein_adiabaticity_sweep(
    adiabatic_ratios;
    hopping = 1.0,
    coupling,
    temperature::Real,
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(adiabatic_ratios, :adiabatic_ratio; positive = true)
    return _continued_model_sweep(
        values;
        build_problem = ratio -> holstein_poisson_problem(hopping = hopping, phonon_frequency = hopping / ratio, coupling = coupling, dimension = dimension),
        beta_from_value = _ -> inv(temperature),
        row_temperature = _ -> temperature,
        sweep_variable = :adiabatic_ratio,
        options = options,
    )
end

"""
    continued_peierls_temperature_sweep(temperatures; hopping=1, phonon_frequency=1, coupling, dimension=1, options=OptimizerOptions())

Run forward and backward Peierls temperature continuations and select the
lower-free-energy branch at each point.
"""
function continued_peierls_temperature_sweep(
    temperatures;
    hopping = 1.0,
    phonon_frequency = 1.0,
    coupling,
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(temperatures, :temperature; positive = true)
    return _continued_model_sweep(
        values;
        build_problem = _ -> peierls_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension),
        beta_from_value = temperature -> inv(temperature),
        row_temperature = identity,
        sweep_variable = :temperature,
        options = options,
    )
end

"""
    continued_peierls_coupling_sweep(couplings; hopping=1, phonon_frequency=1, temperature, dimension=1, options=OptimizerOptions())

Run a Peierls coupling continuation at fixed reduced temperature.
"""
function continued_peierls_coupling_sweep(
    couplings;
    hopping = 1.0,
    phonon_frequency = 1.0,
    temperature::Real,
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(couplings, :coupling; nonnegative = true)
    return _continued_model_sweep(
        values;
        build_problem = coupling -> peierls_poisson_problem(hopping = hopping, phonon_frequency = phonon_frequency, coupling = coupling, dimension = dimension),
        beta_from_value = _ -> inv(temperature),
        row_temperature = _ -> temperature,
        sweep_variable = :coupling,
        options = options,
    )
end

"""
    continued_peierls_adiabaticity_sweep(adiabatic_ratios; hopping=1, coupling, temperature, dimension=1, options=OptimizerOptions())

Run a Peierls adiabaticity continuation where
`adiabatic_ratio = hopping / phonon_frequency`.
"""
function continued_peierls_adiabaticity_sweep(
    adiabatic_ratios;
    hopping = 1.0,
    coupling,
    temperature::Real,
    dimension::Integer = 1,
    options::OptimizerOptions = OptimizerOptions(),
)
    values = _validate_sweep_values(adiabatic_ratios, :adiabatic_ratio; positive = true)
    return _continued_model_sweep(
        values;
        build_problem = ratio -> peierls_poisson_problem(hopping = hopping, phonon_frequency = hopping / ratio, coupling = coupling, dimension = dimension),
        beta_from_value = _ -> inv(temperature),
        row_temperature = _ -> temperature,
        sweep_variable = :adiabatic_ratio,
        options = options,
    )
end

function _validate_sweep_values(values, name::Symbol; positive::Bool = false, nonnegative::Bool = false)
    collected = Float64.(collect(values))
    !isempty(collected) || throw(ArgumentError("$(name) sweep values must not be empty."))
    positive && !all(>(0), collected) && throw(DomainError(collected, "$(name) sweep values must be positive."))
    nonnegative && !all(>=(0), collected) && throw(DomainError(collected, "$(name) sweep values must be non-negative."))
    return collected
end

_transport_solve_grid(model::HolsteinModel, temperatures, frequencies) = _holstein_solve_grid(model, temperatures, frequencies)
_transport_solve_grid(model::PeierlsModel, temperatures, frequencies) = _lattice_solve_grid(model, temperatures, frequencies)
_transport_solve_grid(model::CompositePolaronModel, temperatures, frequencies) = _lattice_solve_grid(model, temperatures, frequencies)

function _transport_kernel_cache(sidebands, rate::Real, dimension::Integer, frequencies; broadening::Real, laguerre_points::Integer)
    cache = Dict{Float64,ComplexF64}()
    for base_frequency in frequencies, sideband in sidebands
        shift = Float64(base_frequency) + sideband.frequency
        haskey(cache, shift) && continue
        cache[shift] = first_return_laplace_d(
            ComplexF64(broadening, -shift),
            rate,
            dimension;
            laguerre_points = laguerre_points,
        )
    end
    return cache
end

function _cached_transport_sum(cache, sidebands, frequency::Real)
    accumulator = ComplexF64(0.0, 0.0)
    for sideband in sidebands
        accumulator += sideband.weight * cache[Float64(frequency) + sideband.frequency]
    end
    return accumulator
end

function _continued_model_sweep(values; build_problem, beta_from_value, row_temperature, sweep_variable::Symbol, options::OptimizerOptions)
    solve_point(value, initial) = _solve_sweep_point(
        build_problem(value),
        beta_from_value(value),
        row_temperature(value),
        initial,
        options,
    )
    return _continued_sweep(values; solve_point = solve_point, objective_key = :free_energy, sweep_variable = sweep_variable)
end

function _solve_sweep_point(problem, beta, temperature, initial, options)
    result = solve_variational(
        problem,
        beta;
        options = options,
        initial_parameters_override = initial,
        use_multistart = initial === nothing && options.multistart,
    )
    solution = solution_result(problem, temperature, result)
    mobility = mobility_result(problem, solution, options)
    row = merge(
        _variational_result_row(problem, result; temperature = temperature),
        _mobility_row(problem, mobility),
    )
    return row, result.parameters
end

function _continued_sweep(values; solve_point, objective_key::Symbol, sweep_variable::Symbol)
    forward_rows = _solve_sweep_branch(values, solve_point, :forward)
    backward_rows = _solve_sweep_branch(reverse(values), solve_point, :backward)
    backward_by_value = Dict(row[sweep_variable] => row for row in backward_rows)
    selected = NamedTuple[]
    for forward_row in forward_rows
        backward_row = backward_by_value[forward_row[sweep_variable]]
        best = forward_row[objective_key] <= backward_row[objective_key] ? forward_row : backward_row
        push!(selected, merge(
            best,
            (;
                selected_branch = best.sweep_direction,
                Symbol(:forward_, objective_key) => forward_row[objective_key],
                Symbol(:backward_, objective_key) => backward_row[objective_key],
            ),
        ))
    end
    return selected
end

function _solve_sweep_branch(values, solve_point, branch::Symbol)
    rows = NamedTuple[]
    previous_parameters = nothing
    for value in values
        row, previous_parameters = solve_point(value, previous_parameters)
        push!(rows, merge(row, (; sweep_direction = branch)))
    end
    branch == :backward && reverse!(rows)
    return rows
end

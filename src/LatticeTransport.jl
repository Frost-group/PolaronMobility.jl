const default_lattice_response_hopping = :bare
const default_lattice_broadening = 1e-6
const default_lattice_sideband_tolerance = 1e-12
const default_low_temperature_threshold = 40.0
const default_holstein_low_temperature_occupancy_threshold = 1e-8

"""
    lattice_impedance(mobility)

Reduced lattice impedance corresponding to a reduced per-carrier mobility or
conductivity. In the reduced convention `σ(Ω) = μ(Ω)` and `Z(Ω) = 1 / μ(Ω)`.
"""
function lattice_impedance(mobility::Number)
    value = ComplexF64(mobility)
    iszero(abs(value)) && return ComplexF64(Inf, 0.0)
    return ComplexF64(inv(value))
end

"""
    lattice_conductivity(mobility)

Reduced per-carrier lattice conductivity. In the reduced convention it is
identical to the reduced mobility.
"""
lattice_conductivity(mobility::Number) = ComplexF64(mobility)

"""
    lattice_holstein_phonon_factor(model, beta, time)

Exact-cloud Holstein blip factor
`exp[-S Cβ (1-cos(ω₀ t)) + i S sin(ω₀ t)]`, where `S = 2(g/ω₀)^2`.
"""
function lattice_holstein_phonon_factor(model::HolsteinModel, beta::Real, time::Real)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    time >= 0 || throw(DomainError(time, "time must be non-negative."))
    iszero(model.coupling) && return ComplexF64(1.0, 0.0)
    s = 2 * (model.coupling / model.phonon_frequency)^2
    angle = model.phonon_frequency * Float64(time)
    thermal_factor = beta == Inf ? 1.0 : Float64(coth(0.5 * Float64(beta) * model.phonon_frequency))
    return ComplexF64(exp(-s * thermal_factor * (1 - cos(angle)) + im * s * sin(angle)))
end

"""
    lattice_peierls_phonon_factor(model, rate, beta, time; response_hopping=:bare)

Normalized Peierls current-vertex cloud
`(J² + α² D⁾(t)) / (J² + α² D⁾(0))`, where `J` is the selected response
hopping scale.
"""
function lattice_peierls_phonon_factor(
    model::PeierlsModel,
    rate::Real,
    beta::Real,
    time::Real;
    response_hopping = default_lattice_response_hopping,
)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    time >= 0 || throw(DomainError(time, "time must be non-negative."))
    hopping_scale = _resolve_response_hopping(model, rate, response_hopping)
    iszero(model.coupling) && return ComplexF64(1.0, 0.0)
    nb = beta == Inf ? 0.0 : inv(exp(Float64(beta) * model.phonon_frequency) - 1)
    angle = model.phonon_frequency * Float64(time)
    dgreater = ComplexF64((nb + 1) * cis(-angle) + nb * cis(angle))
    denominator = hopping_scale^2 + model.coupling^2 * (2 * nb + 1)
    return ComplexF64((hopping_scale^2 + model.coupling^2 * dgreater) / denominator)
end

"""
    holstein_transport_sidebands(model, beta; tolerance=1e-12)

Build the exact-cloud Holstein transport sidebands used by the general-`d`
CTMC conductivity and mobility kernels.
"""
function holstein_transport_sidebands(
    model::HolsteinModel,
    beta::Real;
    tolerance::Real = default_lattice_sideband_tolerance,
)
    s = 2 * (model.coupling / model.phonon_frequency)^2
    return _holstein_sidebands(s, model.phonon_frequency, beta; tolerance = tolerance)
end

"""
    peierls_transport_sidebands(model, beta; response_hopping=:bare, tolerance=1e-12)

Build the normalized Peierls current-vertex sidebands: one zero-phonon channel
and finite-temperature `±ωP` assisted channels.
"""
function peierls_transport_sidebands(
    model::PeierlsModel,
    beta::Real;
    response_hopping = default_lattice_response_hopping,
    tolerance::Real = default_lattice_sideband_tolerance,
)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    tolerance > 0 || throw(DomainError(tolerance, "tolerance must be positive."))
    hopping_scale = _resolve_response_hopping(model, model.hopping, response_hopping)
    nb = beta == Inf ? 0.0 : inv(exp(Float64(beta) * model.phonon_frequency) - 1)
    denominator = hopping_scale^2 + model.coupling^2 * (2 * nb + 1)
    weights = [
        (; frequency = 0.0, weight = ComplexF64(hopping_scale^2 / denominator, 0.0)),
        (; frequency = -model.phonon_frequency, weight = ComplexF64(model.coupling^2 * (nb + 1) / denominator, 0.0)),
    ]
    nb > 0 && push!(weights, (; frequency = model.phonon_frequency, weight = ComplexF64(model.coupling^2 * nb / denominator, 0.0)))
    return _normalize_sidebands(filter(sideband -> abs(sideband.weight) > tolerance, weights))
end

"""
    holstein_peierls_transport_sidebands(model, beta; response_hopping=:bare, tolerance=1e-12)

Build the composite Holstein-Peierls transport sidebands by convolving the
independent Holstein cloud ladder with the Peierls current-vertex channels.
"""
function holstein_peierls_transport_sidebands(
    model::CompositePolaronModel,
    beta::Real;
    response_hopping = default_lattice_response_hopping,
    tolerance::Real = default_lattice_sideband_tolerance,
)
    sidebands = [(; frequency = 0.0, weight = ComplexF64(1.0, 0.0))]
    for component in model.models
        component_sidebands = _transport_sidebands(
            component,
            beta;
            response_hopping = response_hopping,
            tolerance = tolerance,
        )
        sidebands = _convolve_sidebands(sidebands, component_sidebands; tolerance = tolerance)
    end
    return _normalize_sidebands(sidebands)
end

"""
    lattice_mobility_factor(model, rate, beta, frequency; kwargs...)

Dimensionless lattice transport kernel `μ(Ω) / μE`, evaluated from the CTMC
first-return kernel and model-specific transport sidebands.
"""
function lattice_mobility_factor(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    response_hopping = default_lattice_response_hopping,
    time_cutoff = nothing,
    rtol::Real = 1e-6,
    broadening::Real = default_lattice_broadening,
    sideband_tolerance::Real = default_lattice_sideband_tolerance,
    laguerre_points::Integer = default_lattice_laguerre_points,
)
    rate > 0 || throw(DomainError(rate, "rate must be positive."))
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    broadening > 0 || throw(DomainError(broadening, "broadening must be positive."))
    sideband_tolerance > 0 || throw(DomainError(sideband_tolerance, "sideband_tolerance must be positive."))
    _ = time_cutoff
    _ = rtol
    return _lattice_transport_kernel(
        model,
        rate,
        beta,
        frequency;
        response_hopping = response_hopping,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
    )
end

function lattice_mobility_factor(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequencies::AbstractVector{<:Real};
    kwargs...,
)
    return [lattice_mobility_factor(model, rate, beta, frequency; kwargs...) for frequency in frequencies]
end

function lattice_mobility_factor(
    model::AbstractPolaronModel,
    trial::PoissonTrial,
    solution;
    frequencies = [0.0],
    kwargs...,
)
    rate = hasproperty(solution, :rate) ? getproperty(solution, :rate) : poisson_rate(trial, solution.parameters)
    beta = hasproperty(solution, :beta) ? getproperty(solution, :beta) : solution.beta
    return lattice_mobility_factor(model, rate, beta, _as_vector(frequencies); kwargs...)
end

"""
    lattice_mobility(model, rate, beta, frequency; kwargs...)

Reduced per-carrier lattice mobility for Holstein, Peierls, or compatible
composite lattice models evaluated at Poisson rate `rate`.
"""
function lattice_mobility(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    kwargs...,
)
    factor = lattice_mobility_factor(model, rate, beta, frequency; kwargs...)
    mobility_einstein = _einstein_mobility(rate, beta)
    return _scale_lattice_response(factor, mobility_einstein, frequency)
end

function lattice_mobility(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequencies::AbstractVector{<:Real};
    kwargs...,
)
    return [lattice_mobility(model, rate, beta, frequency; kwargs...) for frequency in frequencies]
end

function lattice_mobility(
    model::AbstractPolaronModel,
    trial::PoissonTrial,
    solution;
    frequencies = [0.0],
    kwargs...,
)
    rate = hasproperty(solution, :rate) ? getproperty(solution, :rate) : poisson_rate(trial, solution.parameters)
    beta = hasproperty(solution, :beta) ? getproperty(solution, :beta) : solution.beta
    return lattice_mobility(model, rate, beta, _as_vector(frequencies); kwargs...)
end

"""
    lattice_conductivity(model, rate, beta, frequency; kwargs...)

Reduced per-carrier lattice conductivity. In the reduced convention this is
equal to `lattice_mobility(...)`.
"""
function lattice_conductivity(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    kwargs...,
)
    return lattice_conductivity(lattice_mobility(model, rate, beta, frequency; kwargs...))
end

"""
    lattice_impedance(model, rate, beta, frequency; kwargs...)

Reduced lattice impedance corresponding to the CTMC conductivity kernel.
"""
function lattice_impedance(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    kwargs...,
)
    return lattice_impedance(lattice_mobility(model, rate, beta, frequency; kwargs...))
end

function _lattice_dc_mobility_result(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    temperature::Real;
    response_hopping = default_lattice_response_hopping,
    time_cutoff = nothing,
    rtol::Real = 1e-6,
    broadening::Real = default_lattice_broadening,
    sideband_tolerance::Real = default_lattice_sideband_tolerance,
    laguerre_points::Integer = default_lattice_laguerre_points,
    kappa_source::Symbol = :per_temperature,
)
    mobility_einstein = _einstein_mobility(rate, beta)
    mobility_factor = lattice_mobility_factor(
        model,
        rate,
        beta,
        0.0;
        response_hopping = response_hopping,
        time_cutoff = time_cutoff,
        rtol = rtol,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
    )
    mobility = _scale_lattice_dc_response(mobility_factor, mobility_einstein)
    total_jump_rate = 2 * _lattice_dimension(model) * Float64(rate)
    hopping_scale = _resolve_response_hopping(model, rate, response_hopping)
    components = _component_mobilities(
        model,
        rate,
        beta,
        0.0;
        response_hopping = response_hopping,
        time_cutoff = time_cutoff,
        rtol = rtol,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
    )
    diagnostics = _transport_diagnostics(
        model,
        rate,
        beta,
        0.0,
        mobility_factor;
        response_hopping = response_hopping,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
        kappa_source = kappa_source,
    )
    return LatticeMobilityResult(
        Float64(temperature),
        Float64(beta),
        Float64(mobility),
        Float64(mobility_einstein),
        Float64(real(mobility_factor)),
        Float64(rate),
        Float64(inv(max(total_jump_rate, 1e-15))),
        Float64(total_jump_rate),
        Float64(hopping_scale),
        components,
        diagnostics,
    )
end

function _lattice_response_result(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    temperature::Real,
    frequency::Real;
    response_hopping = default_lattice_response_hopping,
    time_cutoff = nothing,
    rtol::Real = 1e-6,
    broadening::Real = default_lattice_broadening,
    sideband_tolerance::Real = default_lattice_sideband_tolerance,
    laguerre_points::Integer = default_lattice_laguerre_points,
    kappa_source::Symbol = :per_temperature,
)
    mobility_einstein = _einstein_mobility(rate, beta)
    mobility_factor = lattice_mobility_factor(
        model,
        rate,
        beta,
        frequency;
        response_hopping = response_hopping,
        time_cutoff = time_cutoff,
        rtol = rtol,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
    )
    mobility = _scale_lattice_response(mobility_factor, mobility_einstein, frequency)
    conductivity = lattice_conductivity(mobility)
    hopping_scale = _resolve_response_hopping(model, rate, response_hopping)
    components = _component_mobilities(
        model,
        rate,
        beta,
        frequency;
        response_hopping = response_hopping,
        time_cutoff = time_cutoff,
        rtol = rtol,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
    )
    diagnostics = _transport_diagnostics(
        model,
        rate,
        beta,
        frequency,
        mobility_factor;
        response_hopping = response_hopping,
        broadening = broadening,
        sideband_tolerance = sideband_tolerance,
        laguerre_points = laguerre_points,
        kappa_source = kappa_source,
    )
    return LatticeResponseResult(
        Float64(temperature),
        Float64(beta),
        Float64(frequency),
        mobility,
        mobility_factor,
        conductivity,
        lattice_impedance(conductivity),
        Float64(hopping_scale),
        components,
        diagnostics,
    )
end

_lattice_dimension(model::HolsteinModel) = model.dimension
_lattice_dimension(model::PeierlsModel) = model.dimension
_lattice_dimension(model::CompositePolaronModel) = model.dimension

function _component_mobilities(
    model::HolsteinModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    kwargs...,
)
    return (; holstein = lattice_mobility(model, rate, beta, frequency; kwargs...))
end

function _component_mobilities(
    model::PeierlsModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    kwargs...,
)
    return (; peierls = lattice_mobility(model, rate, beta, frequency; kwargs...))
end

function _component_mobilities(
    model::CompositePolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    kwargs...,
)
    pairs = Pair{Symbol,ComplexF64}[]
    counts = Dict{Symbol,Int}()
    for component in model.models
        name = _component_name(component)
        count = get(counts, name, 0) + 1
        counts[name] = count
        label = count == 1 ? name : Symbol(name, count)
        value = lattice_mobility(component, rate, beta, frequency; kwargs...)
        push!(pairs, label => value)
    end
    return NamedTuple(pairs)
end

_component_name(::HolsteinModel) = :holstein
_component_name(::PeierlsModel) = :peierls

function _lattice_transport_kernel(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real;
    response_hopping,
    broadening::Real,
    sideband_tolerance::Real,
    laguerre_points::Integer,
)
    sidebands = _transport_sidebands(
        model,
        beta;
        response_hopping = response_hopping,
        tolerance = sideband_tolerance,
    )
    accumulator = ComplexF64(0.0, 0.0)
    for sideband in sidebands
        shift = Float64(frequency) + sideband.frequency
        accumulator += sideband.weight * first_return_laplace_d(
            ComplexF64(broadening, -shift),
            rate,
            _lattice_dimension(model);
            laguerre_points = laguerre_points,
        )
    end
    return accumulator
end

function lattice_transport_memory_factor(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real;
    kwargs...,
)
    return Float64(real(lattice_mobility_factor(model, rate, beta, 0.0; kwargs...)))
end

function _transport_sidebands(
    model::HolsteinModel,
    beta::Real;
    response_hopping = default_lattice_response_hopping,
    tolerance::Real,
)
    _ = response_hopping
    return holstein_transport_sidebands(model, beta; tolerance = tolerance)
end

function _transport_sidebands(
    model::PeierlsModel,
    beta::Real;
    response_hopping = default_lattice_response_hopping,
    tolerance::Real,
)
    return peierls_transport_sidebands(model, beta; response_hopping = response_hopping, tolerance = tolerance)
end

function _transport_sidebands(
    model::CompositePolaronModel,
    beta::Real;
    response_hopping = default_lattice_response_hopping,
    tolerance::Real,
)
    return holstein_peierls_transport_sidebands(model, beta; response_hopping = response_hopping, tolerance = tolerance)
end

function _holstein_sidebands(
    s::Real,
    phonon_frequency::Real,
    beta::Real;
    tolerance::Real,
)
    if iszero(s)
        return [(; frequency = 0.0, weight = ComplexF64(1.0, 0.0))]
    end
    thermal_scale = beta == Inf ? Inf : Float64(beta) * phonon_frequency
    if beta == Inf || thermal_scale > default_low_temperature_threshold
        return _zero_temperature_sidebands(s, phonon_frequency; tolerance = tolerance)
    end
    nb = inv(exp(thermal_scale) - 1)
    if nb <= default_holstein_low_temperature_occupancy_threshold
        return _zero_temperature_sidebands(s, phonon_frequency; tolerance = tolerance)
    end
    return _finite_temperature_sidebands(s, beta, phonon_frequency; tolerance = tolerance)
end

function _holstein_sidebands(model::HolsteinModel, beta::Real; tolerance::Real)
    return holstein_transport_sidebands(model, beta; tolerance = tolerance)
end

function _zero_temperature_sidebands(s::Real, phonon_frequency::Real; tolerance::Real)
    weights = NamedTuple[]
    weight = exp(-Float64(s))
    push!(weights, (; frequency = 0.0, weight = ComplexF64(weight, 0.0)))
    lmax = Int(ceil(s + 12 * sqrt(s + 1) + 40))
    current = weight
    for ell in 1:lmax
        current *= s / ell
        current <= tolerance && ell > s + 8 * sqrt(s + 1) && continue
        push!(weights, (; frequency = ell * Float64(phonon_frequency), weight = ComplexF64(0.5 * current, 0.0)))
        push!(weights, (; frequency = -ell * Float64(phonon_frequency), weight = ComplexF64(0.5 * current, 0.0)))
    end
    return _normalize_sidebands(weights)
end

function _finite_temperature_sidebands(s::Real, beta::Real, phonon_frequency::Real; tolerance::Real)
    x = Float64(beta) * phonon_frequency
    nb = inv(exp(x) - 1)
    lambda_plus = Float64(s) * (nb + 1)
    lambda_minus = Float64(s) * nb
    pmf_plus = _truncated_poisson_pmf(lambda_plus; tolerance = tolerance)
    pmf_minus = _truncated_poisson_pmf(lambda_minus; tolerance = tolerance)
    ell_weights = Dict{Int,Float64}()
    for (n, pn) in enumerate(pmf_plus), (m, pm) in enumerate(pmf_minus)
        weight = pn * pm
        weight <= tolerance && continue
        ell = (n - 1) - (m - 1)
        ell_weights[ell] = get(ell_weights, ell, 0.0) + weight
    end
    sidebands = NamedTuple[]
    for ell in sort!(collect(keys(ell_weights)))
        weight = ell_weights[ell]
        weight <= tolerance && continue
        push!(sidebands, (; frequency = ell * Float64(phonon_frequency), weight = ComplexF64(weight, 0.0)))
    end
    return _normalize_sidebands(sidebands)
end

function _truncated_poisson_pmf(lambda::Float64; tolerance::Real)
    lambda >= 0 || throw(DomainError(lambda, "Poisson mean must be non-negative."))
    lambda == 0 && return [1.0]
    nmax = Int(ceil(lambda + 12 * sqrt(lambda + 1) + 40))
    pmf = Vector{Float64}(undef, nmax + 1)
    mode = floor(Int, lambda)
    pmf[mode + 1] = exp(-lambda + mode * log(lambda) - loggamma(mode + 1))
    for n in mode:-1:1
        pmf[n] = pmf[n + 1] * n / lambda
    end
    for n in (mode + 1):nmax
        pmf[n + 1] = pmf[n] * lambda / n
    end
    total = sum(pmf)
    total > 0 || return [1.0]
    pmf ./= total
    while length(pmf) > 1 && pmf[end] <= tolerance
        pop!(pmf)
    end
    return pmf
end

function _convolve_sidebands(left::Vector, right::Vector; tolerance::Real)
    combined = NamedTuple[]
    for left_sideband in left, right_sideband in right
        _accumulate_sideband!(
            combined,
            left_sideband.frequency + right_sideband.frequency,
            left_sideband.weight * right_sideband.weight;
            tolerance = tolerance,
        )
    end
    return combined
end

function _normalize_sidebands(sidebands)
    isempty(sidebands) && return [(; frequency = 0.0, weight = ComplexF64(1.0, 0.0))]
    total_weight = sum(real(sideband.weight) for sideband in sidebands)
    abs(total_weight) > 0 || return [(; frequency = 0.0, weight = ComplexF64(1.0, 0.0))]
    return [(
        ; frequency = Float64(sideband.frequency),
          weight = ComplexF64(sideband.weight / total_weight),
    ) for sideband in sidebands]
end

function _accumulate_sideband!(sidebands, frequency::Real, weight::ComplexF64; tolerance::Real)
    for index in eachindex(sidebands)
        if isapprox(sidebands[index].frequency, frequency; atol = tolerance, rtol = 0.0)
            current = sidebands[index]
            sidebands[index] = (; frequency = current.frequency, weight = current.weight + weight)
            return sidebands
        end
    end
    push!(sidebands, (; frequency = Float64(frequency), weight = weight))
    return sidebands
end

function _transport_diagnostics(
    model::AbstractPolaronModel,
    rate::Real,
    beta::Real,
    frequency::Real,
    mobility_factor::ComplexF64;
    response_hopping,
    broadening::Real,
    sideband_tolerance::Real,
    laguerre_points::Integer,
    kappa_source::Symbol,
)
    sidebands = _transport_sidebands(
        model,
        beta;
        response_hopping = response_hopping,
        tolerance = sideband_tolerance,
    )
    return (
        broadening = Float64(broadening),
        sideband_tolerance = Float64(sideband_tolerance),
        laguerre_points = Int(laguerre_points),
        sideband_count = length(sidebands),
        sideband_weight_sum = Float64(sum(real(sideband.weight) for sideband in sidebands)),
        transport_kernel_real = Float64(real(mobility_factor)),
        transport_kernel_imag = Float64(imag(mobility_factor)),
        transport_frequency = Float64(frequency),
        kappa_source = kappa_source,
        rate = Float64(rate),
    )
end

function _resolve_response_hopping(model::AbstractPolaronModel, rate::Real, response_hopping)
    value = if response_hopping === :variational
        Float64(rate)
    elseif response_hopping === :bare
        hasproperty(model, :hopping) || throw(ArgumentError("response_hopping = :bare requires a model with a hopping field."))
        Float64(getproperty(model, :hopping))
    elseif response_hopping isa Real
        Float64(response_hopping)
    else
        throw(ArgumentError("response_hopping must be :variational, :bare, or a positive real hopping scale."))
    end
    value > 0 || throw(DomainError(value, "response_hopping must be positive."))
    return value
end

_einstein_mobility(rate::Real, beta::Real) = beta == Inf ? Inf : Float64(beta) * Float64(rate)

function _scale_lattice_dc_response(factor::ComplexF64, mobility_einstein::Real)
    mobility_einstein == Inf && return Inf
    return Float64(mobility_einstein) * Float64(real(factor))
end

function _scale_lattice_response(factor::ComplexF64, mobility_einstein::Real, frequency::Real)
    if mobility_einstein == Inf
        return iszero(frequency) ? ComplexF64(Inf, 0.0) : factor
    end
    return ComplexF64(Float64(mobility_einstein) * factor)
end

_first_return_laplace_1d(s::ComplexF64, rate::Real) = first_return_laplace_d(s, rate, 1)

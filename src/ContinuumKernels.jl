ball_surface(dimension::Integer) = 2 * π^(dimension / 2) / gamma(dimension / 2)

function phonon_propagator(τ, ω)
    return exp(-τ * ω)
end

function phonon_propagator(τ, ω, β)
    if β == Inf
        return phonon_propagator(τ, ω)
    end
    n = inv(exp(β * ω) - 1)
    return n > eps(Float64) ? n * exp(τ * ω) + (1 + n) * exp(-τ * ω) : exp(-τ * ω)
end

function polaron_propagator(τ, v, w)
    c = (v^2 - w^2) / v^3
    return (1 - v * c) * τ + c * (phonon_propagator(0, v) - phonon_propagator(τ, v)) + sqrt(eps(Float64))
end

function polaron_propagator(τ, v, w, β)
    β == Inf && return polaron_propagator(τ, v, w)
    c = (v^2 - w^2) / v^3
    return c * (phonon_propagator(0, v, β) - phonon_propagator(τ, v, β)) +
           (1 - c * v) * τ * (1 - τ / β) +
           sqrt(eps(Float64))
end

"""
    mean_square_displacement(trial, parameters, tau, beta)

Evaluate the Gaussian trial mean-square-displacement kernel used by Fröhlich
interaction integrals. Scalar Feynman, finite-mode Gaussian, and experimental
nonlocal Gaussian trials share this extension point.
"""
function mean_square_displacement(trial::GaussianFeynmanTrial, parameters, τ, β)
    v, w = unpack(trial, parameters)
    return beta_kernel(τ, v, w, β)
end

function beta_kernel(τ, v, w, β)
    return β == Inf ? polaron_propagator(τ, v, w) : polaron_propagator(τ, v, w, β)
end

function multi_gaussian_kappa(index::Integer, v::AbstractVector, w::AbstractVector)
    value = v[index]^2 - w[index]^2
    for other in eachindex(v)
        other == index && continue
        value *= (v[other]^2 - w[index]^2) / (w[other]^2 - w[index]^2)
    end
    return value
end

function multi_gaussian_h(index::Integer, v::AbstractVector, w::AbstractVector)
    value = v[index]^2 - w[index]^2
    for other in eachindex(v)
        other == index && continue
        value *= (w[other]^2 - v[index]^2) / (v[other]^2 - v[index]^2)
    end
    return value
end

function multi_gaussian_coupling(row::Integer, column::Integer, v::AbstractVector, w::AbstractVector)
    return w[row] * multi_gaussian_kappa(row, v, w) * multi_gaussian_h(column, v, w) /
           (4 * (v[column]^2 - w[row]^2))
end

function polaron_propagator(τ, v::AbstractVector, w::AbstractVector)
    length(v) == length(w) || throw(ArgumentError("v and w must have the same length."))
    value = τ
    for index in eachindex(v)
        value += multi_gaussian_h(index, v, w) / v[index]^2 *
                 ((1 - exp(-v[index] * τ)) / v[index] - τ)
    end
    return value + sqrt(eps(Float64))
end

function polaron_propagator(τ, v::AbstractVector, w::AbstractVector, β)
    β == Inf && return polaron_propagator(τ, v, w)
    length(v) == length(w) || throw(ArgumentError("v and w must have the same length."))
    value = τ * (1 - τ / β)
    for index in eachindex(v)
        numerator = 1 + exp(-v[index] * β) - exp(-v[index] * τ) - exp(v[index] * (τ - β))
        oscillator = numerator / (v[index] * (1 - exp(-v[index] * β)))
        value += multi_gaussian_h(index, v, w) / v[index]^2 * (oscillator - τ * (1 - τ / β))
    end
    return value + sqrt(eps(Float64))
end

function mean_square_displacement(trial::MultiGaussianTrial, parameters, τ, β)
    v, w = unpack(trial, parameters)
    return beta_kernel(τ, v, w, β)
end

function mean_square_displacement(trial::NonlocalGaussianTrial, parameters, τ, β)
    amplitudes = nonlocal_amplitudes(trial, parameters)
    base = β == Inf ? τ : τ * (1 - τ / β)
    value = base
    for (amplitude, frequency) in zip(amplitudes, trial.basis_frequencies)
        oscillator = if β == Inf
            (1 - exp(-frequency * τ)) / frequency
        else
            numerator = 1 + exp(-frequency * β) - exp(-frequency * τ) - exp(frequency * (τ - β))
            numerator / (frequency * (1 - exp(-frequency * β)))
        end
        value += amplitude * (oscillator - base)
    end
    return value + sqrt(eps(Float64))
end

"""
    profile_function(trial::ProfileGaussianTrial, parameters, omega)

Evaluate the positive Gaussian memory profile `Γ(ω)` for a
`ProfileGaussianTrial`.
"""
function profile_function(trial::ProfileGaussianTrial, parameters, omega)
    amplitudes = profile_amplitudes(trial, parameters)
    return sum(
        amplitude * frequency^2 / (omega^2 + frequency^2)
        for (amplitude, frequency) in zip(amplitudes, trial.basis_frequencies)
    )
end

function mean_square_displacement(trial::ProfileGaussianTrial, parameters, τ, β)
    decomposition = _profile_decomposition(trial, parameters)
    return _profile_mean_square_displacement(decomposition, τ, β)
end

function _profile_mean_square_displacement(decomposition, τ, β)
    residues, poles = decomposition
    base = β == Inf ? τ : τ * (1 - τ / β)
    isempty(residues) && return base + sqrt(eps(Float64))
    value = base
    for (residue, pole) in zip(residues, poles)
        oscillator = _profile_oscillator_kernel(pole, τ, β)
        value += residue / pole^2 * (base - oscillator)
    end
    return value + sqrt(eps(Float64))
end

function _profile_oscillator_kernel(pole, τ, β)
    if β == Inf
        return (1 - exp(-pole * τ)) / pole
    end
    numerator = 1 + exp(-pole * β) - exp(-pole * τ) - exp(pole * (τ - β))
    return numerator / (pole * (1 - exp(-pole * β)))
end

function _profile_decomposition(trial::ProfileGaussianTrial, parameters)
    amplitudes = profile_amplitudes(trial, parameters)
    iszero(sum(amplitudes)) && return Float64[], Float64[]

    frequencies2 = abs2.(trial.basis_frequencies)
    numerator = _profile_numerator_polynomial(frequencies2)
    denominator = copy(numerator)
    for index in eachindex(frequencies2)
        contribution = _profile_numerator_polynomial(frequencies2, skip = index)
        scale = amplitudes[index] * frequencies2[index]
        for coefficient in eachindex(contribution)
            denominator[coefficient] += scale * contribution[coefficient]
        end
    end

    roots = _polynomial_roots(denominator)
    residues = Float64[]
    poles = Float64[]
    derivative = _polynomial_derivative(denominator)
    for root in roots
        abs(imag(root)) < 1e-7 || throw(DomainError(root, "profile denominator roots must be real and negative."))
        pole2 = -real(root)
        pole2 > 0 || throw(DomainError(root, "profile denominator roots must be negative."))
        residue = real(_polynomial_eval(numerator, root) / _polynomial_eval(derivative, root))
        push!(poles, sqrt(pole2))
        push!(residues, residue)
    end
    order = sortperm(poles)
    return residues[order], poles[order]
end

function _profile_numerator_polynomial(frequencies2; skip::Union{Nothing,Int} = nothing)
    coefficients = [1.0]
    for index in eachindex(frequencies2)
        index == skip && continue
        coefficients = _poly_mul_x_plus(coefficients, frequencies2[index])
    end
    return coefficients
end

function _poly_mul_x_plus(coefficients, constant)
    result = zeros(Float64, length(coefficients) + 1)
    for index in eachindex(coefficients)
        result[index] += constant * coefficients[index]
        result[index + 1] += coefficients[index]
    end
    return result
end

function _polynomial_derivative(coefficients)
    length(coefficients) == 1 && return [0.0]
    return [index * coefficients[index + 1] for index in 1:(length(coefficients) - 1)]
end

function _polynomial_eval(coefficients, x)
    value = zero(x) + last(coefficients)
    for index in (length(coefficients) - 1):-1:1
        value = value * x + coefficients[index]
    end
    return value
end

function _polynomial_roots(coefficients)
    degree = length(coefficients) - 1
    degree == 0 && return ComplexF64[]
    leading = coefficients[end]
    abs(leading) > eps(Float64) || throw(ArgumentError("profile denominator polynomial has zero leading coefficient."))
    companion = zeros(Float64, degree, degree)
    for index in 2:degree
        companion[index, index - 1] = 1.0
    end
    normalized = coefficients[1:degree] ./ leading
    for index in 1:degree
        companion[index, degree] = -normalized[index]
    end
    return ComplexF64.(eigvals(companion))
end

function _profile_integral_0inf(f; points::Integer = 256)
    first_x = (0.5) * (π / 2) / points
    total = zero(f(tan(first_x)))
    step = (π / 2) / points
    for index in 1:points
        x = (index - 0.5) * step
        ω = tan(x)
        total += f(ω) / cos(x)^2
    end
    return step * total
end

trial_free_energy(v, w) = -3 * (v - w) / 2

function trial_free_energy(v, w, β)
    β == Inf && return trial_free_energy(v, w)
    return 3 / β * (
        log(v / w) -
        (v - w) * β / 2 -
        log(1 - exp(-v * β)) +
        log(1 - exp(-w * β))
    )
end

trial_correction_energy(v, w) = 3 * (v^2 - w^2) / (4v)

function trial_correction_energy(v, w, β)
    β == Inf && return trial_correction_energy(v, w)
    return 3 / 4 * (v^2 - w^2) / v * (coth(v * β / 2) - 2 / (v * β))
end

function trial_energy(v, w; dimension::Integer = 3)
    return trial_free_energy(v, w) * dimension / 3, trial_correction_energy(v, w) * dimension / 3
end

function trial_energy(v, w, β; dimension::Integer = 3)
    return trial_free_energy(v, w, β) * dimension / 3, trial_correction_energy(v, w, β) * dimension / 3
end

function trial_free_energy(v::AbstractVector, w::AbstractVector)
    length(v) == length(w) || throw(ArgumentError("v and w must have the same length."))
    return -3 * sum(v .- w) / 2
end

function trial_free_energy(v::AbstractVector, w::AbstractVector, β)
    β == Inf && return trial_free_energy(v, w)
    length(v) == length(w) || throw(ArgumentError("v and w must have the same length."))
    total = zero(promote_type(eltype(v), eltype(w), typeof(β)))
    for index in eachindex(v)
        total += log(v[index] / w[index]) -
                 (v[index] - w[index]) * β / 2 -
                 log(1 - exp(-v[index] * β)) +
                 log(1 - exp(-w[index] * β))
    end
    return 3 * total / β
end

function trial_correction_energy(v::AbstractVector, w::AbstractVector)
    length(v) == length(w) || throw(ArgumentError("v and w must have the same length."))
    total = zero(promote_type(eltype(v), eltype(w)))
    for row in eachindex(v), column in eachindex(v)
        total += multi_gaussian_coupling(row, column, v, w) / (v[column] * w[row])
    end
    return 3 * total
end

function trial_correction_energy(v::AbstractVector, w::AbstractVector, β)
    β == Inf && return trial_correction_energy(v, w)
    length(v) == length(w) || throw(ArgumentError("v and w must have the same length."))
    total = zero(promote_type(eltype(v), eltype(w), typeof(β)))
    for row in eachindex(v), column in eachindex(v)
        total += multi_gaussian_coupling(row, column, v, w) / (v[column] * w[row]) *
                 (coth(β * v[column] / 2) - 2 / (β * v[column]))
    end
    return 3 * total
end

function trial_energy(v::AbstractVector, w::AbstractVector; dimension::Integer = 3)
    return trial_free_energy(v, w) * dimension / 3, trial_correction_energy(v, w) * dimension / 3
end

function trial_energy(v::AbstractVector, w::AbstractVector, β; dimension::Integer = 3)
    return trial_free_energy(v, w, β) * dimension / 3, trial_correction_energy(v, w, β) * dimension / 3
end

"""
    frohlich_alpha(epsilon_optic, epsilon_static, freq_THz, band_mass)

Fröhlich coupling for a single LO phonon branch. `freq_THz` is in THz.
Follows the Feynman 1955 convention with the historical `4πϵ₀` normalization.
"""
function frohlich_alpha(epsilon_optic::Real, epsilon_static::Real, freq_THz::Real, band_mass::Real)
    ω = Float64(freq_THz) * 2π * 1e12
    return 1 / (2 * 4π * epsilon0) *
           (1 / Float64(epsilon_optic) - 1 / Float64(epsilon_static)) *
           (q^2 / (hbar * ω)) *
           sqrt(2 * Float64(band_mass) * me * ω / hbar)
end

"""
    frohlich_alpha(epsilon_optic, epsilon_ionic, epsilon_ionic_total, freq_THz, band_mass)

Partial Fröhlich coupling for one mode in a multi-phonon decomposition.
"""
function frohlich_alpha(
    epsilon_optic::Real,
    epsilon_ionic::Real,
    epsilon_ionic_total::Real,
    freq_THz::Real,
    band_mass::Real,
)
    ry = q^4 * me / (2 * hbar^2)
    ω = Float64(freq_THz) * 2π * 1e12
    epsilon_static = Float64(epsilon_ionic_total) + Float64(epsilon_optic)
    return sqrt(Float64(band_mass) * ry / (hbar * ω)) *
           Float64(epsilon_ionic) /
           (4π * epsilon0 * Float64(epsilon_optic) * epsilon_static)
end

"""
    dielectric_ionic_mode(freq_THz, infrared_activity, volume_m3)

Ionic dielectric contribution of one infrared-active phonon mode.
"""
function dielectric_ionic_mode(freq_THz::Real, infrared_activity::Real, volume_m3::Real)
    ω = Float64(freq_THz) * 2π * 1e12
    dielectric = q^2 * Float64(infrared_activity) / (3 * Float64(volume_m3) * ω^2 * amu)
    return dielectric / epsilon0
end

"""
    dielectric_ionic_total(modes, volume_m3)

Sum ionic dielectric contributions for a mode matrix whose first column is
frequency in THz and second column is infrared activity.
"""
function dielectric_ionic_total(modes::AbstractMatrix{<:Real}, volume_m3::Real)
    size(modes, 2) >= 2 || throw(ArgumentError("modes must have frequency and infrared-activity columns."))
    return sum(dielectric_ionic_mode(row[1], row[2], volume_m3) for row in eachrow(modes))
end

function frohlich_coupling(k, alpha, omega; dimension::Integer = 3)
    r_p = inv(sqrt(2))
    return omega^2 * alpha * r_p * gamma((dimension - 1) / 2) * (2 * sqrt(π) / k)^(dimension - 1)
end

function frohlich_interaction_energy(v, w, alpha::Real, omega::Real; dimension::Integer = 3, rtol::Real = 1e-4)
    coupling = frohlich_coupling(1, alpha, omega; dimension = dimension)
    integrand(τ) = phonon_propagator(τ, omega) / sqrt(polaron_propagator(τ, v, w) * omega)
    integral = quadgk(integrand, 0, Inf; rtol = rtol)[1]
    return coupling * ball_surface(dimension) / (2π)^dimension * sqrt(π / 2) * integral
end

function frohlich_interaction_energy(v, w, alpha::Real, omega::Real, beta::Real; dimension::Integer = 3, rtol::Real = 1e-4)
    beta == Inf && return frohlich_interaction_energy(v, w, alpha, omega; dimension = dimension, rtol = rtol)
    coupling = frohlich_coupling(1, alpha, omega; dimension = dimension)
    integrand(τ) = phonon_propagator(τ, omega, beta) / sqrt(polaron_propagator(τ, v, w, beta) * omega)
    integral = quadgk(integrand, 0, beta / 2; rtol = rtol)[1]
    return coupling * ball_surface(dimension) / (2π)^dimension * sqrt(π / 2) * integral
end

function frohlich_interaction_energy(
    v,
    w,
    alpha::AbstractVector{<:Real},
    omega::AbstractVector{<:Real};
    dimension::Integer = 3,
    rtol::Real = 1e-4,
)
    length(alpha) == length(omega) || throw(ArgumentError("alpha and omega must have the same length."))
    return sum(frohlich_interaction_energy(v, w, alpha[j], omega[j]; dimension = dimension, rtol = rtol) for j in eachindex(alpha))
end

function frohlich_interaction_energy(
    v,
    w,
    alpha::AbstractVector{<:Real},
    omega::AbstractVector{<:Real},
    beta::Real;
    dimension::Integer = 3,
    rtol::Real = 1e-4,
)
    length(alpha) == length(omega) || throw(ArgumentError("alpha and omega must have the same length."))
    return sum(frohlich_interaction_energy(v, w, alpha[j], omega[j], beta; dimension = dimension, rtol = rtol) for j in eachindex(alpha))
end

"""
    frohlich_energy(v, w, alpha, omega[, beta]; dimension=3)

Feynman/Osaka/Hellwarth variational free energy components. The returned
`EnergyComponents.total` is the variational energy/free energy to minimize.
"""
function frohlich_energy(v, w, alpha, omega; dimension::Integer = 3, rtol::Real = 1e-4)
    A, C = trial_energy(v, w; dimension = dimension)
    B = frohlich_interaction_energy(v, w, alpha, omega; dimension = dimension, rtol = rtol)
    return EnergyComponents(-(A + B + C), A, B, C)
end

function frohlich_energy(v, w, alpha, omega, beta::Real; dimension::Integer = 3, rtol::Real = 1e-4)
    beta == Inf && return frohlich_energy(v, w, alpha, omega; dimension = dimension, rtol = rtol)
    A, C = trial_energy(v, w, beta; dimension = dimension)
    B = frohlich_interaction_energy(v, w, alpha, omega, beta; dimension = dimension, rtol = rtol)
    return EnergyComponents(-(A + B + C), A, B, C)
end

function frohlich_structure_factor(t, v, w, alpha::Real, omega::Real; dimension::Integer = 3)
    coupling = frohlich_coupling(1, alpha, omega; dimension = dimension) * omega
    propagator = polaron_propagator(im * t, v, w) * omega / 2
    integral = ball_surface(dimension) / (2π)^dimension * sqrt(π) / 4 / propagator^(3 / 2)
    return 2 / dimension * coupling * integral * phonon_propagator(im * t, omega)
end

function frohlich_structure_factor(t, v, w, alpha::Real, omega::Real, beta::Real; dimension::Integer = 3)
    beta == Inf && return frohlich_structure_factor(t, v, w, alpha, omega; dimension = dimension)
    coupling = frohlich_coupling(1, alpha, omega; dimension = dimension) * omega
    propagator = polaron_propagator(im * t, v, w, beta) * omega / 2
    integral = ball_surface(dimension) / (2π)^dimension * sqrt(π) / 4 / propagator^(3 / 2)
    return 2 / dimension * coupling * integral * phonon_propagator(im * t, omega, beta)
end

function memory_integral(frequency::Real, structure_factor; cutoff::Real = 1e4, rtol::Real = 1e-4)
    if iszero(frequency)
        return quadgk(t -> -im * t * imag(structure_factor(t)), 0, cutoff; rtol = rtol)[1]
    end
    return quadgk(t -> (1 - exp(im * frequency * t)) / frequency * imag(structure_factor(t)), 0, cutoff; rtol = rtol)[1]
end

"""
    frohlich_memory_function(frequency, v, w, alpha, omega[, beta]; dimension=3, cutoff=1e4, rtol=1e-4)

Complex Fröhlich memory function for Feynman variational parameters `v,w`.
Scalar and multimode `alpha, omega` inputs are supported.
"""
function frohlich_memory_function(
    frequency::Real,
    v,
    w,
    alpha::Real,
    omega::Real,
    beta::Real = Inf;
    dimension::Integer = 3,
    cutoff::Real = 1e4,
    rtol::Real = 1e-4,
)
    structure_factor(t) = frohlich_structure_factor(t, v, w, alpha, omega, beta; dimension = dimension)
    return ComplexF64(memory_integral(frequency, structure_factor; cutoff = cutoff, rtol = rtol))
end

function frohlich_memory_function(
    frequency::Real,
    v,
    w,
    alpha::AbstractVector{<:Real},
    omega::AbstractVector{<:Real},
    beta::Real = Inf;
    dimension::Integer = 3,
    cutoff::Real = 1e4,
    rtol::Real = 1e-4,
)
    length(alpha) == length(omega) || throw(ArgumentError("alpha and omega must have the same length."))
    return sum(
        frohlich_memory_function(frequency, v, w, alpha[j], omega[j], beta; dimension = dimension, cutoff = cutoff, rtol = rtol)
        for j in eachindex(alpha)
    )
end

"""
    frohlich_complex_impedance(frequency, v, w, alpha, omega[, beta]; kwargs...)

Complex reduced impedance `Z(Ω) = -im * (Ω + Σ(Ω))` from the Fröhlich memory
function.
"""
function frohlich_complex_impedance(frequency, v, w, alpha, omega, beta = Inf; dimension::Integer = 3, cutoff::Real = 1e4, rtol::Real = 1e-4)
    return -im * (frequency + frohlich_memory_function(frequency, v, w, alpha, omega, beta; dimension = dimension, cutoff = cutoff, rtol = rtol))
end

"""
    frohlich_complex_conductivity(frequency, v, w, alpha, omega[, beta]; kwargs...)

Complex reduced conductivity, defined as the inverse of
`frohlich_complex_impedance`.
"""
function frohlich_complex_conductivity(frequency, v, w, alpha, omega, beta = Inf; dimension::Integer = 3, cutoff::Real = 1e4, rtol::Real = 1e-4)
    return inv(frohlich_complex_impedance(frequency, v, w, alpha, omega, beta; dimension = dimension, cutoff = cutoff, rtol = rtol))
end

function inverse_frohlich_mobility(v, w, alpha, omega, beta; dimension::Integer = 3, cutoff::Real = 1e4, rtol::Real = 1e-4)
    beta == Inf && return 0.0
    structure_factor(t) = frohlich_structure_factor(t, v, w, alpha, omega, beta; dimension = dimension)
    return abs(imag(memory_integral(0.0, structure_factor; cutoff = cutoff, rtol = rtol)))
end

function inverse_frohlich_mobility(
    v,
    w,
    alpha::AbstractVector{<:Real},
    omega::AbstractVector{<:Real},
    beta;
    dimension::Integer = 3,
    cutoff::Real = 1e4,
    rtol::Real = 1e-4,
)
    return sum(inverse_frohlich_mobility(v, w, alpha[j], omega[j], beta; dimension = dimension, cutoff = cutoff, rtol = rtol) for j in eachindex(alpha))
end

"""
    frohlich_mobility(v, w, alpha, omega, beta; dimension=3, cutoff=1e4, rtol=1e-4)

Reduced DC mobility from the zero-frequency Fröhlich memory function.
"""
function frohlich_mobility(v, w, alpha, omega, beta; dimension::Integer = 3, cutoff::Real = 1e4, rtol::Real = 1e-4)
    invμ = inverse_frohlich_mobility(v, w, alpha, omega, beta; dimension = dimension, cutoff = cutoff, rtol = rtol)
    return iszero(invμ) ? Inf : inv(invμ)
end

function inverse_fhip_low_temperature_mobility(v, w, alpha::Real, omega::Real, beta::Real)
    beta == Inf && return 0.0
    μ = (w / v)^3 * 3 / (4 * omega^2 * alpha * beta) * exp(omega * beta) *
        exp((v^2 - w^2) * omega / (w^2 * v))
    return inv(μ)
end

inverse_fhip_low_temperature_mobility(v, w, alpha::AbstractVector{<:Real}, omega::AbstractVector{<:Real}, beta) =
    sum(inverse_fhip_low_temperature_mobility(v, w, alpha[j], omega[j], beta) for j in eachindex(alpha))

"""
    fhip_low_temperature_mobility(v, w, alpha, omega, beta)

Low-temperature FHIP reduced mobility reference for Feynman variational
parameters.
"""
function fhip_low_temperature_mobility(v, w, alpha, omega, beta)
    invμ = inverse_fhip_low_temperature_mobility(v, w, alpha, omega, beta)
    return iszero(invμ) ? Inf : inv(invμ)
end

function inverse_kadanoff_low_temperature_mobility(v, w, alpha::Real, omega::Real, beta::Real)
    beta == Inf && return (0.0, 0.0, Inf)
    μ_devreese = (w / v)^3 / (2 * omega^2 * alpha) * exp(omega * beta) *
                 exp((v^2 - w^2) * omega / (w^2 * v))
    phonon_occupation = exp(-beta * omega)
    fictitious_mass = (v^2 - w^2) / w^2
    gamma0 = 2 * alpha * phonon_occupation * sqrt(fictitious_mass + 1) *
             exp(-fictitious_mass * omega / v) * omega^2
    μ_kadanoff = inv((fictitious_mass + 1) * gamma0)
    return inv(μ_devreese), inv(μ_kadanoff), inv(gamma0)
end

function inverse_kadanoff_low_temperature_mobility(v, w, alpha::AbstractVector{<:Real}, omega::AbstractVector{<:Real}, beta)
    inv_devreese = 0.0
    inv_kadanoff = 0.0
    inv_tau = 0.0
    for j in eachindex(alpha)
        d, k, τ = inverse_kadanoff_low_temperature_mobility(v, w, alpha[j], omega[j], beta)
        inv_devreese += d
        inv_kadanoff += k
        inv_tau += isinf(τ) ? 0.0 : inv(τ)
    end
    return inv_devreese, inv_kadanoff, iszero(inv_tau) ? Inf : inv(inv_tau)
end

"""
    kadanoff_low_temperature_mobility(v, w, alpha, omega, beta)

Low-temperature Kadanoff mobility references. Returns
`(devreese_mobility, kadanoff_mobility, relaxation_time)`.
"""
function kadanoff_low_temperature_mobility(v, w, alpha, omega, beta)
    inv_devreese, inv_kadanoff, relaxation_time = inverse_kadanoff_low_temperature_mobility(v, w, alpha, omega, beta)
    return iszero(inv_devreese) ? Inf : inv(inv_devreese),
           iszero(inv_kadanoff) ? Inf : inv(inv_kadanoff),
           relaxation_time
end

function inverse_hellwarth_mobility(v, w, alpha::Real, omega::Real, beta::Real; rtol::Real = 1e-3)
    beta == Inf && return (0.0, 0.0)
    R = (v^2 - w^2) / (w^2 * v)
    b = R * beta / sinh(beta * v / 2)
    a = sqrt((beta / 2)^2 + R * beta * coth(beta * v / 2))
    integrand(u, bvalue) = (u^2 + a^2 - bvalue * cos(v * u) + eps(Float64))^(-3 / 2) * cos(omega * u)
    K = quadgk(u -> integrand(u, b), 0, Inf; rtol = rtol)[1]
    K0 = quadgk(u -> integrand(u, 0), 0, Inf; rtol = rtol)[1]
    scale = alpha / (3 * sqrt(π)) * beta^(5 / 2) / sinh(omega * beta / 2) * (v^3 / w^3) * omega^(3 / 2)
    return scale * K, scale * K0
end

function inverse_hellwarth_mobility(v, w, alpha::AbstractVector{<:Real}, omega::AbstractVector{<:Real}, beta; rtol::Real = 1e-3)
    inv_full = 0.0
    inv_b0 = 0.0
    for j in eachindex(alpha)
        full, b0 = inverse_hellwarth_mobility(v, w, alpha[j], omega[j], beta; rtol = rtol)
        inv_full += full
        inv_b0 += b0
    end
    return inv_full, inv_b0
end

"""
    hellwarth_mobility(v, w, alpha, omega, beta; rtol=1e-3)

Hellwarth finite-temperature reduced mobility. Returns `(full, b0)` where `b0`
is the simplified `b = 0` comparison value.
"""
function hellwarth_mobility(v, w, alpha, omega, beta; rtol::Real = 1e-3)
    inv_full, inv_b0 = inverse_hellwarth_mobility(v, w, alpha, omega, beta; rtol = rtol)
    return iszero(inv_full) ? Inf : inv(inv_full),
           iszero(inv_b0) ? Inf : inv(inv_b0)
end

const lattice_bessel_asymptotic_threshold = 700.0
const default_lattice_laguerre_points = 120
const _laguerre_rule_cache = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const lattice_ellipk_tolerance = 1e-14
const lattice_ellipk_maxiter = 256

function _scaled_besseli0(argument::Real)
    x = argument
    x < lattice_bessel_asymptotic_threshold && return besselix(0, x)
    prefactor = inv(sqrt(2π * x))
    invx = inv(x)
    return prefactor * (1 + invx / 8 + 9 * invx^2 / 128 + 225 * invx^3 / 3072)
end

function _scaled_besseli1(argument::Real)
    x = argument
    x < lattice_bessel_asymptotic_threshold && return besselix(1, x)
    prefactor = inv(sqrt(2π * x))
    invx = inv(x)
    return prefactor * (1 - 3 * invx / 8 - 15 * invx^2 / 128 - 315 * invx^3 / 3072)
end

function _laguerre_rule(points::Integer)
    points > 0 || throw(ArgumentError("Laguerre quadrature points must be positive."))
    rule = get(_laguerre_rule_cache, Int(points), nothing)
    rule !== nothing && return rule
    diagonal = Float64[2 * index - 1 for index in 1:points]
    off_diagonal = Float64[index for index in 1:(points - 1)]
    decomposition = eigen(SymTridiagonal(diagonal, off_diagonal))
    nodes = Float64.(decomposition.values)
    weights = Float64.(decomposition.vectors[1, :] .^ 2)
    _laguerre_rule_cache[Int(points)] = (nodes, weights)
    return nodes, weights
end

function _complex_ellipk(parameter::ComplexF64)
    a = ComplexF64(1.0, 0.0)
    b = sqrt(ComplexF64(1.0, 0.0) - parameter)
    for _ in 1:lattice_ellipk_maxiter
        an = (a + b) / 2
        bn = sqrt(a * b)
        abs(an - bn) <= lattice_ellipk_tolerance * max(abs(an), 1.0) && return π / (2 * an)
        a = an
        b = bn
    end
    throw(ArgumentError("complex elliptic K AGM iteration did not converge."))
end

function _square_lattice_green_function(s::ComplexF64, rate::Real)
    a = s + ComplexF64(4 * Float64(rate), 0.0)
    modulus_parameter = (4 * Float64(rate) / a)^2
    return ComplexF64(2 * _complex_ellipk(modulus_parameter) / (π * a))
end

function _square_lattice_green_function_real(omega::Real, rate::Real)
    a = omega + 4 * rate
    modulus_parameter = (4 * rate / a)^2
    return 2 * ellipk(modulus_parameter) / (π * a)
end

function _cubic_lattice_green_function(
    s::ComplexF64,
    rate::Real;
    rtol::Real = 1e-8,
)
    integrand(angle) = _square_lattice_green_function(
        s + ComplexF64(2 * Float64(rate) * (1 - cos(angle)), 0.0),
        rate,
    )
    return ComplexF64(quadgk(integrand, 0, π; rtol = rtol)[1] / π)
end

function _cubic_lattice_green_function_real(
    omega::Real,
    rate::Real;
    rtol::Real = 1e-8,
)
    integrand(angle) = _square_lattice_green_function_real(omega + 2 * rate * (1 - cos(angle)), rate)
    return quadgk(integrand, 0, π; rtol = rtol)[1] / π
end

"""
    periodic_phonon_kernel(u, beta, omega)

Periodic imaginary-time Einstein phonon kernel
`Dβ(u;ω) = (exp(-ωu) + exp(-ω(β-u))) / (1 - exp(-βω))` for
`0 <= u <= beta`. For `beta == Inf`, this reduces to `exp(-ωu)`.
"""
function periodic_phonon_kernel(u::Real, beta::Real, omega::Real)
    omega > 0 || throw(DomainError(omega, "omega must be positive."))
    u >= 0 || throw(DomainError(u, "u must be non-negative."))
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    beta == Inf && return exp(-Float64(omega) * Float64(u))
    reduced_u = mod(Float64(u), Float64(beta))
    denominator = -expm1(-Float64(omega) * Float64(beta))
    return (exp(-Float64(omega) * reduced_u) + exp(-Float64(omega) * (Float64(beta) - reduced_u))) / denominator
end

"""
    lattice_q0(rate, dimension, tau)

Return the origin-return probability for a symmetric continuous-time Poisson
walk on a `dimension`-dimensional hypercubic lattice at imaginary time `tau`.
The scaled Bessel representation evaluates
`[exp(-2κτ) I₀(2κτ)]^d` stably.
"""
function lattice_q0(rate::Real, dimension::Integer, tau::Real)
    rate >= 0 || throw(DomainError(rate, "rate must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    tau >= 0 || throw(DomainError(tau, "tau must be non-negative."))
    argument = 2 * rate * tau
    return _scaled_besseli0(argument)^Int(dimension)
end

"""
    lattice_q1(rate, dimension, tau)

Nearest-neighbor displacement probability along one oriented lattice direction
for the symmetric CTMC walk, evaluated as
`exp(-2κτ) I₁(2κτ) [exp(-2κτ) I₀(2κτ)]^(d-1)`.
"""
function lattice_q1(rate::Real, dimension::Integer, tau::Real)
    rate >= 0 || throw(DomainError(rate, "rate must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    tau >= 0 || throw(DomainError(tau, "tau must be non-negative."))
    iszero(tau) && return 0.0
    argument = 2 * rate * tau
    return _scaled_besseli1(argument) * _scaled_besseli0(argument)^(Int(dimension) - 1)
end

"""
    lattice_current_kernel(rate, dimension, time)

Normalized open-path CTMC bond-current kernel used by the lattice-FHIP
mobility formulas. In reduced units
`Mκ(t) = 2d κ [q0(t) - q1(t)]`, whose integral is exactly one.
"""
function lattice_current_kernel(rate::Real, dimension::Integer, time::Real)
    rate > 0 || throw(DomainError(rate, "rate must be positive."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    time >= 0 || throw(DomainError(time, "time must be non-negative."))
    d = Int(dimension)
    return 2 * d * Float64(rate) * (lattice_q0(rate, d, time) - lattice_q1(rate, d, time))
end

"""
    site_return_bridge(rate, dimension, u, beta)

CTMC bridge site autocorrelation on a periodic imaginary-time loop. In the
large-`beta` limit it approaches `lattice_q0(rate, dimension, u)`.
"""
function site_return_bridge(rate::Real, dimension::Integer, u::Real, beta::Real)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    beta == Inf && return lattice_q0(rate, dimension, u)
    reduced_u = mod(Float64(u), Float64(beta))
    denominator = lattice_q0(rate, dimension, beta)
    return lattice_q0(rate, dimension, reduced_u) *
           lattice_q0(rate, dimension, Float64(beta) - reduced_u) /
           denominator
end

"""
    bond_order_bridge(rate, dimension, u, beta)

Periodic CTMC bridge approximation to a bond-order autocorrelation. Its
infinite-`beta` limit is `2d * (q0(u) + q1(u))`, the local bond-order weight
used by the full-periodic Peierls free energy.
"""
function bond_order_bridge(rate::Real, dimension::Integer, u::Real, beta::Real)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    d = Int(dimension)
    if beta == Inf
        return 2 * d * (lattice_q0(rate, d, u) + lattice_q1(rate, d, u))
    end
    reduced_u = mod(Float64(u), Float64(beta))
    mirrored = Float64(beta) - reduced_u
    denominator = lattice_q0(rate, d, beta)
    return 2 * d * (
        lattice_q0(rate, d, reduced_u) * lattice_q0(rate, d, mirrored) +
        lattice_q1(rate, d, reduced_u) * lattice_q1(rate, d, mirrored)
    ) / denominator
end

"""
    bond_current_bridge(rate, dimension, u, beta)

Periodic CTMC bridge approximation to a directed bond-current autocorrelation.
The infinite-`beta` limit is `2d * (q0(u) - q1(u))`.
"""
function bond_current_bridge(rate::Real, dimension::Integer, u::Real, beta::Real)
    beta > 0 || beta == Inf || throw(DomainError(beta, "beta must be positive or Inf."))
    d = Int(dimension)
    if beta == Inf
        return 2 * d * (lattice_q0(rate, d, u) - lattice_q1(rate, d, u))
    end
    reduced_u = mod(Float64(u), Float64(beta))
    mirrored = Float64(beta) - reduced_u
    denominator = lattice_q0(rate, d, beta)
    return 2 * d * (
        lattice_q0(rate, d, reduced_u) * lattice_q0(rate, d, mirrored) -
        lattice_q1(rate, d, reduced_u) * lattice_q1(rate, d, mirrored)
    ) / denominator
end

"""
    holstein_integral_d(rate, dimension, omega; laguerre_points=120)

Zero-temperature Holstein CTMC influence integral
`∫₀∞ exp(-ω t) ive0(2κt)^d dt`. For `d = 1` this reduces to
`1 / sqrt(ω(ω + 4κ))`; for `d = 2` the square-lattice Green's function gives
an elliptic-integral closed form; for `d = 3` the cubic-lattice Green's
function is reduced to a single finite-interval quadrature over that `d = 2`
kernel.
"""
function holstein_integral_d(
    rate::Real,
    dimension::Integer,
    omega::Real;
    laguerre_points::Integer = default_lattice_laguerre_points,
)
    rate >= 0 || throw(DomainError(rate, "rate must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    omega > 0 || throw(DomainError(omega, "omega must be positive."))
    d = Int(dimension)
    if d == 1
        return inv(sqrt(omega * (omega + 4 * rate)))
    end
    if d == 2
        return _square_lattice_green_function_real(omega, rate)
    end
    if d == 3
        return _cubic_lattice_green_function_real(omega, rate)
    end
    nodes, weights = _laguerre_rule(laguerre_points)
    total = zero(rate / omega)
    scale = inv(omega)
    for (node, weight) in zip(nodes, weights)
        argument = 2 * rate * node * scale
        total += weight * _scaled_besseli0(argument)^d
    end
    return scale * total
end

"""
    peierls_integral_d(rate, dimension, omega; laguerre_points=120)

Zero-temperature Peierls bond integral
`∫₀∞ exp(-ω t) 2d [ive0(2κt) + ive1(2κt)] ive0(2κt)^(d-1) dt`.
Rather than integrating this directly, the implementation rewrites it as
`2d (G₀ + G_a)` using the origin and nearest-neighbor lattice Green's
functions.
"""
function peierls_integral_d(
    rate::Real,
    dimension::Integer,
    omega::Real;
    laguerre_points::Integer = default_lattice_laguerre_points,
)
    rate >= 0 || throw(DomainError(rate, "rate must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    omega > 0 || throw(DomainError(omega, "omega must be positive."))
    d = Int(dimension)
    green = holstein_integral_d(rate, d, omega; laguerre_points = laguerre_points)
    first_return = if d == 1
        (omega + 2 * rate - sqrt(omega * (omega + 4 * rate))) / (2 * rate)
    else
        ((omega + 2 * d * rate) - inv(green)) / (2 * d * rate)
    end
    neighbor = green * first_return
    return 2 * d * real(green + neighbor)
end

"""
    lattice_green_function_d(s, rate, dimension; laguerre_points=120)

General-`d` CTMC return Green's function
`G₀^(d)(s) = ∫₀∞ exp(-s t) ive0(2κt)^d dt`, evaluated on the retarded sheet
for `real(s) > 0`. For `d = 1` an analytic expression is used, for `d = 2`
the square-lattice Green's function is expressed through the complete elliptic
integral `K`, and for `d = 3` the cubic-lattice Green's function is reduced to
a single finite-interval quadrature over the `d = 2` kernel.
"""
function lattice_green_function_d(
    s::ComplexF64,
    rate::Real,
    dimension::Integer;
    laguerre_points::Integer = default_lattice_laguerre_points,
)
    rate >= 0 || throw(DomainError(rate, "rate must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    real(s) > 0 || throw(DomainError(s, "real(s) must be positive for the retarded Laplace transform."))
    d = Int(dimension)
    if iszero(rate)
        return inv(s)
    end
    if d == 1
        return inv(sqrt(s * (s + ComplexF64(4 * Float64(rate), 0.0))))
    end
    if d == 2
        return _square_lattice_green_function(s, rate)
    end
    if d == 3
        return _cubic_lattice_green_function(s, rate)
    end
    nodes, weights = _laguerre_rule(laguerre_points)
    epsilon = real(s)
    phase_scale = imag(s) / epsilon
    spatial_scale = 2 * Float64(rate) / epsilon
    total = ComplexF64(0.0, 0.0)
    for (node, weight) in zip(nodes, weights)
        total += weight * cis(-phase_scale * node) * _scaled_besseli0(spatial_scale * node)^d
    end
    return total / epsilon
end

"""
    first_return_laplace_d(s, rate, dimension; laguerre_points=120)

General-`d` CTMC first-return kernel
`(s + 2dκ - 1/G₀^(d)(s)) / (2dκ)`. For `d = 1`, this reduces to the closed
analytic first-return Laplace transform.
"""
function first_return_laplace_d(
    s::ComplexF64,
    rate::Real,
    dimension::Integer;
    laguerre_points::Integer = default_lattice_laguerre_points,
)
    rate > 0 || throw(DomainError(rate, "rate must be positive."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    d = Int(dimension)
    if d == 1
        root = sqrt(s * (s + ComplexF64(4 * Float64(rate), 0.0)))
        return ComplexF64((s + 2 * Float64(rate) - root) / (2 * Float64(rate)))
    end
    green = lattice_green_function_d(s, rate, d; laguerre_points = laguerre_points)
    return ComplexF64((s + 2 * d * Float64(rate) - inv(green)) / (2 * d * Float64(rate)))
end

const RealVector = Vector{Float64}

"""
    AbstractPolaronModel

Supertype for dimensionless physical polaron models. Concrete models provide
couplings, dispersions, lattice parameters, or material metadata.
"""
abstract type AbstractPolaronModel end

"""
    AbstractContinuumModel

Supertype for continuum/path-integral polaron models with Gaussian trial
processes.
"""
abstract type AbstractContinuumModel <: AbstractPolaronModel end

"""
    AbstractLatticeModel

Supertype for lattice/worldline polaron models with jump-process trial spaces.
"""
abstract type AbstractLatticeModel <: AbstractPolaronModel end

"""
    AbstractMaterial

Supertype for material metadata objects that live at API/unit-conversion
boundaries. Concrete materials derive dimensionless model parameters but are
not used inside low-level physics kernels.
"""
abstract type AbstractMaterial end

"""
    AbstractTrialProcess

Supertype for variational trial processes. A concrete trial defines its
parameter names, initial parameters, bounds, and free/entropy objective terms.
"""
abstract type AbstractTrialProcess end

"""
    AbstractJumpTrial

Supertype for jump-process trial families used by lattice polaron models.
"""
abstract type AbstractJumpTrial <: AbstractTrialProcess end

"""
    AbstractGaussianTrial

Supertype for Gaussian trial processes used with continuum Fröhlich models.
Concrete trials define a mean-square-displacement kernel plus trial entropy
terms.
"""
abstract type AbstractGaussianTrial <: AbstractTrialProcess end

"""
    VariationalProblem(model, trial)

Pair a concrete polaron model with a compatible trial process. Generic solvers
minimize `free_energy + entropy_cost + interaction_free_energy`.
"""
struct VariationalProblem{M<:AbstractPolaronModel,T<:AbstractTrialProcess}
    model::M
    trial::T
end

"""
    EnergyComponents

Fröhlich energy/free-energy breakdown: total, trial oscillator contribution,
interaction contribution, and trial correction.
"""
struct EnergyComponents{T<:Number}
    total::T
    trial_free::T
    interaction::T
    trial_correction::T
end

"""
    OptimizerOptions(; lower=nothing, upper=nothing, initial_parameters=nothing, ...)

Numerical controls shared by all variational solvers. Bounds and initial
parameters use the parameter order reported by `parameter_names(trial)`.
Trial-specific defaults live on the trial type, not in these global options.
"""
struct OptimizerOptions
    lower::RealVector
    upper::RealVector
    initial_parameters::Union{Nothing,RealVector}
    gradient_tolerance::Float64
    quadrature_rtol::Float64
    memory_cutoff::Float64
    warn_on_nonconvergence::Bool
    warm_start::Bool
    multistart::Bool
    adaptive_bounds::Bool
    max_bound_expansions::Int
    bound_hit_fraction::Float64
end

function OptimizerOptions(;
    lower = nothing,
    upper = nothing,
    initial_parameters = nothing,
    gradient_tolerance::Real = sqrt(eps(Float64)),
    quadrature_rtol::Real = 1e-4,
    memory_cutoff::Real = 1e4,
    warn_on_nonconvergence::Bool = true,
    warm_start::Bool = true,
    multistart::Bool = true,
    adaptive_bounds::Bool = true,
    max_bound_expansions::Integer = 6,
    bound_hit_fraction::Real = 0.01,
)
    lower_values = lower === nothing ? Float64[] : _as_vector(lower)
    upper_values = upper === nothing ? Float64[] : _as_vector(upper)
    if !isempty(lower_values) && !isempty(upper_values)
        length(lower_values) == length(upper_values) ||
            throw(ArgumentError("lower and upper bounds must have the same length."))
        all(upper_values .> lower_values) ||
            throw(ArgumentError("all upper bounds must be greater than lower bounds."))
    end
    initial = initial_parameters === nothing ? nothing : _as_vector(initial_parameters)
    reference_length = !isempty(lower_values) ? length(lower_values) : length(upper_values)
    if initial !== nothing && reference_length > 0 && length(initial) != reference_length
        throw(ArgumentError("initial_parameters length must match supplied bounds."))
    end
    return OptimizerOptions(
        lower_values,
        upper_values,
        initial,
        Float64(gradient_tolerance),
        Float64(quadrature_rtol),
        Float64(memory_cutoff),
        warn_on_nonconvergence,
        warm_start,
        multistart,
        adaptive_bounds,
        Int(max_bound_expansions),
        Float64(bound_hit_fraction),
    )
end

"""
    GaussianFeynmanTrial(; dimension=3, initial_v=3.11, initial_w=2.87, ...)

Feynman's two-parameter Gaussian trial process for the continuum Fröhlich
polaron. The optimized coordinates are `(w, delta)` with `v = w + delta`.
"""
struct GaussianFeynmanTrial <: AbstractGaussianTrial
    dimension::Int
    initial_w::Float64
    initial_delta::Float64
    lower_w::Float64
    lower_delta::Float64
    upper_w::Float64
    upper_delta::Float64
end

function GaussianFeynmanTrial(;
    dimension::Integer = 3,
    initial_v::Real = 3.11,
    initial_w::Real = 2.87,
    lower_w::Real = 1e-8,
    lower_delta::Real = 0.0,
    upper_w::Real = 60.0,
    upper_delta::Real = 60.0,
)
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    initial_w > 0 || throw(ArgumentError("initial_w must be positive."))
    initial_v > initial_w || throw(ArgumentError("initial_v must be greater than initial_w."))
    upper_w > lower_w || throw(ArgumentError("upper_w must be greater than lower_w."))
    upper_delta > lower_delta || throw(ArgumentError("upper_delta must be greater than lower_delta."))
    return GaussianFeynmanTrial(
        Int(dimension),
        Float64(initial_w),
        Float64(initial_v - initial_w),
        Float64(lower_w),
        Float64(lower_delta),
        Float64(upper_w),
        Float64(upper_delta),
    )
end

"""
    MultiGaussianTrial(; modes, dimension=3, initial_v=nothing, initial_w=nothing, ...)

Finite-mode Gaussian trial process for Fröhlich variational calculations. The
solver parameter order is `w1, delta1, w2, delta2, ...`, with
`v_i = w_i + delta_i`. For `modes = 1` this reduces to Feynman's two-parameter
Gaussian trial.
"""
struct MultiGaussianTrial <: AbstractGaussianTrial
    modes::Int
    dimension::Int
    initial_w::RealVector
    initial_delta::RealVector
    lower_w::RealVector
    lower_delta::RealVector
    upper_w::RealVector
    upper_delta::RealVector
    min_separation::Float64
end

function MultiGaussianTrial(;
    modes::Integer,
    dimension::Integer = 3,
    initial_v = nothing,
    initial_w = nothing,
    lower_w::Real = 1e-8,
    lower_delta::Real = 1e-8,
    upper_w::Real = 60.0,
    upper_delta::Real = 60.0,
    min_separation::Real = 1e-7,
)
    modes > 0 || throw(ArgumentError("modes must be positive."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    lower_w > 0 || throw(ArgumentError("lower_w must be positive."))
    lower_delta >= 0 || throw(ArgumentError("lower_delta must be non-negative."))
    upper_w > lower_w || throw(ArgumentError("upper_w must be greater than lower_w."))
    upper_delta > lower_delta || throw(ArgumentError("upper_delta must be greater than lower_delta."))
    min_separation > 0 || throw(ArgumentError("min_separation must be positive."))

    default_w = [2.87 * (1 + 0.12 * (index - 1)) for index in 1:modes]
    default_v = [default_w[index] + (0.24 / index) for index in 1:modes]
    w = initial_w === nothing ? default_w : _as_vector(initial_w)
    v = initial_v === nothing ? default_v : _as_vector(initial_v)
    length(w) == modes || throw(ArgumentError("initial_w length must match modes."))
    length(v) == modes || throw(ArgumentError("initial_v length must match modes."))
    all(>(0), w) || throw(ArgumentError("all initial_w values must be positive."))
    all(v .> w) || throw(ArgumentError("all initial_v values must be greater than initial_w."))

    return MultiGaussianTrial(
        Int(modes),
        Int(dimension),
        Float64.(w),
        Float64.(v .- w),
        fill(Float64(lower_w), modes),
        fill(Float64(lower_delta), modes),
        fill(Float64(upper_w), modes),
        fill(Float64(upper_delta), modes),
        Float64(min_separation),
    )
end

"""
    NonlocalGaussianTrial(; basis_frequencies, dimension=3, ...)

Experimental Gaussian trial specified directly by a finite nonlocal memory
kernel basis. Parameters are positive basis amplitudes `a1, a2, ...`, and the
mean-square displacement is built from a sum of oscillator kernels at fixed
`basis_frequencies`.

This type is intended as an extensible kernel scaffold for general functional
approaches in the spirit of Adamowski, Gerlach, and Leschke. Its entropy term is
a configurable quadratic regularizer, not a literature-pinned closed-form
Fröhlich action.
"""
struct NonlocalGaussianTrial <: AbstractGaussianTrial
    basis_frequencies::RealVector
    dimension::Int
    initial_amplitudes::RealVector
    lower_amplitudes::RealVector
    upper_amplitudes::RealVector
    regularization::Float64
end

function NonlocalGaussianTrial(;
    basis_frequencies,
    dimension::Integer = 3,
    initial_amplitudes = nothing,
    lower_amplitude::Real = 0.0,
    upper_amplitude::Real = 10.0,
    regularization::Real = 1e-3,
)
    frequencies = _validate_positive("basis_frequencies", _as_vector(basis_frequencies))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    lower_amplitude >= 0 || throw(ArgumentError("lower_amplitude must be non-negative."))
    upper_amplitude > lower_amplitude || throw(ArgumentError("upper_amplitude must be greater than lower_amplitude."))
    regularization >= 0 || throw(ArgumentError("regularization must be non-negative."))
    amplitudes = initial_amplitudes === nothing ? fill(0.1, length(frequencies)) : _as_vector(initial_amplitudes)
    length(amplitudes) == length(frequencies) ||
        throw(ArgumentError("initial_amplitudes length must match basis_frequencies."))
    all(>=(0), amplitudes) || throw(ArgumentError("initial_amplitudes must be non-negative."))
    return NonlocalGaussianTrial(
        frequencies,
        Int(dimension),
        Float64.(amplitudes),
        fill(Float64(lower_amplitude), length(frequencies)),
        fill(Float64(upper_amplitude), length(frequencies)),
        Float64(regularization),
    )
end

"""
    ProfileGaussianTrial(; basis_frequencies, dimension=3, ...)

General Gaussian profile-function trial for Fröhlich polarons. The quadratic
trial action is represented by a positive memory profile
`Γ(ω) = sum(aᵢ νᵢ² / (ω² + νᵢ²))`, with optimized amplitudes `a1, a2, ...`
and fixed positive basis frequencies `νᵢ`.

The zero-temperature functional uses the standard Gaussian log-determinant
entropy term. The mean-square displacement is evaluated by analytically
decomposing the rational profile into oscillator kernels, so the one-pole
choice `Γ(ω) = (v² - w²)/(ω² + w²)` reproduces Feynman's one-mode trial
exactly.
"""
struct ProfileGaussianTrial <: AbstractGaussianTrial
    basis_frequencies::RealVector
    dimension::Int
    initial_amplitudes::RealVector
    lower_amplitudes::RealVector
    upper_amplitudes::RealVector
    matsubara_terms::Int
end

function ProfileGaussianTrial(;
    basis_frequencies,
    dimension::Integer = 3,
    initial_amplitudes = nothing,
    lower_amplitude::Real = 0.0,
    upper_amplitude::Real = 20.0,
    matsubara_terms::Integer = 4096,
)
    frequencies = _validate_positive("basis_frequencies", _as_vector(basis_frequencies))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    lower_amplitude >= 0 || throw(ArgumentError("lower_amplitude must be non-negative."))
    upper_amplitude > lower_amplitude || throw(ArgumentError("upper_amplitude must be greater than lower_amplitude."))
    matsubara_terms > 0 || throw(ArgumentError("matsubara_terms must be positive."))
    amplitudes = initial_amplitudes === nothing ? fill(0.1, length(frequencies)) : _as_vector(initial_amplitudes)
    length(amplitudes) == length(frequencies) ||
        throw(ArgumentError("initial_amplitudes length must match basis_frequencies."))
    all(>=(0), amplitudes) || throw(ArgumentError("initial_amplitudes must be non-negative."))
    return ProfileGaussianTrial(
        frequencies,
        Int(dimension),
        Float64.(amplitudes),
        fill(Float64(lower_amplitude), length(frequencies)),
        fill(Float64(upper_amplitude), length(frequencies)),
        Int(matsubara_terms),
    )
end

"""
    PoissonTrial(; dimension=1, bare_hopping=1, lower_rate=1e-12, upper_rate=max(10bare_hopping, 1))

Continuous-time symmetric nearest-neighbor hopping trial process for a Holstein
polaron on a hypercubic lattice. The single variational parameter is `rate`.
"""
struct PoissonTrial <: AbstractJumpTrial
    dimension::Int
    bare_hopping::Float64
    lower_rate::Float64
    upper_rate::Float64
end

function PoissonTrial(;
    dimension::Integer = 1,
    bare_hopping::Real = 1.0,
    lower_rate::Real = 1e-12,
    upper_rate::Real = max(10 * Float64(bare_hopping), 1.0),
)
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    bare_hopping > 0 || throw(ArgumentError("bare_hopping must be positive."))
    lower_rate >= 0 || throw(ArgumentError("lower_rate must be non-negative."))
    upper_rate > lower_rate || throw(ArgumentError("upper_rate must be greater than lower_rate."))
    return PoissonTrial(Int(dimension), Float64(bare_hopping), Float64(lower_rate), Float64(upper_rate))
end

"""
    FrohlichMaterial(...)

Polar material metadata and dimensionless Fröhlich coupling data. Constructors
derive mode couplings from dielectric data at the Unitful boundary.
"""
struct FrohlichMaterial <: AbstractMaterial
    alpha::RealVector
    band_mass::Float64
    phonon_frequencies::RealVector
    effective_frequency::Float64
    optical_dielectric::Float64
    static_dielectric::Float64
    ionic_dielectric::RealVector
    infrared_activity::RealVector
    cell_volume::Union{Nothing,Float64}
end

"""
    HolsteinMaterial(...)

Organic/lattice material metadata for a reduced Holstein model. Constructors
derive dimensionless hopping and local coupling from transfer integrals,
phonon wavenumbers, and local relaxation energies at the API boundary.
"""
struct HolsteinMaterial <: AbstractMaterial
    hopping_meV::Float64
    phonon_frequency_cm1::Float64
    phonon_frequency_THz::Float64
    phonon_energy_meV::Float64
    holstein_energy_meV::Float64
    hopping::Float64
    phonon_frequency::Float64
    coupling::Float64
    dimension::Int
    lattice_constant_angstrom::Union{Nothing,Float64}
    label::String
end

"""
    PeierlsMaterial(...)

Organic/lattice material metadata for a reduced Peierls bond-coupled model.
Constructors derive reduced hopping, phonon frequency, and bond-modulation
coupling from transfer integrals, phonon wavenumbers, and Peierls relaxation
energies at the API boundary.
"""
struct PeierlsMaterial <: AbstractMaterial
    hopping_meV::Float64
    phonon_frequency_cm1::Float64
    phonon_frequency_THz::Float64
    phonon_energy_meV::Float64
    peierls_energy_meV::Float64
    hopping::Float64
    phonon_frequency::Float64
    coupling::Float64
    dimension::Int
    lattice_constant_angstrom::Union{Nothing,Float64}
    label::String
end

"""
    FrohlichModel(alpha, phonon_frequencies=1; dimension=3, band_mass=1, ...)

Dimensionless continuum Fröhlich model. `alpha` and `phonon_frequencies` may be
scalars or equal-length vectors for multimode materials.
"""
struct FrohlichModel <: AbstractContinuumModel
    alpha::RealVector
    phonon_frequencies::RealVector
    dimension::Int
    band_mass::Float64
    effective_frequency::Float64
    material::Union{Nothing,FrohlichMaterial}
end

function FrohlichModel(
    alpha,
    phonon_frequencies = 1.0;
    dimension::Integer = 3,
    band_mass::Real = 1.0,
    effective_frequency::Real = first(_as_vector(phonon_frequencies)),
    material::Union{Nothing,FrohlichMaterial} = nothing,
)
    α = _validate_positive("alpha", _as_vector(alpha))
    ω = _validate_positive("phonon_frequencies", _as_vector(phonon_frequencies))
    length(α) == length(ω) ||
        throw(ArgumentError("alpha and phonon_frequencies must have the same length."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    band_mass > 0 || throw(ArgumentError("band_mass must be positive."))
    effective_frequency > 0 || throw(ArgumentError("effective_frequency must be positive."))
    return FrohlichModel(α, ω, Int(dimension), Float64(band_mass), Float64(effective_frequency), material)
end

"""
    HolsteinModel(; dimension=1, hopping=1, phonon_frequency=1, coupling)

Dimensionless Holstein model with local coupling, bare nearest-neighbor hopping,
and Einstein phonon frequency. `effective_frequency` stores the physical
reference frequency in THz when the model is derived from a `HolsteinMaterial`.
"""
struct HolsteinModel <: AbstractLatticeModel
    dimension::Int
    hopping::Float64
    phonon_frequency::Float64
    coupling::Float64
    effective_frequency::Float64
    material::Union{Nothing,HolsteinMaterial}
end

function HolsteinModel(;
    dimension::Integer = 1,
    hopping::Real = 1.0,
    phonon_frequency::Real = 1.0,
    coupling::Real,
    effective_frequency::Real = phonon_frequency,
    material::Union{Nothing,HolsteinMaterial} = nothing,
)
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    hopping > 0 || throw(ArgumentError("hopping must be positive."))
    phonon_frequency > 0 || throw(ArgumentError("phonon_frequency must be positive."))
    coupling >= 0 || throw(ArgumentError("coupling must be non-negative."))
    effective_frequency > 0 || throw(ArgumentError("effective_frequency must be positive."))
    return HolsteinModel(
        Int(dimension),
        Float64(hopping),
        Float64(phonon_frequency),
        Float64(coupling),
        Float64(effective_frequency),
        material,
    )
end

"""
    PeierlsModel(; dimension=1, hopping=1, phonon_frequency=1, coupling)

Dimensionless Peierls lattice model with bond-modulation coupling, bare
nearest-neighbor hopping, and Einstein bond phonon frequency. The Peierls
coupling acts on hopping/bond operators rather than site density.
"""
struct PeierlsModel <: AbstractLatticeModel
    dimension::Int
    hopping::Float64
    phonon_frequency::Float64
    coupling::Float64
    effective_frequency::Float64
    material::Union{Nothing,PeierlsMaterial}
end

function PeierlsModel(;
    dimension::Integer = 1,
    hopping::Real = 1.0,
    phonon_frequency::Real = 1.0,
    coupling::Real,
    effective_frequency::Real = phonon_frequency,
    material::Union{Nothing,PeierlsMaterial} = nothing,
)
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    hopping > 0 || throw(ArgumentError("hopping must be positive."))
    phonon_frequency > 0 || throw(ArgumentError("phonon_frequency must be positive."))
    coupling >= 0 || throw(ArgumentError("coupling must be non-negative."))
    effective_frequency > 0 || throw(ArgumentError("effective_frequency must be positive."))
    return PeierlsModel(
        Int(dimension),
        Float64(hopping),
        Float64(phonon_frequency),
        Float64(coupling),
        Float64(effective_frequency),
        material,
    )
end

"""
    CompositePolaronModel(models...)

Compatible polaron models that share one path/trial space. In this pass,
composition is intended for lattice CTMC models such as Holstein plus Peierls
with a common `PoissonTrial`.
"""
struct CompositePolaronModel{M<:Tuple} <: AbstractLatticeModel
    models::M
    dimension::Int
    hopping::Float64
    effective_frequency::Float64
end

function CompositePolaronModel(models::M) where {M<:Tuple}
    !isempty(models) || throw(ArgumentError("CompositePolaronModel requires at least one component model."))
    all(model -> model isa AbstractLatticeModel, models) ||
        throw(ArgumentError("all composite components must be AbstractLatticeModel instances."))
    first_model = first(models)
    hasproperty(first_model, :dimension) && hasproperty(first_model, :hopping) ||
        throw(ArgumentError("composite models must expose dimension and hopping fields."))
    dimension = getproperty(first_model, :dimension)
    hopping = getproperty(first_model, :hopping)
    effective_frequency = hasproperty(first_model, :effective_frequency) ? getproperty(first_model, :effective_frequency) : 1.0
    for model in models
        hasproperty(model, :dimension) && hasproperty(model, :hopping) ||
            throw(ArgumentError("composite models must expose dimension and hopping fields."))
        getproperty(model, :dimension) == dimension ||
            throw(ArgumentError("all composite models must have the same dimension."))
        isapprox(getproperty(model, :hopping), hopping; rtol = 1e-8, atol = 1e-12) ||
            throw(ArgumentError("all composite models must share the same reduced hopping."))
    end
    return CompositePolaronModel{M}(models, Int(dimension), Float64(hopping), Float64(effective_frequency))
end

"""
    VariationalResult

Optimized variational parameters and objective decomposition at one inverse
temperature `beta`.
"""
struct VariationalResult
    parameters::RealVector
    parameter_names::Vector{Symbol}
    beta::Float64
    free_energy::Float64
    reference_free_energy::Float64
    entropy_cost::Float64
    interaction_energy::Float64
    diagnostics::NamedTuple
end

"""
    ObservableResult

Generic observable container for extension APIs that only need temperature,
inverse temperature, and diagnostics.
"""
struct ObservableResult
    beta::Float64
    temperature::Float64
    diagnostics::NamedTuple
end

"""
    FrohlichSolution

Fröhlich-specific optimized solution at one temperature, including Feynman
`v,w`, energy components, effective-mass estimates, and radius diagnostics.
"""
struct FrohlichSolution
    temperature::Float64
    beta::Float64
    dimension::Int
    variational::VariationalResult
    parameters::NamedTuple
    v::Float64
    w::Float64
    energy::EnergyComponents{Float64}
    spring_constant::Float64
    fictitious_mass::Float64
    asymptotic_mass::Float64
    reduced_mass::Float64
    radius::Float64
end

"""
    FrohlichMobilityResult

Fröhlich DC mobility estimates at one temperature, including memory-function,
FHIP, Kadanoff, and Hellwarth references.
"""
struct FrohlichMobilityResult
    temperature::Float64
    beta::Float64
    mobility::Float64
    fhip_low_temperature::Float64
    kadanoff_devreese_low_temperature::Float64
    kadanoff_low_temperature::Float64
    relaxation_time::Float64
    hellwarth::Float64
    hellwarth_b0::Float64
end

"""
    FrohlichResponseResult

Fröhlich frequency-dependent memory function, impedance, and conductivity at one
temperature/frequency point.
"""
struct FrohlichResponseResult
    temperature::Float64
    beta::Float64
    frequency::Float64
    memory_function::ComplexF64
    impedance::ComplexF64
    conductivity::ComplexF64
end

"""
    LatticeSolution

Optimized Poisson variational solution for lattice polaron models at one
temperature.
"""
struct LatticeSolution
    temperature::Float64
    beta::Float64
    variational::VariationalResult
    rate::Float64
    free_energy::Float64
    entropy_cost::Float64
    interaction_energy::Float64
    diagnostics::NamedTuple
end

"""
    LatticeMobilityResult

Reduced DC transport result for optimized lattice Poisson trials. `mobility`
stores the transport-corrected lattice mobility, `mobility_einstein` keeps the
bare Einstein estimate, and `mobility_factor = mobility / mobility_einstein`
stores the reduced CTMC transport kernel at zero frequency for finite
temperature. `diagnostics`
records transport-specific metadata such as broadening, sideband
normalization, and rate provenance.
"""
struct LatticeMobilityResult{C<:NamedTuple,D<:NamedTuple}
    temperature::Float64
    beta::Float64
    mobility::Float64
    mobility_einstein::Float64
    mobility_factor::Float64
    diffusion_constant::Float64
    mean_waiting_time::Float64
    total_jump_rate::Float64
    response_hopping::Float64
    component_mobilities::C
    diagnostics::D
end

"""
    LatticeResponseResult

Reduced lattice response for Holstein, Peierls, and compatible composite
lattice Poisson models. `mobility_factor` stores the finite-frequency reduced
CTMC transport kernel. At finite temperature the package reports
`mobility = conductivity = βκ * mobility_factor`; at zero temperature and
finite frequency it instead reports the reduced optical kernel directly,
`mobility = conductivity = mobility_factor`, avoiding the divergent Einstein
prefactor. `impedance` is always the inverse reduced conductivity. `diagnostics`
records transport-specific metadata such as broadening, sideband
normalization, and rate provenance.
"""
struct LatticeResponseResult{C<:NamedTuple,D<:NamedTuple}
    temperature::Float64
    beta::Float64
    frequency::Float64
    mobility::ComplexF64
    mobility_factor::ComplexF64
    conductivity::ComplexF64
    impedance::ComplexF64
    response_hopping::Float64
    component_mobilities::C
    diagnostics::D
end

"""
    PolaronResult

Full solve result over temperature and frequency grids for a concrete
model-trial pair. `zero_temperature` stores the `beta = Inf` solution used as
the warm-start anchor even for lattice solves.
"""
struct PolaronResult{P<:VariationalProblem,S,MO,RO}
    problem::P
    temperatures::RealVector
    frequencies::RealVector
    zero_temperature::S
    solutions::Vector{S}
    mobilities::Vector{MO}
    responses::Matrix{RO}
end

_as_vector(x::Real) = [Float64(x)]
_as_vector(xs::AbstractVector{<:Real}) = Float64.(collect(xs))
_as_vector(xs::Tuple) = Float64.(collect(xs))

function _validate_positive(name::AbstractString, xs::AbstractVector{<:Real})
    all(>(0), xs) || throw(DomainError(xs, "$name must contain only positive values."))
    return xs
end

function Base.show(io::IO, problem::VariationalProblem)
    print(io, "VariationalProblem(model=$(typeof(problem.model)), trial=$(typeof(problem.trial)))")
end

function Base.show(io::IO, model::FrohlichModel)
    print(io, "FrohlichModel(alpha=$(model.alpha), phonon_frequencies=$(model.phonon_frequencies), dimension=$(model.dimension))")
end

function Base.show(io::IO, model::HolsteinModel)
    print(io, "HolsteinModel(coupling=$(model.coupling), hopping=$(model.hopping), phonon_frequency=$(model.phonon_frequency), dimension=$(model.dimension))")
end

function Base.show(io::IO, model::PeierlsModel)
    print(io, "PeierlsModel(coupling=$(model.coupling), hopping=$(model.hopping), phonon_frequency=$(model.phonon_frequency), dimension=$(model.dimension))")
end

function Base.show(io::IO, model::CompositePolaronModel)
    print(io, "CompositePolaronModel($(length(model.models)) components, dimension=$(model.dimension), hopping=$(model.hopping))")
end

function Base.show(io::IO, result::PolaronResult)
    print(
        io,
        "PolaronResult(model=$(nameof(typeof(result.problem.model))), trial=$(nameof(typeof(result.problem.trial))), ",
        "$(length(result.temperatures)) temperatures, $(length(result.frequencies)) frequencies)",
    )
end

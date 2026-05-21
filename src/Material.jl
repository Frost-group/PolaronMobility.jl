"""
    FrohlichMaterial(epsilon_optic, epsilon_static, band_mass, phonon_frequency)

Construct a single-mode polar material and derive its dimensionless Fröhlich
coupling from dielectric constants, band mass, and phonon frequency.
"""
function FrohlichMaterial(
    epsilon_optic::Real,
    epsilon_static::Real,
    band_mass::Real,
    phonon_frequency::Real,
)
    alpha = frohlich_alpha(epsilon_optic, epsilon_static, phonon_frequency, band_mass)
    ionic = Float64(epsilon_static - epsilon_optic)
    return FrohlichMaterial(
        [Float64(alpha)],
        Float64(band_mass),
        [Float64(phonon_frequency)],
        Float64(phonon_frequency),
        Float64(epsilon_optic),
        Float64(epsilon_static),
        [ionic],
        Float64[],
        nothing,
    )
end

"""
    FrohlichMaterial(epsilon_optic, epsilon_static, band_mass, phonon_frequencies, infrared_activity, volume)

Construct a multi-mode polar material. The Hellwarth B-scheme effective
frequency is used for reduced-unit conversion, while per-mode couplings are
derived from the infrared activities and mode frequencies.
"""
function FrohlichMaterial(
    epsilon_optic::Real,
    epsilon_static::Real,
    band_mass::Real,
    phonon_frequencies::AbstractVector{<:Real},
    infrared_activity::AbstractVector{<:Real},
    volume::Real,
)
    length(phonon_frequencies) == length(infrared_activity) ||
        throw(ArgumentError("phonon_frequencies and infrared_activity must have the same length."))
    frequencies = Float64.(collect(phonon_frequencies))
    activities = Float64.(collect(infrared_activity))
    effective_frequency = hellwarth_b_scheme(hcat(frequencies, activities))
    ionic = dielectric_ionic_mode.(frequencies, activities, Float64(volume))
    alpha = frohlich_alpha.(epsilon_optic, ionic, sum(ionic), frequencies, band_mass)
    return FrohlichMaterial(
        Float64.(alpha),
        Float64(band_mass),
        frequencies,
        Float64(effective_frequency),
        Float64(epsilon_optic),
        Float64(epsilon_static),
        Float64.(ionic),
        activities,
        Float64(volume),
    )
end

"""
    material_to_problem(material::FrohlichMaterial; dimension=3, trial=:feynman, kwargs...)

Convert a `FrohlichMaterial` into a reduced-unit Fröhlich
`VariationalProblem`. Pass Kelvin temperatures and THz frequencies directly to
`solve`.

`trial` may be `:feynman`, `:multi_gaussian`, `:profile_gaussian`, or
`:nonlocal_gaussian`. Profile and nonlocal trials use the material's reduced
phonon frequencies as their default basis frequencies.
"""
function material_to_problem(
    material::FrohlichMaterial;
    dimension::Integer = 3,
    trial::Symbol = :feynman,
    modes::Integer = max(1, min(length(material.alpha), 2)),
    basis_frequencies = nothing,
    initial_v = nothing,
    initial_w = nothing,
    initial_amplitudes = nothing,
    matsubara_terms::Integer = 4096,
    regularization::Real = 1e-3,
)
    reduced_ω = material.phonon_frequencies ./ material.effective_frequency
    basis = basis_frequencies === nothing ? reduced_ω : basis_frequencies
    common = (;
        coupling = material.alpha,
        phonon_frequency = reduced_ω,
        dimension = dimension,
        band_mass = material.band_mass,
        effective_frequency = material.effective_frequency,
        material = material,
    )
    if trial == :feynman
        return frohlich_feynman_problem(; common..., initial_v = initial_v, initial_w = initial_w)
    elseif trial == :multi_gaussian
        return frohlich_multi_gaussian_problem(; common..., modes = modes, initial_v = initial_v, initial_w = initial_w)
    elseif trial == :profile_gaussian
        return frohlich_profile_gaussian_problem(;
            common...,
            basis_frequencies = basis,
            initial_amplitudes = initial_amplitudes,
            matsubara_terms = matsubara_terms,
        )
    elseif trial == :nonlocal_gaussian
        return frohlich_nonlocal_gaussian_problem(;
            common...,
            basis_frequencies = basis,
            initial_amplitudes = initial_amplitudes,
            regularization = regularization,
        )
    end
    throw(ArgumentError("trial must be :feynman, :multi_gaussian, :profile_gaussian, or :nonlocal_gaussian."))
end

"""
    HolsteinMaterial(; hopping_meV, phonon_frequency_cm1, holstein_energy_meV, dimension=1, lattice_constant_angstrom=nothing, label="")

Construct a lattice material for a reduced Holstein model. `hopping_meV` is the
magnitude of the nearest-neighbor transfer integral, `phonon_frequency_cm1` is
the local effective phonon wavenumber, and `holstein_energy_meV` is the local
polaron binding/relaxation energy. The reduced model uses
`phonon_frequency = 1`, `hopping = J/(hν)`, and
`coupling = sqrt(E_H/(hν))`.
"""
function HolsteinMaterial(;
    hopping_meV::Real,
    phonon_frequency_cm1::Real,
    holstein_energy_meV::Real,
    dimension::Integer = 1,
    lattice_constant_angstrom = nothing,
    label::AbstractString = "",
)
    hopping_meV > 0 || throw(ArgumentError("hopping_meV must be positive; pass the magnitude of the transfer integral."))
    phonon_frequency_cm1 > 0 || throw(ArgumentError("phonon_frequency_cm1 must be positive."))
    holstein_energy_meV >= 0 || throw(ArgumentError("holstein_energy_meV must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    if lattice_constant_angstrom !== nothing
        lattice_constant_angstrom > 0 || throw(ArgumentError("lattice_constant_angstrom must be positive when supplied."))
    end

    frequency_THz = wavenumber_THz(phonon_frequency_cm1)
    phonon_energy = wavenumber_meV(phonon_frequency_cm1)
    return HolsteinMaterial(
        Float64(hopping_meV),
        Float64(phonon_frequency_cm1),
        frequency_THz,
        phonon_energy,
        Float64(holstein_energy_meV),
        Float64(hopping_meV) / phonon_energy,
        1.0,
        sqrt(Float64(holstein_energy_meV) / phonon_energy),
        Int(dimension),
        lattice_constant_angstrom === nothing ? nothing : Float64(lattice_constant_angstrom),
        String(label),
    )
end

"""
    PeierlsMaterial(; hopping_meV, phonon_frequency_cm1, peierls_energy_meV, dimension=1, lattice_constant_angstrom=nothing, label="")

Construct a lattice material for a reduced Peierls bond-coupled model.
`hopping_meV` is the magnitude of the nearest-neighbor transfer integral,
`phonon_frequency_cm1` is the bond phonon wavenumber, and
`peierls_energy_meV` sets the bond-modulation coupling through
`coupling = sqrt(E_P/(hν_P))`.
"""
function PeierlsMaterial(;
    hopping_meV::Real,
    phonon_frequency_cm1::Real,
    peierls_energy_meV::Real,
    dimension::Integer = 1,
    lattice_constant_angstrom = nothing,
    label::AbstractString = "",
)
    hopping_meV > 0 || throw(ArgumentError("hopping_meV must be positive; pass the magnitude of the transfer integral."))
    phonon_frequency_cm1 > 0 || throw(ArgumentError("phonon_frequency_cm1 must be positive."))
    peierls_energy_meV >= 0 || throw(ArgumentError("peierls_energy_meV must be non-negative."))
    dimension > 0 || throw(ArgumentError("dimension must be positive."))
    if lattice_constant_angstrom !== nothing
        lattice_constant_angstrom > 0 || throw(ArgumentError("lattice_constant_angstrom must be positive when supplied."))
    end

    frequency_THz = wavenumber_THz(phonon_frequency_cm1)
    phonon_energy = wavenumber_meV(phonon_frequency_cm1)
    return PeierlsMaterial(
        Float64(hopping_meV),
        Float64(phonon_frequency_cm1),
        frequency_THz,
        phonon_energy,
        Float64(peierls_energy_meV),
        Float64(hopping_meV) / phonon_energy,
        1.0,
        sqrt(Float64(peierls_energy_meV) / phonon_energy),
        Int(dimension),
        lattice_constant_angstrom === nothing ? nothing : Float64(lattice_constant_angstrom),
        String(label),
    )
end

"""
    rubrene_holstein_material(; direction=:high_mobility, dimension=1, lattice_constant_angstrom=nothing)

Return a one-mode Holstein material parameterization for rubrene from
Ordejón et al., Phys. Rev. B 96, 035202 (2017), Table II. The local-coupling
defaults are `E_H = 106.8 meV` and `ω_H = 1208.9 cm^-1`.

`direction` selects the transfer integral magnitude: `:high_mobility` or
`:AA_plus_b` uses `134.0 meV`, `:AC` uses `28.9 meV`, `:AB` uses `4.1 meV`,
and `:AA_plus_2b` uses `10.7 meV`.
"""
function rubrene_holstein_material(;
    direction::Symbol = :high_mobility,
    dimension::Integer = 1,
    lattice_constant_angstrom = nothing,
)
    hopping = _rubrene_hopping(direction)
    return HolsteinMaterial(
        hopping_meV = hopping,
        phonon_frequency_cm1 = 1208.9,
        holstein_energy_meV = 106.8,
        dimension = dimension,
        lattice_constant_angstrom = lattice_constant_angstrom,
        label = "rubrene Ordejon2017 Holstein $(direction)",
    )
end

"""
    rubrene_peierls_material(; direction=:high_mobility, dimension=1, lattice_constant_angstrom=nothing)

Return a one-mode Peierls material parameterization for rubrene from Ordejón
et al., Phys. Rev. B 96, 035202 (2017), Table II defaults
`J = 134.0 meV`, `E_P = 21.9 meV`, and `ω_P = 117.9 cm^-1`.
"""
function rubrene_peierls_material(;
    direction::Symbol = :high_mobility,
    dimension::Integer = 1,
    lattice_constant_angstrom = nothing,
)
    hopping = _rubrene_hopping(direction)
    return PeierlsMaterial(
        hopping_meV = hopping,
        phonon_frequency_cm1 = 117.9,
        peierls_energy_meV = 21.9,
        dimension = dimension,
        lattice_constant_angstrom = lattice_constant_angstrom,
        label = "rubrene Ordejon2017 Peierls $(direction)",
    )
end

function _rubrene_hopping(direction::Symbol)
    hoppings = Dict(
        :high_mobility => 134.0,
        :AA_plus_b => 134.0,
        :AC => 28.9,
        :AB => 4.1,
        :AA_plus_2b => 10.7,
    )
    haskey(hoppings, direction) ||
        throw(ArgumentError("direction must be :high_mobility, :AA_plus_b, :AC, :AB, or :AA_plus_2b."))
    return hoppings[direction]
end

"""
    lambda_holstein(material)

Dimensionless Holstein coupling ratio `λ = E_H / (2 d J)` for a
`HolsteinMaterial`.
"""
lambda_holstein(material::HolsteinMaterial) =
    material.holstein_energy_meV / (2 * material.dimension * material.hopping_meV)

"""
    lambda_peierls(material)

Dimensionless Peierls coupling ratio `λ_P = E_P / (2 d J)` for a
`PeierlsMaterial`.
"""
lambda_peierls(material::PeierlsMaterial) =
    material.peierls_energy_meV / (2 * material.dimension * material.hopping_meV)

"""
    material_to_problem(material::HolsteinMaterial; trial=:poisson, dimension=material.dimension)

Convert a `HolsteinMaterial` into a reduced-unit Holstein-Poisson
`VariationalProblem`. Pass Kelvin temperatures and THz frequencies directly to
`solve`; they are reduced using the material phonon frequency.
"""
function material_to_problem(
    material::HolsteinMaterial;
    trial::Symbol = :poisson,
    dimension::Integer = material.dimension,
)
    trial == :poisson || throw(ArgumentError("HolsteinMaterial currently supports only trial = :poisson."))
    return holstein_poisson_problem(
        hopping = material.hopping,
        phonon_frequency = material.phonon_frequency,
        coupling = material.coupling,
        dimension = dimension,
        effective_frequency = material.phonon_frequency_THz,
        material = material,
    )
end

"""
    material_to_problem(material::PeierlsMaterial; trial=:poisson, dimension=material.dimension)

Convert a `PeierlsMaterial` into a reduced-unit Peierls-Poisson
`VariationalProblem`. Pass Kelvin temperatures and THz frequencies directly to
`solve`; they are reduced using the material Peierls phonon frequency.
"""
function material_to_problem(
    material::PeierlsMaterial;
    trial::Symbol = :poisson,
    dimension::Integer = material.dimension,
)
    trial == :poisson || throw(ArgumentError("PeierlsMaterial currently supports only trial = :poisson."))
    return peierls_poisson_problem(
        hopping = material.hopping,
        phonon_frequency = material.phonon_frequency,
        coupling = material.coupling,
        dimension = dimension,
        effective_frequency = material.phonon_frequency_THz,
        material = material,
    )
end

"""
    rubrene_holstein_peierls_problem(; direction=:high_mobility, dimension=1, lattice_constant_angstrom=nothing)

Build a combined rubrene Holstein + Peierls Poisson variational problem using
Ordejón Table II defaults. The component models are reduced to one shared
Holstein phonon reference frequency before composition, so Kelvin and THz
inputs to `solve` remain unambiguous.
"""
function rubrene_holstein_peierls_problem(;
    direction::Symbol = :high_mobility,
    dimension::Integer = 1,
    lattice_constant_angstrom = nothing,
)
    holstein_material = rubrene_holstein_material(
        direction = direction,
        dimension = dimension,
        lattice_constant_angstrom = lattice_constant_angstrom,
    )
    peierls_material = rubrene_peierls_material(
        direction = direction,
        dimension = dimension,
        lattice_constant_angstrom = lattice_constant_angstrom,
    )
    reference_frequency = holstein_material.phonon_frequency_THz
    reference_energy = energy_meV(1.0, reference_frequency)
    holstein = HolsteinModel(
        dimension = dimension,
        hopping = holstein_material.hopping_meV / reference_energy,
        phonon_frequency = holstein_material.phonon_energy_meV / reference_energy,
        coupling = sqrt(holstein_material.holstein_energy_meV / reference_energy),
        effective_frequency = reference_frequency,
        material = holstein_material,
    )
    peierls = PeierlsModel(
        dimension = dimension,
        hopping = peierls_material.hopping_meV / reference_energy,
        phonon_frequency = peierls_material.phonon_energy_meV / reference_energy,
        coupling =
            sqrt(peierls_material.peierls_energy_meV *
                peierls_material.phonon_energy_meV) /
            reference_energy,
        effective_frequency = reference_frequency,
        material = peierls_material,
    )
    return VariationalProblem(
        combine_models(holstein, peierls),
        PoissonTrial(dimension = dimension, bare_hopping = holstein.hopping),
    )
end

function Base.show(io::IO, material::FrohlichMaterial)
    print(io, "FrohlichMaterial(alpha=$(material.alpha), band_mass=$(material.band_mass), phonon_frequencies=$(material.phonon_frequencies), effective_frequency=$(material.effective_frequency))")
end

function Base.show(io::IO, material::HolsteinMaterial)
    print(io, "HolsteinMaterial(label=\"$(material.label)\", hopping=$(material.hopping), coupling=$(material.coupling), phonon_frequency=$(material.phonon_frequency))")
end

function Base.show(io::IO, material::PeierlsMaterial)
    print(io, "PeierlsMaterial(label=\"$(material.label)\", hopping=$(material.hopping), coupling=$(material.coupling), phonon_frequency=$(material.phonon_frequency))")
end

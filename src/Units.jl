const hbar = 1.054571817e-34
const ħ = hbar
const q = 1.602176634e-19
const eV = q
const me = 9.1093837015e-31
const kB = 1.380649e-23
const epsilon0 = 8.85418682e-12
const ϵ_0 = epsilon0
const amu = 1.660_539_066_60e-27
const speed_of_light_cm_s = 2.99792458e10

"""
    MaterialUnits(;)

Default Unitful output units for material-derived Fröhlich results: THz,
`cm^2/V/s`, meV, and Å.
"""
struct MaterialUnits
    frequency_unit
    mobility_unit
    energy_unit
    radius_unit
end

"""
    FrohlichResultUnits

Unitful arrays converted from a dimensionless continuum `PolaronResult` by
`material_units`.
"""
struct FrohlichResultUnits
    mobility::Vector{typeof(1.0u"cm^2/V/s")}
    energy::Vector{typeof(1.0u"meV")}
    radius::Vector{typeof(1.0u"Å")}
    frequency::Vector{typeof(1.0u"THz")}
end

"""
    LatticeResultUnits

Unitful arrays converted from a material-derived lattice `PolaronResult` by
`material_units`. Mobilities require a lattice constant because the lattice
mobility scale is `e a^2 / ħ`.
"""
struct LatticeResultUnits
    mobility::Vector{typeof(1.0u"cm^2/V/s")}
    mobility_einstein::Vector{typeof(1.0u"cm^2/V/s")}
    energy::Vector{typeof(1.0u"meV")}
    frequency::Vector{typeof(1.0u"THz")}
end

MaterialUnits() = MaterialUnits(u"THz", u"cm^2/V/s", u"meV", u"Å")

reduced_temperature(temperature_K::Real, reference_frequency_THz::Real) =
    kB * Float64(temperature_K) / (hbar * 2π * Float64(reference_frequency_THz) * 1e12)

reduced_frequency(frequency_THz::Real, reference_frequency_THz::Real) =
    Float64(frequency_THz) / Float64(reference_frequency_THz)

"""
    wavenumber_THz(frequency_cm1)

Convert a phonon wavenumber in `cm^-1` to cycles-per-second THz.
"""
wavenumber_THz(frequency_cm1::Real) = Float64(frequency_cm1) * speed_of_light_cm_s / 1e12

"""
    wavenumber_meV(frequency_cm1)

Convert a phonon wavenumber in `cm^-1` to meV.
"""
wavenumber_meV(frequency_cm1::Real) = energy_meV(1.0, wavenumber_THz(frequency_cm1))

"""
    mobility_cm2_per_v_s(mobility, frequency_THz, band_mass)

Convert reduced mobility to `cm^2 V^-1 s^-1` using the reference phonon
frequency and carrier band mass in electron-mass units.
"""
mobility_cm2_per_v_s(mobility::Real, frequency_THz::Real, band_mass::Real) =
    Float64(mobility) * q / (2π * Float64(frequency_THz) * 1e12 * Float64(band_mass) * me) * 1e4

"""
    energy_meV(energy, frequency_THz)

Convert reduced energy to meV using a reference phonon frequency in THz.
"""
energy_meV(energy::Real, frequency_THz::Real) =
    Float64(energy) * hbar * 2π * Float64(frequency_THz) * 1e12 / q * 1e3

"""
    frequency_THz(frequency, reference_frequency_THz)

Convert a reduced frequency to THz.
"""
frequency_THz(frequency::Real, reference_frequency_THz::Real) =
    Float64(frequency) * Float64(reference_frequency_THz)

"""
    radius_angstrom(radius, frequency_THz, band_mass)

Convert a reduced polaron radius to Å using the reference phonon frequency and
carrier band mass in electron-mass units.
"""
radius_angstrom(radius::Real, frequency_THz::Real, band_mass::Real) =
    Float64(radius) * sqrt(hbar / (2π * Float64(frequency_THz) * 1e12 * Float64(band_mass) * me)) * 1e10

"""
    lattice_mobility_cm2_per_v_s(mobility, lattice_constant_angstrom)

Convert reduced lattice mobility to `cm^2 V^-1 s^-1` for a lattice spacing in
Å. The scale follows from `μ = μ̃ e a^2 / ħ`.
"""
lattice_mobility_cm2_per_v_s(mobility::Real, lattice_constant_angstrom::Real) =
    Float64(mobility) * q * (Float64(lattice_constant_angstrom) * 1e-10)^2 / hbar * 1e4

"""
    material_units(result::PolaronResult; units=MaterialUnits())

Return Unitful mobility, energy, radius, and frequency arrays for a material
derived Fröhlich result without mutating the original dimensionless result.
"""
function material_units(
    result::PolaronResult{P};
    units::MaterialUnits = MaterialUnits(),
) where {P<:VariationalProblem{<:FrohlichModel,<:AbstractGaussianTrial}}
    p = result.problem.model
    return FrohlichResultUnits(
        [mobility_cm2_per_v_s(m.mobility, p.effective_frequency, p.band_mass) * units.mobility_unit for m in result.mobilities],
        [energy_meV(s.energy.total, p.effective_frequency) * units.energy_unit for s in result.solutions],
        [radius_angstrom(s.radius, p.effective_frequency, p.band_mass) * units.radius_unit for s in result.solutions],
        [frequency_THz(Ω, p.effective_frequency) * units.frequency_unit for Ω in result.frequencies],
    )
end

"""
    material_units(result::PolaronResult; units=MaterialUnits())

Return Unitful lattice mobilities, solution energies, and response frequencies
for a material-derived Holstein result. A lattice constant is required because
the mobility conversion uses the lattice scale `e a^2 / ħ`.
"""
function material_units(
    result::PolaronResult{P};
    units::MaterialUnits = MaterialUnits(),
) where {P<:VariationalProblem{<:HolsteinModel,PoissonTrial}}
    model = result.problem.model
    material = model.material
    material === nothing && throw(ArgumentError("Holstein material unit conversion requires a HolsteinMaterial-derived problem."))
    lattice_constant = material.lattice_constant_angstrom
    lattice_constant === nothing && throw(ArgumentError("Holstein material mobility conversion requires lattice_constant_angstrom."))
    return LatticeResultUnits(
        [lattice_mobility_cm2_per_v_s(m.mobility, lattice_constant) * units.mobility_unit for m in result.mobilities],
        [lattice_mobility_cm2_per_v_s(m.mobility_einstein, lattice_constant) * units.mobility_unit for m in result.mobilities],
        [energy_meV(s.free_energy, model.effective_frequency) * units.energy_unit for s in result.solutions],
        [frequency_THz(Ω, model.effective_frequency) * units.frequency_unit for Ω in result.frequencies],
    )
end

"""
    material_units(result::PolaronResult; units=MaterialUnits())

Return Unitful lattice mobilities, solution energies, and response frequencies
for material-derived Peierls or compatible composite lattice results. A lattice
constant is required for mobility conversion.
"""
function material_units(
    result::PolaronResult{P};
    units::MaterialUnits = MaterialUnits(),
) where {M<:Union{PeierlsModel,CompositePolaronModel},P<:VariationalProblem{M,PoissonTrial}}
    model = result.problem.model
    material = _lattice_unit_material(model)
    material === nothing && throw(ArgumentError("Peierls material unit conversion requires a material-derived problem."))
    lattice_constant = material.lattice_constant_angstrom
    lattice_constant === nothing && throw(ArgumentError("Peierls material mobility conversion requires lattice_constant_angstrom."))
    effective_frequency = _lattice_unit_frequency(model)
    return LatticeResultUnits(
        [lattice_mobility_cm2_per_v_s(m.mobility, lattice_constant) * units.mobility_unit for m in result.mobilities],
        [lattice_mobility_cm2_per_v_s(m.mobility_einstein, lattice_constant) * units.mobility_unit for m in result.mobilities],
        [energy_meV(s.free_energy, effective_frequency) * units.energy_unit for s in result.solutions],
        [frequency_THz(Ω, effective_frequency) * units.frequency_unit for Ω in result.frequencies],
    )
end

_lattice_unit_material(model::PeierlsModel) = model.material
_lattice_unit_material(model::HolsteinModel) = model.material

function _lattice_unit_material(model::CompositePolaronModel)
    for component in model.models
        material = _lattice_unit_material(component)
        material !== nothing && return material
    end
    return nothing
end

_lattice_unit_frequency(model::PeierlsModel) = model.effective_frequency
_lattice_unit_frequency(model::HolsteinModel) = model.effective_frequency
_lattice_unit_frequency(model::CompositePolaronModel) = model.effective_frequency

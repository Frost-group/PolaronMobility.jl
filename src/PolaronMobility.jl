"""
PolaronMobility.jl - https://github.com/jarvist/PolaronMobility.jl
  Codes by Jarvist Moore Frost, 2017-2024; Bradley A.A. Martin 2020-2024
  PolaronMobility.jl is a Julia package which calculates the temperature-dependent polaron
  mobility for a material. We implement the Feynman variational approach with
  Osaka/Hellwarth free-energy for the polaron system.  
"""
module PolaronMobility

 # Type to hold the Frohlich polaron data.
export FrohlichPolaron, frohlichpolaron                                  

 # Type to hold material specific data.
export Material, material, save_material, load_material            

# Type to hold Holstein polaron data.
export HolsteinPolaron, holsteinpolaron

# Frohlich alpha, energy and minimisation functions.   
export frohlichalpha, ϵ_ionic_mode, frohlich_energy, A, B, C, feynmanvw

# Hellwarth effective mode schemes.
export HellwarthBScheme, HellwarthAScheme                 

# Polaron memory functions
export frohlich_memory_function, holstein_memory_function                    

# Complex Response functions
export frohlich_complex_impedence, frohlich_complex_conductivity, holstein_complex_conductivity, holstein_complex_impedence

# Mobility functions
export frohlich_mobility, Hellwarth_mobility, Kadanoff_mobility_lowT, FHIP_mobility_lowT 

# Functions to save and load Polaron types to/from .jld format
export save_frohlich_polaron, load_frohlich_polaron, save_holstein_polaron, load_holstein_polaron

# Polaron effective mass
export polaron_effective_mass

# Generalised polaron functions
export polaron_memory_function, vw_variation, polaron_propagator, phonon_propagator, trial_energy

# Holstein functions
export holstein_interaction_energy, holstein_energy, holstein_structure_factor, holstein_memory_function, holstein_mobility

export frohlich_structure_factor

# K-space functions
export spherical_k_integral, holstein_coupling, frohlich_coupling,holstein_interaction_energy_k_space, frohlich_interaction_energy_k_space, holstein_energy_k_space, frohlich_energy_k_space, holstein_structure_factor_k_space, frohlich_structure_factor_k_space, holstein_memory_function_k_space, holstein_mobility_k_space, frohlich_memory_function_k_space, frohlich_mobility_k_space

##### load in library routines... #####
# stdlib
using LinearAlgebra, Printf, JLD

#Key numerics used:

# One-dimensional numerical integration in Julia using adaptive Gauss-Kronrod quadrature
import QuadGK.quadgk, QuadGK.quadgk_count

# Special function arising from cutoff k-space integrals (c.f. Holstein model).
import SpecialFunctions.erf, SpecialFunctions.gamma, SpecialFunctions.gamma_inc

# Using the powerful Julia Optim package to optimise the variational parameters
using Optim

# Import required Unitful units, for mutability (generating new polaron units)
using Unitful
using Unitful: @unit, Dimension, Dimensions, NoDims, NoUnits, Units, dimension, uconvert, ustrip

# Unitful functions
export puconvert, punit, pustrip, m0_pu, e_pu, ħ_pu, k_pu, ω0_pu, a0_pu, E0_pu, β0_pu, T0_pu, μ0_pu, t0_pu, addunits!, auconvert!, suconvert!, phuconvert, phunit, phustrip

include("PolaronFunctions.jl")  # General use functions for polaron model.
include("LegacyFunctions.jl")   # Outdated functions that still have some use.
include("TrialPolaron.jl")      # Feynman's Trial Spring-Mass polaron model.
include("FrohlichTheory.jl")    # Frohlich polaron theory.
include("HolsteinTheory.jl")    # Holstein poalron theory.
include("HellwarthTheory.jl")   # Multimode -> equivalent effective mode.
include("Material.jl")          # Material type and constructors.
include("FrohlichPolaron.jl")   # Polaron type and constructors.
include("HolsteinPolaron.jl")   # Holstein type and constructors.
include("PolaronUnits.jl")      # Implements internal Feynman 'polaron units', and SI conversions

# Register newly defined units with Unitful
Unitful.register(PolaronMobility)
__init__() = Unitful.register(PolaronMobility)

# Utility functions
export reduce_array, P, P_plus_one, ball_surface

end # module


# HalidePerovskites
# Uses the newly forked off FeynmanKadanoffOsakaHellwarth.jl package
# Codes by Jarvist Moore Frost, 2017
# These codes were developed with Julia 0.5.0, and requires the Optim and Plots packages.

# This file, is based on doing further calculations beyond:
# https://arxiv.org/abs/1704.05404
# Polaron mobility in halide perovskites
# Jarvist Moore Frost
# (Submitted on 18 Apr 2017 [v1])

push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility 

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

polaronmobility(x...) = PolaronMobility.polaronmobility(x..., verbose=true)

for T in [300, 93]

#####
# Call simulation

# CsSnX3 X={Cl,Br,I}
# L. Huang, W.Lambrecht - PRB 88, 165203 (2013)
# Dielectric consts, from TABLE VII
# Effective masses from TABLE VI, mh*
cm1=2.997e10 # cm-1 to Herz

#Ts,Kμs, Hμs, FHIPμs, ks, Ms, As, Bs, Cs, Fs, Taus
#effectivemass=0.12 # the bare-electron band effective-mass. 
# --> 0.12 for electrons and 0.15 for holes, in MAPI. See 2014 PRB.
# MAPI  4.5, 24.1, 2.25THz - 75 cm^-1 ; α=

println()
println("#"^75)
println("    MAPI @ T = ",T)
MAPIe=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.12)
#MAPIh=polaronmobility(T, 4.5, 24.1, 2.25E12, 0.15)

println()
println("#"^75)
println("    CsPbBr3 @ T = ",T)
# CsPbBr3 — ϵ_optic=4.96, ϵ_static=16.5, phonon_freq=3.08 THz
# ϵ from Huang & Lambrecht PRB 90 195201 (2014), Table VII (cubic phase)
CsPbBr3e_mat = polaronmobility(T, 4.96, 16.5, 3.08E12, 0.12)
#CsPbBr3h_mat = polaronmobllity(T, 4.96, 16.5,  3.08E12, 0.15)

end

println("That's me!")


# Titania - polaron mobility 
println("   My Oberon! what visions have I seen!")
println("   Methought I was enamour'd of an ass.")

push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility 

##### load in library routines... #####
# Plot figures with Plots, which defaults to Pyplot backend
#using Plots
#default(grid=false) # No silly dotted grid lines
#default(size=(400,300)) # A good small size for two-column EPS output
#default(size=(800,600)) # Nice size for small-ish PNGs for slides

# Physical constants
const hbar = const ħ = 1.05457162825e-34;          # kg m2 / s 
const eV = const q = const ElectronVolt = 1.602176487e-19;                         # kg m2 / s2 
const me=MassElectron = 9.10938188e-31;                          # kg
const Boltzmann = const kB =  1.3806504e-23;                  # kg m2 / K s2 
const ε_0 = 8.854E-12 #Units: C2N−1m−2, permittivity of free space

#####
# Call simulation

# Rutile 
# Hendry et al. PRB 69, 081101 (2004)
# DOI: https://doi.org/10.1103/PhysRevB.69.081101
f=24E12 # p.3 LHS, 2nd para
meff=1.2 # ¬  - perpendicular # p.3 LHS, 4th para
#meff=0.6 # || - parallel

# F. A. GRANT
# Rev. Mod. Phys. 31, 646 – Published 1 July 1959
# DOI: https://doi.org/10.1103/RevModPhys.31.646
# Early paper --> some of this static dielectric will be free charges in the defective Titania
ϵ_static=173.0 # p. 650 LHS 
ϵ_optic=8.4

# Optical dielectric constants via:
# https://www.azom.com/article.aspx?ArticleID=1179 - optical properties
ϵ_optic=2.49^2 # Anatase: 6.2

anatase=polaronmobility(10:20:1000, ϵ_optic, ϵ_static, f, meff)
savepolaron("anatase",anatase)

ϵ_optic=2.903^2 # Rutile: 8.427

rutile=polaronmobility(10:20:1000, ϵ_optic, ϵ_static, f, meff)
savepolaron("rutile",rutile)

# Mobility of Electronic Charge Carriers in Titanium Dioxide
# T. Bak, M. K. Nowotny, L. R. Sheppard, and J. Nowotny*
# J. Phys. Chem. C 2008, 112, 12981–12987
# DOI: https://doi.org/10.1021/jp801028j
# Lots of nice comparison; single crystal mobility data; Arrhenius plots etc.

println("That's me!")


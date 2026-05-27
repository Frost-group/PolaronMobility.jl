# ASE-thresholds.jl
# Polaron parameters for halide perovskites: MAPI, CsPbI3, CsPbBr3
# Uses the modern PolaronMobility.jl v2.x API
#   --> patched to Brad's rewrite 2026-05-27
#   I think this now agrees with what I calculated with my old code!
#
# MAPI parameters from Frost 2017 PRB: https://arxiv.org/abs/1704.05404
# CsPbI3 parameters from Frost 2017 (DFT, PBEsol 9×9×9 k-mesh)
# CsPbBr3 parameters from Huang & Lambrecht PRB 90 195201 (2014) and
#   Yettapu et al. Nano Lett. 16, 4838 (2016) for effective masses
#
# Mott criterion: n_Mott = (0.26 / r_p)³
#   Edwards & Sienko, Phys. Rev. B 17, 2575 (1978)
#   where r_p is the polaron radius (Å → cm for density in cm⁻³)
#
# Jarvist Moore Frost, 2026

using PolaronMobility
using Unitful
using Printf

Temps = [300.0, 93.0]  # K

# --- Material parameters ---
# FrohlichMaterial(ϵ_optic, ϵ_static, m_eff, phonon_freq_THz)

# MAPbI3 — Frost 2017 PRB
MAPIe_mat = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
MAPIh_mat = FrohlichMaterial(4.5, 24.1, 0.15, 2.25)

# CsPbI3 — Frost 2017 DFT; ϵ_optic=6.1, ϵ_ionic=12.0 → ϵ_static=18.1
CsPbI3e_mat = FrohlichMaterial(6.1, 18.1, 0.12, 2.57)
CsPbI3h_mat = FrohlichMaterial(6.1, 18.1, 0.15, 2.57)

# CsPbBr3 — ϵ_optic=4.96, ϵ_static=16.5, phonon_freq=3.08 THz
# ϵ from Huang & Lambrecht PRB 90 195201 (2014), Table VII (cubic phase)
# m_eff from Yettapu et al. Nano Lett. 16 4838 (2016); me*=0.22, mh*=0.24
CsPbBr3e_mat = FrohlichMaterial(4.96, 16.5, 0.22, 3.08)
CsPbBr3h_mat = FrohlichMaterial(4.96, 16.5, 0.24, 3.08)

# Order: MAPI, CsPbI3, then CsPbBr3 last
materials = [
    ("MAPbI₃ electron", MAPIe_mat),
    ("MAPbI₃ hole", MAPIh_mat),
    ("CsPbI₃ electron", CsPbI3e_mat),
    ("CsPbI₃ hole", CsPbI3h_mat),
    ("CsPbBr₃ electron", CsPbBr3e_mat),
    ("CsPbBr₃ hole", CsPbBr3h_mat),
]

"""
    mott_density(r_Å)
Mott criterion carrier density (cm⁻³) from polaron radius r_Å (in Å).
n_Mott = (0.26 / r_p)³; Edwards & Sienko PRB 17 2575 (1978).
"""
mott_density(r_Å) = 1e24 * (2 * (2 * r_Å)^3)^-1  # Å → cm, then cube

for T in Temps
    println("="^72)
    println("  Halide Perovskite Polaron Parameters @ T = $(T) K")
    println("="^72)

    for (name, mat) in materials
        prob = material_to_problem(mat)
        res = solve(
            prob;
            temperatures = [T],
            frequencies = [0.0],
            options = OptimizerOptions(multistart = false),
        )
        units_res = material_units(res)
        sol = res.solutions[1]

        μ_val = ustrip(units_res.mobility[1])
        R_val = ustrip(units_res.radius[1])
        n_Mott = mott_density(R_val)

        println("\n--- $name ---")
        println("  α     = ", round.(mat.alpha, digits=3), " v = ", round(sol.v, digits=3), " w = ", round(sol.w, digits=3))
        println("  μ     = ", round(μ_val, digits=1), " cm²/Vs")
        println("  R_p   = ", round(R_val, digits=2), " Å")
        println("  n_Mott= ", @sprintf("%.2e", n_Mott), " cm⁻³")
    end
    println()
end

println("="^72)
println("  Done.")

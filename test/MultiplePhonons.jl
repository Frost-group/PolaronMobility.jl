# Check values of the new multimode polaron codes

@testset verbose = true "MultiplePhonons" begin

# Physical constants
ħ = 1.05457162825e-34;          # kg m2 / s
eV = 1.602176487e-19;           # kg m2 / s2
me = 9.10938188e-31;            # kg
kB =  1.3806504e-23;            # kg m2 / K s2
ε_0 = 8.854E-12                 # Units: C2N−1m−2, permittivity of free space
c = 3e8

MAPI= [
# 96.20813558773261 0.4996300522819191
# 93.13630357703363 1.7139631746083817
# 92.87834578121567 0.60108592692181
# 92.4847918585963 0.0058228799414729
# 92.26701437594754 0.100590086574602
# 89.43972834606603 0.006278895133832249
# 46.89209141511332 0.2460894564364346
# 46.420949316788 0.14174282581124137
# 44.0380222871706 0.1987196948553428
# 42.89702947649343 0.011159939465770681
# 42.67180170168193 0.02557751102757614
# 41.46971205834201 0.012555230726601503
# 37.08982543385215 0.00107488277468418
# 36.53555265689563 0.02126940080871224
# 30.20608114002676 0.009019481779712388
# 27.374810898415028 0.03994453721421388
# 26.363055017011728 0.05011922682554448
# 9.522966890022039 0.00075631870522737
4.016471586720514 0.08168931020200264
3.887605410774121 0.006311654262282101
3.5313112232401513 0.05353548710183397
2.755392921480459 0.021303020776321225
2.4380741812443247 0.23162784335484837
2.2490917637719408 0.2622203718355982
2.079632190634424 0.23382298607799906
2.0336707697261187 0.0623239656843172
1.5673011873879714 0.0367465760261409
1.0188379384951798 0.0126328938653956
1.0022960504442775 0.006817361620021601
0.9970130778462072 0.0103757951973341
0.9201781906386209 0.01095811116040592
0.800604081794174 0.0016830270365341532
0.5738689505255512 0.00646428491253749
#0.022939578929507105 8.355742795827834e-05   # Acoustic modes!
#0.04882611767873102 8.309858592685e-06
#0.07575149723846182 2.778248540373041e-05
]

volume = (6.29e-10)^3
ϵ_optic = 4.5
m_eff = 0.12
phonon_freq = MAPI[:, 1]
ir_activity = MAPI[:, 2]
ϵ_static = 24.1

# Ionic dielectric and decomposed alpha

ω = 2π .* phonon_freq
ϵ_ionic = [ϵ_ionic_mode(i, j, volume) for (i, j) in zip(phonon_freq, ir_activity)]

α = [multi_frohlichalpha(ϵ_optic, i, sum(ϵ_ionic), j, m_eff) for (i, j) in zip(ϵ_ionic, phonon_freq)]

@testset "Ionic dielectric and decomposed alphas" begin

    @test ϵ_ionic ≈ [0.2999680470756664, 0.0247387647244569, 0.2543132184018061, 0.16621617310133838, 2.3083204422506296, 3.0707979601813267, 3.2026782087486407, 0.892674135624958, 0.8861579096771846, 0.7209278375756829, 0.40199759805819046, 0.6183279038315278, 0.7666391525823296, 0.1555444994147378, 1.1627710200840813] rtol = 1e-3

    @test α ≈ [0.03401013445306177, 0.002850969846883158, 0.03075081562607006, 0.02275292381052607, 0.33591418423943553, 0.46526818717696034, 0.5046331098089347, 0.14223560721522646, 0.16083871312929882, 0.1622911897190622, 0.09123913334086006, 0.14070972715961402, 0.18159786148348117, 0.039500396977008016, 0.3487735470060312] rtol = 1e-3

    @test sum(α) ≈ 2.663366500992453 rtol = 1e-3

    println('\n', "Ionic dielectric and alphas:")
    println("-------------------------------")
    println("  ω (THz)    ϵ_ionic     αs")
    println("-------------------------------")
    display(hcat(ω, ϵ_ionic, α))
    println("\n-------------------------------")
end

# Variations

v_0, w_0 = var_params(α; v = 0.0, w = 0.0, ω = ω) # Athermal

β = [i .* ħ / (kB * 300) * 1e12 for i in ω]
v, w = var_params(α, β; v = 0.0, w = 0.0, ω = ω) # Thermal

@testset "Multiple mode variations" begin
    
    @test v_0 ≈ 3.292283619446986 rtol = 1e-3
    @test w_0 ≈ 2.679188425097246 rtol = 1e-3
    
    @test v ≈ 35.19211042393129 rtol = 1e-3
    @test w ≈ 32.454157668863225 rtol = 1e-3

    println("\nVariational parameters:")
    println("Athermal v = $(v_0[1]) | athermal w = $(w_0[1])")
    println("300K v = $(v[1]) | 300K w = $(w[1])")
end

# Energies

E = multi_F(v_0, w_0, α; ω = ω) * 1000 * ħ / eV * 1e12 # Enthalpy

F = multi_F(v, w, α, β; ω = ω) * 1000 * ħ / eV * 1e12 # Free energy

@testset "Multiple mode energies" begin
    
    @test E ≈ -19.50612170650821 rtol = 1e-3

    @test F ≈ -42.79764110613318 rtol = 1e-3

    println("\nPolaron energies")
    println("0K energy: $E meV")
    println("300K energy: $F meV")
end

# Mobility

μ = polaron_mobility(β, α, v, w; ω = ω) * eV / (1e12 * me * m_eff) * 100^2

@testset "Multiple mode mobility" begin
    
    @test μ ≈ 160.50844330430286 rtol = 1e-3

    println("\nMobility at 300K: $μ")
end

# Impedence 

Z_dc = polaron_complex_impedence(0, β, α, v, w; ω = ω) / sqrt(ε_0 * ϵ_optic) / c # DC limit

Z = polaron_complex_impedence(50, β, α, v, w; ω = ω) / sqrt(ε_0 * ϵ_optic) / c # AC limit

@testset "Multiple mode impedence" begin

    @test Z_dc ≈ 0.04822198080982958 + 0.0im rtol = 1e-3

    @test Z ≈ 0.015597562999157006 - 0.15169951660903833im rtol = 1e-3

    println("\nComplex Impedence:")
    println("DC limit: $Z_dc")
    println("At 50 THz: $Z")
end

# Conductivity

σ_dc = polaron_complex_conductivity(0, β, α, v, w; ω = ω) / sqrt(ε_0 * ϵ_optic) / c # DC limit

σ = polaron_complex_conductivity(50, β, α, v, w; ω = ω) / sqrt(ε_0 * ϵ_optic) / c # AC limit

@testset "Multiple mode conductivity" begin
    
    @test σ_dc ≈ 5.783096153975994e-6 - 0.0im rtol = 1e-3

    @test σ ≈ 1.8703663429355958e-7 + 1.8190897521649982e-6im rtol = 1e-3

    println("\nComplex Conductivity:")
    println("DC limit: $σ_dc")
    println("At 50 THz: $σ\n")
end

# Single Mode Tests

Hellwarth_B_freq = HellwarthBScheme(hcat(phonon_freq, ir_activity))

@testset "Hellwarth effective frequency" begin
    
    @test Hellwarth_B_freq ≈ 2.25 rtol = 1e-3
end 

@testset "Single effective mode MAPI" begin 

    singlemode_polaron = make_polaron(ϵ_optic, ϵ_static, Hellwarth_B_freq, m_eff, [0.0, 300.0], [0.0, 3.0]; volume = volume, ir_activity = nothing, rtol = 1e-4, verbose = true, threads = true)

    println('\n', singlemode_polaron)

    @test singlemode_polaron.α ≈ 2.393 rtol = 1e-3 
    @test singlemode_polaron.v ≈ [3.308644142915268; 19.847591395925644] rtol = 1e-3
    @test singlemode_polaron.w ≈ [2.6633969095604466; 16.948206590039813] rtol = 1e-3
    @test singlemode_polaron.F ≈ [-23.02903831886734, -35.46521250753788] rtol = 1e-3
    @test singlemode_polaron.Z ≈ [0.0 + 0.0im 0.11385135083260453 + 0.0im; 0.018788147630145906 - 0.0320853372512867im 0.11055356479444363 - 0.0023745543208784164im] rtol = 1e-3
    @test singlemode_polaron.σ ≈ [Inf + 0.0im 8.783382829337691 - 0.0im; 13.590340404155473 + 23.208815675232696im 9.04121795743805 + 0.1941942189449604im]
    @test singlemode_polaron.μ ≈ [Inf, 136.42796295778325] rtol = 1e-3
end

# Multimode Tests

@testset "Multiple mode (15) MAPI" begin

    multimode_polaron = make_polaron(ϵ_optic, ϵ_static, phonon_freq, m_eff, [0.0, 300.0], [0.0, 3.0]; volume = volume, ir_activity = ir_activity, rtol = 1e-4, verbose = true, threads = true)

    println('\n', multimode_polaron)

    @test multimode_polaron.α ≈ [0.03401013445306177, 0.002850969846883158, 0.03075081562607006, 0.02275292381052607, 0.33591418423943553, 0.46526818717696034, 0.5046331098089347, 0.14223560721522646, 0.16083871312929882, 0.1622911897190622, 0.09123913334086006, 0.14070972715961402, 0.18159786148348117, 0.039500396977008016, 0.3487735470060312] rtol = 1e-3
    @test multimode_polaron.v ≈ [3.292288128236545; 35.19197149158517] rtol = 1e-3
    @test multimode_polaron.w ≈ [2.679192692280829; 32.45402163343302] rtol = 1e-3
    @test multimode_polaron.F ≈ [-19.506121706507198, -42.79764112300015] rtol = 1e-3
    @test multimode_polaron.Z ≈ [0.0 + 0.0im 0.09684805201838559 + 0.0im; 0.013916089370749993 - 0.02698207327877506im 0.09373671239869472 - 0.00586253350488456im] rtol = 1e-3
    @test multimode_polaron.σ ≈ [Inf + 0.0im 10.325452904413197 - 0.0im; 15.098451823667599 + 29.274570078468706im 10.626611965355764 + 0.6646154649133272im] rtol = 1e-3
    @test multimode_polaron.μ ≈ [Inf, 160.50598091483337] rtol = 1e-3
end

print('\n')

end
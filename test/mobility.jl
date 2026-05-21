@testset "Mobility" begin
    material = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
    problem = material_to_problem(material)
    result = solve(problem; temperatures = [0, 300])

    @test result.mobilities[1].mobility == Inf
    @test result.mobilities[2].mobility ≈ 0.1315855458846205 rtol = 1e-3
    @test result.mobilities[2].hellwarth ≈ 0.13161105783014404 rtol = 1e-3
    @test mobility_cm2_per_v_s(result.mobilities[2].mobility, material.effective_frequency, material.band_mass) ≈ 136.42 rtol = 0.02

    @test result.mobilities[2].fhip_low_temperature > result.mobilities[2].mobility
    @test result.mobilities[2].relaxation_time > 0
end

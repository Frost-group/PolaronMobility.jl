@testset "Response" begin
    material = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
    problem = material_to_problem(material)
    result = solve(problem; temperatures = [0, 300], frequencies = [0, 3])

    zero_dc = result.responses[1, 1]
    @test isinf(real(zero_dc.memory_function))
    @test isinf(imag(zero_dc.impedance))
    @test zero_dc.conductivity == 0 + 0im

    ac = result.responses[2, 2]
    @test ac.memory_function ≈ -1.1750219539716098 + 7.379516901672299im rtol = 1e-3
    @test ac.impedance ≈ 7.379516901672299 - 0.15831137936172346im rtol = 1e-3
    @test ac.conductivity ≈ 0.13544788933177754 + 0.0029057379334531375im rtol = 1e-3

    @test frohlich_complex_impedance(
        ac.frequency,
        result.solutions[2].v,
        result.solutions[2].w,
        problem.model.alpha,
        problem.model.phonon_frequencies,
        result.solutions[2].beta,
    ) ≈ ac.impedance rtol = 1e-3
end

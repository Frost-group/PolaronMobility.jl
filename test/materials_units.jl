@testset "Materials And Units" begin
    material = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
    problem = material_to_problem(material)
    result = solve(problem; temperatures = [0, 300], frequencies = [3])
    units = material_units(result)

    @test material.alpha[1] ≈ 2.3939410167951287 rtol = 1e-8
    @test result.temperatures ≈ [0.0, 2.7782158848126093] rtol = 1e-8
    @test result.frequencies ≈ [1.3333333333333333]
    @test units.mobility[2] ≈ 136.42248372298718u"cm^2/V/s" rtol = 1e-6
    @test units.frequency[1] ≈ 3.0u"THz"

    before = result.mobilities[2].mobility
    material_units(result)
    @test result.mobilities[2].mobility == before

    multimode = FrohlichMaterial(4.5, 24.1, 0.12, [4.0, 2.0], [0.1, 0.2], (6.3e-10)^3)
    multimode_problem = material_to_problem(multimode)
    @test length(multimode_problem.model.alpha) == 2
    @test length(multimode_problem.model.phonon_frequencies) == 2
    @test multimode_problem.model.phonon_frequencies ≈ multimode.phonon_frequencies ./ multimode.effective_frequency

    multimode_result = solve(
        multimode_problem;
        temperatures = [0, 300],
        frequencies = [0, 3],
        options = OptimizerOptions(multistart = false),
    )
    @test multimode_result.temperatures[1] == 0.0
    @test multimode_result.temperatures[2] > 0
    @test multimode_result.frequencies[2] > 0
    @test size(multimode_result.responses) == (2, 2)
    @test all(isfinite, [multimode_result.solutions[1].energy.total, multimode_result.solutions[2].energy.total])

    material_trials = (
        material_to_problem(multimode; trial = :feynman),
        material_to_problem(multimode; trial = :multi_gaussian, modes = 2),
        material_to_problem(multimode; trial = :profile_gaussian, matsubara_terms = 64),
        material_to_problem(multimode; trial = :nonlocal_gaussian),
    )
    @test material_trials[1] isa VariationalProblem{FrohlichModel,GaussianFeynmanTrial}
    @test material_trials[2] isa VariationalProblem{FrohlichModel,MultiGaussianTrial}
    @test material_trials[3] isa VariationalProblem{FrohlichModel,ProfileGaussianTrial}
    @test material_trials[4] isa VariationalProblem{FrohlichModel,NonlocalGaussianTrial}
    for problem_variant in material_trials
        small = solve(
            problem_variant;
            temperatures = [0],
            frequencies = [0],
            options = OptimizerOptions(multistart = false, adaptive_bounds = false, quadrature_rtol = 1e-3),
        )
        @test small.temperatures == [0.0]
        @test size(small.responses) == (1, 1)
        @test isfinite(small.zero_temperature.energy.total)
    end
    @test_throws ArgumentError material_to_problem(multimode; trial = :unknown_trial)

    rubrene = rubrene_holstein_material(lattice_constant_angstrom = 7.2)
    rubrene_problem = material_to_problem(rubrene)
    rubrene_result = solve(
        rubrene_problem;
        temperatures = [300],
        frequencies = [0, 10],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    rubrene_units = material_units(rubrene_result)

    @test rubrene_problem isa VariationalProblem{HolsteinModel,PoissonTrial}
    @test rubrene.hopping ≈ 134.0 / wavenumber_meV(1208.9)
    @test rubrene.coupling ≈ sqrt(106.8 / wavenumber_meV(1208.9))
    @test lambda_holstein(rubrene) ≈ 106.8 / (2 * 134.0)
    @test rubrene_result.temperatures[1] ≈ PolaronMobility.reduced_temperature(300, rubrene.phonon_frequency_THz)
    @test rubrene_result.frequencies[2] ≈ PolaronMobility.reduced_frequency(10, rubrene.phonon_frequency_THz)
    @test size(rubrene_result.responses) == (2, 1)
    @test isfinite(rubrene_result.solutions[1].free_energy)
    @test isfinite(ustrip(rubrene_units.mobility_einstein[1]))
    @test rubrene_units.frequency[2] ≈ 10u"THz"

    @test_throws ArgumentError material_to_problem(rubrene; trial = :feynman)
    no_lattice = rubrene_holstein_material()
    @test_throws ArgumentError material_units(solve(material_to_problem(no_lattice); temperatures = [300]))

    peierls_material = rubrene_peierls_material(lattice_constant_angstrom = 7.2)
    peierls_problem = material_to_problem(peierls_material)
    peierls_result = solve(
        peierls_problem;
        temperatures = [300],
        frequencies = [0, 10],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    peierls_units = material_units(peierls_result)

    @test peierls_material isa PeierlsMaterial
    @test peierls_problem isa VariationalProblem{PeierlsModel,PoissonTrial}
    @test peierls_material.hopping ≈ 134.0 / wavenumber_meV(117.9)
    @test peierls_material.coupling ≈ sqrt(21.9 / wavenumber_meV(117.9))
    @test lambda_peierls(peierls_material) ≈ 21.9 / (2 * 134.0)
    @test peierls_result.temperatures[1] ≈ PolaronMobility.reduced_temperature(300, peierls_material.phonon_frequency_THz)
    @test peierls_result.frequencies[2] ≈ PolaronMobility.reduced_frequency(10, peierls_material.phonon_frequency_THz)
    @test size(peierls_result.responses) == (2, 1)
    @test isfinite(peierls_result.solutions[1].free_energy)
    @test isfinite(ustrip(peierls_units.mobility_einstein[1]))
    @test peierls_units.frequency[2] ≈ 10u"THz"
    @test_throws ArgumentError material_to_problem(peierls_material; trial = :feynman)

    rubrene_composite = rubrene_holstein_peierls_problem(lattice_constant_angstrom = 7.2)
    composite_result = solve(
        rubrene_composite;
        temperatures = [300],
        frequencies = [0, 10],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test rubrene_composite.model isa CompositePolaronModel
    @test length(rubrene_composite.model.models) == 2
    @test composite_result isa PolaronResult
    @test composite_result.temperatures[1] ≈ PolaronMobility.reduced_temperature(300, rubrene_composite.model.effective_frequency)
    @test size(composite_result.responses) == (2, 1)

    composite_units = material_units(composite_result)
    @test isfinite(ustrip(composite_units.mobility[1]))
    @test composite_result.mobilities[1].mobility <= composite_result.mobilities[1].mobility_einstein
    @test composite_result.mobilities[1].mobility ≈
          composite_result.mobilities[1].mobility_einstein * composite_result.mobilities[1].mobility_factor rtol = 1e-12
end

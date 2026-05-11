@testset "Peierls And Composite Lattice Models" begin
    zero_problem = peierls_poisson_problem(coupling = 0.0, hopping = 1.4)
    zero_result = solve(
        zero_problem;
        temperatures = [0.5],
        frequencies = [0.0, 0.25],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test zero_result isa PolaronResult
    @test zero_result.solutions[1].rate ≈ 1.4 rtol = 1e-5
    @test zero_result.solutions[1].free_energy ≈ -2 * 1.4 rtol = 1e-5
    @test size(zero_result.responses) == (2, 1)
    @test response_table(zero_result)[1].model == :peierls

    weak = peierls_poisson_problem(coupling = 0.25, hopping = 1.0)
    strong = peierls_poisson_problem(coupling = 0.5, hopping = 1.0)
    params = [1.0]
    weak_shift = interaction_free_energy(weak.model, weak.trial, params, Inf)
    strong_shift = interaction_free_energy(strong.model, strong.trial, params, Inf)
    @test weak_shift < 0
    @test strong_shift / weak_shift ≈ 4.0 rtol = 1e-6

    weak_result = solve(weak; temperatures = [1.0], frequencies = [0.0], options = OptimizerOptions(multistart = false, adaptive_bounds = false))
    strong_result = solve(strong; temperatures = [1.0], frequencies = [0.0], options = OptimizerOptions(multistart = false, adaptive_bounds = false))
    @test strong_result.solutions[1].free_energy < weak_result.solutions[1].free_energy
    @test strong_result.solutions[1].rate <= weak_result.solutions[1].rate
    @test strong_result.mobilities[1].component_mobilities.peierls != 0

    holstein = holstein_poisson_problem(coupling = 0.4, hopping = 1.0)
    peierls = peierls_poisson_problem(coupling = 0.3, hopping = 1.0)
    composite_model = combine_models(holstein.model, peierls.model)
    composite_problem = VariationalProblem(composite_model, PoissonTrial(bare_hopping = 1.0))
    beta = 2.0
    @test interaction_free_energy(composite_model, composite_problem.trial, params, beta) ≈
          interaction_free_energy(holstein.model, holstein.trial, params, beta) +
          interaction_free_energy(peierls.model, peierls.trial, params, beta)

    zero_peierls = peierls_poisson_problem(coupling = 0.0, hopping = 1.0)
    holstein_plus_zero = VariationalProblem(combine_models(holstein.model, zero_peierls.model), PoissonTrial(bare_hopping = 1.0))
    pure_holstein = solve(holstein; temperatures = [1.0], frequencies = [0.0], options = OptimizerOptions(multistart = false, adaptive_bounds = false))
    composite_holstein = solve(holstein_plus_zero; temperatures = [1.0], frequencies = [0.0], options = OptimizerOptions(multistart = false, adaptive_bounds = false))
    @test composite_holstein.solutions[1].free_energy ≈ pure_holstein.solutions[1].free_energy rtol = 1e-5
    @test composite_holstein.mobilities[1].mobility ≈ pure_holstein.mobilities[1].mobility rtol = 1e-5
    expected_peierls_component = lattice_mobility(
        zero_peierls.model,
        composite_holstein.solutions[1].rate,
        composite_holstein.solutions[1].beta,
        0.0,
    )
    @test composite_holstein.mobilities[1].component_mobilities.peierls ≈
          expected_peierls_component rtol = 1e-10

    zero_holstein = holstein_poisson_problem(coupling = 0.0, hopping = 1.0)
    peierls_plus_zero = VariationalProblem(combine_models(zero_holstein.model, peierls.model), PoissonTrial(bare_hopping = 1.0))
    pure_peierls = solve(peierls; temperatures = [1.0], frequencies = [0.0], options = OptimizerOptions(multistart = false, adaptive_bounds = false))
    composite_peierls = solve(peierls_plus_zero; temperatures = [1.0], frequencies = [0.0], options = OptimizerOptions(multistart = false, adaptive_bounds = false))
    @test composite_peierls.solutions[1].free_energy ≈ pure_peierls.solutions[1].free_energy rtol = 1e-5

    @test_throws ArgumentError PeierlsMaterial(hopping_meV = 1.0, phonon_frequency_cm1 = -1.0, peierls_energy_meV = 1.0)
    @test_throws ArgumentError PeierlsMaterial(hopping_meV = 1.0, phonon_frequency_cm1 = 100.0, peierls_energy_meV = -1.0)
    @test_throws ArgumentError combine_models(frohlich_feynman_problem(coupling = 1.0).model, holstein.model)
end

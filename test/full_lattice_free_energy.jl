@testset "Full Periodic Lattice Free Energy" begin
    options = OptimizerOptions(multistart = false, adaptive_bounds = false, quadrature_rtol = 1e-4)

    for builder in (holstein_poisson_problem, peierls_poisson_problem)
        problem = builder(coupling = 0.0, hopping = 1.3, dimension = 1)
        result = solve_variational(problem, 10.0; options = options)
        @test result.parameters[1] ≈ 1.3 rtol = 1e-5
        @test result.free_energy ≈ -2.6 rtol = 1e-5
    end

    holstein = holstein_poisson_problem(coupling = 0.7, hopping = 1.0, phonon_frequency = 1.0)
    holstein_zero = interaction_free_energy(holstein.model, holstein.trial, [0.9], Inf; rtol = 1e-5)
    holstein_cold = interaction_free_energy(holstein.model, holstein.trial, [0.9], 80.0; rtol = 1e-5)
    @test holstein_zero < 0
    @test holstein_cold ≈ holstein_zero rtol = 2e-2

    peierls_zero = peierls_poisson_problem(coupling = 0.0, hopping = 1.0)
    peierls_weak = peierls_poisson_problem(coupling = 0.25, hopping = 1.0)
    peierls_strong = peierls_poisson_problem(coupling = 0.5, hopping = 1.0)
    params = [1.0]
    @test interaction_free_energy(peierls_zero.model, peierls_zero.trial, params, Inf) == 0.0
    weak_shift = interaction_free_energy(peierls_weak.model, peierls_weak.trial, params, Inf)
    strong_shift = interaction_free_energy(peierls_strong.model, peierls_strong.trial, params, Inf)
    @test weak_shift < 0
    @test strong_shift / weak_shift ≈ 4.0 rtol = 1e-6

    old_closure = 2 * peierls_weak.model.dimension * params[1]^2 * return_probability(peierls_weak.trial, 0.5, params)
    @test PolaronMobility.peierls_bond_correlation(peierls_weak.model, peierls_weak.trial, params, 0.5) != old_closure

    composite = combine_models(holstein.model, peierls_weak.model)
    trial = PoissonTrial(bare_hopping = 1.0)
    @test interaction_free_energy(composite, trial, params, 3.0) ≈
          interaction_free_energy(holstein.model, trial, params, 3.0) +
          interaction_free_energy(peierls_weak.model, trial, params, 3.0)
end

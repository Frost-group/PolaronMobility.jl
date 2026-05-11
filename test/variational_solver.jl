@testset "Variational Solver" begin
    schultz = [
        3.00 3.44 2.55 -3.1333
        5.00 4.02 2.13 -5.4401
        7.00 5.81 1.60 -8.1127
        9.00 9.85 1.28 -11.486
        11.0 15.5 1.15 -15.710
    ]

    for row in eachrow(schultz)
        alpha, expected_v, expected_w, expected_energy = row
        solution = solve_frohlich(alpha).zero_temperature
        @test solution.v ≈ expected_v atol = 0.1
        @test solution.w ≈ expected_w atol = 0.1
        @test solution.energy.total ≈ expected_energy atol = 0.01
    end

    result = solve_frohlich(2.393156008589176; temperatures = [0.0, 2.7782158848126093])
    @test [s.v for s in result.solutions] ≈ [3.3081654923629096, 19.856578327838974] rtol = 1e-3
    @test [s.w for s in result.solutions] ≈ [2.663611744517749, 16.958531394655713] rtol = 1e-3

    problem = frohlich_feynman_problem(coupling = 3.0)
    variational = solve_variational(problem, Inf)
    params = variational.parameters
    @test objective(problem, params, Inf) ≈
          free_energy(problem.trial, params, Inf) +
          entropy_cost(problem.trial, params, Inf) +
          interaction_free_energy(problem.model, problem.trial, params, Inf)
    @test variational.diagnostics.optimizer_success
    @test variational.diagnostics.optimizer_hit_upper_bound_1 == false
    @test feynman_v(variational) ≈ 3.4212966497238306 rtol = 1e-3
    @test feynman_w(variational) ≈ 2.5603045739035077 rtol = 1e-3

    tight = OptimizerOptions(upper = [2.55, 0.8], multistart = false, adaptive_bounds = true)
    expanded = solve_variational(problem, Inf; options = tight)
    @test expanded.diagnostics.optimizer_bound_expansions >= 1
    @test expanded.diagnostics.optimizer_bound_upper_1 > 2.55

    warm_started = solve_variational(problem, Inf; initial_parameters_override = [5.0, 1.0], use_multistart = false)
    @test warm_started.free_energy ≈ variational.free_energy atol = 1e-6

    single_mode = frohlich_multi_gaussian_problem(coupling = 3.0, modes = 1)
    single_solution = solve(single_mode; temperatures = 0.0, frequencies = 0.0, options = OptimizerOptions(multistart = false)).zero_temperature
    @test single_solution.v ≈ feynman_v(variational) rtol = 1e-5
    @test single_solution.w ≈ feynman_w(variational) rtol = 1e-5
    @test single_solution.energy.total ≈ variational.free_energy rtol = 1e-5

    two_mode = frohlich_multi_gaussian_problem(coupling = 3.0, modes = 2)
    two_mode_result = solve_variational(two_mode, Inf; options = OptimizerOptions(multistart = false))
    @test length(two_mode_result.parameters) == 4
    @test two_mode_result.parameter_names == [:w1, :delta1, :w2, :delta2]
    @test two_mode_result.free_energy <= single_solution.energy.total + 1e-5
    @test length(multi_gaussian_v(two_mode_result)) == 2
    @test length(multi_gaussian_w(two_mode_result)) == 2

    @test_throws ArgumentError solve_variational(two_mode, Inf; initial_parameters_override = [1.0, 0.1], use_multistart = false)
    @test_throws DomainError objective(two_mode, [1.0, 0.1, 1.0, 0.2], Inf)

    two_mode_multistart = solve(
        two_mode;
        temperatures = [0.0, 0.5],
        frequencies = [0.0, 1.0],
        options = OptimizerOptions(multistart = true),
    )
    @test size(two_mode_multistart.responses) == (2, 2)
    @test isfinite(two_mode_multistart.zero_temperature.energy.total)

    nonlocal = frohlich_nonlocal_gaussian_problem(coupling = 1.0, basis_frequencies = [1.0, 2.0])
    @test nonlocal isa VariationalProblem{FrohlichModel,NonlocalGaussianTrial}
    @test parameter_names(nonlocal.trial) == [:a1, :a2]
    nonlocal_result = solve(
        nonlocal;
        temperatures = [0.0, 0.5],
        frequencies = [0.0, 1.0],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false, upper = [0.4, 0.4]),
    )
    @test size(nonlocal_result.responses) == (2, 2)
    @test hasproperty(nonlocal_result.zero_temperature.parameters, :a1)
    @test isfinite(nonlocal_result.zero_temperature.energy.total)
    @test isnan(nonlocal_result.mobilities[2].fhip_low_temperature)
    @test mean_square_displacement(nonlocal.trial, [0.1, 0.2], 0.5, Inf) > 0

    profile = frohlich_profile_gaussian_problem(
        coupling = 1.0,
        basis_frequencies = [1.0, 2.0],
        matsubara_terms = 128,
    )
    @test profile isa VariationalProblem{FrohlichModel,ProfileGaussianTrial}
    @test parameter_names(profile.trial) == [:a1, :a2]
    @test profile_function(profile.trial, [0.2, 0.1], 0.5) > 0
    @test entropy_cost(profile.trial, [0.0, 0.0], Inf) ≈ 0.0 atol = 1e-10
    @test mean_square_displacement(profile.trial, [0.0, 0.0], 0.5, Inf) ≈ 0.5 rtol = 1e-4
    @test mean_square_displacement(profile.trial, [0.0, 0.0], 0.25, 1.0) ≈ 0.25 * 0.75 rtol = 1e-2

    feynman_v_value = feynman_v(variational)
    feynman_w_value = feynman_w(variational)
    profile_amplitude = (feynman_v_value^2 - feynman_w_value^2) / feynman_w_value^2
    profile_embedding = frohlich_profile_gaussian_problem(
        coupling = 3.0,
        basis_frequencies = [feynman_w_value],
        initial_amplitudes = [profile_amplitude],
    )
    @test objective(profile_embedding, [profile_amplitude], Inf) ≈ variational.free_energy rtol = 1e-8
    @test entropy_cost(profile_embedding.trial, [profile_amplitude], Inf) ≈ entropy_cost(problem.trial, params, Inf) rtol = 1e-8
    for tau in (1e-4, 1e-2, 0.1, 0.5, 2.0)
        @test mean_square_displacement(profile_embedding.trial, [profile_amplitude], tau, Inf) ≈
              mean_square_displacement(problem.trial, params, tau, Inf) rtol = 1e-8
    end

    profile_result = solve(
        profile;
        temperatures = [0.0, 0.5],
        frequencies = [0.0, 1.0],
        options = OptimizerOptions(
            multistart = false,
            adaptive_bounds = false,
            upper = [0.4, 0.4],
            quadrature_rtol = 1e-3,
        ),
    )
    @test size(profile_result.responses) == (2, 2)
    @test hasproperty(profile_result.zero_temperature.parameters, :a1)
    @test isfinite(profile_result.zero_temperature.energy.total)
    @test isnan(profile_result.mobilities[2].hellwarth)
end

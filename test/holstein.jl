@testset "Holstein Poisson" begin
    problem = holstein_poisson_problem(coupling = 0.0, hopping = 1.0, phonon_frequency = 1.0, dimension = 1)
    @test problem isa VariationalProblem{HolsteinModel,PoissonTrial}

    variational = solve_variational(problem, 20.0)
    @test variational.parameters[1] ≈ 1.0 atol = 1e-5
    @test variational.free_energy ≈ -2.0 atol = 1e-5
    @test objective(problem, variational.parameters, 20.0) ≈
          free_energy(problem.trial, variational.parameters, 20.0) +
          entropy_cost(problem.trial, variational.parameters, 20.0) +
          interaction_free_energy(problem.model, problem.trial, variational.parameters, 20.0)

    weak = solve_variational(holstein_poisson_problem(coupling = 0.5), 20.0)
    strong = solve_variational(holstein_poisson_problem(coupling = 2.0), 20.0)
    @test strong.parameters[1] < weak.parameters[1]
    @test strong.free_energy < weak.free_energy

    result = solve_holstein(coupling = 1.5, temperatures = [0.05, 0.5], frequencies = [0.0, 1.0])
    @test result isa PolaronResult
    @test length(result.solutions) == 2
    @test size(result.responses) == (2, 2)
    @test result.mobilities[1].mobility > 0
    @test result.mobilities[1].mobility_einstein > 0
    @test 0 < result.mobilities[1].mobility_factor <= 1
    @test real(result.responses[1, 1].mobility) ≈ result.mobilities[1].mobility rtol = 1e-5
    @test real(result.responses[1, 1].mobility_factor) ≈ result.mobilities[1].mobility_factor rtol = 1e-5
    @test imag(result.responses[1, 1].mobility) ≈ 0 atol = 1e-8
    @test result.responses[1, 1] isa LatticeResponseResult

    rates = range(1e-4, 10.0; length = 300)
    scan_problem = holstein_poisson_problem(coupling = 1.5, phonon_frequency = 0.75)
    scan_result = solve_variational(scan_problem, 25.0)
    dense_min = minimum(objective(scan_problem, [rate], 25.0) for rate in rates)
    @test scan_result.free_energy <= dense_min + 5e-4
end

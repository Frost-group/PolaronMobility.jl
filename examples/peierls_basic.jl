using PolaronMobility

println("== Peierls Poisson trial ==")

problem = peierls_poisson_problem(
    hopping = 1.0,
    phonon_frequency = 0.8,
    coupling = 0.35,
    dimension = 1,
)

result = solve(
    problem;
    temperatures = [0.5, 1.0, 2.0],
    frequencies = [0.0, 0.25, 0.5],
    options = OptimizerOptions(multistart = false, adaptive_bounds = false),
)

println("optimized rates = ", [solution.rate for solution in result.solutions])
println("free energies = ", [solution.free_energy for solution in result.solutions])
println("response rows = ", length(response_table(result)))

println("\n== Holstein + Peierls composition ==")

holstein = holstein_poisson_problem(hopping = 1.0, phonon_frequency = 1.0, coupling = 0.4)
peierls = peierls_poisson_problem(hopping = 1.0, phonon_frequency = 0.8, coupling = 0.35)
combined = VariationalProblem(
    combine_models(holstein.model, peierls.model),
    PoissonTrial(bare_hopping = 1.0),
)
combined_result = solve(
    combined;
    temperatures = [1.0],
    frequencies = [0.0, 0.5],
    options = OptimizerOptions(multistart = false, adaptive_bounds = false),
)

println("combined optimized rate = ", combined_result.solutions[1].rate)
println("combined interaction = ", combined_result.solutions[1].interaction_energy)

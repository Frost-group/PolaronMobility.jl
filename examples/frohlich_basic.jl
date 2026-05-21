using PolaronMobility

println("== Fröhlich one-mode Feynman trial ==")

problem = frohlich_feynman_problem(coupling = 3.0, phonon_frequency = 1.0)
result = solve(
    problem;
    temperatures = [0.0, 0.5, 1.0],
    frequencies = [0.0, 1.0],
    options = OptimizerOptions(multistart = false),
)

println("zero-temperature energy = ", result.zero_temperature.energy.total)
println("zero-temperature v,w = ", (result.zero_temperature.v, result.zero_temperature.w))
println("T=0.5 mobility = ", result.mobilities[2].mobility)
println("Ω=1 conductivity at T=0.5 = ", result.responses[2, 2].conductivity)

variational = solve_variational(problem, Inf; options = OptimizerOptions(multistart = false))
println("single variational solve = ", (; v = feynman_v(variational), w = feynman_w(variational), energy = variational.free_energy))

params = variational.parameters
parts = (
    free = free_energy(problem.trial, params, Inf),
    entropy = entropy_cost(problem.trial, params, Inf),
    interaction = interaction_free_energy(problem.model, problem.trial, params, Inf),
)
println("objective decomposition = ", parts)

println("direct kernel energy = ", frohlich_energy(result.zero_temperature.v, result.zero_temperature.w, [3.0], [1.0]).total)
println("direct memory function Ω=1 = ", frohlich_memory_function(1.0, result.zero_temperature.v, result.zero_temperature.w, [3.0], [1.0]))
println("direct impedance Ω=1 = ", frohlich_complex_impedance(1.0, result.zero_temperature.v, result.zero_temperature.w, [3.0], [1.0]))
println("direct conductivity Ω=1 = ", frohlich_complex_conductivity(1.0, result.zero_temperature.v, result.zero_temperature.w, [3.0], [1.0]))

println("solution rows = ", length(solution_table(result)))
println("mobility rows = ", length(mobility_table(result)))
println("response rows = ", length(response_table(result)))

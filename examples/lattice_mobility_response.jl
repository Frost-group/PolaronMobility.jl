using PolaronMobility

println("== General-d lattice CTMC mobility response ==")

options = OptimizerOptions(multistart = false, adaptive_bounds = false)

holstein = holstein_poisson_problem(coupling = 0.5, hopping = 1.0, dimension = 1)
peierls = peierls_poisson_problem(coupling = 0.35, hopping = 1.0, dimension = 1)

for problem in (holstein, peierls)
    local lattice_result = solve(problem; temperatures = [1.0], frequencies = [0.0, 0.5, 1.0], options = options)
    println(typeof(problem.model), " rate = ", lattice_result.solutions[1].rate)
    println(typeof(problem.model), " DC mobility factor = ", lattice_result.mobilities[1].mobility_factor)
    println(typeof(problem.model), " mobility(Ω=1) = ", lattice_result.responses[3, 1].mobility)
end

composite = VariationalProblem(
    combine_models(holstein.model, peierls.model),
    PoissonTrial(dimension = 1, bare_hopping = 1.0),
)
composite_result = solve(composite; temperatures = [1.0], frequencies = [0.0, 0.5, 1.0], options = options)

println("composite component mobilities at Ω=0 = ", composite_result.responses[1, 1].component_mobilities)
println("composite conductivity at Ω=1 = ", composite_result.responses[3, 1].conductivity)

guide_rows = holstein_transport_sweep(
    coupling = 0.5,
    dimension = 2,
    temperatures = [0.5, 1.0],
    frequencies = [0.0],
    broadening = 0.02,
    options = options,
)
println("guide-style frozen-kappa rows = ", length(guide_rows))
println("first guide row sideband count = ", first(guide_rows).sideband_count)

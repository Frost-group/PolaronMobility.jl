using PolaronMobility

println("== Full-periodic lattice free energies ==")

options = OptimizerOptions(multistart = false, adaptive_bounds = false, quadrature_rtol = 1e-3)

holstein = holstein_poisson_problem(coupling = 0.7, hopping = 1.0, phonon_frequency = 1.0)
peierls = peierls_poisson_problem(coupling = 0.4, hopping = 1.0, phonon_frequency = 1.0)

rate = [1.0]
println("Holstein T=0 interaction = ", interaction_free_energy(holstein.model, holstein.trial, rate, Inf; rtol = 1e-3))
println("Holstein beta=10 interaction = ", interaction_free_energy(holstein.model, holstein.trial, rate, 10.0; rtol = 1e-3))
println("Peierls T=0 interaction = ", interaction_free_energy(peierls.model, peierls.trial, rate, Inf; rtol = 1e-3))
println("Peierls beta=10 interaction = ", interaction_free_energy(peierls.model, peierls.trial, rate, 10.0; rtol = 1e-3))

println("site bridge C(τ=0.5,β=4) = ", site_return_bridge(1.0, 1, 0.5, 4.0))
println("bond-order bridge C(τ=0.5,β=4) = ", bond_order_bridge(1.0, 1, 0.5, 4.0))
println("bond-current bridge C(τ=0.5,β=4) = ", bond_current_bridge(1.0, 1, 0.5, 4.0))

combined = VariationalProblem(combine_models(holstein.model, peierls.model), PoissonTrial(bare_hopping = 1.0))
result = solve(combined; temperatures = [1.0], frequencies = [0.0], options = options)

println("combined optimized rate = ", result.solutions[1].rate)
println("combined full-periodic free energy = ", result.solutions[1].free_energy)

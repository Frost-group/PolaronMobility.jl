using PolaronMobility

println("== Coded model limits ==")

options = OptimizerOptions(multistart = false, adaptive_bounds = false, quadrature_rtol = 1e-3)

println("\nFröhlich weak and stronger coupling")
weak_frohlich = solve_frohlich(0.1; temperatures = [0.0], frequencies = [0.0], options = options)
stronger_frohlich = solve_frohlich(3.0; temperatures = [0.0], frequencies = [0.0], options = options)
println("weak E + alpha = ", weak_frohlich.zero_temperature.energy.total + 0.1)
println("stronger v/w = ", stronger_frohlich.zero_temperature.v / stronger_frohlich.zero_temperature.w)

println("\nHolstein zero, weak, and high-temperature trends")
bare_holstein = solve_holstein(coupling = 0.0, hopping = 1.0; temperatures = [1.0], frequencies = [0.0], options = options)
weak_holstein = solve_holstein(coupling = 0.3, hopping = 1.0; temperatures = [1.0], frequencies = [0.0], options = options)
hot_holstein = solve_holstein(coupling = 0.3, hopping = 1.0; temperatures = [5.0], frequencies = [0.0], options = options)
println("zero-coupling rate = ", bare_holstein.solutions[1].rate)
println("weak interaction energy = ", weak_holstein.solutions[1].interaction_energy)
println("high-temperature mobility factor = ", hot_holstein.mobilities[1].mobility_factor)

println("\nPeierls zero and weak-coupling scaling")
bare_peierls = solve_peierls(coupling = 0.0, hopping = 1.0; temperatures = [1.0], frequencies = [0.0], options = options)
weak_peierls = peierls_poisson_problem(coupling = 0.25, hopping = 1.0)
stronger_peierls = peierls_poisson_problem(coupling = 0.5, hopping = 1.0)
weak_shift = interaction_free_energy(weak_peierls.model, weak_peierls.trial, [1.0], Inf; rtol = 1e-3)
stronger_shift = interaction_free_energy(stronger_peierls.model, stronger_peierls.trial, [1.0], Inf; rtol = 1e-3)
println("zero-coupling rate = ", bare_peierls.solutions[1].rate)
println("Peierls g^2 scaling ratio = ", stronger_shift / weak_shift)

println("\nAdiabaticity markers")
adiabatic = holstein_poisson_problem(coupling = 0.5, hopping = 2.0, phonon_frequency = 0.5)
antiadiabatic = holstein_poisson_problem(coupling = 0.5, hopping = 0.5, phonon_frequency = 2.0)
println("adiabatic J/omega = ", adiabatic.model.hopping / adiabatic.model.phonon_frequency)
println("antiadiabatic J/omega = ", antiadiabatic.model.hopping / antiadiabatic.model.phonon_frequency)

using PolaronMobility

println("== Holstein and Peierls mobility-factor comparison ==")

temperature = 1.0
beta = inv(temperature)
rate = 1.0
frequencies = [0.0, 0.5, 1.0, 2.0]

holstein = holstein_poisson_problem(coupling = 0.5, hopping = 1.0)
peierls = peierls_poisson_problem(coupling = 0.5, hopping = 1.0)
composite = combine_models(holstein.model, peierls.model)

for frequency in frequencies
    holstein_factor = lattice_mobility_factor(holstein.model, rate, beta, frequency)
    peierls_factor = lattice_mobility_factor(peierls.model, rate, beta, frequency)
    composite_factor = lattice_mobility_factor(composite, rate, beta, frequency)
    println("Ω=", frequency, " FH=", holstein_factor, " FP=", peierls_factor, " Ftotal=", composite_factor)
end

println("zero-coupling lattice mobility = ", lattice_mobility(holstein_poisson_problem(coupling = 0.0).model, rate, beta, 0.0))

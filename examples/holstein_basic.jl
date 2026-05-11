using PolaronMobility

println("== Holstein Poisson CTMC trial ==")

problem = holstein_poisson_problem(
    hopping = 1.0,
    phonon_frequency = 1.0,
    coupling = 1.5,
    dimension = 1,
)

result = solve(
    problem;
    temperatures = [0.25, 0.5, 1.0],
    frequencies = [0.0, 1.0],
    options = OptimizerOptions(multistart = false),
)

println("optimized rate at T=0.25 = ", result.solutions[1].rate)
println("Einstein mobility at T=0.25 = ", result.mobilities[1].mobility_einstein)
println("lattice-FHIP mobility at T=0.25 = ", result.mobilities[1].mobility)
println("mobility factor at T=0.25 = ", result.mobilities[1].mobility_factor)
println("complex mobility at Ω=1,T=0.25 = ", result.responses[2, 1].mobility)
println("mobility factor at Ω=1,T=0.25 = ", result.responses[2, 1].mobility_factor)
println("impedance at Ω=1,T=0.25 = ", result.responses[2, 1].impedance)

trial = problem.trial
rate = [result.solutions[1].rate]
println("return probability τ=0 = ", return_probability(trial, 0.0, rate))
println("return probability τ=1 = ", return_probability(trial, 1.0, rate))
println("DC component mobility = ", result.mobilities[1].component_mobilities.holstein)

rows = holstein_frequency_sweep(
    [0.0, 0.5, 1.0];
    coupling = 1.5,
    temperatures = [0.5],
    options = OptimizerOptions(multistart = false),
)
println("frequency sweep rows = ", length(rows))

using PolaronMobility

println("== Fröhlich Gaussian trial families ==")

options = OptimizerOptions(multistart = false, adaptive_bounds = false, quadrature_rtol = 1e-3)

feynman = frohlich_feynman_problem(coupling = 2.0)
multi = frohlich_multi_gaussian_problem(coupling = 2.0, modes = 2)
profile = frohlich_profile_gaussian_problem(
    coupling = 2.0,
    basis_frequencies = [0.75, 1.5],
    matsubara_terms = 128,
)
nonlocal = frohlich_nonlocal_gaussian_problem(
    coupling = 2.0,
    basis_frequencies = [0.75, 1.5],
)

feynman_result = solve(feynman; temperatures = [0.0, 0.5], frequencies = [0.0], options = options)
multi_result = solve(multi; temperatures = [0.0, 0.5], frequencies = [0.0], options = options)
profile_result = solve(profile; temperatures = [0.0, 0.5], frequencies = [0.0], options = options)

for trial_result in (feynman_result, multi_result, profile_result)
    println(typeof(trial_result.problem.trial), " variational energy = ", trial_result.zero_temperature.energy.total)
    println(typeof(trial_result.problem.trial), " parameters = ", trial_result.zero_temperature.parameters)
end

v = feynman_result.zero_temperature.v
w = feynman_result.zero_temperature.w
feynman_profile_amplitude = (v^2 - w^2) / w^2
embedded_profile = frohlich_profile_gaussian_problem(
    coupling = 2.0,
    basis_frequencies = [w],
    initial_amplitudes = [feynman_profile_amplitude],
)
embedded_result = solve(embedded_profile; temperatures = [0.0], frequencies = [0.0], options = options)
println("Profile embedding of Feynman energy = ", embedded_result.zero_temperature.energy.total)

nonlocal_result = solve(nonlocal; temperatures = [0.0, 0.5], frequencies = [0.0], options = options)
println(typeof(nonlocal.trial), " pseudo-objective = ", nonlocal_result.zero_temperature.energy.total)
println("NonlocalGaussianTrial is an experimental kernel scaffold, not a validated variational energy bound.")

profile_trial = profile.trial
profile_parameters = [0.1, 0.2]
println("Γ(ω=1) = ", profile_function(profile_trial, profile_parameters, 1.0))
println("D(τ=0.5,T=0) = ", mean_square_displacement(profile_trial, profile_parameters, 0.5, Inf))
println("D(τ=0.25,β=1) = ", mean_square_displacement(profile_trial, profile_parameters, 0.25, 1.0))

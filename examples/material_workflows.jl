using PolaronMobility
using Unitful

println("== Material-derived Fröhlich workflows ==")

single_mode = FrohlichMaterial(4.5, 24.1, 0.12, 2.25)
single_problem = material_to_problem(single_mode)
single_result = solve(
    single_problem;
    temperatures = [0, 300],
    frequencies = [0, 3],
    options = OptimizerOptions(multistart = false),
)
single_units = material_units(single_result)

println("single-mode alpha = ", single_mode.alpha)
println("300 K mobility = ", single_units.mobility[2])
println("3 THz reduced conductivity = ", single_result.responses[2, 2].conductivity)

multimode = FrohlichMaterial(
    4.5,
    24.1,
    0.12,
    [4.0, 2.0],
    [0.1, 0.2],
    (6.3e-10)^3,
)

println("multimode reduced frequencies = ", multimode.phonon_frequencies ./ multimode.effective_frequency)
println("multimode couplings = ", multimode.alpha)

for trial in (:feynman, :multi_gaussian, :profile_gaussian, :nonlocal_gaussian)
    local trial_problem = material_to_problem(
        multimode;
        trial = trial,
        modes = 2,
        matsubara_terms = 64,
    )
    local trial_result = solve(
        trial_problem;
        temperatures = [0],
        frequencies = [0],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false, quadrature_rtol = 1e-3),
    )
    println(trial, " zero-temperature energy = ", trial_result.zero_temperature.energy.total)
end

println("Hellwarth B effective frequency = ", hellwarth_b_scheme([4.0 0.1; 2.0 0.2]))
println("ionic dielectric total = ", dielectric_ionic_total([4.0 0.1; 2.0 0.2], (6.3e-10)^3))

println("\n== Material-derived Rubrene Holstein workflow ==")

rubrene = rubrene_holstein_material(lattice_constant_angstrom = 7.2)
rubrene_problem = material_to_problem(rubrene)
rubrene_result = solve(
    rubrene_problem;
    temperatures = [300],
    frequencies = [0, 10],
    options = OptimizerOptions(multistart = false, adaptive_bounds = false),
)
rubrene_units = material_units(rubrene_result)

println("rubrene label = ", rubrene.label)
println("rubrene reduced hopping = ", rubrene.hopping)
println("rubrene reduced coupling = ", rubrene.coupling)
println("rubrene lambda = ", lambda_holstein(rubrene))
println("300 K Einstein mobility = ", rubrene_units.mobility_einstein[1])
println("300 K lattice-FHIP mobility = ", rubrene_units.mobility[1])

println("\n== Material-derived Rubrene Peierls workflow ==")

rubrene_peierls = rubrene_peierls_material(lattice_constant_angstrom = 7.2)
peierls_problem = material_to_problem(rubrene_peierls)
peierls_result = solve(
    peierls_problem;
    temperatures = [300],
    frequencies = [0, 10],
    options = OptimizerOptions(multistart = false, adaptive_bounds = false),
)
peierls_units = material_units(peierls_result)

println("rubrene Peierls reduced hopping = ", rubrene_peierls.hopping)
println("rubrene Peierls reduced coupling = ", rubrene_peierls.coupling)
println("rubrene Peierls lambda = ", lambda_peierls(rubrene_peierls))
println("300 K Peierls Einstein mobility = ", peierls_units.mobility_einstein[1])

println("\n== Combined Rubrene Holstein + Peierls workflow ==")

combined_problem = rubrene_holstein_peierls_problem(lattice_constant_angstrom = 7.2)
combined_result = solve(
    combined_problem;
    temperatures = [300],
    frequencies = [0, 10],
    options = OptimizerOptions(multistart = false, adaptive_bounds = false),
)

println("combined components = ", length(combined_problem.model.models))
println("combined optimized rate = ", combined_result.solutions[1].rate)
println("combined response rows = ", length(response_table(combined_result)))

@testset "Lattice Transport" begin
    rate = 0.9
    beta = 2.0
    zero_model = holstein_poisson_problem(coupling = 0.0).model
    expected_zero_factor = first_return_laplace_d(
        ComplexF64(PolaronMobility.default_lattice_broadening, 0.0),
        rate,
        1,
    )

    @test lattice_current_kernel(rate, 1, 0.0) ≈ 2 * rate
    @test lattice_holstein_phonon_factor(zero_model, beta, 1.0) == 1.0 + 0.0im
    @test lattice_mobility_factor(zero_model, rate, beta, 0.0) ≈ expected_zero_factor atol = 1e-12
    @test lattice_mobility(zero_model, rate, beta, 0.0) ≈ beta * rate * expected_zero_factor atol = 1e-12
    @test lattice_impedance(beta * rate) ≈ inv(beta * rate)

    zero_result = solve(
        holstein_poisson_problem(coupling = 0.0, hopping = 1.0);
        temperatures = [0.5],
        frequencies = [0.0],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test zero_result.mobilities[1].mobility_factor ≈ real(PolaronMobility._first_return_laplace_1d(
        ComplexF64(PolaronMobility.default_lattice_broadening, 0.0),
        zero_result.mobilities[1].diffusion_constant,
    )) atol = 1e-12
    @test zero_result.mobilities[1].mobility ≈ zero_result.mobilities[1].mobility_einstein * zero_result.mobilities[1].mobility_factor atol = 1e-12
    @test zero_result.mobilities[1].diagnostics.kappa_source == :per_temperature

    guided_model = holstein_poisson_problem(coupling = 0.55, hopping = 1.0, phonon_frequency = 1.2, dimension = 1).model
    guided_frequency = 0.37
    guided_broadening = 1e-6
    guide_sidebands = holstein_transport_sidebands(guided_model, beta; tolerance = PolaronMobility.default_lattice_sideband_tolerance)
    guide_expected = sum(
        sideband.weight * first_return_laplace_d(
            ComplexF64(guided_broadening, -(guided_frequency + sideband.frequency)),
            rate,
            1,
        ) for sideband in guide_sidebands
    )
    @test lattice_mobility_factor(
        guided_model,
        rate,
        beta,
        guided_frequency;
        broadening = guided_broadening,
    ) ≈ guide_expected atol = 1e-10

    t0_model = holstein_poisson_problem(coupling = 0.55, hopping = 1.0, phonon_frequency = 1.2, dimension = 1).model
    t0_sidebands = holstein_transport_sidebands(t0_model, Inf; tolerance = 1e-12)
    @test abs(sum(real(sideband.weight) for sideband in t0_sidebands) - 1) ≤ 1e-12

    crossover_model = holstein_poisson_problem(coupling = 1.8, hopping = 1.0, phonon_frequency = 1.0, dimension = 1).model
    crossover_sidebands = holstein_transport_sidebands(crossover_model, inv(0.03162277660168379); tolerance = 1e-12)
    @test length(crossover_sidebands) > 1
    @test abs(sum(real(sideband.weight) for sideband in crossover_sidebands) - 1) ≤ 1e-12

    strong_model = holstein_poisson_problem(coupling = 4.0, hopping = 1.0, phonon_frequency = 1.0, dimension = 1).model
    strong_highT_sidebands = holstein_transport_sidebands(strong_model, inv(30.0); tolerance = 1e-12)
    @test length(strong_highT_sidebands) > 100
    @test abs(sum(real(sideband.weight) for sideband in strong_highT_sidebands) - 1) ≤ 1e-12

    stronger_model = holstein_poisson_problem(coupling = 8.0, hopping = 1.0, phonon_frequency = 1.0, dimension = 1).model
    stronger_highT_sidebands = holstein_transport_sidebands(stronger_model, inv(10.0); tolerance = 1e-12)
    @test length(stronger_highT_sidebands) > 100
    @test abs(sum(real(sideband.weight) for sideband in stronger_highT_sidebands) - 1) ≤ 1e-12

    general_d_factor = lattice_mobility_factor(
        holstein_poisson_problem(coupling = 0.25, hopping = 1.0, phonon_frequency = 1.1, dimension = 2).model,
        rate,
        beta,
        0.2;
        laguerre_points = 80,
    )
    @test isfinite(real(general_d_factor))
    @test isfinite(imag(general_d_factor))

    holstein = holstein_poisson_problem(coupling = 0.4, hopping = 1.0)
    peierls = peierls_poisson_problem(coupling = 0.3, hopping = 1.0)
    FH = lattice_mobility_factor(holstein.model, rate, beta, 0.0)
    FP = lattice_mobility_factor(peierls.model, rate, beta, 0.0)
    @test 0 < real(FH) <= 1
    @test 0 < real(FP) <= 1
    @test imag(FH) ≥ 0
    @test isfinite(imag(FP))

    composite = combine_models(holstein.model, peierls.model)
    FC = lattice_mobility_factor(composite, rate, beta, 0.0)
    @test 0 < real(FC) <= 1
    @test abs(sum(real(sideband.weight) for sideband in holstein_peierls_transport_sidebands(composite, beta)) - 1) ≤ 1e-12

    result = solve(
        VariationalProblem(composite, PoissonTrial(bare_hopping = 1.0));
        temperatures = 1.0,
        frequencies = [0.0, 0.5],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test size(result.responses) == (2, 1)
    @test result.responses[1, 1] isa LatticeResponseResult
    @test result.mobilities[1].mobility > 0
    @test hasproperty(result.responses[1, 1].component_mobilities, :holstein)
    @test hasproperty(result.responses[1, 1].component_mobilities, :peierls)
    @test result.responses[1, 1].diagnostics.sideband_weight_sum ≈ 1.0 atol = 1e-12

    zero_temperature_result = solve(
        holstein_poisson_problem(coupling = 1.0);
        temperatures = [0.0],
        frequencies = [0.5, 1.0],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test zero_temperature_result.responses[1, 1].mobility ≈ zero_temperature_result.responses[1, 1].mobility_factor atol = 1e-12
    @test zero_temperature_result.responses[2, 1].conductivity ≈ zero_temperature_result.responses[2, 1].mobility_factor atol = 1e-12
    @test isfinite(real(zero_temperature_result.responses[1, 1].conductivity))
    @test isfinite(imag(zero_temperature_result.responses[1, 1].conductivity))
    @test zero_temperature_result.responses[1, 1].impedance ≈ inv(zero_temperature_result.responses[1, 1].conductivity) atol = 1e-12

    zero_temperature_peierls = solve(
        peierls_poisson_problem(coupling = 0.8);
        temperatures = [0.0],
        frequencies = [0.5],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test zero_temperature_peierls.responses[1, 1].mobility ≈ zero_temperature_peierls.responses[1, 1].mobility_factor atol = 1e-12
    @test isfinite(real(zero_temperature_peierls.responses[1, 1].conductivity))
end

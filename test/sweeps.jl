@testset "Continuation Sweeps" begin
    rows = continued_frohlich_temperature_sweep(
        [0.25, 0.5];
        coupling = 1.0,
        phonon_frequency = 1.0,
        options = OptimizerOptions(multistart = false),
    )

    @test length(rows) == 2
    @test rows[1].temperature == 0.25
    @test rows[1].beta == 4.0
    @test rows[2].temperature == 0.5
    @test rows[2].beta == 2.0
    @test rows[1].selected_branch in (:forward, :backward)
    @test hasproperty(rows[1], :forward_free_energy)
    @test hasproperty(rows[1], :backward_free_energy)
    @test hasproperty(rows[1], :w)
    @test hasproperty(rows[1], :delta)
    @test hasproperty(rows[1], :optimizer_success)

    coupling_rows = continued_frohlich_coupling_sweep(
        [1.0, 2.0];
        phonon_frequency = 1.0,
        temperature = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test [row.coupling for row in coupling_rows] == [1.0, 2.0]

    adiabatic_rows = continued_frohlich_adiabaticity_sweep(
        [0.5, 2.0];
        coupling = 1.0,
        temperature = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test adiabatic_rows[1].adiabatic_ratio == 0.5
    @test adiabatic_rows[1].phonon_frequency == 2.0
    @test adiabatic_rows[2].adiabatic_ratio == 2.0
    @test adiabatic_rows[2].phonon_frequency == 0.5

    holstein_rows = continued_holstein_temperature_sweep(
        [0.25, 0.5];
        coupling = 1.0,
        options = OptimizerOptions(multistart = false),
    )
    @test length(holstein_rows) == 2
    @test holstein_rows[1].temperature == 0.25
    @test holstein_rows[1].beta == 4.0
    @test holstein_rows[1].selected_branch in (:forward, :backward)
    @test hasproperty(holstein_rows[1], :rate)
    @test hasproperty(holstein_rows[1], :mobility)
    @test hasproperty(holstein_rows[1], :mobility_einstein)
    @test hasproperty(holstein_rows[1], :mobility_factor)
    @test hasproperty(holstein_rows[1], :lambda_holstein)

    holstein_low_temperature_rows = continued_holstein_temperature_sweep(
        [10.0^x for x in -2.0:0.1:-1.0];
        coupling = 1.8,
        options = OptimizerOptions(multistart = false),
    )
    low_temperature_factors = [row.mobility_factor for row in holstein_low_temperature_rows]
    @test maximum(low_temperature_factors) < 0.1

    holstein_strong_temperature_rows = continued_holstein_temperature_sweep(
        [10.0^x for x in -1.2:0.1:-0.8];
        coupling = 4.0,
        options = OptimizerOptions(multistart = false),
    )
    strong_temperature_factors = [row.mobility_factor for row in holstein_strong_temperature_rows]
    @test maximum(strong_temperature_factors) < 0.01

    holstein_couplings = continued_holstein_coupling_sweep(
        [0.0, 1.0];
        temperature = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test [row.coupling for row in holstein_couplings] == [0.0, 1.0]

    holstein_adiabaticity = continued_holstein_adiabaticity_sweep(
        [0.5, 2.0];
        coupling = 1.0,
        temperature = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test holstein_adiabaticity[1].adiabatic_ratio == 0.5
    @test holstein_adiabaticity[1].phonon_frequency == 2.0
    @test holstein_adiabaticity[2].adiabatic_ratio == 2.0
    @test holstein_adiabaticity[2].phonon_frequency == 0.5

    peierls_rows = continued_peierls_temperature_sweep(
        [0.25, 0.5];
        coupling = 1.0,
        options = OptimizerOptions(multistart = false),
    )
    @test length(peierls_rows) == 2
    @test peierls_rows[1].temperature == 0.25
    @test peierls_rows[1].beta == 4.0
    @test peierls_rows[1].selected_branch in (:forward, :backward)
    @test hasproperty(peierls_rows[1], :rate)
    @test hasproperty(peierls_rows[1], :mobility)
    @test hasproperty(peierls_rows[1], :mobility_einstein)
    @test hasproperty(peierls_rows[1], :mobility_factor)
    @test hasproperty(peierls_rows[1], :lambda_peierls)

    peierls_couplings = continued_peierls_coupling_sweep(
        [0.0, 1.0];
        temperature = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test [row.coupling for row in peierls_couplings] == [0.0, 1.0]

    peierls_adiabaticity = continued_peierls_adiabaticity_sweep(
        [0.5, 2.0];
        coupling = 1.0,
        temperature = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test peierls_adiabaticity[1].adiabatic_ratio == 0.5
    @test peierls_adiabaticity[1].phonon_frequency == 2.0
    @test peierls_adiabaticity[2].adiabatic_ratio == 2.0
    @test peierls_adiabaticity[2].phonon_frequency == 0.5

    frohlich_result = solve(
        frohlich_feynman_problem(coupling = 1.0);
        temperatures = [0.0, 0.5],
        frequencies = [0.0, 1.0],
        options = OptimizerOptions(multistart = false),
    )
    @test length(solution_table(frohlich_result)) == 2
    @test length(mobility_table(frohlich_result)) == 2
    @test length(response_table(frohlich_result)) == 4
    @test length(sweep_table(frohlich_result)) == 4
    @test hasproperty(first(response_table(frohlich_result)), :conductivity_abs)
    @test hasproperty(first(response_table(frohlich_result)), :w)

    multi_result = solve(
        frohlich_multi_gaussian_problem(coupling = 1.0, modes = 2);
        temperatures = 0.0,
        frequencies = [0.0, 1.0],
        options = OptimizerOptions(multistart = false),
    )
    multi_rows = response_table(multi_result)
    @test length(multi_rows) == 2
    @test hasproperty(first(multi_rows), :w1)
    @test hasproperty(first(multi_rows), :delta2)

    frequency_rows = frequency_sweep(
        frohlich_feynman_problem(coupling = 1.0);
        temperatures = 0.5,
        frequencies = [0.0, 1.0],
        options = OptimizerOptions(multistart = false),
    )
    @test length(frequency_rows) == 2
    @test [row.frequency for row in frequency_rows] == [0.0, 1.0]

    holstein_frequency_rows = holstein_frequency_sweep(
        [0.0, 1.0];
        coupling = 1.0,
        temperatures = 0.5,
        options = OptimizerOptions(multistart = false),
    )
    @test length(holstein_frequency_rows) == 2
    @test hasproperty(first(holstein_frequency_rows), :mobility_abs)
    @test hasproperty(first(holstein_frequency_rows), :mobility_factor_abs)
    @test hasproperty(first(holstein_frequency_rows), :conductivity_abs)

    guide_rows = holstein_transport_sweep(
        coupling = 0.6,
        temperatures = [0.5, 1.0],
        frequencies = [0.0, 0.5],
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test length(guide_rows) == 4
    @test first(guide_rows).kappa_source == :zero_temperature
    @test hasproperty(first(guide_rows), :cached_frequency_shifts)
    @test hasproperty(first(guide_rows), :sideband_weight_sum)

    per_temperature_rows = holstein_transport_sweep(
        coupling = 0.6,
        temperatures = [0.5],
        frequencies = [0.0],
        kappa_source = :per_temperature,
        options = OptimizerOptions(multistart = false, adaptive_bounds = false),
    )
    @test first(per_temperature_rows).kappa_source == :per_temperature

    mktempdir() do dir
        path = write_sweep_csv(joinpath(dir, "frequency.csv"), frequency_rows)
        @test isfile(path)
        text = read(path, String)
        @test occursin("frequency", text)
        @test occursin("conductivity_abs", text)
    end

    @test_throws ArgumentError plot_frequency_sweep(frequency_rows)

    if Base.find_package("Plots") !== nothing
        ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
        @eval using Plots
        figure = plot_frequency_sweep(frequency_rows)
        mktempdir() do dir
            output = joinpath(dir, "frequency.png")
            savefig(figure, output)
            @test isfile(output)
        end
    end
end

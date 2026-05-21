@testset "API" begin
    @test !isdefined(PolaronMobility, :FrohlichProblem)
    @test isdefined(PolaronMobility, :AbstractMaterial)
    @test isdefined(PolaronMobility, :AbstractContinuumModel)
    @test isdefined(PolaronMobility, :AbstractLatticeModel)
    @test isdefined(PolaronMobility, :VariationalProblem)
    @test isdefined(PolaronMobility, :FrohlichModel)
    @test isdefined(PolaronMobility, :PolaronResult)
    @test isdefined(PolaronMobility, :GaussianFeynmanTrial)
    @test isdefined(PolaronMobility, :MultiGaussianTrial)
    @test isdefined(PolaronMobility, :NonlocalGaussianTrial)
    @test isdefined(PolaronMobility, :ProfileGaussianTrial)
    @test isdefined(PolaronMobility, :AbstractJumpTrial)
    @test isdefined(PolaronMobility, :FrohlichSolution)
    @test isdefined(PolaronMobility, :LatticeSolution)
    @test isdefined(PolaronMobility, :FrohlichMobilityResult)
    @test isdefined(PolaronMobility, :FrohlichResponseResult)
    @test isdefined(PolaronMobility, :LatticeMobilityResult)
    @test isdefined(PolaronMobility, :LatticeResponseResult)
    @test isdefined(PolaronMobility, :FrohlichMaterial)
    @test isdefined(PolaronMobility, :HolsteinMaterial)
    @test isdefined(PolaronMobility, :PeierlsMaterial)
    @test isdefined(PolaronMobility, :PeierlsModel)
    @test isdefined(PolaronMobility, :CompositePolaronModel)
    @test isdefined(PolaronMobility, :FrohlichResultUnits)
    @test isdefined(PolaronMobility, :LatticeResultUnits)
    @test isdefined(PolaronMobility, :frohlich_feynman_problem)
    @test isdefined(PolaronMobility, :solve_frohlich)
    @test isdefined(PolaronMobility, :frohlich_multi_gaussian_problem)
    @test isdefined(PolaronMobility, :frohlich_nonlocal_gaussian_problem)
    @test isdefined(PolaronMobility, :frohlich_profile_gaussian_problem)
    @test isdefined(PolaronMobility, :peierls_poisson_problem)
    @test isdefined(PolaronMobility, :solve_peierls)
    @test isdefined(PolaronMobility, :combine_models)
    @test isdefined(PolaronMobility, :frequency_sweep)
    @test isdefined(PolaronMobility, :solution_table)
    @test isdefined(PolaronMobility, :mobility_table)
    @test isdefined(PolaronMobility, :response_table)
    @test isdefined(PolaronMobility, :write_sweep_csv)
    @test isdefined(PolaronMobility, :continued_peierls_temperature_sweep)
    @test isdefined(PolaronMobility, :continued_peierls_coupling_sweep)
    @test isdefined(PolaronMobility, :continued_peierls_adiabaticity_sweep)
    @test isdefined(PolaronMobility, :periodic_phonon_kernel)
    @test isdefined(PolaronMobility, :lattice_q0)
    @test isdefined(PolaronMobility, :lattice_q1)
    @test isdefined(PolaronMobility, :lattice_current_kernel)
    @test isdefined(PolaronMobility, :site_return_bridge)
    @test isdefined(PolaronMobility, :bond_order_bridge)
    @test isdefined(PolaronMobility, :bond_current_bridge)
    @test isdefined(PolaronMobility, :holstein_integral_d)
    @test isdefined(PolaronMobility, :peierls_integral_d)
    @test isdefined(PolaronMobility, :lattice_green_function_d)
    @test isdefined(PolaronMobility, :first_return_laplace_d)
    @test isdefined(PolaronMobility, :lattice_holstein_phonon_factor)
    @test isdefined(PolaronMobility, :lattice_peierls_phonon_factor)
    @test isdefined(PolaronMobility, :holstein_transport_sidebands)
    @test isdefined(PolaronMobility, :peierls_transport_sidebands)
    @test isdefined(PolaronMobility, :holstein_peierls_transport_sidebands)
    @test isdefined(PolaronMobility, :lattice_mobility_factor)
    @test isdefined(PolaronMobility, :lattice_mobility)
    @test isdefined(PolaronMobility, :lattice_conductivity)
    @test isdefined(PolaronMobility, :lattice_impedance)
    @test isdefined(PolaronMobility, :lattice_transport_sweep)
    @test isdefined(PolaronMobility, :holstein_transport_sweep)
    @test isdefined(PolaronMobility, :peierls_transport_sweep)
    @test isdefined(PolaronMobility, :holstein_peierls_transport_sweep)
    @test isdefined(PolaronMobility, :frohlich_alpha)
    @test isdefined(PolaronMobility, :fhip_low_temperature_mobility)
    @test isdefined(PolaronMobility, :kadanoff_low_temperature_mobility)
    @test isdefined(PolaronMobility, :mobility_cm2_per_v_s)
    @test !isdefined(PolaronMobility, :SolverOptions)
    @test !isdefined(PolaronMobility, :FrohlichPolaron)
    @test !isdefined(PolaronMobility, :VariationalSolution)
    @test !isdefined(PolaronMobility, :MobilityResult)
    @test !isdefined(PolaronMobility, :ResponseResult)
    @test !isdefined(PolaronMobility, :frohlichpolaron)
    @test !isdefined(PolaronMobility, :feynmanvw)
    @test !isdefined(PolaronMobility, :frohlichalpha)
    @test !isdefined(PolaronMobility, :fhip_mobility_lowT)
    @test !isdefined(PolaronMobility, :kadanoff_mobility_lowT)
    @test !isdefined(PolaronMobility, :frohlich_complex_impedence)
    @test !isdefined(PolaronMobility, :FrohlichResult)
    @test !isdefined(PolaronMobility, :HolsteinResult)
    @test !isdefined(PolaronMobility, :PeierlsResult)
    @test !isdefined(PolaronMobility, :HolsteinSolution)
    @test !isdefined(PolaronMobility, :PeierlsSolution)
    @test !isdefined(PolaronMobility, :HolsteinMobilityResult)
    @test !isdefined(PolaronMobility, :HolsteinResponseResult)
    @test !isdefined(PolaronMobility, :PeierlsMobilityResult)
    @test !isdefined(PolaronMobility, :PeierlsResponseResult)
    @test !isdefined(PolaronMobility, :ResultUnits)
    @test !isdefined(PolaronMobility, :HolsteinResultUnits)
    @test !isdefined(PolaronMobility, :continued_sweep)
    @test !isdefined(PolaronMobility, :variational_result_row)
    @test !isdefined(PolaronMobility, :Material)

    variational_problem = frohlich_feynman_problem(coupling = 1.0)
    @test variational_problem isa VariationalProblem{FrohlichModel,GaussianFeynmanTrial}
    @test variational_problem.model isa AbstractContinuumModel
    @test variational_problem.trial isa AbstractGaussianTrial

    multi_problem = frohlich_multi_gaussian_problem(coupling = 1.0, modes = 2)
    @test multi_problem isa VariationalProblem{FrohlichModel,MultiGaussianTrial}
    @test parameter_names(multi_problem.trial) == [:w1, :delta1, :w2, :delta2]
    @test_throws ArgumentError frohlich_multi_gaussian_problem(coupling = 1.0, modes = 0)

    result = solve(variational_problem; temperatures = 0.0, frequencies = 0.0)

    @test result isa PolaronResult
    @test result.problem isa VariationalProblem{FrohlichModel,GaussianFeynmanTrial}
    @test result.problem.model.alpha == [1.0]
    @test result.temperatures == [0.0]
    @test result.frequencies == [0.0]
    @test length(result.solutions) == 1
    @test size(result.responses) == (1, 1)
    @test result.zero_temperature.variational isa VariationalResult
    @test result.zero_temperature.parameters.w == result.zero_temperature.w
    @test feynman_v(result.zero_temperature) == result.zero_temperature.v
    @test fieldtype(typeof(result.zero_temperature), :v) === Float64
    @test fieldtype(typeof(result.zero_temperature), :energy) === EnergyComponents{Float64}
    @test eltype(result.solutions) === FrohlichSolution
    @test eltype(result.mobilities) === FrohlichMobilityResult
    @test eltype(result.responses) === FrohlichResponseResult

    holstein_problem = holstein_poisson_problem(coupling = 1.0)
    holstein = solve(holstein_problem; temperatures = 1.0, frequencies = 0.0)
    @test holstein isa PolaronResult
    @test holstein.problem.model isa AbstractLatticeModel
    @test holstein.problem.trial isa AbstractJumpTrial
    @test holstein.temperatures == [1.0]
    @test holstein.frequencies == [0.0]
    @test length(holstein.solutions) == 1
    @test size(holstein.responses) == (1, 1)
    @test eltype(holstein.mobilities) <: LatticeMobilityResult
    @test eltype(holstein.responses) <: LatticeResponseResult
    @test hasproperty(holstein.mobilities[1], :diagnostics)
    @test hasproperty(holstein.responses[1, 1], :diagnostics)

    peierls_problem = peierls_poisson_problem(coupling = 0.5)
    peierls = solve(peierls_problem; temperatures = 1.0, frequencies = 0.0)
    @test peierls isa PolaronResult
    @test peierls.temperatures == [1.0]
    @test peierls.frequencies == [0.0]
    @test length(peierls.solutions) == 1
    @test size(peierls.responses) == (1, 1)
    @test eltype(peierls.mobilities) <: LatticeMobilityResult
    @test eltype(peierls.responses) <: LatticeResponseResult

    composite = VariationalProblem(combine_models(holstein_problem.model, peierls_problem.model), PoissonTrial())
    composite_result = solve(composite; temperatures = 1.0, frequencies = 0.0)
    @test composite_result isa PolaronResult
    @test composite_result.problem.model isa CompositePolaronModel

    @test_throws ArgumentError FrohlichModel([1.0, 2.0], [1.0])
    @test_throws DomainError FrohlichModel(-1.0)
    @test_throws DomainError solve(variational_problem; temperatures = [-1.0])
    @test_throws MethodError OptimizerOptions(initial_v = 3.11)
end

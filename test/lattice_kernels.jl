@testset "Lattice Kernels" begin
    rate = 0.8
    dimension = 2
    beta = 4.0
    u = 0.75

    @test lattice_q0(rate, dimension, 0.0) == 1.0
    @test lattice_q1(rate, dimension, 0.0) == 0.0
    @test periodic_phonon_kernel(u, beta, 1.2) ≈ periodic_phonon_kernel(beta - u, beta, 1.2)

    @test site_return_bridge(rate, dimension, u, beta) ≈ site_return_bridge(rate, dimension, beta - u, beta)
    @test bond_order_bridge(rate, dimension, u, beta) ≈ bond_order_bridge(rate, dimension, beta - u, beta)
    @test bond_current_bridge(rate, dimension, u, beta) ≈ bond_current_bridge(rate, dimension, beta - u, beta)

    @test site_return_bridge(rate, dimension, u, Inf) ≈ lattice_q0(rate, dimension, u)
    @test bond_order_bridge(rate, dimension, u, Inf) ≈
          2 * dimension * (lattice_q0(rate, dimension, u) + lattice_q1(rate, dimension, u))
    @test bond_current_bridge(rate, dimension, u, Inf) ≈
          2 * dimension * (lattice_q0(rate, dimension, u) - lattice_q1(rate, dimension, u))

    s = ComplexF64(0.05, -0.3)
    @test first_return_laplace_d(s, rate, 1) ≈ PolaronMobility._first_return_laplace_1d(s, rate)
    @test isfinite(real(lattice_green_function_d(s, rate, 2; laguerre_points = 80)))
    @test isfinite(imag(lattice_green_function_d(s, rate, 2; laguerre_points = 80)))
    direct_green(d) = QuadGK.quadgk(t -> exp(-s * t) * lattice_q0(rate, d, t), 0, Inf; rtol = 1e-8)[1]
    @test lattice_green_function_d(s, rate, 2; laguerre_points = 80) ≈ direct_green(2) rtol = 1e-8
    @test lattice_green_function_d(s, rate, 3; laguerre_points = 80) ≈ direct_green(3) rtol = 1e-8
    @test holstein_integral_d(rate, 1, 1.2) ≈ inv(sqrt(1.2 * (1.2 + 4 * rate))) atol = 1e-12
    @test peierls_integral_d(rate, 2, 1.2; laguerre_points = 80) ≈
          2 * 2 * real(
              lattice_green_function_d(ComplexF64(1.2, 0.0), rate, 2; laguerre_points = 80) *
              (1 + first_return_laplace_d(ComplexF64(1.2, 0.0), rate, 2; laguerre_points = 80))
          ) atol = 1e-12
    @test peierls_integral_d(rate, 2, 1.2; laguerre_points = 80) > 0

    @test_throws DomainError lattice_q0(-1.0, dimension, u)
    @test_throws ArgumentError lattice_q0(rate, 0, u)
    @test_throws DomainError periodic_phonon_kernel(u, beta, -1.0)
end

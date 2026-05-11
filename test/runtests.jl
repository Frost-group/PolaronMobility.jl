using PolaronMobility
using QuadGK
using Test
using Unitful

@testset "PolaronMobility" begin
    include("api.jl")
    include("kernels.jl")
    include("lattice_kernels.jl")
    include("variational_solver.jl")
    include("mobility.jl")
    include("response.jl")
    include("materials_units.jl")
    include("sweeps.jl")
    include("holstein.jl")
    include("peierls.jl")
    include("full_lattice_free_energy.jl")
    include("lattice_transport.jl")
    include("docs_smoke.jl")
end

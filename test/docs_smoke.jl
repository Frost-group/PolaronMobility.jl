@testset "Docs And Examples" begin
    missing_docstrings = Symbol[]
    for name in names(PolaronMobility)
        name === :PolaronMobility && continue
        Base.Docs.hasdoc(PolaronMobility, name) || push!(missing_docstrings, name)
    end
    @test isempty(missing_docstrings)

    examples_dir = joinpath(@__DIR__, "..", "examples")
    example_files = sort(filter(endswith(".jl"), readdir(examples_dir; join = true)))
    @test !isempty(example_files)
    for example in example_files
        @test_nowarn redirect_stdout(devnull) do
            include(example)
        end
    end
end

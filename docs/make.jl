pushfirst!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using PolaronMobility, Documenter

makedocs(
    sitename = "PolaronMobility.jl",
    modules = [PolaronMobility],
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Scientific Discussion" => "scientific_discussion.md",
        "Lattice Transport" => "lattice_transport.md",
        "API Reference" => "functions.md",
    ],
)

deploydocs(repo = "github.com/Frost-group/PolaronMobility.jl.git")

push!(LOAD_PATH,"../src/") # load module from local directory

using PolaronMobility, Documenter

makedocs(sitename="PolaronMobility.jl documentation")

deploydocs(repo="github.com/Frost-group/PolaronMobility.jl.git",)


using PolaronMobility

println("== General-d CTMC transport sweep ==")

options = OptimizerOptions(multistart = false, adaptive_bounds = false)

holstein_rows = holstein_transport_sweep(
    coupling = 0.45,
    hopping = 1.0,
    phonon_frequency = 1.1,
    dimension = 2,
    temperatures = [0.5, 1.0],
    frequencies = [0.0],
    broadening = 0.02,
    options = options,
)

println("rows = ", length(holstein_rows))
println("first row mobility = ", first(holstein_rows).mobility)
println("first row kernel = ", first(holstein_rows).mobility_factor_real, " + ", first(holstein_rows).mobility_factor_imag, "im")

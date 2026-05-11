using PolaronMobility

println("== Sweeps, tables, CSV, and optional plotting ==")

options = OptimizerOptions(multistart = false, quadrature_rtol = 1e-3)

coupling_rows = continued_frohlich_coupling_sweep(
    [0.5, 1.0];
    temperature = 0.5,
    options = options,
)
temperature_rows = continued_holstein_temperature_sweep(
    [0.25, 0.5];
    coupling = 1.0,
    options = options,
)
frequency_rows = frohlich_frequency_sweep(
    [0.0, 0.5, 1.0];
    coupling = 1.0,
    temperatures = [0.5],
    options = options,
)

println("Fröhlich coupling sweep rows = ", length(coupling_rows))
println("Holstein temperature sweep rows = ", length(temperature_rows))
println("Fröhlich frequency sweep rows = ", length(frequency_rows))

mktempdir() do dir
    csv_path = write_sweep_csv(joinpath(dir, "frequency.csv"), frequency_rows)
    println("wrote temporary CSV = ", isfile(csv_path))
    println("CSV path = ", csv_path)

    try
        ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
        @eval using Plots
        figure = plot_frequency_sweep(frequency_rows)
        png_path = joinpath(dir, "frequency.png")
        savefig(figure, png_path)
        println("wrote temporary plot = ", isfile(png_path))
    catch error
        println("plotting skipped: ", typeof(error))
    end
end

module PolaronMobilityPlotsExt

using Plots
import PolaronMobility:
    plot_adiabaticity_sweep,
    plot_coupling_sweep,
    plot_frequency_sweep,
    plot_response_components,
    plot_temperature_sweep

function _values(rows, column::Symbol)
    return [getproperty(row, column) for row in rows if hasproperty(row, column)]
end

function _plot_rows(rows, x::Symbol, y::Symbol; xlabel = string(x), ylabel = string(y), title = "", kwargs...)
    xs = _values(rows, x)
    ys = _values(rows, y)
    length(xs) == length(ys) || throw(ArgumentError("rows must contain both $x and $y columns."))
    return Plots.plot(xs, ys; marker = :circle, xlabel = xlabel, ylabel = ylabel, title = title, kwargs...)
end

function plot_coupling_sweep(rows::AbstractVector{<:NamedTuple}; y::Symbol = :free_energy, kwargs...)
    return _plot_rows(rows, :coupling, y; xlabel = "coupling", ylabel = string(y), title = "Coupling sweep", kwargs...)
end

function plot_temperature_sweep(rows::AbstractVector{<:NamedTuple}; y::Symbol = :free_energy, kwargs...)
    return _plot_rows(rows, :temperature, y; xlabel = "temperature", ylabel = string(y), title = "Temperature sweep", kwargs...)
end

function plot_adiabaticity_sweep(rows::AbstractVector{<:NamedTuple}; y::Symbol = :free_energy, kwargs...)
    return _plot_rows(rows, :adiabatic_ratio, y; xlabel = "adiabatic ratio", ylabel = string(y), title = "Adiabaticity sweep", kwargs...)
end

function plot_frequency_sweep(rows::AbstractVector{<:NamedTuple}; y::Symbol = :conductivity_abs, kwargs...)
    default_y = any(hasproperty(row, y) for row in rows) ? y : :mobility_abs
    return _plot_rows(rows, :frequency, default_y; xlabel = "frequency", ylabel = string(default_y), title = "Frequency sweep", kwargs...)
end

function plot_response_components(
    rows::AbstractVector{<:NamedTuple};
    components::Tuple{Vararg{Symbol}} = (:conductivity_real, :conductivity_imag, :conductivity_abs),
    kwargs...,
)
    available = [component for component in components if any(hasproperty(row, component) for row in rows)]
    isempty(available) && (available = [component for component in (:mobility_real, :mobility_imag, :mobility_abs) if any(hasproperty(row, component) for row in rows)])
    isempty(available) && throw(ArgumentError("rows do not contain recognized response component columns."))
    plot = Plots.plot(; xlabel = "frequency", ylabel = "response", title = "Response components", kwargs...)
    frequencies = _values(rows, :frequency)
    for component in available
        Plots.plot!(plot, frequencies, _values(rows, component); label = string(component))
    end
    return plot
end

end

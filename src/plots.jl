
plot(xdata, ydata::Array{uwreal}; kwargs...) = plot(xdata, value.(ydata), yerr=err.(ydata); kwargs...)
function plot!(pl::Plots.Plot, xdata, ydata::Array{uwreal}; kwargs...) 
    Plots.plot!(pl, xdata, value.(ydata); yerr=err.(ydata), kwargs...)
end

plot(ydata::Array{uwreal}; kwargs...) = plot(1:length(ydata), ydata; kwargs...)

function plot(corrws::AbstractCorrelatorAnalysis; kwargs...)
	pl = plot(corrws.xdata, corrws.ydata; kwargs...)
	display(pl)
    return nothing
end


function plot(hws::CorrelatorFitHistories; column::Int64 = 2, kwargs...)
    plot(hws.tmin, hws.fitp[column, :]; kwargs...)
end


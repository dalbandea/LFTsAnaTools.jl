
function Plots.plot(xdata, ydata::Array{uwreal}; kwargs...) 
    if ydata[1].err == zero(eltype(ydata[1].err))
        uwerr.(ydata)
    end
    Plots.plot(xdata, value.(ydata), yerr=err.(ydata); kwargs...)
end
function Plots.plot!(pl::Plots.Plot, xdata, ydata::Array{uwreal}; kwargs...) 
    if ydata[1].err == zero(eltype(ydata[1].err))
        uwerr.(ydata)
    end
    Plots.plot!(pl, xdata, value.(ydata); yerr=err.(ydata), kwargs...)
end

Plots.plot(ydata::Array{uwreal}; kwargs...) = Plots.plot(1:length(ydata), ydata; kwargs...)

function Plots.plot(corrws::AbstractCorrelatorAnalysis; kwargs...)
	pl = Plots.plot(corrws.xdata, corrws.ydata; kwargs...)
	display(pl)
    return nothing
end


function Plots.plot(hws::CorrelatorFitHistories; column::Int64 = 2, kwargs...)
    Plots.plot(hws.tmin, hws.fitp[column, :]; kwargs...)
end


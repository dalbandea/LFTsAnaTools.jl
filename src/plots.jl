
RecipesBase.@recipe function plot(uwy::Array{uwreal})
    uwerr.(uwy)
    yerror := ADerrors.err.(uwy)
    seriestype --> :scatter
    ADerrors.value.(uwy)
end

RecipesBase.@recipe function plot(xdata, uwy::Array{uwreal})
    uwerr.(uwy)
    yerror := ADerrors.err.(uwy)
    seriestype --> :scatter
    xdata, ADerrors.value.(uwy)
end

RecipesBase.@recipe function plot(corrws::AbstractCorrelatorAnalysis)
    uwerr.(corrws.ydata)
    yerror := ADerrors.err.(corrws.ydata)
    seriestype --> :scatter
    corrws.xdata, ADerrors.value.(corrws.ydata)
end

RecipesBase.@recipe function plot(emass::EffectiveMass)
    uwerr.(emass.ydata)
    yerror := ADerrors.err.(emass.ydata)
    seriestype --> :scatter
    emass.xdata, ADerrors.value.(emass.ydata)
end

RecipesBase.@recipe function plot(xdata, f::UtilsFunc, params)
    uwys = [f(x,params) for x in xdata]
    uwerr.(uwys)
    ribbons := ADerrors.err.(uwys)
    xdata, ADerrors.value.(uwys)
end


# function Plots.plot(xdata, ydata::Array{uwreal}; kwargs...) 
#     if ydata[1].err == zero(eltype(ydata[1].err))
#         uwerr.(ydata)
#     end
#     Plots.plot(xdata, value.(ydata), yerr=err.(ydata); kwargs...)
# end
# function Plots.plot!(pl::Plots.Plot, xdata, ydata::Array{uwreal}; kwargs...) 
#     if ydata[1].err == zero(eltype(ydata[1].err))
#         uwerr.(ydata)
#     end
#     Plots.plot!(pl, xdata, value.(ydata); yerr=err.(ydata), kwargs...)
# end

# Plots.plot(ydata::Array{uwreal}; kwargs...) = Plots.plot(1:length(ydata), ydata; kwargs...)

# function Plots.plot(corrws::AbstractCorrelatorAnalysis; kwargs...)
# 	pl = Plots.plot(corrws.xdata, corrws.ydata; kwargs...)
# 	display(pl)
#     return nothing
# end


# function Plots.plot(hws::CorrelatorFitHistories; column::Int64 = 2, kwargs...)
#     Plots.plot(hws.tmin, hws.fitp[column, :]; kwargs...)
# end


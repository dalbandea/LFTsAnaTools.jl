abstract type AbstractObservableAnalysis end
abstract type AbstractCorrelatorAnalysis <: AbstractObservableAnalysis end

abstract type CorrelatorSymmetry end
abstract type SymmetricCorrelator <: CorrelatorSymmetry end
abstract type AntisymmetricCorrelator <: CorrelatorSymmetry end


mutable struct CorrelatorFitHistories <: AbstractCorrelatorAnalysis
    fitp
    cs
    cse
    tmin
    function CorrelatorFitHistories()
        return new([],[],[],[])
    end
end
export CorrelatorFitHistories

mutable struct CorrelatorAnalysis <: AbstractCorrelatorAnalysis
    filepath::String
    type
    history
    T
    xdata::AbstractArray
    ydata::Vector{uwreal}
    tmin::Int64
    tmax::Int64
    histories::CorrelatorFitHistories
    ID::String
    burnout::Int64
    function CorrelatorAnalysis(filepath::String; type = SymmetricCorrelator, ensemble_ID::String = filepath, burnout::Int64 = 1)
        history = read(type, filepath)[burnout:end, :]
        T = length(history[1,:])
        x = []
        ydata = uwreal[]
        tmin = 1
        tmax = T
        histories = CorrelatorFitHistories()
        return new(filepath, type, history, T, x, ydata, tmin, tmax, histories, ensemble_ID, burnout)
    end
end

function reset_histories!(corrws::AbstractCorrelatorAnalysis)
    corrws.tmin = 1
    corrws.histories = CorrelatorFitHistories()
    return nothing
end
export reset_histories!

function uwreal(corrws::AbstractCorrelatorAnalysis)
	# Load correlation data
    corrs = corrws.history

	corrws.xdata = range(0, length=corrws.T)
	corrws.ydata = Vector{uwreal}(undef, corrws.T);

	for i in eachindex(corrws.ydata)
		corrws.ydata[i] = uwreal(corrs[:,i], corrws.ID)
	end

	return nothing
end
export uwreal

function uwrealsym(corrws::AbstractCorrelatorAnalysis)
	# Load correlation data
    corrs = corrws.history

	T2p1 = convert(Int64, corrws.T/2+1)
	corrws.xdata = range(0, length=T2p1)
	corrws.ydata = Vector{uwreal}(undef, T2p1);

	for i in eachindex(corrws.ydata)
        corrws.ydata[i] = uwreal((corrs[:,i] .+ sign(corrws.type) * corrs[:,1+(corrws.T-i+1)%corrws.T])/2, corrws.ID)
	end
	corrws.ydata[1] = uwreal(corrs[:, 1], corrws.ID )
	corrws.ydata[end] = uwreal(corrs[:, T2p1], corrws.ID )

	return nothing
end
export uwrealsym

sign(::Type{SymmetricCorrelator}) = one(eltype(Float64))
sign(::Type{AntisymmetricCorrelator}) = -one(eltype(Float64))

# function correlator_fit_function(corrws::AbstractCorrelatorAnalysis) 
#     error("correlator_fit_function not defined for $(typeof(corrws))")
# end

correlator_fit_function(corrws::AbstractCorrelatorAnalysis) = correlator_fit_function(corrws.type, corrws)

function tmin_fit(corrws::AbstractCorrelatorAnalysis, prms0::Vector{Float64}; Tmax::Int64=length(corrws.ydata)-1, plot_tmin_fit::Bool=false, yscale::Symbol=:identity)

    ydata = corrws.ydata
    xdata = corrws.xdata

	if (Tmax > xdata[end][1])
		throw("Tmax value greater than greatest value in xdata")
	end

	pos_tmin = findfirst(x->x[1]==corrws.tmin, xdata)
	pos_Tmax = findfirst(x->x[1]==Tmax, xdata)
	fit_region = pos_tmin:pos_Tmax

    f = correlator_fit_function(corrws)

	# Fit
	fitp, cse, cs = fit_data(f, xdata[fit_region], ydata[fit_region], prms0)
    uwerr.(fitp)

	# # Plot fit
	# if(plot_tmin_fit==true)
	# 	if(size(xdata[1],1)==1)
	# 		pl = plot_fit(xdata, ydata, f, fitp)
	# 		plot!(pl, xdata, value.(ydata), yerr=err.(ydata), seriestype=:scatter, title="tmin = "*string(tmin), yscale=yscale)
	# 		# plot!(ylim=(9000,10000))
	# 		display(pl)
	# 	else
	# 		res = f.(xdata, [fitp for i in 1:length(xdata)])
	# 		uwerr.(res)
	# 		pl_x = (hcat(xdata...) |> permutedims)[:,1]
	# 		pl = plot( pl_x, value.(res), ribbons=err.(res), reuse=false, title="tmin="*string(tmin)*", 1")
	# 		plot!(pl, pl_x, value.(ydata), yerr=err.(ydata), seriestype=:scatter)
	# 		display(pl)
	# 	end
	# end
	
	return fitp, cse, cs

end
export tmin_fit

function tmin_loop(corrws::AbstractCorrelatorAnalysis, prms0::Vector{Float64}; update_prms::Bool=true, Tmax::Int64=length(corrws.ydata)-1, plot_tmin_fit::Bool=false, plot_column::Int64=0)

	tminvalues = corrws.tmin:corrws.tmax

	prms = prms0

    first_tmin = corrws.tmin

	for itmin in tminvalues
        corrws.tmin = itmin
		# Fit
		fitp, cse, cs = tmin_fit(corrws, prms)
		display(cs)

        if itmin == first_tmin
            append!(corrws.histories.fitp, fitp)
            corrws.histories.fitp = convert(Array{uwreal}, corrws.histories.fitp)
        else
            corrws.histories.fitp = hcat(corrws.histories.fitp, fitp)
        end
        push!(corrws.histories.cs, cs)
        push!(corrws.histories.cse, cse)
        push!(corrws.histories.tmin, itmin)

		if(update_prms)
			prms = value.(fitp)
		end
	end

	# if(plot_column > 0)
	# 	pl = plot(fitps[:,1], value.(fitps[:,plot_column]), yerr=err.(fitps[:,plot_column]), seriestype=:scatter)
	# 	display(pl)
	# end

	return nothing
end
export tmin_loop


function derivate_sym_correlator!(corrws::AbstractCorrelatorAnalysis)

    corr_x = corrws.xdata
    corr_y = corrws.ydata

	T2p1 = length(corrws.ydata)
	dcorr_x = corr_x[2:end-1]
	dcorr_y = Vector{uwreal}(undef, T2p1-2 )

	for i in eachindex(dcorr_y)
		dcorr_y[i] = (corr_y[i+2]-corr_y[i])/2
	end

    corrws.xdata = dcorr_x
    corrws.ydata = dcorr_y

	return nothing
end
export derivate_sym_correlator!


mutable struct StraightLine <: UtilsFunc end
export StraightLine

fit(f::UtilsFunc, xdata, ydata, prms0) = fit_data(f::UtilsFunc, xdata, ydata, prms0)
export fit

(s::StraightLine)(x,p) = p[1] * x + p[2]

import Utils: nparameters
nparameters(::StraightLine) = 2

function Plots.plot(xdata, ydata::Array{uwreal}, f::UtilsFunc, fitp::Array{uwreal}; kwargs...)
    pl = plot_fit(xdata, ydata, f, fitp)
    plot!(pl, xdata, ydata; kwargs...)
end

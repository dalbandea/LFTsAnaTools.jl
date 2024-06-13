
using LFTSampling
using LFTU1
import ADerrors: uwreal
import Base: read
using DelimitedFiles


LFTsAnaTools.sign(::Type{T}) where T <: AbstractCorrelator = sign(get_symmetry(T))

get_symmetry(::Type{U1PionCorrelator}) = SymmetricCorrelator
correlator_fit_function(::Type{U1PionCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)


get_symmetry(::Type{U1PCACCorrelator}) = AntisymmetricCorrelator
correlator_fit_function(::Type{U1PCACCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)

export get_symmetry

abstract type NeutralPion end
export NeutralPion
abstract type EtaPrime end
export EtaPrime
abstract type ChargedPion <: LFTSampling.AbstractCorrelator end
export ChargedPion
abstract type QEtaPrime <: LFTSampling.AbstractCorrelator end
export QEtaPrime


"""
Reads histories conn-ifljfl.txt, disc-ifljfl.txt and delta-ifl.txt for ifl,jfl=1,2, returns
them as a vector of dimension (3, 2, 2, ncfgs, T) standing for (conn/disc/delta, ifl,
jfl, ncfgs, T). In the case of delta, notice that only Delta[1,:,:,:] is filled,
i.e. Delta[2,:,:,:] = 0.0
"""
function Base.read(type::Type{NeutralPion}, dirfilespath; burnout = 0, prefix="")
    ncfgs, T = size(readdlm(joinpath(dirfilespath, "$(prefix)conn-11.txt"),','))
    corrs = (
             ID = dirfilespath,
             conn = zeros(Float64, 2, 2, ncfgs-burnout, T),
             disc = zeros(Float64, 2, 2, ncfgs-burnout, T),
             Delta = zeros(Float64, 2, 2, ncfgs-burnout, T)
            )
    for ifl in 1:2
        deltafile = joinpath(dirfilespath, "$(prefix)delta-$ifl.txt")
        # corrs.Delta[1, ifl, :, :] .= readdlm(deltafile, ',')[burnout+1:end, :]
        for jfl in ifl:2
            connfile = joinpath(dirfilespath, "$(prefix)conn-$ifl$jfl.txt")
            discfile = joinpath(dirfilespath, "$(prefix)disc-$ifl$jfl.txt")
            corrs.conn[ifl, jfl, :, :] .= readdlm(connfile, ',')[burnout+1:end, :]
            corrs.disc[ifl, jfl, :, :] .= readdlm(discfile, ',')[burnout+1:end, :]
        end
    end
    return stack([corrs.conn, corrs.disc, corrs.Delta], dims=1), ncfgs - burnout, T
end

function Base.read(type::Type{EtaPrime}, dirfilespath; burnout = 0, prefix="")
    ncfgs, T = size(readdlm(joinpath(dirfilespath, "$(prefix)conn-11.txt"),','))
    corrs = (
             ID = dirfilespath,
             conn = zeros(Float64, 2, 2, ncfgs-burnout, T),
             disc = zeros(Float64, 2, 2, ncfgs-burnout, T),
             Delta = zeros(Float64, 2, 2, ncfgs-burnout, T)
            )
    for ifl in 1:2
        deltafile = joinpath(dirfilespath, "$(prefix)delta-$ifl.txt")
        # corrs.Delta[1, ifl, :, :] .= readdlm(deltafile, ',')[burnout+1:end, :]
        for jfl in ifl:2
            connfile = joinpath(dirfilespath, "$(prefix)conn-$ifl$jfl.txt")
            discfile = joinpath(dirfilespath, "$(prefix)disc-$ifl$jfl.txt")
            corrs.conn[ifl, jfl, :, :] .= readdlm(connfile, ',')[burnout+1:end, :]
            corrs.disc[ifl, jfl, :, :] .= readdlm(discfile, ',')[burnout+1:end, :]
        end
    end
    return stack([corrs.conn, corrs.disc, corrs.Delta], dims=1), ncfgs - burnout, T
end


"""
Reads history in conn_11.txt
"""
function Base.read(type::Type{ChargedPion}, dirfilespath::String; burnout = 0, prefix = "")
    ncfgs, T = size(readdlm(joinpath(dirfilespath, "$(prefix)conn-11.txt"),','))
    connfile = joinpath(dirfilespath, "$(prefix)conn-12.txt")
    corr = readdlm(connfile, ',')[burnout+1:end, :]
    return corr, ncfgs - burnout, T
end

LFTsAnaTools.get_symmetry(::Type{ChargedPion}) = SymmetricCorrelator

function Base.read(type::Type{QEtaPrime}, dirfilespath::String; burnout = 0, prefix = "")
    ncfgs, T = size(readdlm(joinpath(dirfilespath, "eta-prime-qq.txt"),','))
    connfile = joinpath(dirfilespath, "eta-prime-qq.txt")
    corr = readdlm(connfile, ',')[burnout+1:end, :]
    return corr, ncfgs - burnout, T
end

LFTsAnaTools.get_symmetry(::Type{QEtaPrime}) = SymmetricCorrelator


"""
Builds neutral pion into corrws.ydata
"""
function ADerrors.uwreal(corrws::AbstractCorrelatorAnalysis, type::Type{NeutralPion})
    uwhist = Array{uwreal}(undef, 2, 2, 2, corrws.T)
    for is in 1:2, ifl in 1:2, jfl in ifl:2, t in 1:corrws.T
        uwhist[is, ifl, jfl, t] = uwreal(corrws.history[is, ifl, jfl, :, t], corrws.ID)
    end
    ypi0 = (uwhist[1,1,1,:] + uwhist[1,2,2,:]) - (uwhist[2,1,1,:] + uwhist[2,2,2,:]) + 2*uwhist[2,1,2,:]
    corrws.ydata = ypi0
    corrws.xdata = 0:corrws.T-1
    return nothing
end

"""
Builds eta prime into corrws.ydata
"""
function ADerrors.uwreal(corrws::AbstractCorrelatorAnalysis, type::Type{EtaPrime})
    uwhist = Array{uwreal}(undef, 2, 2, 2, corrws.T)
    for is in 1:2, ifl in 1:2, jfl in ifl:2, t in 1:corrws.T
        uwhist[is, ifl, jfl, t] = uwreal(corrws.history[is, ifl, jfl, :, t], corrws.ID)
    end
    ypi0 = (uwhist[1,1,1,:] + uwhist[1,2,2,:]) - (uwhist[2,1,1,:] + uwhist[2,2,2,:]) - 2*uwhist[2,1,2,:]
    corrws.ydata = ypi0
    corrws.xdata = 0:corrws.T-1
    return nothing
end

"""
Symmetrizes a history (matrix of dimension nconfgsxT) returning uwreal over all
the configs for each timeslice
"""
function LFTsAnaTools.uwrealsym(hist, ID)
    T = size(hist)[2]
	T2p1 = convert(Int64, T/2+1)
	xdata = range(0, length=T2p1)
	ydata = Vector{uwreal}(undef, T2p1);
	for i in eachindex(ydata)
        ydata[i] = uwreal((hist[:,i] .+ hist[:,1+(T-i+1)%T])/2, ID)
	end
	ydata[1] = uwreal(hist[:, 1], ID)
	ydata[end] = uwreal(hist[:, T2p1], ID)
	return ydata
end

"""
Symmetrizes a history (matrix of dimension nconfgsxT) and performs reweighting
given the weights `weights` returning uwreal over all the configs for each
timeslice
"""
function uwrealsymrw(hist, weights, ID)
    T = size(hist)[2]
	T2p1 = convert(Int64, T/2+1)
	xdata = range(0, length=T2p1)
	ydata = Vector{uwreal}(undef, T2p1);
	for i in eachindex(ydata)
        ydata[i] = reweight((hist[:,i] .+ hist[:,1+(T-i+1)%T])/2, weights, ID)
	end
	ydata[1] = reweight(hist[:, 1], weights, ID)
	ydata[end] = reweight(hist[:, T2p1], weights, ID)
	return ydata
end
export uwrealsymrw


"""
Symmetrizes an uwreal vector
"""
function LFTsAnaTools.uwrealsym(uwa::Vector{uwreal})
    T = length(uwa)
	T2p1 = convert(Int64, T/2+1)
	xdata = range(0, length=T2p1)
	ydata = Vector{uwreal}(undef, T2p1);
	for i in eachindex(ydata)
        ydata[i] = (uwa[i] + uwa[1+(T-i+1)%T])/2
	end
	ydata[1] = uwa[1]
	ydata[end] = uwa[T2p1]
	return ydata
end

"""
Builds symmetric neutral pion into corrws.ydata. If corrws.reweights != nothing, it builds the reweighted correlator
"""
function LFTsAnaTools.uwrealsym(corrws::AbstractCorrelatorAnalysis, type::Type{NeutralPion})
	T2p1 = convert(Int64, corrws.T/2+1)
	corrws.xdata = range(0, length=T2p1)
    uwhist = Array{uwreal}(undef, 2, 2, 2, T2p1)
    for is in 1:2, ifl in 1:2, jfl in ifl:2
        if corrws.reweights == nothing
            uwhist[is, ifl, jfl, :] .= LFTsAnaTools.uwrealsym(corrws.history[is,ifl, jfl, :, :], corrws.ID)
        else
            uwhist[is, ifl, jfl, :] .= uwrealsymrw(corrws.history[is,ifl, jfl, :, :], corrws.reweights, corrws.ID)
        end
    end
    ypi0 = 1/2 * (uwhist[1,1,1,:] .+ uwhist[1,2,2,:]) .- 1/2 * (uwhist[2,1,1,:] .+ uwhist[2,2,2,:]) .+ uwhist[2,1,2,:]
    # uwDD = make_disconnected_pieces(corrws) # this does not seem to affect at all...
    # ypi0 = ypi0 + (1/2*uwrealsym(uwDD[1,1,:]) + 1/2*uwrealsym(uwDD[2,2,:]) - uwrealsym(uwDD[1,2,:]))
    corrws.ydata = ypi0
    return nothing
end

"""
Builds symmetric neutral pion into corrws.ydata. If corrws.reweights != nothing, it builds the reweighted correlator
"""
function LFTsAnaTools.uwrealsym(corrws::AbstractCorrelatorAnalysis, type::Type{EtaPrime})
	T2p1 = convert(Int64, corrws.T/2+1)
	corrws.xdata = range(0, length=T2p1)
    uwhist = Array{uwreal}(undef, 2, 2, 2, T2p1)
    for is in 1:2, ifl in 1:2, jfl in ifl:2
        if corrws.reweights == nothing
            uwhist[is, ifl, jfl, :] .= LFTsAnaTools.uwrealsym(corrws.history[is,ifl, jfl, :, :], corrws.ID)
        else
            uwhist[is, ifl, jfl, :] .= uwrealsymrw(corrws.history[is,ifl, jfl, :, :], corrws.reweights, corrws.ID)
        end
    end
    ypi0 = 1/2 * (uwhist[1,1,1,:] .+ uwhist[1,2,2,:]) .- 1/2 * (uwhist[2,1,1,:] .+ uwhist[2,2,2,:]) .- uwhist[2,1,2,:]
    # uwDD = make_disconnected_pieces(corrws) # this does not seem to affect at all...
    # ypi0 = ypi0 + (1/2*uwrealsym(uwDD[1,1,:]) + 1/2*uwrealsym(uwDD[2,2,:]) - uwrealsym(uwDD[1,2,:]))
    corrws.ydata = ypi0
    return nothing
end



"""
Returns disconnected traces expectation values, <tr[γD_ifl]><tr[γD_jfl]>,
storing it in uwDD[ifl, jfl, t]
"""
function make_disconnected_pieces(corrws::AbstractCorrelatorAnalysis)
    uwDelta = Array{uwreal}(undef, 2, corrws.T)
    for ifl in 1:2, t in 1:corrws.T
        uwDelta[ifl, t] = uwreal(corrws.history[3,1,ifl,:,t], corrws.ID)
    end
    uwDD = Array{uwreal}(undef, 2, 2, corrws.T)
    uwDD .= 0.0
    for ifl in 1:2, jfl in 1:2
        for t in 1:corrws.T, tt in 1:corrws.T
            uwDD[ifl, jfl, t] += uwDelta[ifl, tt] * uwDelta[jfl, (tt+t-1-1)%corrws.T+1] / corrws.T
        end
    end
    return uwDD
end



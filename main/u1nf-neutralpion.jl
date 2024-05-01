using Revise
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTU1
using LFTsAnaTools
using ADerrors
using LaTeXStrings
using Plots
using TOML
using DelimitedFiles
using Utils
import Base: read
import ADerrors: uwreal
import LFTsAnaTools: uwrealsym, correlator_fit_function

plotlyjs()


# fdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/Nfsim-b4.0-L24-m[0.02, 0.02]_D2024-04-18-18-42-12.272"

# ncfgs = 12009

# fdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/Nfsim-b4.0-L24-m[-0.02, 0.06]_D2024-04-18-18-51-24.285/"

# ncfgs = 12022


length(ARGS) == 1 || error("Only one argument is expected! (Path to input directory)")
isdir(ARGS[1]) || error("Path provided is not a directory")
wdir = ARGS[1]


abstract type NeutralPion end

function Base.read(type::Type{NeutralPion}, dirfilespath; burnout = 0)
    ncfgs, T = size(readdlm(joinpath(dirfilespath, "conn-11.txt"),','))
    corrs = (
             ID = wdir,
             conn = zeros(Float64, 2, 2, ncfgs-burnout, T),
             disc = zeros(Float64, 2, 2, ncfgs-burnout, T)
            )
    for ifl in 1:2, jfl in ifl:2
        connfile = joinpath(wdir, "conn-$ifl$jfl.txt")
        discfile = joinpath(wdir, "disc-$ifl$jfl.txt")
        corrs.conn[ifl, jfl, :, :] .= readdlm(connfile, ',')[burnout+1:end, :]
        corrs.disc[ifl, jfl, :, :] .= readdlm(discfile, ',')[burnout+1:end, :]
    end
    return stack([corrs.conn, corrs.disc], dims=1), ncfgs - burnout, T
end

function ADerrors.uwreal(corrws::AbstractCorrelatorAnalysis, type::Type{NeutralPion})
    uwhist = Array{uwreal}(undef, 2, 2, 2, corrws.T)
    for is in 1:2, ifl in 1:2, jfl in ifl:2, t in 1:corrws.T
        uwhist[is, ifl, jfl, t] = uwreal(corrws.history[is, ifl, jfl, :, t], corrws.ID)
    end
    ypi0 = 2*(uwhist[1,1,1,:] + uwhist[1,2,2,:]) - (uwhist[2,1,1,:] + uwhist[2,2,2,:]) + 2*uwhist[2,1,2,:]
    corrws.ydata = ypi0
    corrws.xdata = 0:corrws.T-1
    return nothing
end

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

function LFTsAnaTools.uwrealsym(corrws::AbstractCorrelatorAnalysis, type::Type{NeutralPion})
	T2p1 = convert(Int64, corrws.T/2+1)
	corrws.xdata = range(0, length=T2p1)
    uwhist = Array{uwreal}(undef, 2, 2, 2, T2p1)
    for is in 1:2, ifl in 1:2, jfl in ifl:2
        uwhist[is, ifl, jfl, :] .= LFTsAnaTools.uwrealsym(corrws.history[is,ifl, jfl, :, :], corrws.ID)
    end
    ypi0 = 2*(uwhist[1,1,1,:] + uwhist[1,2,2,:]) - (uwhist[2,1,1,:] + uwhist[2,2,2,:]) + 2*uwhist[2,1,2,:]
    corrws.ydata = ypi0
    return nothing
end

correlator_fit_function(type::Type{NeutralPion}, corrws) = SymCorrelator(T=corrws.T, s=1, c=0)

function parse_mass(filepath)
    masses = replace(filter(!isspace, filepath), r".*\-m(\[.*\]).*" => s"\1")
    return Meta.eval(Meta.parse(masses))
end


function plot_and_save(emass::EffectiveMass)
    title = emass.title
    txtpath = joinpath(wdir, title)*"-m$(emass.m0).txt"
    pdfpath = joinpath(wdir, title)*"-m$(emass.m0).pdf"
    htmlpath = joinpath(wdir, title)*"-m$(emass.m0).html"
    mywrite(txtpath, emass)
    pl = plot(emass, title=title*" for m0=$(emass.m0)", seriestype=:scatter)
    savefig(pl, pdfpath)
    savefig(pl, htmlpath)
    return nothing
end

a = CorrelatorAnalysis(wdir, type = NeutralPion, ensemble_ID = string(parse_mass(wdir)))

LFTsAnaTools.uwrealsym(a)

mpi = pion_effective_mass(a)
mpifit = pion_fit_effective_mass(a)

uwerr(mpi)
uwerr(mpifit)

plot_and_save(mpi)
plot_and_save(mpifit)

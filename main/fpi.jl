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
import LFTsAnaTools: correlator_fit_function
import BDIO
using Utils

plotlyjs()

length(ARGS) == 1 || error("Only one argument is expected! (Path to input file)")
isdir(ARGS[1]) || error("Path provided is not a directory")

wdir = ARGS[1]

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/Nf2sim-b5.0-L24-m0.02_D2024-03-12-15-54-02.07/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta4.0/Nf2sim-b4.0-L64-m0.02_D2024-03-14-18-42-47.905/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta4.0/Nf2sim-b4.0-L64-m-0.01_D2024-03-14-18-42-47.839/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta4.0/Nf2sim-b4.0-L64-m-0.04_D2024-03-14-18-42-47.751/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta4.0/Nf2sim-b4.0-L64-m-0.07_D2024-03-14-18-42-47.667/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta5.0/Nf2sim-b5.0-L64-m-0.03_D2024-03-26-11-31-44.306/measurements" 

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta5.0/Nf2sim-b5.0-L64-m-0.06_D2024-03-26-11-31-44.479/measurements" 

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta5.0/Nf2sim-b5.0-L64-m0.005_D2024-03-26-11-31-44.427/measurements" 

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta5.0/Nf2sim-b5.0-L64-m0.025_D2024-03-26-11-31-44.367/measurements" 

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta6.0/Nf2sim-b6.0-L64-m-0.025_D2024-03-26-11-28-33.148/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta6.0/Nf2sim-b6.0-L64-m-0.05_D2024-03-26-11-28-33.085/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta6.0/Nf2sim-b6.0-L64-m0.025_D2024-03-26-11-28-33.197/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/1-Smilga/L64/beta6.0/Nf2sim-b6.0-L64-m0.0_D2024-03-26-11-28-33.261/measurements"

# wdir = "/home/david/"


"""
Returns vector of Strings with occurrences in the directory using regex, but must escape \$ with \\.
"""
function filterdir(dir::String, texts::Array{String})
    dirfiles = readdir(dir)
    occurrences = filter(s -> occursin(Regex("$(texts[1])"), s), dirfiles)
    for text in texts
        occurrences = filter(s -> occursin(Regex("$text"), s), occurrences)
    end
    return joinpath.(dir, occurrences)
end
filterdir(dir::String, text::String) = filterdir(dir, [text])

function parse_mass(filepath)
    mass = replace(filepath, r".*\-m(\-?[0-9]+\.[0-9]+).*" => s"\1")
    return f64(mass)
end

function extract_correlators(dict, type, burnout)
    wss = []
    for k in keys(dict)
        ws = CorrelatorAnalysis(dict[k], type = type, burnout = burnout, ensemble_ID = k)
        push!(wss, ws)
    end
    return wss
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

import Base: read

function Base.read(type, dirfilespath::String; burnout = 0, prefix = "")
    ncfgs, T = size(readdlm(dirfilespath,','))
    connfile = dirfilespath
    corr = readdlm(connfile, ',')[burnout+1:end, :]
    return corr, ncfgs - burnout, T
end

PPpath = filterdir(wdir, "pion.txt")[1]
APpath = filterdir(wdir, "pcac.txt")[1]

ppcorr = CorrelatorAnalysis(PPpath, type = SymmetricCorrelator)
apcorr = CorrelatorAnalysis(APpath, type = AntisymmetricCorrelator)

uwrealsym(ppcorr)
uwrealsym(apcorr)


LFTsAnaTools.correlator_fit_function(::Type{SymmetricCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)

LFTsAnaTools.correlator_fit_function(::Type{AntisymmetricCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=-1, c=0)

# Create analysis directory

anapath = joinpath(wdir, "analysis-$(size(ppcorr.history)[1])")
mkpath(anapath)

title = "L = $(ppcorr.T)"

# Fit

reset_histories!(ppcorr)
ppcorr.tmin = 0
ppcorr.tmax = 23
ppcorr.Tmax = 25
ppmeff, APP = pion_fit_effective_mass_and_constant(ppcorr)

reset_histories!(apcorr)
apcorr.tmin = 1
apcorr.tmax = 23
apcorr.Tmax = 25
apmeff, AAP = pion_fit_effective_mass_and_constant(apcorr)


# Plot fits pi0

fitpath = joinpath(anapath, "fits-PP")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(SymmetricCorrelator, ppcorr)
for i in eachindex(ppcorr.histories.tmin)
    tmin = ppcorr.histories.tmin[i]
    prms = ppcorr.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:ppcorr.Tmax, f, prms)
    plot!(pl, ppcorr.xdata[1:ppcorr.Tmax], ppcorr.ydata[1:ppcorr.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end

fitpath = joinpath(anapath, "fits-AP")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(AntisymmetricCorrelator, apcorr)
for i in eachindex(apcorr.histories.tmin)
    tmin = apcorr.histories.tmin[i]
    prms = apcorr.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:apcorr.Tmax, f, prms)
    plot!(pl, apcorr.xdata[1:apcorr.Tmax], apcorr.ydata[1:apcorr.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end

# Fit ∂AP

derivate_sym_correlator!(apcorr)
apcorr.ydata .= -apcorr.ydata
LFTsAnaTools.correlator_fit_function(::Type{AntisymmetricCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)

reset_histories!(apcorr)
apcorr.tmin = 2
apcorr.tmax = 23
apcorr.Tmax = 25
dapmeff, dAAP = pion_fit_effective_mass_and_constant(apcorr)

fitpath = joinpath(anapath, "fits-dAP")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(AntisymmetricCorrelator, apcorr)
for i in eachindex(apcorr.histories.tmin)
    tmin = apcorr.histories.tmin[i]
    prms = apcorr.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:apcorr.Tmax, f, prms)
    plot!(pl, apcorr.xdata[2:apcorr.Tmax], apcorr.ydata[2:apcorr.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end

# Save to BDIO

bdpath = joinpath(anapath, "BDIO-mass-amplitude")
mkpath(bdpath)

write(joinpath(bdpath, "apmeff.bdio"), apmeff)
write(joinpath(bdpath, "ppmeff.bdio"), ppmeff)
write(joinpath(bdpath, "app.bdio"), APP)
write(joinpath(bdpath, "aap.bdio"), AAP)

# Plot 

meffpath = joinpath(anapath, "plots-meffs-amplitudes")
mkpath(meffpath)

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meffs")
plot!(pl, apmeff, label = "AP")
plot!(pl, ppmeff, label = "PP")
savefig(pl, joinpath(meffpath, "meffs-ap-pp.pdf"))
savefig(pl, joinpath(meffpath, "meffs-ap-pp.html"))

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "AAP")
plot!(pl, AAP, label = "AP")
savefig(pl, joinpath(meffpath, "aap.pdf"))
savefig(pl, joinpath(meffpath, "aap.html"))

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "APP")
plot!(pl, APP, label = "PP")
savefig(pl, joinpath(meffpath, "app.pdf"))
savefig(pl, joinpath(meffpath, "app.html"))


# Compute Fpi

fpipath = joinpath(anapath, "fpi")
mkpath(fpipath)

mpi = ppmeff.ydata[15]
T = ppcorr.T

Fpi = AAP.ydata[15]/sqrt(APP.ydata[15]) * sqrt(1/mpi*sinh(mpi*T/2))
uwerr(Fpi)

writedlm(joinpath(fpipath, "fpi.txt"), hcat(Fpi.mean, Fpi.err), ',')
write(joinpath(fpipath, "fpi.bdio"), Fpi)

Gpi = sqrt(APP.ydata[15]) * sqrt(2*mpi*sinh(mpi*T/2))
uwerr(Gpi)

writedlm(joinpath(fpipath, "gpi.txt"), hcat(Gpi.mean, Gpi.err), ',')
write(joinpath(fpipath, "fpi.bdio"), Gpi)

Fpi2 = dAAP.ydata[15]/sqrt(APP.ydata[15]) * sqrt(sinh(mpi*T/2))/mpi^(3/2)
uwerr(Fpi2)

writedlm(joinpath(fpipath, "fpi2.txt"), hcat(Fpi2.mean, Fpi2.err), ',')
write(joinpath(fpipath, "fpi2.bdio"), Fpi2)


# I think that if I overkill with statistics I don't need more sophisticated
# methods for fitting... just save mpip, APP and AAP to bdio files with plots,
# select best fit, then in a separate script compute Fpi -> save to BDIO


# derivate_sym_correlator!(apcorr)
# LFTsAnaTools.correlator_fit_function(::Type{AntisymmetricCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)

# plot(apcorr)

# LFTsAnaTools.correlator_fit_function(::Type{AntisymmetricCorrelator}, corrws::AbstractCorrelatorAnalysis) = APCorrelator2(T = corrws.T)

# import Utils: nparameters
# @Base.kwdef mutable struct APCorrelator2 <: UtilsFunc 
#     T
# end
# Utils.nparameters(s::APCorrelator2) = 2
# (s::APCorrelator2)(x,p) = p[1] * sinh( ppmeff.ydata[end].mean*(x-s.T/2) ) + p[2]*0.0

# reset_histories!(apcorr)
# apcorr.tmin = 1
# apcorr.tmax = 10
# apcorr.Tmax = 12
# apmeff, AAP = pion_fit_effective_mass_and_constant(apcorr)

# plot(AAP)

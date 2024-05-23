using Revise
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTU1
using LFTsAnaTools
using Plots
using DelimitedFiles
using ADerrors

plotlyjs()

wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/L64/Nfsim-b4.0-L64-m[-0.02,0.06]_D2024-04-30-13-10-42.216/measurements"

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/L24/Nfsim-b4.0-L24-m[-0.02, 0.06]_D2024-04-18-18-51-24.285/measurements"

# wdir = "/home/david/test/piontrash/"

length(ARGS) == 1 || error("Only one argument is expected! (Path to input directory)")
isdir(ARGS[1]) || error("Path provided is not a directory")
wdir = ARGS[1]


function parse_mass(filepath)
    masses = replace(filter(!isspace, filepath), r".*\-m(\[.*\]).*" => s"\1")
    return Meta.eval(Meta.parse(masses))
end

pi0 = CorrelatorAnalysis(wdir, type = NeutralPion, ensemble_ID = string(parse_mass(wdir)), burnout = 5, prefix = "ex")
pip = CorrelatorAnalysis(wdir, type = ChargedPion, ensemble_ID = string(parse_mass(wdir)), burnout = 5, prefix = "ex")

pi0.reweights = vec(readdlm(joinpath(wdir, "reweighting.txt")))[pi0.burnout+1:end]
pip.reweights = pi0.reweights
# pi0.reweights = nothing
# reset_histories!(pi0)

LFTsAnaTools.uwrealsym(pi0)
LFTsAnaTools.uwrealsym(pip)

# Compute effective masses
mpi0 = pion_effective_mass(pi0)
mpip = pion_effective_mass(pip)

reset_histories!(pi0)
pi0.tmin = 0
pi0.Tmax = 12
pi0.tmax = 10
mpi0fit = pion_fit_effective_mass(pi0)

reset_histories!(pip)
pip.tmin = 0
pip.Tmax = convert(Int64, pip.T/2)
pip.tmax = convert(Int64, pip.T/2-2)
mpipfit = pion_fit_effective_mass(pip)

uwerr.(pi0.ydata)
uwerr.(pip.ydata)
uwerr(mpi0)
uwerr(mpi0fit)
uwerr(mpip)
uwerr(mpipfit)


### Output analysis

anapath = joinpath(wdir, "analysis")
mkpath(anapath)

title = "L = $(pip.T)"

# Plot pi0 meffs


meffpi0path = joinpath(anapath, "meffs-pi0")
mkpath(meffpi0path)

mywrite(joinpath(meffpi0path, "meff-pi0.txt"), mpi0)
mywrite(joinpath(meffpi0path, "mefffit-pi0.txt"), mpi0fit)

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meffs pi0")
plot!(pl, mpi0, label = "cosh")
plot!(pl, mpi0fit, label = "fit")
savefig(pl, joinpath(meffpi0path, "meffs-pi0.pdf"))
savefig(pl, joinpath(meffpi0path, "meffs-pi0.html"))


# Plot pi+ meffs

meffpippath = joinpath(anapath, "meffs-pip")
mkpath(meffpippath)

mywrite(joinpath(meffpippath, "meff-pip.txt"), mpip)
mywrite(joinpath(meffpippath, "mefffit-pip.txt"), mpipfit)

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meffs pi+")
plot!(pl, mpip, label = "cosh")
plot!(pl, mpipfit, label = "fit")
savefig(pl, joinpath(meffpippath, "meffs-pip.pdf"))
savefig(pl, joinpath(meffpippath, "meffs-pip.html"))


# Plot pi+ pi0 fit meff

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meff")
plot!(pl, mpipfit, label="Mpi+")
plot!(pl, mpi0fit, label="Mpi0")
savefig(pl, joinpath(meffpi0path, "mefffits-pi0pip.pdf"))
savefig(pl, joinpath(meffpi0path, "mefffits-pi0pip.html"))


# Plot correlators

corrpath = joinpath(anapath, "correlators")
mkpath(corrpath)

basefname = "corr-pip-pi0"
txtpathpip = joinpath(corrpath,"corr-pip.txt")
txtpathpi0 = joinpath(corrpath,"corr-pi0.txt")
mywrite(txtpathpip, pip)
mywrite(txtpathpi0, pi0)

pl = plot()
plot!(pl, xlabel = "t", ylabel = "C(t)", title = title)
plot!(pl, pi0, label = "pi+")
plot!(pl, pip, label = "pi0")
savefig(pl, joinpath(corrpath, basefname*".pdf"))
savefig(pl, joinpath(corrpath, basefname*".html"))
plot!(pl, yscale = :log10)
savefig(pl, joinpath(corrpath, basefname*"-log10.pdf"))
savefig(pl, joinpath(corrpath, basefname*"-log10.html"))


# Plot fits

fitpath = joinpath(anapath, "fits")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(NeutralPion, pi0)
for tmin in pi0.histories.tmin
    prms = pi0.histories.fitp[:,tmin+1]
    pl = plot()
    plot!(pl, 0:pi0.Tmax, f, prms)
    plot!(pl, pi0.xdata[1:21], pi0.ydata[1:21], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end


# Plot connected and disconnected correlators

uwhist = Array{uwreal}(undef, 2, 2, 2, pi0.T)
for is in 1:2, ifl in 1:2, jfl in ifl:2, t in 1:pi0.T
    uwhist[is, ifl, jfl, t] = uwreal(pi0.history[is, ifl, jfl, :, t], pi0.ID)
end
ypi0conn = 1/2*(uwhist[1,1,1,:] + uwhist[1,2,2,:]) 
ypi0disc =  1/2*(uwhist[2,1,1,:] + uwhist[2,2,2,:]) - uwhist[2,1,2,:]

pl = plot(ypi0conn)
savefig(pl, joinpath(corrpath, "conn.pdf"))

pl = plot(ypi0conn, yaxis=:log10)
savefig(pl, joinpath(corrpath, "conn-log10.pdf"))

pl = plot(ypi0disc)
savefig(pl, joinpath(corrpath, "disc.pdf"))


# Plot ΔM vs t


dmpath = joinpath(anapath, "deltaM")
mkpath(dmpath)

function get_deltam_mbar(piws)
    ms = parse_mass(piws.ID)
    return abs(ms[1] - ms[2]), (ms[1] + ms[2])/2
end

deltam, mbar = get_deltam_mbar(pi0)

dMs = [mpip.ydata[i] - mpi0.ydata[i] for i in eachindex(mpi0.ydata)]
deltaMs = EffectiveMass(mpi0.xdata, dMs)

mywrite(joinpath(dmpath,"deltaM-eff.txt"), deltaMs)

pl = plot()
plot!(pl, deltaMs)
plot!(title = "dm=$(string(deltam))")
savefig(pl, joinpath(dmpath, "deltaM-eff.pdf"))
savefig(pl, joinpath(dmpath, "deltaM-eff.html"))



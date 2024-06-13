using Revise
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTU1
using LFTsAnaTools
using Plots
using DelimitedFiles
using ADerrors
using Utils
using BDIO
import LFTsAnaTools: correlator_fit_function

plotlyjs()

# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/L64/Nfsim-b4.0-L64-m[0.0,0.04]_D2024-04-30-13-07-06.72/measurements"
# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/L64/Nfsim-b4.0-L64-m[-0.04,0.08]_D2024-04-30-13-16-23.079/measurements"
# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/L64/Nfsim-b4.0-L64-m[-0.06,0.1]_D2024-04-30-13-18-05.106/measurements"
# wdir = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTU1.jl/results/2-non-degenerate/L64/Nfsim-b4.0-L64-m[-0.08,0.12]_D2024-04-30-13-20-42.039/measurements"


length(ARGS) == 1 || error("Only one argument is expected! (Path to input directory)")
isdir(ARGS[1]) || error("Path provided is not a directory")
wdir = ARGS[1]

function parse_mass(filepath)
    masses = replace(filter(!isspace, filepath), r".*\-m(\[.*\]).*" => s"\1")
    return Meta.eval(Meta.parse(masses))
end

LFTsAnaTools.correlator_fit_function(::Type{ChargedPion}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)
LFTsAnaTools.correlator_fit_function(type::Type{NeutralPion}, corrws) = SymCorrelator(T=corrws.T, s=1, c=0)
LFTsAnaTools.correlator_fit_function(type::Type{EtaPrime}, corrws) = SymCorrelator(T=corrws.T, s=1, c=0)
# LFTsAnaTools.correlator_fit_function(type::Type{QEtaPrime}, corrws) = SymCorrelator(T=corrws.T, s=1, c=0)


# Extract correlators

pi0 = CorrelatorAnalysis(wdir, type = NeutralPion, ensemble_ID = string(parse_mass(wdir)), burnout = 0, prefix = "ex")
pip = CorrelatorAnalysis(wdir, type = ChargedPion, ensemble_ID = string(parse_mass(wdir)), burnout = 0, prefix = "ex")
eta = CorrelatorAnalysis(wdir, type = EtaPrime, ensemble_ID = wdir, burnout = 0, prefix = "ex")

# Output analysis directory


anapath = joinpath(wdir, "analysis-$(size(pip.history)[1])")
mkpath(anapath)

title = "L = $(eta.T)"

# Reweight

pi0.reweights = vec(readdlm(joinpath(wdir, "reweighting.txt")))[pi0.burnout+1:end]
pip.reweights = pi0.reweights
eta.reweights = pi0.reweights

# Build uw correlators

LFTsAnaTools.uwrealsym(pi0)
LFTsAnaTools.uwrealsym(pip)
LFTsAnaTools.uwrealsym(eta)

uwerr.(pip.ydata)
uwerr.(pi0.ydata)
uwerr.(eta.ydata)

# Plot uw correlators

corrpath = joinpath(anapath, "correlators")
mkpath(corrpath)

basefname = "corr-pip-pi0"
txtpathpip = joinpath(corrpath,"corr-pip.txt")
txtpathpi0 = joinpath(corrpath,"corr-pi0.txt")
txtpatheta = joinpath(corrpath,"corr-eta.txt")
mywrite(txtpathpip, pip)
mywrite(txtpathpi0, pi0)
mywrite(txtpatheta, eta)

pl = plot()
plot!(pl, xlabel = "t", ylabel = "C(t)", title = title)
plot!(pl, pip, label = "pi+")
plot!(pl, pi0, label = "pi0")

savefig(pl, joinpath(corrpath, basefname*".pdf"))
savefig(pl, joinpath(corrpath, basefname*".html"))

plot!(pl, yscale = :log10)

savefig(pl, joinpath(corrpath, basefname*"-log10.pdf"))
savefig(pl, joinpath(corrpath, basefname*"-log10.html"))

basefname = "corr-eta"
pl = plot()
plot!(pl, xlabel = "t", ylabel = "C(t)", title = title)
plot!(pl, eta, label = "eta")
plot!(pl, yscale=:log10)

savefig(pl, joinpath(corrpath, basefname*".pdf"))
savefig(pl, joinpath(corrpath, basefname*".html"))


# Compute effective masses
mpi0 = pion_effective_mass(pi0)
mpip = pion_effective_mass(pip)
meta = pion_effective_mass(eta)
mqeta = pion_effective_mass(qeta)


# Fits

reset_histories!(pi0)
pi0.tmin = 0
pi0.Tmax = 16
pi0.tmax = 14
mpi0fit = pion_fit_effective_mass(pi0)

reset_histories!(pip)
pip.tmin = 0
pip.Tmax = convert(Int64, pip.T/2)
pip.tmax = convert(Int64, pip.T/2-2)
mpipfit = pion_fit_effective_mass(pip)

reset_histories!(eta)
eta.tmin = 0
eta.Tmax = 13
eta.tmax = 11
metafit = pion_fit_effective_mass(eta)

uwerr(mpip)
uwerr(mpipfit)

uwerr(mpi0)
uwerr(mpi0fit)

uwerr(meta)
uwerr(metafit)


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



# Plot fits pi0

fitpath = joinpath(anapath, "fits-pi0")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(NeutralPion, pi0)
for i in eachindex(pi0.histories.tmin)
    tmin = pi0.histories.tmin[i]
    prms = pi0.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:pi0.Tmax, f, prms)
    plot!(pl, pi0.xdata[1:pi0.Tmax], pi0.ydata[1:pi0.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end


# Plot fits eta

fitpath = joinpath(anapath, "fits-eta")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(EtaPrime, eta)
for i in eachindex(eta.histories.tmin)
    tmin = eta.histories.tmin[i]
    prms = eta.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:eta.Tmax, f, prms)
    plot!(pl, eta.xdata[1:eta.Tmax], eta.ydata[1:eta.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end

# Save masses to BDIO

bdpath = joinpath(anapath, "BDIO-masses")
mkpath(bdpath)

write(joinpath(bdpath, "mpip.bdio"), mpip)
write(joinpath(bdpath, "mpi0.bdio"), mpi0)
write(joinpath(bdpath, "meta.bdio"), meta)


# Compute effective masses of the derivative pi0

derivate_sym_correlator!(pi0)
pi0.ydata .= -pi0.ydata


LFTsAnaTools.correlator_fit_function(type::Type{NeutralPion}, corrws) = SymCorrelator(T=corrws.T, s=-1, c=0)

mdpi0 = pion_effective_mass(pi0)

reset_histories!(pi0)
pi0.tmin = 1
pi0.Tmax = 17
pi0.tmax = 15
mdpi0fit = pion_fit_effective_mass(pi0)

uwerr.(pi0.ydata)
uwerr.(pip.ydata)


# Plot fits dpi0

fitpath = joinpath(anapath, "fits-dpi0")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(NeutralPion, pi0)
for i in eachindex(pi0.histories.tmin)
    tmin = pi0.histories.tmin[i]
    prms = pi0.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:pi0.Tmax, f, prms)
    plot!(pl, pi0.xdata[1:pi0.Tmax], pi0.ydata[1:pi0.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end


# Plot pi0 meffs

if mpi0fit.ydata[3].mean < 0
    mpi0fit.ydata .= -mpi0fit.ydata
end

if mdpi0fit.ydata[3].mean < 0
    mdpi0fit.ydata .= -mdpi0fit.ydata
end

meffpi0path = joinpath(anapath, "meffs-pi0")
mkpath(meffpi0path)

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meffs pi0")
plot!(pl, mpi0, label = "cosh")
plot!(pl, mpi0fit, label = "fit")
plot!(pl, mdpi0, label = "cosh der.")
plot!(pl, mdpi0fit, label = "fit der.")

savefig(pl, joinpath(meffpi0path, "meffs-pi0.pdf"))
savefig(pl, joinpath(meffpi0path, "meffs-pi0.html"))

mywrite(joinpath(meffpi0path, "meff-pi0.txt"), mpi0)
mywrite(joinpath(meffpi0path, "mefffit-pi0.txt"), mpi0fit)
mywrite(joinpath(meffpi0path, "meff-dpi0.txt"), mdpi0)
mywrite(joinpath(meffpi0path, "mefffit-dpi0.txt"), mdpi0fit)



# Plot pi+ pi0 fit meff

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meff")
plot!(pl, mpipfit, label="Mpi+")
plot!(pl, mpi0fit, label="Mpi0")
savefig(pl, joinpath(meffpi0path, "mefffits-pi0pip.pdf"))
savefig(pl, joinpath(meffpi0path, "mefffits-pi0pip.html"))


# Plot ΔM vs t

dmpath = joinpath(anapath, "deltaM")
mkpath(dmpath)

function get_deltam_mbar(piws)
    ms = parse_mass(piws.ID)
    return abs(ms[1] - ms[2]), (ms[1] + ms[2])/2
end

deltam, mbar = get_deltam_mbar(pi0)

dMs = [mpip.ydata[22] - mpi0.ydata[i] for i in eachindex(mpi0.ydata)]
deltaMs = EffectiveMass(mpi0.xdata, dMs)

ddMs = [mpip.ydata[22] - mdpi0.ydata[i] for i in eachindex(mdpi0.ydata)]
deltadMs = EffectiveMass(mdpi0.xdata, ddMs)

dMsf = [mpip.ydata[22] - mpi0fit.ydata[i] for i in eachindex(mpi0fit.ydata)]
deltaMsf = EffectiveMass(mpi0fit.xdata, dMsf)

ddMsf = [mpip.ydata[22] - mdpi0fit.ydata[i] for i in eachindex(mdpi0fit.ydata)]
deltadMsf = EffectiveMass(mdpi0fit.xdata, ddMsf)

pl = plot()
plot!(pl, deltaMs, label="cosh")
plot!(pl, deltaMsf, label="fit")
plot!(pl, deltadMs, label="cosh der.")
plot!(pl,deltadMsf, label="fit der.")
plot!(title = "dm=$(string(deltam))")

savefig(pl, joinpath(dmpath, "deltaM-eff.pdf"))
savefig(pl, joinpath(dmpath, "deltaM-eff.html"))

mywrite(joinpath(dmpath,"deltaM-eff.txt"), deltaMs)
mywrite(joinpath(dmpath,"deltadM-eff.txt"), deltadMs)
mywrite(joinpath(dmpath,"deltaMf-eff.txt"), deltaMsf)
mywrite(joinpath(dmpath,"deltadMf-eff.txt"), deltadMsf)



# Eta derivative

derivate_sym_correlator!(eta)
eta.ydata .= .-eta.ydata

LFTsAnaTools.correlator_fit_function(type::Type{EtaPrime}, corrws) = SymCorrelator(T=corrws.T, s=-1, c=0)

mdeta = pion_effective_mass(eta)

reset_histories!(eta)
eta.tmin = 1
eta.Tmax = 14
eta.tmax = 12
mdetafit = pion_fit_effective_mass(eta)


# Eta derivative fit plots

fitpath = joinpath(anapath, "fits-deta")
mkpath(fitpath)

f = LFTsAnaTools.correlator_fit_function(EtaPrime, eta)
for i in eachindex(eta.histories.tmin)
    tmin = eta.histories.tmin[i]
    prms = eta.histories.fitp[:,i]
    pl = plot()
    plot!(pl, 0:eta.Tmax, f, prms)
    plot!(pl, eta.xdata[1:eta.Tmax], eta.ydata[1:eta.Tmax], yaxis=:log10)
    plot!(pl, title = "tmin = $tmin", xlabel = "t", ylabel = "C(t)")
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.html"))
    savefig(pl, joinpath(fitpath, "fit-tmin$tmin-log10.pdf"))
end


# Etamass plots

if metafit.ydata[3].mean < 0
    metafit.ydata .= -metafit.ydata
end

if mdetafit.ydata[3].mean < 0
    mdetafit.ydata .= -mdetafit.ydata
end

meffetapath = joinpath(anapath, "meffs-eta")
mkpath(meffetapath)

pl = plot(title = title)
plot!(pl, xlabel = "tmin", ylabel = "Meffs eta")
plot!(pl,meta, label="eta cosh")
plot!(pl,metafit, label="eta fit")
plot!(pl,mdeta, label="eta cosh der.")
plot!(pl,mdetafit, label="eta fit der.")
# plot!(pl,mqetafit, label="qeta fit")

savefig(pl, joinpath(meffetapath, "meffs-eta.pdf"))
savefig(pl, joinpath(meffetapath, "meffs-eta.html"))

mywrite(joinpath(meffetapath, "meff-eta.txt"), meta)
mywrite(joinpath(meffetapath, "mefffit-eta.txt"), metafit)
mywrite(joinpath(meffetapath, "meff-deta.txt"), mdeta)
mywrite(joinpath(meffetapath, "mefffit-qeta.txt"), mqetafit)


# Save masses to BDIO

bdpath = joinpath(anapath, "BDIO-masses")
mkpath(bdpath)

write(joinpath(bdpath, "mpipfit.bdio"), mpipfit)
write(joinpath(bdpath, "mpi0fit.bdio"), mpi0fit)
write(joinpath(bdpath, "metafit.bdio"), metafit)

write(joinpath(bdpath, "mdpi0.bdio"), mdpi0)
write(joinpath(bdpath, "mdeta.bdio"), mdeta)

write(joinpath(bdpath, "mdpi0fit.bdio"), mdpi0fit)
write(joinpath(bdpath, "mdetafit.bdio"), mdetafit)



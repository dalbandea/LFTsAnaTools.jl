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
LFTsAnaTools.correlator_fit_function(type::Type{EtaPrime}, corrws) = SymCorrelator(T=corrws.T, s=1, c=0)
# LFTsAnaTools.correlator_fit_function(type::Type{QEtaPrime}, corrws) = SymCorrelator(T=corrws.T, s=1, c=0)


# Extract correlators

pip = CorrelatorAnalysis(wdir, type = ChargedPion, ensemble_ID = string(parse_mass(wdir)), burnout = 0, prefix = "")
eta = CorrelatorAnalysis(wdir, type = EtaPrime, ensemble_ID = wdir, burnout = 0, prefix = "")

# Output analysis directory


anapath = joinpath(wdir, "analysis-$(size(pip.history)[1])")
mkpath(anapath)

title = "L = $(eta.T)"

# Build uw correlators

LFTsAnaTools.uwrealsym(pip)
LFTsAnaTools.uwrealsym(eta)

uwerr.(pip.ydata)
uwerr.(eta.ydata)

# Plot uw correlators

corrpath = joinpath(anapath, "correlators")
mkpath(corrpath)

basefname = "corrs"
txtpathpip = joinpath(corrpath,"corr-pip.txt")
txtpatheta = joinpath(corrpath,"corr-eta.txt")
mywrite(txtpathpip, pip)
mywrite(txtpatheta, eta)

pl = plot()
plot!(pl, xlabel = "t", ylabel = "C(t)", title = title)
plot!(pl, pip, label = "pi+")

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
mpip = pion_effective_mass(pip)
meta = pion_effective_mass(eta)
mqeta = pion_effective_mass(qeta)


# Fits

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
write(joinpath(bdpath, "meta.bdio"), meta)







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
write(joinpath(bdpath, "metafit.bdio"), metafit)

write(joinpath(bdpath, "mdeta.bdio"), mdeta)

write(joinpath(bdpath, "mdetafit.bdio"), mdetafit)




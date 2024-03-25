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
plotlyjs()

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


length(ARGS) == 1 || error("Only one argument is expected! (Path to input file)")
isfile(ARGS[1]) || error("Path provided is not a file")

if length(ARGS) == 1
    infile = ARGS[1]
else
    infile = "main/infile.in"
end
pdata = TOML.parsefile(infile)

burnout = pdata["Analysis parameters"]["burnout"]
safetmin = pdata["Analysis parameters"]["safetmin"]
wdir = pdata["Paths"]["wdir"]
ppfiles = readlines(pdata["Paths"]["pppath"])
apfiles = readlines(pdata["Paths"]["appath"])

m0 = parse_mass.(ppfiles)

ppdict = Dict(string.(m0) .=> ppfiles)
apdict = Dict(string.(m0) .=> apfiles)

ppwss = extract_correlators(ppdict, U1PionCorrelator, burnout)
apwss = extract_correlators(apdict, U1PCACCorrelator, burnout)

uwrealsym.(ppwss)
uwrealsym.(apwss)

mPCACs = compute_mpcac.(ppwss, apwss)
mpis = pion_effective_mass.(ppwss)
fit_mpis = pion_fit_effective_mass.(ppwss)

uwerr.(mPCACs)
uwerr.(mpis)
uwerr.(fit_mpis)


# Plot and save results

plot_and_save.(mPCACs)
plot_and_save.(mpis)
plot_and_save.(fit_mpis)

for i in eachindex(m0)
    htmlpath = joinpath(wdir, "jointeffmass")*"-m$(fit_mpis[i].m0).html"
    pl = plot(mpis[i], title = "m0=$(fit_mpis[i].m0)")
    plot!(pl, fit_mpis[i])
    savefig(pl, htmlpath)
end


# Perform analysis with safe tmin

defaultdir = joinpath(wdir, "defaults")
mkpath(defaultdir)

pickedmpcacs = [emass.ydata[safetmin] for emass in mPCACs]
pickedmpis = [emass.ydata[safetmin] for emass in fit_mpis]
pickedm0s = [emass.m0 for emass in mPCACs]

s = StraightLine()
fitps, cs, cse = fit(s, pickedm0s, pickedmpcacs, [1.0, 1.0])
uwerr.(fitps)
mc = -fitps[2] / fitps[1]
uwerr(mc)
fitx = -0.3:0.01:0.15
fity = [s(x, fitps) for x in fitx]
uwerr.(fity)
pl = plot(pickedm0s, pickedmpcacs, s, fitps, seriestype=:scatter, title="m_c = $(ADerrors.value(mc)) +/- $(ADerrors.err(mc))", xlabel=L"m_0", ylabel=L"m_{PCAC}", legend=false)
savefig(pl, joinpath(defaultdir, "critical-mass.pdf"))

mRs = ADerrors.value.([s(m, fitps) for m in pickedm0s])

open(joinpath(defaultdir,"mR-vs-mpi.txt"), "w") do io
    for i in eachindex(mRs)
        mywrite(io, mRs[i])
        write(io, ",")
        mywrite(io, pickedmpis[i])
        write(io, "\n")
    end
end


logtext = 
"""
Default tmin = $safetmin

Fit results:
- cs = $cs
- cse = $cse
- fitps[1] = $(fitps[1])
- fitps[2] = $(fitps[2])
- mc = $mc

Contents of mR-vs-mpi.txt: renormalized mass, value(pion mass), error(pion mass)
"""

mywrite(joinpath(defaultdir, "defaultlog.txt"), logtext)

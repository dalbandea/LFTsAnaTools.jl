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



length(ARGS) == 1 || error("Only one argument is expected! (Path to input file)")
isfile(ARGS[1]) || error("Path provided is not a file")

if length(ARGS) == 1
    infile = ARGS[1]
else
    infile = "main/infile.in"
end
pdata = TOML.parsefile(infile)

burnout = pdata["Analysis parameters"]["burnout"]
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

plotlyjs()

for i in eachindex(m0)
    htmlpath = joinpath(wdir, "jointeffmass")*"-m$(m0[i]).html"
    pl = plot(mpis[i], title = "m0=$(m0[i])")
    plot!(pl, fit_mpis[i])
    savefig(pl, htmlpath)
end

using Revise
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTU1
using LFTsAnaTools
using ADerrors
using LaTeXStrings
using Plots
using ArgParse

f64(s::String) = parse(Float64, s)

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

function compute_mpcac(ppws, apws)
    dapws = deepcopy(apws)
    derivate_sym_correlator!(dapws)
    mpcac = -1/2 * (dapws.ydata ./ ppws.ydata[2:end-1]) # not sure about the sign
    return EffectiveMass(ppws.xdata[2:end-1], mpcac, f64(ppws.ID), title="mpcac")
end

function pion_effective_mass(ppws::CorrelatorAnalysis)
    mpi = log.(ppws.ydata[1:end-1] ./ ppws.ydata[2:end])
    return EffectiveMass(ppws.xdata[1:end-1], mpi, f64(ppws.ID), title="Mpi-eff")
end

function pion_fit_effective_mass!(ppws::CorrelatorAnalysis)
    try
        tmin_loop(ppws, [1.0, 1.0])
    catch
    end
    return nothing
end

function pion_fit_effective_mass(ppws::CorrelatorAnalysis)
    pion_fit_effective_mass!(ppws)
    fit_mpis = ppws.histories.fitp[2,1:end-1]
    T = length(fit_mpis)
    return EffectiveMass(ppws.xdata[1:T], fit_mpis, f64(ppws.ID), title="Mpi-fiteff")
end

mutable struct EffectiveMass
    xdata
    ydata
    m0
    title
    function EffectiveMass(xdata, ydata, m0; title)
        return new(xdata, ydata, m0, title)
    end
end

ADerrors.uwerr(emass::EffectiveMass) = ADerrors.uwerr.(emass.ydata)
Plots.plot(emass::EffectiveMass; kwargs...) = Plots.plot(emass.xdata, emass.ydata; kwargs...)

mywrite(io, obs::T) where T <: Real  = write(io, "$(obs)")
mywrite(io, obs::uwreal)  = write(io, "$(ADerrors.value(obs)),$(ADerrors.err(obs))")
function mywrite(file, emass::EffectiveMass)
    open(file, "w") do io
        for i in eachindex(emass.ydata)
            mywrite(io, emass.xdata[i])
            write(io, ",")
            mywrite(io, emass.ydata[i])
            write(io, "\n")
        end
    end
end

function plot_and_save(emass::EffectiveMass)
    title = emass.title
    txtpath = joinpath(wdir, title)*"-m$(emass.m0).txt"
    pdfpath = joinpath(wdir, title)*"-m$(emass.m0).pdf"
    mywrite(txtpath, emass)
    pl = plot(emass, title=title*" for m0=$(emass.m0)", seriestype=:scatter)
    savefig(pl, pdfpath)
    return nothing
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--burnout"
        help = "number of configurations to skip in the analysis"
        required = false
        arg_type = Int
        default = 1

        "--pppath"
        help = "path to directory to save configurations and logs"
        required = true
        arg_type = String

        "--appath"
        help = "path to directory to save configurations and logs"
        required = true
        arg_type = String

        "--wdir"
        help = "path to directory to save configurations and logs"
        required = true
        arg_type = String
    end

    return parse_args(s)
end


parsed_args = parse_commandline()

burnout = parsed_args["burnout"]
wdir = parsed_args["wdir"]
ppfiles = readlines(parsed_args["pppath"])
apfiles = readlines(parsed_args["appath"])

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


mutable struct EffectiveMass
    xdata
    ydata
    m0
    title
    function EffectiveMass(xdata, ydata, m0 = nothing; title = "")
        return new(xdata, ydata, m0, title)
    end
end
export EffectiveMass


f64(s::String) = parse(Float64, s)
export f64


ADerrors.uwerr(emass::EffectiveMass) = ADerrors.uwerr.(emass.ydata)
Plots.plot(emass::EffectiveMass; kwargs...) = Plots.plot(emass.xdata, emass.ydata; seriestype=:scatter, label = emass.title, xlabel="t",  kwargs...)
Plots.plot!(pl::Plots.Plot, emass::EffectiveMass; kwargs...) = Plots.plot!(pl, emass.xdata, emass.ydata; seriestype=:scatter, label = emass.title, kwargs...)

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
function mywrite(file, message::String)
    open(file, "w") do io
        write(io, message)
    end
end

export mywrite


function compute_mpcac(ppws, apws)
    dapws = deepcopy(apws)
    derivate_sym_correlator!(dapws)
    mpcac = -1/2 * (dapws.ydata ./ ppws.ydata[2:end-1]) # not sure about the sign
    return EffectiveMass(ppws.xdata[2:end-1], mpcac, f64(ppws.ID), title="mpcac")
end
export compute_mpcac

function pion_effective_mass(ppws::CorrelatorAnalysis)
    mpi = @. acosh((ppws.ydata[1:end-2] + ppws.ydata[3:end])/(2*ppws.ydata[2:end-1]))
    return EffectiveMass(ppws.xdata[2:end-1], mpi, f64(ppws.ID), title="Mpi-eff")
end
export pion_effective_mass

function pion_fit_effective_mass!(ppws::CorrelatorAnalysis)
    try
        tmin_loop(ppws, [1.0, 1.0])
    catch
    end
    return nothing
end
export pion_fit_effective_mass!

function pion_fit_effective_mass(ppws::CorrelatorAnalysis)
    pion_fit_effective_mass!(ppws)
    fit_mpis = ppws.histories.fitp[2,1:end-1]
    T = length(fit_mpis)
    return EffectiveMass(ppws.xdata[1:T], fit_mpis, f64(ppws.ID), title="Mpi-fiteff")
end
export pion_fit_effective_mass

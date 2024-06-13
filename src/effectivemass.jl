
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
# Plots.plot(emass::EffectiveMass; kwargs...) = Plots.plot(emass.xdata, emass.ydata; seriestype=:scatter, label = emass.title, xlabel="t",  kwargs...)
# Plots.plot!(pl::Plots.Plot, emass::EffectiveMass; kwargs...) = Plots.plot!(pl, emass.xdata, emass.ydata; seriestype=:scatter, label = emass.title, kwargs...)

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
function mywrite(file, corrws::CorrelatorAnalysis)
    open(file, "w") do io
        for i in eachindex(corrws.ydata)
            mywrite(io, corrws.xdata[i])
            write(io, ",")
            mywrite(io, corrws.ydata[i])
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


function pion_deriv_effective_mass(ppws::CorrelatorAnalysis)
    mpi = Vector{uwreal}(undef, length(ppws.ydata)-2)
    for i in 1:length(mpi)
        try
            mpi[i] = asinh((ppws.ydata[i] + ppws.ydata[i+2])/(2*ppws.ydata[i+1]))
        catch
            mpi[i] = uwreal([0.0, 0.0], "zero")
        end
    end
    return EffectiveMass(ppws.xdata[2:end-1], mpi, ppws.ID, title="Mpi-eff")
end
export pion_deriv_effective_mass

function pion_effective_mass(ppws::CorrelatorAnalysis)
    mpi = Vector{uwreal}(undef, length(ppws.ydata)-2)
    for i in 1:length(mpi)
        try
            mpi[i] = acosh((ppws.ydata[i] + ppws.ydata[i+2])/(2*ppws.ydata[i+1]))
        catch
            mpi[i] = uwreal([0.0, 0.0], "zero")
        end
    end
    return EffectiveMass(ppws.xdata[2:end-1], mpi, ppws.ID, title="Mpi-eff")
end
export pion_effective_mass

function pion_fit_effective_mass!(ppws::CorrelatorAnalysis)
    try
        prms = ones(nparameters(correlator_fit_function(ppws)))
        show(correlator_fit_function(ppws))
        tmin_loop(ppws, prms)
    catch
    end
    return nothing
end
export pion_fit_effective_mass!

function pion_fit_effective_mass(ppws::CorrelatorAnalysis)
    tmin = ppws.tmin
    tmax = ppws.tmax
    tmin >= ppws.xdata[1] || error("tmin=$tmin outside of range $(ppws.xdata[1]):$(ppws.xdata[end])")
    tmax <= ppws.xdata[end] || error("tmax=$tmax outside of range $(ppws.xdata[1]):$(ppws.xdata[end])")
    pion_fit_effective_mass!(ppws)
    fit_mpis = ppws.histories.fitp[2,:]
    length(fit_mpis) == length(tmin:tmax) || error("Some fits failed")
    return EffectiveMass(ppws.xdata[tmin+1:tmax+1], fit_mpis, ppws.ID, title="Mpi-fiteff")
end
export pion_fit_effective_mass


#####################
# EffectiveMass I/O #
#####################


import BDIO: BDIO_write!, BDIO_read, BDIO_open, BDIO_seek!, BDIO_close!
import Base: read, write

BDIO.BDIO_write!(fb::BDIO.BDIOstream, uwval::uwreal) = write_uwreal(uwval, fb, 8)
function BDIO.BDIO_write!(fb::BDIO.BDIOstream, uwv::Vector{uwreal})
    for i in eachindex(uwv)
        BDIO_write!(fb, uwv[i])
    end
    return nothing
end
function BDIO.BDIO_read(fb::BDIO.BDIOstream, uwv::Vector{uwreal})
    for i in eachindex(uwv)
        uwv[i] = read_uwreal(fb)
        BDIO_seek!(fb)
        BDIO_seek!(fb)
    end
end

function BDIO.BDIO_write!(fb::BDIO.BDIOstream, m::EffectiveMass)
    N = length(m.xdata)
    BDIO.BDIO_start_record!(fb, BDIO.BDIO_BIN_GENERIC, 1)
    BDIO.BDIO_write!(fb, [convert(Int32, N)])
    BDIO.BDIO_write!(fb, collect(m.xdata))
    BDIO.BDIO_write_hash!(fb)
    BDIO.BDIO_write!(fb, m.ydata)
    return nothing
end

function write(bdfile::String, m::EffectiveMass)
    fb = BDIO.BDIO_open(bdfile, "d", "Test file")
    BDIO.BDIO_write!(fb, m)
    BDIO_close!(fb)
end

function read_effmass(bdfile::String)
    fb = BDIO_open(bdfile, "r")
    BDIO_seek!(fb)
    reg = zeros(Int32, 1)
    BDIO_read(fb, reg)
    L = reg[1]
    reg = zeros(Int64, L)
    BDIO_read(fb, reg)
    BDIO_seek!(fb)
    BDIO_seek!(fb)
    xdata = copy(reg)
    reg = Vector{uwreal}(undef, L)
    BDIO_read(fb, reg)
    BDIO_close!(fb)
    return EffectiveMass(xdata, reg)
end
export read_effmass

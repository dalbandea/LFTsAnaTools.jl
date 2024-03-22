module LFTsAnaTools

using ADerrors
import ADerrors: uwreal

using Plots
import Plots: plot
using Utils

include("analysis.jl")
export AbstractCorrelatorAnalysis, AbstractObservableAnalysis, CorrelatorSymmetry, SymmetricCorrelator, AntisymmetricCorrelator
export CorrelatorAnalysis

include("plots.jl")
export plot, plot!  

end # module

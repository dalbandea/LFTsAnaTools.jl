module LFTsAnaTools

using ADerrors
import ADerrors: uwreal

using Plots
using Utils

include("analysis.jl")
export AbstractCorrelatorAnalysis, AbstractObservableAnalysis, CorrelatorSymmetry, SymmetricCorrelator, AntisymmetricCorrelator
export CorrelatorAnalysis

include("effectivemass.jl")

include("extensions/LFTU1ext.jl")

include("plots.jl")

end # module

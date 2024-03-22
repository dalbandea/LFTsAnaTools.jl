module LFTU1Ext

using LFTSampling
using LFTU1
using LFTsAnaTools
using Utils

import LFTsAnaTools: sign, correlator_fit_function

LFTsAnaTools.sign(::Type{T}) where T <: AbstractCorrelator = sign(get_symmetry(T))

get_symmetry(::Type{U1PionCorrelator}) = SymmetricCorrelator
correlator_fit_function(::Type{U1PionCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)


get_symmetry(::Type{U1PCACCorrelator}) = AntisymmetricCorrelator
correlator_fit_function(::Type{U1PCACCorrelator}, corrws::AbstractCorrelatorAnalysis) = SymCorrelator(T=corrws.T, s=1, c=0)

export get_symmetry

end

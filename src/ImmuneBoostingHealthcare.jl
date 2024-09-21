
module ImmuneBoostingHealthcare

using DrWatson
using CategoricalArrays, DataFrames, Distributions, StaticArrays, StatsBase, Turing
import Base: minimum

include("structs.jl")
include("consts.jl")
include("simulations.jl")
include("processdata.jl")
include("analysis.jl")
include("plotting.jl")

export 
    # structs.jl
    AbstractParameters, HCWSEIIRRRp, SEIIRRRSp, WXYYZSEIIRRRSp,
    # consts.jl
    
    # simulations.jl
    betahh, betahp, betaph, betapp, stochasticseiirrrs, stochasticwxyyzseiirrrs,
    # processdata.jl
    insertproportions!,
    # analysis.jl
    calculatebetah, calculatebetahs, calculatebetap, calculatebetaps, calculatebetas, 
    calculatelambdac, calculatelambdacs, countdates, counthospitals, datamatrices, fitmodel, 
    hcwseiirrr, hcwseiirrr_isolating, hcwseiirrr_isolating!, hospitalconditionmatrices, 
    loadchainsdf, predictdiagnoses, predicttotaldiagnoses
    # plotting.jld2
    

end  # module ImmuneBoostingHealthcare

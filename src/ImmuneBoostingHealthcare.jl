
module ImmuneBoostingHealthcare

using DrWatson
using CategoricalArrays, DataFrames, Dates, Distributions, StaticArrays, StatsBase, Turing
import Base: minimum

include("structs.jl")
include("consts.jl")
include("simulations.jl")
include("processdata.jl")
include("analysis.jl")

export 
    # `structs.jl`
    AbstractParameters, HCWSEIIRRRp, SEIIRRRSp, WXYYZSEIIRRRSp,
    # no exports from `consts.jl`
    # `simulations.jl`
    vaccinatestaff,
    # `processdata.jl`
    datetot, insertproportions!, ttodate,
    # `analysis.jl`
    calculatebetah, calculatebetahs, calculatebetap, calculatebetaps, calculatebetas, 
    calculatelambdac, calculatelambdacs, countdates, counthospitals, datamatrices, 
    findobservationssincejuly, fitmodel, fitmodelperhospital, hcwseiirrr, 
    hcwseiirrr_isolating, hcwseiirrr_isolating!, hospitalconditionmatrices, loadchainsdf, 
    loadchainsperhospitaldf, predictdiagnoses, predicttotaldiagnoses, processoutputs, 
    processoutputsdict, processoutputsperhospital, producecounterfactualoutputsdict

end  # module ImmuneBoostingHealthcare

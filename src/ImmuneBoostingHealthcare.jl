
module ImmuneBoostingHealthcare

using DrWatson
using CategoricalArrays, DataFrames, Dates, Distributions, StaticArrays, StatsBase, Turing
import Base: minimum

include("model.jl")

export NegativeCompartmentSizeError

export HCWSEIRRRVCPOutput, HCWSEIRRRVOutput, HCWSEIRRRParameters
export addsusceptible!, advancetime!, becomeimmunefromvaccine!, boostR2!, boostR3!, calculateboosting, calculatelambda, calculatelambdaprime_c, calculatelambdaprime_h, calculatelambdaprime_p, countnewlydiagnosed, diagnose!, gettime, gettotaldiagnosed, gettotalimmune, get_n, get_E, get_I, get_R1, get_R2, get_R3, get_S, get_v, expose!, exposevaccinated!, progress!, recover!, runmodel, runmodel!, vaccinateR2!, vaccinateR3!, vaccinateS!, wane1!, wane2!, wane3!

#=
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
=#
end  # module ImmuneBoostingHealthcare

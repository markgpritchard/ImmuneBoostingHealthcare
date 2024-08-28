
module ImmuneBoostingHealthcare

using DrWatson
using CategoricalArrays, DataFrames, DifferentialEquations, Distributions, StaticArrays
using StatsBase, Turing
import Base: minimum

include("structs.jl")
include("consts.jl")
include("simulations.jl")
include("processdata.jl")
include("analysis.jl")
include("plotting.jl")

export 
    # structs.jl
    SEIIRRRSp, WXYYZSEIIRRRSp,
    modifyp,
    # consts.jl
    
    # simulations.jl
    seiirrrs!, wxyyzseiirrrs!, 
    makeu0!, seiirrrs_u0, simulatehospital, simulatehospitals, wxyyzseiirrrs_u0,
    # processdata.jl
    insertproportions!,
    # analysis.jl
    calculatebetas, countdates, counthospitals, datamatrices, fitmodel, 
    hospitalconditionmatrices, loadchainsdf, predictinfections, predictinfections!, 
    summarizepredictedinfections
    # plotting.jld2
    

end  # module ImmuneBoostingHealthcare

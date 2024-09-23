
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include("analysesimssetup.jl")
include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions

## Unboosted sims 

unboostedoutputs = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsdf("fittedvalues_unboostedsimulation")
    unboostedpsi0df = loadchainsdf("fittedvalues_psi_0_unboostedsimulation", 0)
    processoutputs(
        unboostedsimulation, finaldata, unboosteddf, unboostedpsi0df, vaccinated; 
        dateid=:t
    )    
end

plotchains(unboostedoutputs[:chaindf])
plotchains(unboostedoutputs[:chaindf_psi0])

unboostedtotalsfig = plotoutputs(unboostedoutputs)



## Boosted sims 

boostedoutputs = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsdf("fittedvalues_boostedsimulation")
    boostedpsi0df = loadchainsdf("fittedvalues_psi_0_boostedsimulation", 0)
    processoutputs(
        boostedsimulation, finaldata, boosteddf, boostedpsi0df, vaccinated; 
        dateid=:t
    )    
end

plotchains(boostedoutputs[:chaindf])
plotchains(boostedoutputs[:chaindf_psi0])

boostedtotalsfig = plotoutputs(boostedoutputs)

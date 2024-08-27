
include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare

#using Turing, Pigeons, Random
using Pigeons
using StatsBase

if isfile(datadir("sims", "simulations.jld2"))
    simulations = load(datadir("sims", "simulations.jld2"))
else 
    include("generatesimulations.jl")
end

using CairoMakie 


## Unboosted sims 

## Boosted sims 

## Covid data 

# numbers vaccinated currently simulated 
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for t ∈ axes(vaccinated, 1), j ∈ axes(vaccinated, 2)
        if t ∈ [ 300, 450, 650, 800 ]
            vaccinated[t, j] = 0.8
        end
    end
    vaccinated
end

nhospitals = counthospitals(coviddata)
ndates = countdates(coviddata)

@unpack newstaff, patients, staff = datamatrices(coviddata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(coviddata)

stringency = coviddata.StringencyIndex_Average[1:ndates]
community = coviddata.weeklycases[1:ndates] ./ 56_000_000

coviddf = loadchainsdf("fittedvalues_coviddata")
plotchains(coviddf)

@unpack βp, βh, βc = calculatebetas(coviddf, vpd, psb, stringency)

function summarizepredictedinfections(predictedinfections::Array{<:Real, 3})
    ndates, nhospitals, nsamples = size(predictedinfections)      
    totals = zeros(nsamples, nhospitals)
    means = zeros(nhospitals)
    lcris = zeros(nhospitals)
    ucris = zeros(nhospitals)

    for i ∈ axes(totals, 2), j ∈ axes(totals, 1)
        totals[j, i] = sum(@view predictedinfections[:, i, j])
    end

    for i ∈ 1:nhospitals
        means[i] = mean(totals[:, i])
        lcri, ucri = quantile(predictedtotalinfections[:, i], [ 0.05, 0.95 ])
        lcris[i] = lcri
        ucris[i] = ucri
    end

    return @ntuple totals means lcris ucris
end

function summarizepredictedinfections(args...; kwargs...)
    predictedinfections = predictinfections(args...; kwargs...) 
    return summarizepredictedinfections(predictedinfections)
end

predictedinfections = summarizepredictedinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency
)

predictedinfectionswithoutboost = summarizepredictedinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency; ψ=0
)

fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
scatter!(ax1, totalinfections, predictedinfections.means; color=:blue, markersize=3)
#rangebars!(ax, totalinfections, lcripredictedtotalinfections, ucripredictedtotalinfections)
lines!(ax1, [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(ax2, totalinfections, predictedinfectionswithoutboost.means; color=:blue, markersize=3)
#rangebars!(ax, totalinfections, lcripredictedtotalinfections, ucripredictedtotalinfections)
lines!(ax2, [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])


fig

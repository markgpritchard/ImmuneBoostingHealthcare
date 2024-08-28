
include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare
include(srcdir("PlottingFunctions.jl"))

#using Turing, Pigeons, Random
#using CairoMakie, Pigeons, StatsBase
using CairoMakie, Pigeons
using .PlottingFunctions

if isfile(datadir("sims", "simulations.jld2"))
    simulations = load(datadir("sims", "simulations.jld2"))
else 
    include("generatesimulations.jl")
end

## Unboosted sims 

## Boosted sims 

## Covid data 

nhospitals = counthospitals(coviddata)
ndates = countdates(coviddata)

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

@unpack newstaff, patients, staff = datamatrices(coviddata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(coviddata)

stringency = coviddata.StringencyIndex_Average[1:ndates]
community = coviddata.weeklycases[1:ndates] ./ 56_000_000

coviddf = loadchainsdf("fittedvalues_coviddata")
plotchains(coviddf)

@unpack βp, βh, βc = calculatebetas(coviddf, vpd, psb, stringency)

totalinfections = [ sum(@view newstaff[:, i]) for i ∈ 1:nhospitals ]

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
rangebars!(ax1, totalinfections, predictedinfections.lcris, predictedinfections.ucris; color=( :blue, 0.1 ))
lines!(ax1, [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(ax2, totalinfections, predictedinfectionswithoutboost.means; color=:blue, markersize=3)
#rangebars!(ax, totalinfections, lcripredictedtotalinfections, ucripredictedtotalinfections)
lines!(ax2, [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])


fig



exampleinfections = predictinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency;
    inds=10
)[:, :, 1]


exampleinfections .- newstaff[:, 10]
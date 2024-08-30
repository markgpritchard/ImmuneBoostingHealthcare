
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

#unboostedsimulation = simulations["unboostedsimulation"]
@unpack unboostedsimulation, allbetas = simulations

nhospitals = counthospitals(unboostedsimulation)
ndates = countdates(unboostedsimulation; dateid=:t)

# simulated numbers vaccinated  
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for t ∈ axes(vaccinated, 1), j ∈ axes(vaccinated, 2)
        vaccinated[300, j] = 0.8
    end
    vaccinated
end

@unpack newstaff, patients, staff = datamatrices(unboostedsimulation, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(unboostedsimulation)

stringency = coviddata.StringencyIndex_Average[1:ndates]
community = unboostedsimulation.CommunityCases[1:ndates] ./ 56_000_000

unboosteddf = loadchainsdf("fittedvalues_unboostedsimulation")
plotchains(unboosteddf)

@unpack βp, βh, βc = calculatebetas(unboosteddf, vpd, psb, stringency)

totalinfections = [ sum(@view newstaff[:, i]) for i ∈ 1:nhospitals ]

predictedinfections = summarizepredictedinfections(
    unboosteddf, patients, staff, community, vaccinated, vpd, psb, stringency
)

predictedinfectionswithoutboost = summarizepredictedinfections(
    unboosteddf, patients, staff, community, vaccinated, vpd, psb, stringency; ψ=0
)

medianbetas = let
    @unpack βp, βh, βc = calculatebetas(unboosteddf, vpd, psb, stringency)
    βclcris = zeros(ndates)
    βcmedians = zeros(ndates)
    βcucris = zeros(ndates)
    for i ∈ 1:ndates 
        vals = quantile([ b[i] for b ∈ βc ], [ 0.05, 0.5, 0.95 ])
        βclcris[i] = vals[1]
        βcmedians[i] = vals[2]
        βcucris[i] = vals[3]
    end
    βhlcris = zeros(nhospitals)
    βhmedians = zeros(nhospitals)
    βhucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βh ], [ 0.05, 0.5, 0.95 ])
        βhlcris[i] = vals[1]
        βhmedians[i] = vals[2]
        βhucris[i] = vals[3]
    end
    βplcris = zeros(nhospitals)
    βpmedians = zeros(nhospitals)
    βpucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βp ], [ 0.05, 0.5, 0.95 ])
        βplcris[i] = vals[1]
        βpmedians[i] = vals[2]
        βpucris[i] = vals[3]
    end
    @ntuple βclcris βcmedians βcucris βhlcris βhmedians βhucris βplcris βpmedians βpucris
end

fig = Figure()
axs1 = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
axs2 = [ Axis(fig[i, 2]) for i ∈ 1:3 ]
scatter!(axs1[1], totalinfections, predictedinfections.medians; color=:blue, markersize=3)
rangebars!(
    axs1[1], totalinfections, predictedinfections.lcris, predictedinfections.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[1], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[2], totalinfections, predictedinfectionswithoutboost.medians; 
    color=:blue, markersize=3,
)
rangebars!(
    axs1[2], 
    totalinfections, 
    predictedinfectionswithoutboost.lcris,
    predictedinfectionswithoutboost.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[2], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[3], 
    totalinfections, 
    (predictedinfections.medians .- predictedinfectionswithoutboost.medians) ./ 
        predictedinfections.medians; 
    color=:blue, markersize=3,
)
#scatter!(ax3, totalinfections, predictedinfections.means .- predictedinfectionswithoutboost.means; color=:blue, markersize=3)

scatter!(axs2[2], [ b[1] for b ∈ allbetas ], medianbetas[:βhmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[2], 
    [ b[1] for b ∈ allbetas ], 
    medianbetas[:βhlcris],
    medianbetas[:βhucris]; 
    color=( :blue, 0.1 ),
)
lines!(axs2[2], [ extrema([ b[1] for b ∈ allbetas ])... ], [ extrema([ b[1] for b ∈ allbetas ])... ])

scatter!(axs2[3], [ b[2] for b ∈ allbetas ], medianbetas[:βpmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[3], 
    [ b[2] for b ∈ allbetas ], 
    medianbetas[:βplcris],
    medianbetas[:βpucris]; 
    color=( :blue, 0.1 ),
)
lines!(axs2[3], [ extrema([ b[2] for b ∈ allbetas ])... ], [ extrema([ b[2] for b ∈ allbetas ])... ])


linkxaxes!(axs1...)
linkyaxes!(axs1[1], axs1[2])

fig



## Boosted sims 

boostedsimulation = simulations["boostedsimulation"]

## Covid data 


nhospitals = counthospitals(finaldata)
ndates = countdates(finaldata)



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

@unpack newstaff, patients, staff = datamatrices(finaldata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(finaldata)

stringency = finaldata.StringencyIndex_Average[1:ndates]
community = finaldata.weeklycases[1:ndates] ./ 56_000_000

coviddf = loadchainsdf("fittedvalues_coviddata")
plotchains(coviddf)

@unpack βp, βh, βc = calculatebetas(coviddf, vpd, psb, stringency)

totalinfections = [ sum(@view newstaff[:, i]) for i ∈ 1:nhospitals ]

predictedinfections = summarizepredictedinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency
)

predictedinfectionswithoutboost = summarizepredictedinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency; 
    ψ=0,
)

medianbetas = let
    @unpack βp, βh, βc = calculatebetas(coviddf, vpd, psb, stringency)
    βclcris = zeros(ndates)
    βcmedians = zeros(ndates)
    βcucris = zeros(ndates)
    for i ∈ 1:ndates 
        vals = quantile([ b[i] for b ∈ βc ], [ 0.05, 0.5, 0.95 ])
        βclcris[i] = vals[1]
        βcmedians[i] = vals[2]
        βcucris[i] = vals[3]
    end
    βhlcris = zeros(nhospitals)
    βhmedians = zeros(nhospitals)
    βhucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βh ], [ 0.05, 0.5, 0.95 ])
        βhlcris[i] = vals[1]
        βhmedians[i] = vals[2]
        βhucris[i] = vals[3]
    end
    βplcris = zeros(nhospitals)
    βpmedians = zeros(nhospitals)
    βpucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βp ], [ 0.05, 0.5, 0.95 ])
        βplcris[i] = vals[1]
        βpmedians[i] = vals[2]
        βpucris[i] = vals[3]
    end
    @ntuple βclcris βcmedians βcucris βhlcris βhmedians βhucris βplcris βpmedians βpucris
end

fig = Figure()
axs1 = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
axs2 = [ Axis(fig[i, 2]) for i ∈ 1:3 ]
scatter!(axs1[1], totalinfections, predictedinfections.medians; color=:blue, markersize=3)
rangebars!(
    axs1[1], totalinfections, predictedinfections.lcris, predictedinfections.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[1], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[2], totalinfections, predictedinfectionswithoutboost.medians; 
    color=:blue, markersize=3,
)
rangebars!(
    axs1[2], 
    totalinfections, 
    predictedinfectionswithoutboost.lcris,
    predictedinfectionswithoutboost.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[2], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[3], 
    totalinfections, 
    (predictedinfections.medians .- predictedinfectionswithoutboost.medians) ./ 
        predictedinfections.medians; 
    color=:blue, markersize=3,
)
#scatter!(ax3, totalinfections, predictedinfections.means .- predictedinfectionswithoutboost.means; color=:blue, markersize=3)

scatter!(axs2[1], coviddata.Date[1:832], medianbetas[:βcmedians]; color=:blue, markersize=3)

#scatter!(axs2[2], ordinalrank(medianbetas[:βhmedians]), medianbetas[:βhmedians]; color=:blue, markersize=3)
scatter!(axs2[2], totalinfections, medianbetas[:βhmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[2], 
    totalinfections,
    #ordinalrank(medianbetas[:βhmedians]),
    medianbetas[:βhlcris],
    medianbetas[:βhucris]; 
    color=( :blue, 0.1 ),
)

#scatter!(axs2[3], ordinalrank(medianbetas[:βpmedians]), medianbetas[:βpmedians]; color=:blue, markersize=3)
scatter!(axs2[3], totalinfections, medianbetas[:βpmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[3], 
    totalinfections,
    #ordinalrank(medianbetas[:βpmedians]),
    medianbetas[:βplcris],
    medianbetas[:βpucris]; 
    color=( :blue, 0.1 ),
)


linkxaxes!(axs1...)
linkyaxes!(axs1[1], axs1[2])

fig



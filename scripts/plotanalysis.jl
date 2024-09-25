
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot from simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("analysesimssetup.jl")

## Unboosted sims 

unboostedoutputs = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsdf("fittedvalues_unboostedsimulation")
    processoutputs(unboostedsimulation, finaldata, unboosteddf, vaccinated; dateid=:t)    
end

plotchains(unboostedoutputs[:chaindf])

unboostedtotalsfig = plotoutputs(unboostedoutputs)

#=
plothospitaloutputs(unboostedoutputs)
plothospitaloutputs(unboostedoutputs; firstplot=26)
plothospitaloutputs(unboostedoutputs; firstplot=51)
plothospitaloutputs(unboostedoutputs; firstplot=76)
=#

unboostedoutputs_omega180 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsdf(
        "fittedvalues_unboostedsimulation_omega_0.00556"; 
        omega=0.00556
    )
    processoutputs(unboostedsimulation, finaldata, unboosteddf, vaccinated; dateid=:t)    
end

plotchains(unboostedoutputs_omega180[:chaindf])

unboostedtotalsfig = plotoutputs(unboostedoutputs_omega180)

## Boosted sims 

boostedoutputs = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsdf("fittedvalues_boostedsimulation")
    processoutputs(boostedsimulation, finaldata, boosteddf, vaccinated; dateid=:t)    
end

plotchains(boostedoutputs[:chaindf])

boostedtotalsfig = plotoutputs(boostedoutputs)

plothospitaloutputs(boostedoutputs)
plothospitaloutputs(boostedoutputs; firstplot=26)
plothospitaloutputs(boostedoutputs; firstplot=51)
plothospitaloutputs(boostedoutputs; firstplot=76)

boostedoutputs_omega180 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsdf(
        "fittedvalues_boostedsimulation_omega_0.00556"; 
        omega=0.00556
    )
    processoutputs(boostedsimulation, finaldata, boosteddf, vaccinated; dateid=:t)    
end

plotchains(boostedoutputs_omega180[:chaindf])

boostedtotalsfig = plotoutputs(boostedoutputs_omega180)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots from data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if isfile(datadir("exp_pro", "finaldata.jld2"))
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
else 
    include("loaddata.jl")
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
end

include("analysedatasetup.jl")

dataoutputs = let 
    datadf = loadchainsdf("fittedvalues_coviddata")
    processoutputs(finaldata, datadf, vaccinated)    
end

plotchains(dataoutputs[:chaindf])
#plotchains(dataoutputs[:chaindf_psi0])

scatter(dataoutputs[:totalinfections], vpd)
scatter(dataoutputs[:totalinfections], psb)


datatotalsfig = plotoutputs(dataoutputs)

plothospitaloutputs(dataoutputs)
plothospitaloutputs(dataoutputs; firstplot=26)
plothospitaloutputs(dataoutputs; firstplot=51)
plothospitaloutputs(dataoutputs; firstplot=76)
plothospitaloutputs(dataoutputs; firstplot=101)
plothospitaloutputs(dataoutputs; firstplot=126)

dataoutputs_omega180 = let 
    datadf = loadchainsdf("fittedvalues_coviddata_omega_0.00556"; omega=0.00556)
    processoutputs(finaldata, datadf, vaccinated)  
end

plotchains(dataoutputs_omega180[:chaindf])

datatotals_omega180fig = plotoutputs(dataoutputs_omega180)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prior samples 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

priorsamples = sample(
    fitmodel(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    ), 
    Prior(), 
    1000
)

priorsamples_psi0 = sample(
    fitmodel(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals;
        psiprior=0
    ), 
    Prior(), 
    1000
)

priorsamplefig, priorsampleparameterfig = let 
    priorsamplesdf = DataFrame(priorsamples)
    priorsamplespredicteddiagnoses = predicttotaldiagnoses(
        priorsamplesdf, 
        patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
    priorsamplesbetas = calculatebetas(priorsamplesdf, vpd, psb, nhospitals)
    priorsampleslambdacs = calculatelambdacs(priorsamplesdf, stringency, community) 

    priorsamplesdf_psi0 = DataFrame(priorsamples_psi0)
    priorsamplespredicteddiagnoses_psi0 = predicttotaldiagnoses(
        priorsamplesdf_psi0, 0, 
        patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
    priorsamplesbetas_psi0 = calculatebetas(priorsamplesdf_psi0, vpd, psb, nhospitals)
    priorsampleslambdacs_psi0 = calculatelambdacs(priorsamplesdf_psi0, stringency, community) 

    priorsamplefig = Figure()
    axs = [ Axis(priorsamplefig[i, 1]) for i ∈ 1:3 ]
    scatter!(
        axs[1], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses.mediantotaldiagnoses; 
        color=:blue, markersize=3
    )
    rangebars!(
        axs[1], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses.lcitotaldiagnoses, 
        priorsamplespredicteddiagnoses.ucitotaldiagnoses; 
        color=( :blue, 0.1 ),
    )
    lines!(
        axs[1], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ]
    )

    scatter!(
        axs[2], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses_psi0.mediantotaldiagnoses; 
        color=:blue, markersize=3
    )
    rangebars!(
        axs[2], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses_psi0.lcitotaldiagnoses, 
        priorsamplespredicteddiagnoses_psi0.ucitotaldiagnoses; 
        color=( :blue, 0.1 ),
    )
    lines!(
        axs[2], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ]
    )

    scatter!(
        axs[3], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        (
            priorsamplespredicteddiagnoses.mediantotaldiagnoses .- 
            priorsamplespredicteddiagnoses_psi0.mediantotaldiagnoses
        ) ./ priorsamplespredicteddiagnoses.mediantotaldiagnoses; 
        color=:blue, markersize=3,
    )

    ( priorsamplefig , 1 )
end

## Unboosted sims 

#unboostedsimulation = simulations["unboostedsimulation"]

@unpack unboostedsimulation = simulations

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

stringency = finaldata.StringencyIndex_Average[1:ndates]
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

scatter!(axs2[1], finaldata.Date[1:832], medianbetas[:βcmedians]; color=:blue, markersize=3)

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



##############################################################################################

# hospital RCD #25

using StaticArrays

idnumber = 80
id = unique(finaldata.StringCodes)[idnumber]

h20data = filter(:StringCodes => x -> x == id, finaldata)

@unpack βp, βh, βc = calculatebetas(
    coviddf, 
    h20data.VolumePerBed[1], 
    h20data.ProportionSingleBeds[1], 
    h20data.StringencyIndex_Average    
)

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, h20data.CovidAbsences; color=:black, markersize=3)
for i ∈ axes(coviddf, 1)
    predictedinfections = zeros(ndates) 
    immunevector = SizedVector{3}(zeros(Float64, 3)) 
    for t ∈ 2:ndates 
        foi = βc[i][t-1] * community[t-1] + βp[i][1] * h20data.PatientsProportion[t-1] + βh[i][1] * predictedinfections[t-1]
        v = vaccinated[(t - 1), idnumber]
        immune10 = predictedinfections[t-1] + v * (1 - sum(immunevector) - predictedinfections[t-1])
        # probability of boosting from natural immune boosting plus vaccination
        pb = (1 - exp(-coviddf.ψ[i] * foi)) * (1 - v) + v  
        for x ∈ 1:2
            immune10 += pb * immunevector[x]
            immunevector[x] += -(pb + 2 * coviddf.ω[i]) * immunevector[x] + 
                2 * coviddf.ω[i] * (1 - pb) * immunevector[x+1]
        end
        immunevector[2] += -(2 * coviddf.ω[i] * (1 - pb)) * 
            immunevector[2] + 
            immune10
        predictedinfections[t] = (1 - sum(immunevector)) * (1 - exp(-foi))
    end
    lines!(ax, predictedinfections; color=( :blue, 0.1 ))
end



fig




####

fig, ax = lines([ 1 - exp(-x * 0) for x ∈ 0:0.01:2 ])
lines!(ax, [ (1 - (1 - exp(-x)))^0 for x ∈ 0:0.01:2 ])


fig
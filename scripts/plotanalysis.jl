
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot from simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

include("analysesimssetup.jl")

#=
plothospitaloutputs(unboostedoutputs)
plothospitaloutputs(unboostedoutputs; firstplot=26)
plothospitaloutputs(unboostedoutputs; firstplot=51)
plothospitaloutputs(unboostedoutputs; firstplot=76)
=#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations without natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# no hospital-specific parameters, unboosted immunity lasts 180 days

unboostedoutputs_omega180 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsdf(
        "fittedvalues_unboostedsimulation_omega_0.00556"; 
        omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(unboostedsimulation, finaldata, unboosteddf, vaccinated; dateid=:t)    
end

plotchains(unboostedoutputs_omega180["chaindf"])

unboostedtotals_omega180fig = plotoutputs(unboostedoutputs_omega180)


# no hospital-specific parameters, unboosted immunity lasts 100 days

unboostedoutputs_omega100 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsdf("fittedvalues_unboostedsimulation_omega_0.01"; omega=0.01)
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(unboostedsimulation, finaldata, unboosteddf, vaccinated; dateid=:t)    
end

plotchains(unboostedoutputs_omega100["chaindf"])

unboostedtotals_omega100fig = plotoutputs(unboostedoutputs_omega100)


# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 180 days

jseries = [
    82, 38, 23, 25, 74, 6, 9, 75, 57, 65, 7, 59, 
    81, 31, 3, 50, 10, 94, 15, 42, 19, 8, 47, 21, 27
]

unboostedoutputsperhospital_omega180 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsperhospitaldf(
        "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        unboostedsimulation, finaldata, unboosteddf, vaccinated, jseries; 
        dateid=:t
    )    
end

plotchains(unboostedoutputsperhospital_omega180["chaindf"]; size=( 400, 4800 ))

unboostedtotalsperhospital_omega180fig = plotoutputs(unboostedoutputsperhospital_omega180)


# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 100 days

unboostedoutputsperhospital_omega100 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsperhospitaldf(
        "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01"; 
        jseries, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        unboostedsimulation, finaldata, unboosteddf, vaccinated, jseries; 
        dateid=:t
    )    
end

plotchains(unboostedoutputsperhospital_omega100["chaindf"]; size=( 400, 4800 ))

unboostedtotalsperhospital_omega100fig = plotoutputs(unboostedoutputsperhospital_omega100)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 180 days

jseriesall = 1:nhospitals

unboostedoutputsperhospitalall_omega180 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsperhospitaldf(
        "fittedvalues_unboostedsimulationperhospital_omega_0.00556"; 
        jseries=jseriesall, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        unboostedsimulation, finaldata, unboosteddf, vaccinated, jseriesall; 
        dateid=:t
    )    
end

plotchains(unboostedoutputsperhospitalall_omega180["chaindf"]; size=( 400, 4800 ))

unboostedtotalsperhospitalall_omega180fig = plotoutputs(
    unboostedoutputsperhospitalall_omega180
)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 100 days

unboostedoutputsperhospitalall_omega100 = let 
    @unpack unboostedsimulation = simulations
    unboosteddf = loadchainsperhospitaldf(
        "fittedvalues_unboostedsimulationperhospital_omega_0.01"; 
        jseries=jseriesall, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        unboostedsimulation, finaldata, unboosteddf, vaccinated, jseriesall; 
        dateid=:t
    )    
end

plotchains(unboostedoutputsperhospitalall_omega100["chaindf"]; size=( 400, 4800 ))

unboostedtotalsperhospitalall_omega100fig = plotoutputs(
    unboostedoutputsperhospitalall_omega100
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# no hospital-specific parameters, unboosted immunity lasts 180 days

boostedoutputs_omega180 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsdf(
        "fittedvalues_boostedsimulation_omega_0.00556"; 
        omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(boostedsimulation, finaldata, boosteddf, vaccinated; dateid=:t)    
end

plotchains(boostedoutputs_omega180["chaindf"])

boostedtotals_omega180fig = plotoutputs(boostedoutputs_omega180)


# no hospital-specific parameters, unboosted immunity lasts 100 days

boostedoutputs_omega100 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsdf("fittedvalues_boostedsimulation_omega_0.01"; omega=0.01)
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(boostedsimulation, finaldata, boosteddf, vaccinated; dateid=:t)    
end

plotchains(boostedoutputs_omega100["chaindf"])

boostedtotals_omega100fig = plotoutputs(boostedoutputs_omega100)


# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 180 days

boostedoutputsperhospital_omega180 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsperhospitaldf(
        "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        boostedsimulation, finaldata, boosteddf, vaccinated, jseries; 
        dateid=:t
    )    
end

plotchains(boostedoutputsperhospital_omega180["chaindf"]; size=( 400, 4800 ))

boostedtotalsperhospital_omega180fig = plotoutputs(boostedoutputsperhospital_omega180)


# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 100 days

boostedoutputsperhospital_omega100 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsperhospitaldf(
        "fittedvalues_boostedsimulationperhospital_subset_omega_0.01"; 
        jseries, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        boostedsimulation, finaldata, boosteddf, vaccinated, jseries; 
        dateid=:t
    )    
end

plotchains(boostedoutputsperhospital_omega100["chaindf"]; size=( 400, 4800 ))

boostedtotalsperhospital_omega100fig = plotoutputs(boostedoutputsperhospital_omega100)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 180 days

boostedoutputsperhospitalall_omega180 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsperhospitaldf(
        "fittedvalues_boostedsimulationperhospital_omega_0.00556"; 
        jseries=jseriesall, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        boostedsimulation, finaldata, boosteddf, vaccinated, jseriesall; 
        dateid=:t
    )    
end

plotchains(boostedoutputsperhospitalall_omega180["chaindf"]; size=( 400, 4800 ))

boostedtotalsperhospitalall_omega180fig = plotoutputs(
    boostedoutputsperhospitalall_omega180
)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 100 days

boostedoutputsperhospitalall_omega100 = let 
    @unpack boostedsimulation = simulations
    boosteddf = loadchainsperhospitaldf(
        "fittedvalues_boostedsimulationperhospital_omega_0.01"; 
        jseries=jseriesall, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        boostedsimulation, finaldata, boosteddf, vaccinated, jseriesall; 
        dateid=:t
    )    
end

plotchains(boostedoutputsperhospitalall_omega100["chaindf"]; size=( 400, 4800 ))

boostedtotalsperhospitalall_omega100fig = plotoutputs(
    boostedoutputsperhospitalall_omega100
)


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

# no hospital-specific parameters, unboosted immunity lasts 180 days

dataoutputs_omega180 = let 
    datadf = loadchainsdf("fittedvalues_coviddata_omega_0.00556"; omega=0.00556)
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(finaldata, datadf, vaccinated)    
end

plotchains(dataoutputs_omega180["chaindf"])

dataoutputs_omega180fig = plotoutputs(dataoutputs_omega180)


# no hospital-specific parameters, unboosted immunity lasts 100 days

dataoutputs_omega100 = let 
    datadf = loadchainsdf("fittedvalues_coviddata_omega_0.01"; omega=0.01)
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(finaldata, datadf, vaccinated)    
end

plotchains(dataoutputs_omega100["chaindf"])

dataoutputs_omega100fig = plotoutputs(dataoutputs_omega100)


# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 180 days

dataoutputsperhospital_omega180 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseries)    
end

plotchains(dataoutputsperhospital_omega180["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospital_omega180fig = plotoutputs(dataoutputsperhospital_omega180)


# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 100 days

dataoutputsperhospital_omega100 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.01"; 
        jseries, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(
        finaldata, datadf, vaccinated, jseries; 
        dateid=:t
    )    
end

plotchains(dataoutputsperhospital_omega100["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospital_omega100fig = plotoutputs(dataoutputsperhospital_omega100)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 180 days

dataoutputsperhospitalall_omega180 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_omega_0.00556"; 
        jseries=jseriesall, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseriesall)    
end

plotchains(dataoutputsperhospitalall_omega180["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospitalall_omega180fig = plotoutputs(dataoutputsperhospitalall_omega180)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 100 days

bdataoutputsperhospitalall_omega100 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_omega_0.01"; 
        jseries=jseriesall, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseriesall)    
end

plotchains(dataoutputsperhospitalall_omega100["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospitalall_omega100fig = plotoutputs(dataoutputsperhospitalall_omega100)


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
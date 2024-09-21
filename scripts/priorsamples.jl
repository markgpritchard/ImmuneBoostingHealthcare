
using DrWatson

@quickactivate :ImmuneBoostingHealthcare
using CairoMakie

include("analysedatasetup.jl")

priorsamples = sample(
    fitmodel(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    ), 
    Prior(), 
    1000
)
priorsamplesdf = DataFrame(priorsamples)
priorsamplespredicteddiagnoses = predicttotaldiagnoses(
    priorsamplesdf, 
    patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
)
priorsamplesbetas = calculatebetas(priorsamplesdf, vpd, psb, nhospitals)
priorsampleslambdacs = calculatelambdacs(priorsamplesdf, stringency, community) 

priorsamples_psi0 = sample(
    fitmodel(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals;
        psiprior=0
    ), 
    Prior(), 
    1000
)
priorsamplesdf_psi0 = DataFrame(priorsamples_psi0)
priorsamplespredicteddiagnoses_psi0 = predicttotaldiagnoses(
    priorsamplesdf_psi0, 0, 
    patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
)
priorsamplesbetas_psi0 = calculatebetas(priorsamplesdf_psi0, vpd, psb, nhospitals)
priorsampleslambdacs_psi0 = calculatelambdacs(priorsamplesdf_psi0, stringency, community) 

fig = Figure()
axs1 = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
axs2 = [ Axis(fig[i, 2]) for i ∈ 1:3 ]
axs3 = [ Axis(fig[i, 3]) for i ∈ 1:3 ]
scatter!(
    axs1[1], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplespredicteddiagnoses.mediantotaldiagnoses; 
    color=:blue, markersize=3
)
rangebars!(
    axs1[1], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplespredicteddiagnoses.lcitotaldiagnoses, 
    priorsamplespredicteddiagnoses.ucitotaldiagnoses; 
    color=( :blue, 0.1 ),
)
lines!(
    axs1[1], 
    [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ], 
    [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ]
)

scatter!(
    axs1[2], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplespredicteddiagnoses_psi0.mediantotaldiagnoses; 
    color=:blue, markersize=3
)
rangebars!(
    axs1[2], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplespredicteddiagnoses_psi0.lcitotaldiagnoses, 
    priorsamplespredicteddiagnoses_psi0.ucitotaldiagnoses; 
    color=( :blue, 0.1 ),
)
lines!(
    axs1[2], 
    [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ], 
    [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ]
)

scatter!(
    axs1[3], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    (
        priorsamplespredicteddiagnoses.mediantotaldiagnoses .- 
        priorsamplespredicteddiagnoses_psi0.mediantotaldiagnoses
    ) ./ priorsamplespredicteddiagnoses.mediantotaldiagnoses; 
    color=:blue, markersize=3,
)

scatter!(
    axs2[1], finaldata.Date[1:832], priorsampleslambdacs.medianlambdac; 
    color=:blue, markersize=3
)
rangebars!(
    axs2[1], 
    finaldata.Date[1:832], 
    priorsampleslambdacs.lcilambdac, 
    priorsampleslambdacs.ucilambdac; 
    color=( :blue, 0.1 ),
)

scatter!(
    axs2[2], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplesbetas.medianbetah; 
    color=:blue, markersize=3
)
rangebars!(
    axs2[2], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplesbetas.lcibetah, 
    priorsamplesbetas.ucibetah; 
    color=( :blue, 0.1 ),
)

scatter!(
    axs2[3], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplesbetas.medianbetap; 
    color=:blue, markersize=3
)
rangebars!(
    axs2[3], 
    [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
    priorsamplesbetas.lcibetap, 
    priorsamplesbetas.ucibetap; 
    color=( :blue, 0.1 ),
)

linkxaxes!(axs1...)
linkyaxes!(axs1[1], axs1[2])

fig

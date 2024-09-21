
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include("analysedatasetup.jl")

id = 6 ## temporarily 

#=
priorsamples = sample(
    fitmodel(
        newstaff, patients, staff, vaccinated, community, 
        vpd, psb, stringency, ndates, nhospitals
    ), 
    Prior(), 
    10_000
)
=#
#=
logprobinitialconditions = let
    α1 = 0.1
    α2 = 0.0
    α3 = 0.0
    α4 = 0.1
    α5 = 0.0
    α6 = 0.0
    α7 = 0.1
    α8 = 0.1
    ω = 0.02 
    ψ = 1.0
    sigma2 = 1.0
    
    lp = 0.0 

    T = typeof(α1)

    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community
    
    for j ∈ 1:nhospitals
        predictedinfections = 0.0
        r1 = zero(T)
        r2 = zero(T)
        r3 = zero(T)
        λp = max(zero(T), α1 + α2 * vpd[j] + α3 * psb[j]) .* patients[:, j]
        βh = max(zero(T), α4 + α5 * vpd[j] + α6 * psb[j]) 
        for t ∈ 2:ndates 
            foi = λc[t] + λp[j] + βh * predictedinfections
            ξ = 1 - exp(-ψ * foi)
            predictedinfections = (1 - predictedinfections - r1 - r2 - r3) * (1 - exp(-foi))
            newr1 = predictedinfections + vaccinated[(t - 1), j] + ξ * (r2 + r3) + (1 - (1 - ξ) * 3 * ω) * r1 
            newr2 = (1 - ξ) * 3 * ω * r1 + (1 - (1 - ξ) * 3 * ω - ξ) * r2
            newr3 = (1 - ξ) * 3 * ω * r2 + (1 - (1 - ξ) * 3 * ω - ξ) * r3 
            r1 = newr1 
            r2 = newr2 
            r3 = newr3
            lp += logpdf(Normal(newstaff[t, j], sigma2), predictedinfections)
        end
    end

    lp
end
=#

# This code sets `Pigeons.initialization` to ensure that the first set of parameters used
# give a valid output

function fitdatamodel_target(
    newstaff=newstaff, 
    patients=patients, 
    staff=staff, 
    vaccinated=vaccinated, 
    community=community, 
    vpd=vpd, 
    psb=psb, 
    stringency=stringency, 
    ndates=ndates, 
    nhospitals=nhospitals
)
    return Pigeons.TuringLogPotential(
        fitmodel(
            newstaff, patients, staff, vaccinated, community, 
            vpd, psb, stringency, ndates, nhospitals
        )
    )
end

const FitdatamodelType = typeof(fitdatamodel_target())

function Pigeons.initialization(target::FitdatamodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :α1, 1, 0.1)
    Pigeons.update_state!(result, :α2, 1, 0.0)
    Pigeons.update_state!(result, :α3, 1, 0.0)
    Pigeons.update_state!(result, :α4, 1, 0.1)
    Pigeons.update_state!(result, :α5, 1, 0.0)
    Pigeons.update_state!(result, :α6, 1, 0.0)
    Pigeons.update_state!(result, :α7, 1, 0.1)
    Pigeons.update_state!(result, :α8, 1, 0.1)
    Pigeons.update_state!(result, :ω, 1, 0.02)
    Pigeons.update_state!(result, :ψ, 1, 1.0)
    Pigeons.update_state!(result, :sigma2, 1, 1.0)

    return result
end

fitted_pt = pigeons( ;
    target=fitdatamodel_target(
        newstaff, patients, staff, vaccinated, community, 
        vpd, psb, stringency, ndates, nhospitals
    ),
    n_rounds=0,
    n_chains=5,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(id),
    variational=GaussianReference(),
)

new_pt = fitted_pt

for i ∈ 1:n_rounds
    filename = "fittedvalues_coviddata_id_$(id)_round_$(i).jld2"
    nextfilename = "fittedvalues_coviddata_id_$(id)_round_$(i + 1).jld2"
    isfile(datadir("sims", nextfilename)) && continue
    if isfile(datadir("sims", filename))
        global new_pt = load(datadir("sims", filename))["pt"]
    else
        pt = increment_n_rounds!(new_pt, 1)
        global new_pt = pigeons(pt)
        new_chains = Chains(new_pt)
        resultdict = Dict(
            "chain" => new_chains, 
            "pt" => new_pt, 
            "n_rounds" => i, 
            "n_chains" => 5,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end

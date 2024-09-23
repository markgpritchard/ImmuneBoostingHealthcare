
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include("analysedatasetup.jl")

function fitdatamodel_target(
    patients=patients, 
    staff=staff, 
    vaccinated=vaccinated, 
    community=community, 
    vpd=vpd, 
    psb=psb, 
    stringency=stringency, 
    ndates=ndates, 
    nhospitals=nhospitals;
    psiprior=0,
)
    return Pigeons.TuringLogPotential(
        fitmodel(
            patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals;
            psiprior=0
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
    Pigeons.update_state!(result, :α8, 1, 0.0)
    Pigeons.update_state!(result, :ω, 1, 0.02)
    Pigeons.update_state!(result, :θ, 1, 2/7)
    Pigeons.update_state!(result, :sigma2, 1, 0.05)

    return result
end

fitted_pt = pigeons( ;
    target=fitdatamodel_target(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals;
        psiprior=0,
    ),
    n_rounds=0,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(id),
    variational=GaussianReference(),
)

new_pt = fitted_pt

for i ∈ 1:n_rounds
    filename = "fittedvalues_psi_0_coviddata_id_$(id)_round_$(i).jld2"
    nextfilename = "fittedvalues_psi_0_coviddata_id_$(id)_round_$(i + 1).jld2"
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
            "n_chains" => 4,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end
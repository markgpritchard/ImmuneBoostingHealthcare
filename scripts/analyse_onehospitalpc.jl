
using DrWatson
#=
@quickactivate "ImmuneBoostingHealthcare"

using Pkg 
Pkg.instantiate()
=#
@quickactivate :ImmuneBoostingHealthcare


using Distributions, Pigeons, Random, Turing

include("loaddata.jl")
include("analysesimssetup.jl")

alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))

datarbk = filter(:StringCodes => x -> x == "RBK", finaldata)
insertcols!(datarbk, :Community => community)

λc = 0.02 .* community

patients = let 
    pts = zeros(832)
    prev = 0.0
    for (t, p) ∈ enumerate(datarbk.PatientsProportion) 
        if prev == 0 
            pts[t] = p 
            prev = p
        elseif p > 5 * prev || p < 0.2 * prev
            pts[t] = pts[t-1]
        else 
            pts[t] = p 
            prev = p
        end 
    end
    pts
end

staff = let 
    stf = zeros(832)
    prev = 0.0
    for (t, p) ∈ enumerate(datarbk.StaffProportion) 
        if prev == 0 
            stf[t] = p 
            prev = p
        elseif p > 5 * prev || p < 0.2 * prev
            stf[t] = stf[t-1]
        else 
            stf[t] = p 
            prev = p
        end 
    end
    stf
end

u0 = zeros(16)
u0[1] = 1


@model function fitmodelonehospital( 
    patients, staff, vaccinated, community, stringency, ndates;
    betahprior=truncated(Exponential(1), 0, 10),
    betapprior=truncated(Exponential(1), 0, 10),
    alpha7prior=truncated(Normal(0.2, 1), -1, 10),
    alpha8prior=truncated(Normal(0.1, 0.1), 0, 10),  # require greater stringency leads to less transmission
    omegaprior=truncated(Exponential(0.02), 0, 0.33),
    psiprior=truncated(Exponential(1), 0, 1000),
    thetaprior=Beta(1, 1),
    sigma2prior=Exponential(1),
)
    βh ~ betahprior
    βp ~ betapprior
    α7 ~ alpha7prior
    α8 ~ alpha8prior
     
    if omegaprior isa Number 
        ω = omegaprior
    else
        ω ~ omegaprior
    end

    if psiprior isa Number 
        ψ = psiprior 
    else
        ψ ~ psiprior
    end

    θ ~ thetaprior
    sigma2 ~ sigma2prior

    T = typeof(α7)
   # u0 = zeros(16)
   # u0[1] = 1.0
    isolating = zeros(Float64, ndates)

    λc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ] .* community
    
    p = HCWSEIIRRRp(
        βh, 
        βp, 
        0.5, 
        0.2, 
        ψ, 
        ω, 
        θ
    )
        #println("$p")
    hcwseiirrr_isolating!(isolating, ImmuneBoostingHealthcare.automatic, p, 1:ndates, λc, patients, vaccinated, 1)

    for t ∈ 1:ndates
        if isnan(isolating[t])
            Turing.@addlogprob! -Inf
        else
            Turing.@addlogprob! logpdf(Normal(isolating[t], sigma2), staff[t])
        end
    end
end


function fitdatamodel_target(
    patients=patients, 
    staff=staff, 
    vaccinated=vaccinated, 
    community=community, 
    stringency=stringency, 
    ndates=832;
    omegaprior=truncated(Exponential(0.02), 0, 0.33),
)
    return Pigeons.TuringLogPotential(
        fitmodelonehospital(
            patients, staff, vaccinated, community, stringency, ndates;
            omegaprior
        )
    )
end

const FitdatamodelType = typeof(fitdatamodel_target())

function Pigeons.initialization(target::FitdatamodelType, rng::AbstractRNG, ::Int64)
    result = DynamicPPL.VarInfo(
        rng, target.model, DynamicPPL.SampleFromPrior(), DynamicPPL.PriorContext()
    )
    DynamicPPL.link!!(result, DynamicPPL.SampleFromPrior(), target.model)

    Pigeons.update_state!(result, :βh, 1, 0.075)
    Pigeons.update_state!(result, :βp, 1, 0.075)
    Pigeons.update_state!(result, :α7, 1, 0.1)
    Pigeons.update_state!(result, :α8, 1, 0.0)
    Pigeons.update_state!(result, :ψ, 1, 1.0)
    Pigeons.update_state!(result, :θ, 1, 2/7)
    Pigeons.update_state!(result, :sigma2, 1, 1.0)
    return result
end

fittedvalues_onehospital = load(datadir("sims", "fittedvalues_onehospital_id_1_round_5.jld2"))
fittedvalues_onehospitaldf = DataFrame(fittedvalues_onehospital["chain"])


fitted_pt_omega001 = pigeons( ;
    target=fitdatamodel_target(
        patients, staff, vaccinated, community, stringency, 832; 
        omegaprior=ω,
    ),
    n_rounds=0,
    n_chains=4,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(id),
    variational=GaussianReference(),
)

new_pt = fitted_pt_omega001

for id ∈ 1:4
    for i ∈ 0:12
        filename = "fittedvalues_onehospital_omega_001_id_$(id)_round_$(i)nn.jld2"
        nextfilename = "fittedvalues_onehospital_omega_001_id_$(id)_round_$(i + 1)nn.jld2"
        isfile(datadir("sims", nextfilename)) && continue
        if isfile(datadir("sims", filename))
            global new_pt = load(datadir("sims", filename))["pt"]
        elseif i == 0 
            fitted_pt_omega001 = pigeons( ;
                target=fitdatamodel_target(
                    patients, staff, vaccinated, community, stringency, 832; 
                    omegaprior=ω,
                ),
                n_rounds=0,
                n_chains=4,
                multithreaded=true,
                record=[ traces; record_default() ],
                seed=(id),
                variational=GaussianReference(),
            )
            global new_pt = fitted_pt_omega001
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
end

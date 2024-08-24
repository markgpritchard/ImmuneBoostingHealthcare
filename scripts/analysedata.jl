
include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare

using Turing, Pigeons, Random

if length(ARGS) == 2 
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
else
    id = 1 
    n_rounds = 4
end

nhospitals = counthospitals(coviddata)
ndates = countdates(coviddata)

@unpack newstaff, patients, staff = datamatrices(coviddata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(coviddata)

stringency = coviddata.StringencyIndex_Average[1:ndates]
community = coviddata.weeklycases[1:ndates] ./ 56_000_000

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

function fitdatamodel_target(
    newstaff=newstaff, patients=patients, staff=staff, vaccinated=vaccinated, community=community, 
    vpd=vpd, psb=psb, stringency=stringency, ndates=ndates, nhospitals=nhospitals
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
    Pigeons.update_state!(result, :α2, 1, 0.1)
    Pigeons.update_state!(result, :α3, 1, 0.1)
    Pigeons.update_state!(result, :α4, 1, 0.1)
    Pigeons.update_state!(result, :α5, 1, 0.1)
    Pigeons.update_state!(result, :α6, 1, 0.1)
    Pigeons.update_state!(result, :α7, 1, 0.1)
    Pigeons.update_state!(result, :α8, 1, 0.1)
    Pigeons.update_state!(result, :ω, 1, 0.02)
    Pigeons.update_state!(result, :ψ, 1, 1.0)
    Pigeons.update_state!(result, :sigma2, 1, 0.05)

    return result
end

fitted_pt = pigeons( ;
    target=fitdatamodel_target(),
    n_rounds=0,
    n_chains=10,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(id),
    variational=GaussianReference(),
)

new_pt = fitted_pt

##########################################################################################
#=
using StaticArrays

immunevectorlength = 10

α1 = 0.1
α2 = 0.1
α3 = 0.1
α4 = 0.1
α5 = 0.1
α6 = 0.1
α7 = 0.1
α8 = 0.1
ω = 0.02
ψ = 1.0
sigma2 = 0.05

T = typeof(α1)

    # levels of immunity
    immunevector = SizedVector{immunevectorlength}(zeros(T, immunevectorlength))  

    βp = [ max(zero(T), α1 + α2 * v + α3 * p) for (v, p) ∈ zip(vpd, psb) ]
    βh = [ max(zero(T), α4 + α5 * v + α6 * p) for (v, p) ∈ zip(vpd, psb) ]
    βc = [ max(zero(T), α7 + α8 * (100 - s)) for s ∈ stringency ]

    foi = zeros(T, ndates, nhospitals)
    for t ∈ axes(foi, 1), j ∈ axes(foi, 2)
        foi[t, j] = βp[j] * patients[t, j] + βh[j]* staff[t, j] + βc[t] * community[t]
    end

    predictedinfections = zeros(T, ndates, nhospitals)
    ImmuneBoostingHealthcare._predictinfections!(
        predictedinfections, immunevector, ndates, nhospitals, foi, vaccinated, ψ, ω
    )

    lp = 0
    for t ∈ axes(predictedinfections, 1), j ∈ axes(predictedinfections, 2)
        t <= 14 && continue 
        println(logpdf(Normal(newstaff[t, j], sigma2), predictedinfections[t, j]))
        lp += logpdf(Normal(newstaff[t, j], sigma2), predictedinfections[t, j])
    end

=#
##########################################################################################

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
            "n_chains" => 10,
        )
        safesave(datadir("sims", filename), resultdict)
    end
end


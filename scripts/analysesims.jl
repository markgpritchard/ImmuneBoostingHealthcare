
include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare

using Turing, Pigeons

if isfile(datadir("sims", "simulations.jld2"))
    simulations = load(datadir("sims", "simulations.jld2"))
else 
    include("generatesimulations.jl")
end

if length(ARGS) == 3
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
    sim = ARGS[3]
else
    id = 1 
    n_rounds = 10
    sim = "unboostedsimulation"
end

## no boosting 

simulation = simulations[sim]

nhospitals = counthospitals(simulation)
ndates = countdates(simulation; dateid=:t)

@unpack newstaff, patients, staff = datamatrices(simulation, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(simulation)

stringency = coviddata.StringencyIndex_Average[1:ndates]
community = simulation.CommunityCases[1:ndates] ./ 56_000_000

# numbers vaccinated currently simulated 
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for d ∈ axes(vaccinated, 1), h ∈ axes(vaccinated, 2)
        if d == 300
            vaccinated[d, h] = 0.9
        end
    end
    vaccinated
end

model = fitmodel(
    newstaff, patients, staff, vaccinated, community, 
    vpd, psb, stringency, ndates, nhospitals
)

fitted_pt = pigeons( ;
    target=TuringLogPotential(model),
    n_rounds=0,
    n_chains=10,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(1000 + id),
    variational=GaussianReference(),
)

new_pt = fitted_pt

for i ∈ 1:2#n_rounds
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


include("loaddata.jl")  # includes using DrWatson and @quickactivate :ImmuneBoostingHealthcare

using Turing, Pigeons

if length(ARGS) == 2 
    id = parse(Int, ARGS[1])
    n_rounds = min(10, parse(Int, ARGS[2]))
else
    id = 1 
    n_rounds = 10
end

nhospitals = counthospitals(coviddata)
ndates = countdates(coviddata)

@unpack newstaff, patients, staff = datamatrices(coviddata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(coviddata)

stringency = coviddata.StringencyIndex_Average[1:ndates]
weeklycases = coviddata.weeklycases[1:ndates]

# numbers vaccinated currently simulated 
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for d ∈ axes(vaccinated, 1), h ∈ axes(vaccinated, 2)
        if d ∈ [ 300, 450, 650, 800 ]
            vaccinated[d, h] = 0.8
        end
    end
    vaccinated
end

model = fitmodel(
    newstaff, patients, staff, vaccinated, weeklycases, 
    vpd, psb, stringency, ndates, nhospitals
)

fitted_pt = pigeons( ;
    target=TuringLogPotential(model),
    n_rounds=0,
    n_chains=8,
    multithreaded=true,
    record=[ traces; record_default() ],
    seed=(n_rounds * 10 + id),
    variational=GaussianReference(),
)
fitted_chains = Chains(fitted_pt)

for i ∈ 1:n_rounds 
    fitted_pt = increment_n_rounds!(fitted_pt, 1)
    new_pt = pigeons(fitted_pt)
    new_chains = Chains(new_pt)
    fitted_pt = deepcopy(new_pt)
    output = Dict(
        "chain" => fitted_chains, 
        "pt" => fitted_pt, 
        "n_rounds" => i, 
        "n_chains" => 8,
    )
    safesave(datadir("sims", "fittedvalues_id_$(id)_round_$(i).jld2"), output)
end

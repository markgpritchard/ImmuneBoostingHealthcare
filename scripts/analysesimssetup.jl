
using CategoricalArrays, CSV, DataFrames, Dates, Pigeons, Random, Turing

if isfile(datadir("exp_pro", "finaldata.jld2"))
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
else 
    include("loaddata.jl")
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
end

if length(ARGS) == 3
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
    sim = ARGS[3]
elseif length(ARGS) == 4
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
    sim = ARGS[3]
    ω = parse(Float64, ARGS[4])
else
    id = 1 
    n_rounds = 2
    sim = "unboostedsimulation"
    ω = 0.01
end

if isfile(datadir("sims", "$sim.jld2"))
    simulations = load(datadir("sims", "$sim.jld2"))
else 
    include("generatesimulations.jl")
end

## no boosting 

simulation = simulations[sim]

nhospitals = counthospitals(simulation)
ndates = countdates(simulation; dateid=:t)

@unpack newstaff, patients, staff = datamatrices(simulation, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(simulation)

stringency = simulation.StringencyIndex_Average[1:ndates]
community = simulation.weeklycases[1:ndates] ./ 56_000_000

# numbers of vaccinated healthcare workers, assumed equal in all hospitals 
vaccinated = [ vaccinatestaff(t) for t ∈ 1:ndates ] 

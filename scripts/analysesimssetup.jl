
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

if isfile(datadir("sims", "unboostedsimulation.jld2"))
    unboostedsimulation = load(datadir("sims", "unboostedsimulation.jld2"))
else 
    include("generatesimulations.jl")
end

if isfile(datadir("sims", "midboostedsimulation.jld2"))
    midboostedsimulation = load(datadir("sims", "midboostedsimulation.jld2"))
else 
    include("generatesimulations.jl")
end

if isfile(datadir("sims", "boostedsimulation.jld2"))
    boostedsimulation = load(datadir("sims", "boostedsimulation.jld2"))
else 
    include("generatesimulations.jl")
end

## no boosting 

if sim == "unboostedsimulation"
    simulation = unboostedsimulation["unboostedsimulation"]
elseif sim == "boostedsimulation" 
    simulation = boostedsimulation["boostedsimulation"]
else 
    simulation = midboostedsimulation["midboostedsimulation"]
end

nhospitals = counthospitals(simulation)
ndates = countdates(simulation; dateid=:t)

@unpack newstaff, patients, staff = datamatrices(simulation, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(simulation)

stringency = simulation.StringencyIndex_Average[1:ndates]
community = simulation.weeklycases[1:ndates] ./ 56_000_000

for t ∈ 1:832 
    if t < 200 
        if community[t] > 0.0005 
            community[t] = community[t-1]
        end
    else
        if community[t] > 0.005 
            community[t] = community[t-1]
        end
    end
end


# numbers of vaccinated healthcare workers, assumed equal in all hospitals 
vaccinated = [ vaccinatestaff(t) for t ∈ 1:ndates ] 

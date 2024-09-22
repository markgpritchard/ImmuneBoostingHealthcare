
using CategoricalArrays, CSV, DataFrames, Dates, Pigeons, Random, Turing

if isfile(datadir("exp_pro", "finaldata.jld2"))
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
else 
    include("loaddata.jl")
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
end

if length(ARGS) == 2 
    id = parse(Int, ARGS[1])
    n_rounds = parse(Int, ARGS[2])
else
    id = 1 
    n_rounds = 4
end

nhospitals = counthospitals(finaldata)
ndates = countdates(finaldata)

@unpack newstaff, patients, staff = datamatrices(finaldata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(finaldata)

stringency = finaldata.StringencyIndex_Average[1:ndates]
community = finaldata.weeklycases[1:ndates] ./ 56_000_000
# clean anomalous community values 
for i ∈ [ 106:112; 684:690 ]
    community[i] = (community[i-7] + community[i+7]) / 2
end

# numbers of vaccinated healthcare workers, assumed equal in all hospitals 
vaccinated = [ vaccinatestaff(t) for t ∈ 1:ndates ] 

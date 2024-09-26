
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

using DataFrames, Dates

if length(ARGS) == 1
    option = parse(Int, ARGS[1])
else
    option = 3
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Outputs from simulations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# in case they are needed, omit various combinations of values 

if option == 1  # unboosted simulation 
    include("analysesimssetup.jl")
    @unpack unboostedsimulation = simulations
    data = unboostedsimulation
    filenamekey = "unboostedsimulation"
    shortfilenamekey = "unboosted"
    dateid = :t
elseif option == 2  # boosted simulation 
    include("analysesimssetup.jl")
    @unpack boostedsimulation = simulations
    data = boostedsimulation
    filenamekey = "boostedsimulation"
    shortfilenamekey = "boosted"
    dateid = :t
elseif option == 3  # data
    if isfile(datadir("exp_pro", "finaldata.jld2"))
        finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
    else 
        include("loaddata.jl")
        finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
    end
    
    include("analysedatasetup.jl")
    data = finaldata
    filenamekey = "coviddata"
    shortfilenamekey = "data"
    dateid = :Date
end

df = loadchainsdf("fittedvalues_$filenamekey")
outputs = processoutputs(data, finaldata, df, vaccinated; dateid)
safesave(datadir("sims", "$(shortfilenamekey)outputs.jld2"), outputs)

#=
for i ∈ [ 0, 1], j ∈ [ 0, 1 ], k ∈ [ 0, 1 ], l ∈ [ 0, 1 ]
    i == j == k == l && continue 
    includechains = Int[ ]
    if i == 1 push!(includechains, 1) end
    if j == 1 push!(includechains, 2) end
    if k == 1 push!(includechains, 3) end
    if l == 1 push!(includechains, 4) end
    outputs = processoutputs(
        data, 
        finaldata, 
        filter(:chain => x -> x ∈ includechains, df), 
        vaccinated; 
        dateid
    )
    safesave(
        datadir("sims", "$(shortfilenamekey)outputs_chains_$(includechains).jld2"), outputs
    )
end
=#
df_180 = loadchainsdf("fittedvalues_$(filenamekey)_omega_0.00556"; omega=0.00556)
output_180 = processoutputs(data, finaldata, df_180, vaccinated; dateid)
safesave(datadir("sims", "$(shortfilenamekey)outputs_omega180.jld2"), output_180)
#=
for i ∈ [ 0, 1], j ∈ [ 0, 1 ], k ∈ [ 0, 1 ], l ∈ [ 0, 1 ]
    i == j == k == l && continue 
    includechains = Int[ ]
    if i == 1 push!(includechains, 1) end
    if j == 1 push!(includechains, 2) end
    if k == 1 push!(includechains, 3) end
    if l == 1 push!(includechains, 4) end
    outputs = processoutputs(
        data, 
        finaldata, 
        filter(:chain => x ->  x ∈ includechains, df_180), 
        vaccinated; 
        dateid
    )
    safesave(
        datadir("sims", "$(shortfilenamekey)outputs_chains_omega180_$(includechains).jld2"), 
        outputs
    )
end
=#
df_100 = loadchainsdf("fittedvalues_$(filenamekey)_omega_0.01"; omega=0.01)
output_100 = processoutputs(data, finaldata, df_100, vaccinated; dateid)
safesave(datadir("sims", "$(shortfilenamekey)outputs_omega100.jld2"), output_100)
#=
for i ∈ [ 0, 1], j ∈ [ 0, 1 ], k ∈ [ 0, 1 ], l ∈ [ 0, 1 ]
    i == j == k == l && continue 
    includechains = Int[ ]
    if i == 1 push!(includechains, 1) end
    if j == 1 push!(includechains, 2) end
    if k == 1 push!(includechains, 3) end
    if l == 1 push!(includechains, 4) end
    outputs = processoutputs(
        data, 
        finaldata, 
        filter(:chain => x ->  x ∈ includechains, df_100), 
        vaccinated; 
        dateid=:t
    )
    safesave(
        datadir("sims", "$(shortfilenamekey)outputs_chains_omega100_$(includechains).jld2"), 
        outputs
    )
end
=#
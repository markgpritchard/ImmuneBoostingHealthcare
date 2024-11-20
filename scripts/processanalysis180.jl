
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include("analysesimssetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if isfile(datadir("sims", "alternativevaccinations.jld2"))
    alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))
else 
    include("processanalysisvaccinationdates.jl")
    alternativevaccinations = load(datadir("sims", "alternativevaccinations.jld2"))
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations without natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@info "preparing unboostedoutputs180"
unboostedoutputs180 = processoutputsdict(
    unboostedsimulation["unboostedsimulation"],
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556",
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.00556, #selectchains,
)
safesave(datadir("sims", "unboostedoutputs180.jld2"), unboostedoutputs180)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@info "preparing boostedoutputs180"
boostedoutputs180 = processoutputsdict(
    boostedsimulation["boostedsimulation"],
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.00556, #selectchains,
)
safesave(datadir("sims", "boostedoutputs180.jld2"), boostedoutputs180)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with "mid-level" natural immune boosting (ψ = 0.5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@info "preparing midboostedoutputs180"
midboostedoutputs180 = processoutputsdict(
    midboostedsimulation["midboostedsimulation"],
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.00556, #selectchains,
)
safesave(datadir("sims", "midboostedoutputs180.jld2"), midboostedoutputs180)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots from data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if isfile(datadir("exp_pro", "finaldata.jld2"))
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
else 
    include("loaddata.jl")
    finaldata = load(datadir("exp_pro", "finaldata.jld2"))["finaldata"]
end

include("analysedatasetup.jl")

@info "preparing dataoutputs180"
dataoutputs180 = processoutputsdict(
    newstaff,
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.00556, #selectchains,
)
safesave(datadir("sims", "dataoutputs180.jld2"), dataoutputs180)

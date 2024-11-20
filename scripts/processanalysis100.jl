
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

@info "preparing unboostedoutputs100"
unboostedoutputs100 = processoutputsdict(
    unboostedsimulation["unboostedsimulation"],
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.01, #selectchains,
)
safesave(datadir("sims", "unboostedoutputs100.jld2"), unboostedoutputs100)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@info "preparing boostedoutputs100"
boostedoutputs100 = processoutputsdict(
    boostedsimulation["boostedsimulation"],
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.01, #selectchains,
)
safesave(datadir("sims", "boostedoutputs100.jld2"), boostedoutputs100)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with "mid-level" natural immune boosting (ψ = 0.5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

@info "preparing midboostedoutputs100"
midboostedoutputs100 = processoutputsdict(
    midboostedsimulation["midboostedsimulation"],
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.01, #selectchains,
)
safesave(datadir("sims", "midboostedoutputs100.jld2"), midboostedoutputs100)


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

@info "preparing dataoutputs100"
dataoutputs100 = processoutputsdict(
    newstaff,
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    vaccinated, 
    alternativevaccinations;
    dates=470:831, jseries=1:nhospitals, omega=0.01, #selectchains,
)
safesave(datadir("sims", "dataoutputs100.jld2"), dataoutputs100)

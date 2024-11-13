
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include("analysesimssetup.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

alternativevaccinations = let 
    minus2months = [ vaccinatestaff(t; boostdateoffset=-62) for t ∈ 1:ndates ] 
    minus1month = [ vaccinatestaff(t; boostdateoffset=-31) for t ∈ 1:ndates ] 
    plus1month = [ vaccinatestaff(t; boostdateoffset=30) for t ∈ 1:ndates ] 
    plus2months = [ vaccinatestaff(t; boostdateoffset=61) for t ∈ 1:ndates ] 
    @strdict minus2months minus1month plus1month plus2months
end

safesave(datadir("sims", "alternativevaccinations.jld2"), alternativevaccinations)


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


using DrWatson

@quickactivate :ImmuneBoostingHealthcare

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

altvaccinated_minus2months = [ vaccinatestaff(t; boostdateoffset=-62) for t ∈ 1:ndates ] 
altvaccinated_minus1month = [ vaccinatestaff(t; boostdateoffset=-31) for t ∈ 1:ndates ] 
altvaccinated_plus1month = [ vaccinatestaff(t; boostdateoffset=30) for t ∈ 1:ndates ] 
altvaccinated_plus2months = [ vaccinatestaff(t; boostdateoffset=61) for t ∈ 1:ndates ] 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations without natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unboostedobserveddiagnosesafterjuly = [ 
    (
        d = filter(
            [ :StringCodes, :t ] => (x, y) -> x == c && y ∈ 470:831, 
            unboostedsimulation["unboostedsimulation"]
        );
        sum(d.StaffNewAbsences) ./ d.StaffTotal[1]
    )
    for c ∈ unique(unboostedsimulation["unboostedsimulation"].StringCodes)
]

unboosteddf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x ∈ [ 2, 4 ], unboosteddf_omega180)


## hospital-specific parameters, unboosted immunity lasts 180 days

unboostedoutputsperhospital_omega180 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.00556, #selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_omega180_forcepsi0 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.00556, #selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_omega180_diagnosesafterjuly = predicttotaldiagnoses(
    unboostedoutputsperhospital_omega180["chaindf"], 
    unboostedoutputsperhospital_omega180["patients"], 
    vaccinated, 
    unboostedoutputsperhospital_omega180["community"], 
    unboostedoutputsperhospital_omega180["vpd"], 
    unboostedoutputsperhospital_omega180["psb"], 
    unboostedoutputsperhospital_omega180["stringency"], 
    unboostedoutputsperhospital_omega180["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

unboostedoutputsperhospital_omega180_counterfactuals = producecounterfactualoutputsdict(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556",
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556,# selectchains=[ 2, 4 ],
)


## hospital-specific parameters, unboosted immunity lasts 100 days

unboosteddf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_unboostedsimulationperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x ∈ [ 3, 4 ], unboosteddf_omega100)

unboostedoutputsperhospital_omega100 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.01, #selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_omega100_forcepsi0 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.01, #selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_omega100_diagnosesafterjuly = predicttotaldiagnoses(
    unboostedoutputsperhospital_omega100["chaindf"], 
    unboostedoutputsperhospital_omega100["patients"], 
    vaccinated, 
    unboostedoutputsperhospital_omega100["community"], 
    unboostedoutputsperhospital_omega100["vpd"], 
    unboostedoutputsperhospital_omega100["psb"], 
    unboostedoutputsperhospital_omega100["stringency"], 
    unboostedoutputsperhospital_omega100["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

unboostedoutputsperhospital_omega100_counterfactuals = producecounterfactualoutputsdict(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 3, 4 ],
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

boostedobserveddiagnosesafterjuly = [ 
    (
        d = filter(
            [ :StringCodes, :t ] => (x, y) -> x == c && y ∈ 470:831, 
            boostedsimulation["boostedsimulation"]
        );
        sum(d.StaffNewAbsences) ./ d.StaffTotal[1]
    )
    for c ∈ unique(boostedsimulation["boostedsimulation"].StringCodes)
]

# hospital-specific parameters, unboosted immunity lasts 180 days

boosteddf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_boostedsimulationperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x ∈ [ 2, 3, 4 ], boosteddf_omega180)

boostedoutputsperhospital_omega180 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_omega180_forcepsi0 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_omega180_diagnosesafterjuly = predicttotaldiagnoses(
    boostedoutputsperhospital_omega180["chaindf"], 
    boostedoutputsperhospital_omega180["patients"], 
    vaccinated, 
    boostedoutputsperhospital_omega180["community"], 
    boostedoutputsperhospital_omega180["vpd"], 
    boostedoutputsperhospital_omega180["psb"], 
    boostedoutputsperhospital_omega180["stringency"], 
    boostedoutputsperhospital_omega180["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

boostedoutputsperhospital_omega180_counterfactuals = producecounterfactualoutputsdict(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)


## hospital-specific parameters, unboosted immunity lasts 100 days

boosteddf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_boostedsimulationperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x <= 2, boosteddf_omega100)

boostedoutputsperhospital_omega100 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.01, #selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_forcepsi0 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.01, #selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_diagnosesafterjuly = predicttotaldiagnoses(
    boostedoutputsperhospital_omega100["chaindf"], 
    boostedoutputsperhospital_omega100["patients"], 
    vaccinated, 
    boostedoutputsperhospital_omega100["community"], 
    boostedoutputsperhospital_omega100["vpd"], 
    boostedoutputsperhospital_omega100["psb"], 
    boostedoutputsperhospital_omega100["stringency"], 
    boostedoutputsperhospital_omega100["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

boostedoutputsperhospital_omega100_counterfactuals = producecounterfactualoutputsdict(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with "mid-level" natural immune boosting (ψ = 0.5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

midboostedobserveddiagnosesafterjuly = [ 
    (
        d = filter(
            [ :StringCodes, :t ] => (x, y) -> x == c && y ∈ 470:831, 
            midboostedsimulation["midboostedsimulation"]
        );
        sum(d.StaffNewAbsences) ./ d.StaffTotal[1]
    )
    for c ∈ unique(midboostedsimulation["midboostedsimulation"].StringCodes)
]

# hospital-specific parameters, unboosted immunity lasts 180 days

midboosteddf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x ∈ [ 2, 3, 4 ], midboosteddf_omega180)

midboostedoutputsperhospital_omega180 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

midboostedoutputsperhospital_omega180_forcepsi0 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

midboostedoutputsperhospital_omega180_diagnosesafterjuly = predicttotaldiagnoses(
    midboostedoutputsperhospital_omega180["chaindf"], 
    midboostedoutputsperhospital_omega180["patients"], 
    vaccinated, 
    midboostedoutputsperhospital_omega180["community"], 
    midboostedoutputsperhospital_omega180["vpd"], 
    midboostedoutputsperhospital_omega180["psb"], 
    midboostedoutputsperhospital_omega180["stringency"], 
    midboostedoutputsperhospital_omega180["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

midboostedoutputsperhospital_omega180_counterfactuals = producecounterfactualoutputsdict(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)


## hospital-specific parameters, unboosted immunity lasts 100 days

midboosteddf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_midboostedsimulationperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x <= 2, midboosteddf_omega100)

midboostedoutputsperhospital_omega100 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.01, #selectchains=[ 1, 2 ],
)

midboostedoutputsperhospital_omega100_forcepsi0 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.01, #selectchains=[ 1, 2 ],
)

midboostedoutputsperhospital_omega100_diagnosesafterjuly = predicttotaldiagnoses(
    midboostedoutputsperhospital_omega100["chaindf"], 
    midboostedoutputsperhospital_omega100["patients"], 
    vaccinated, 
    midboostedoutputsperhospital_omega100["community"], 
    midboostedoutputsperhospital_omega100["vpd"], 
    midboostedoutputsperhospital_omega100["psb"], 
    midboostedoutputsperhospital_omega100["stringency"], 
    midboostedoutputsperhospital_omega100["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

midboostedoutputsperhospital_omega100_counterfactuals = producecounterfactualoutputsdict(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)


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

dataobserveddiagnosesafterjuly = [ 
    sum(@view newstaff[470:831, i]) 
    for i ∈ axes(newstaff, 2) 
]

# hospital-specific parameters, unboosted immunity lasts 180 days

datadf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_coviddataperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x == 3, boosteddf_omega180)

dataoutputsperhospital_omega180 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.00556, #selectchains=[ 1, 2 ],
)

dataoutputsperhospital_omega180_forcepsi0 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.00556,# selectchains=3,
)

dataoutputsperhospital_omega180_diagnosesafterjuly = predicttotaldiagnoses(
    dataoutputsperhospital_omega180["chaindf"], 
    dataoutputsperhospital_omega180["patients"], 
    vaccinated, 
    dataoutputsperhospital_omega180["community"], 
    dataoutputsperhospital_omega180["vpd"], 
    dataoutputsperhospital_omega180["psb"], 
    dataoutputsperhospital_omega180["stringency"], 
    dataoutputsperhospital_omega180["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

dataoutputsperhospital_omega180_counterfactuals = producecounterfactualoutputsdict(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=3,
)


## hospital-specific parameters, unboosted immunity lasts 100 days

datadf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_coviddataperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x == 3, boosteddf_omega100)

dataoutputsperhospital_omega100 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.01, #selectchains=[ 1, 2 ],
)

dataoutputsperhospital_omega100_forcepsi0 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    vaccinated, 
    1:nhospitals; 
    forcepsi=0.0, dateid=:t, omega=0.01,# selectchains=3,
)

dataoutputsperhospital_omega100_diagnosesafterjuly = predicttotaldiagnoses(
    dataoutputsperhospital_omega100["chaindf"], 
    dataoutputsperhospital_omega100["patients"], 
    vaccinated, 
    dataoutputsperhospital_omega100["community"], 
    dataoutputsperhospital_omega100["vpd"], 
    dataoutputsperhospital_omega100["psb"], 
    dataoutputsperhospital_omega100["stringency"], 
    dataoutputsperhospital_omega100["ndates"], 
    1:nhospitals;
    daterange=470:831,
)

dataoutputsperhospital_omega100_counterfactuals = producecounterfactualoutputsdict(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    [ 
        altvaccinated_minus2months, 
        altvaccinated_minus1month, 
        altvaccinated_plus1month, 
        altvaccinated_plus2months 
    ], 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=3,
)

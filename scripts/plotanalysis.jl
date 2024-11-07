
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions

include("analysesimssetup.jl")



function subsetdiagnoses(dataset, vaccinated, jseries=1:dataset["nhospitals"]; daterange)
    return predicttotaldiagnoses(
        dataset["chaindf"], 
        dataset["patients"], 
        vaccinated, 
        dataset["community"], 
        dataset["vpd"], 
        dataset["psb"], 
        dataset["stringency"], 
        dataset["ndates"], 
        jseries;
        daterange,
    )
end

function estimateeffectofboosting(boosted, unboosted; nhospitals, CrI=[ 0.05, 0.95 ])
    raweffectofboosting = boosted .- unboosted
   
    medianeffectofboosting = [ quantile(raweffectofboosting[:, j], 0.5) for j ∈ 1:nhospitals ]
    lceffectofboosting = [ quantile(raweffectofboosting[:, j], CrI[1]) for j ∈ 1:nhospitals ]
    uceffectofboosting = [ quantile(raweffectofboosting[:, j], CrI[2]) for j ∈ 1:nhospitals ]
    
    return @ntuple medianeffectofboosting lceffectofboosting uceffectofboosting 
end








function plotcounterfactualvaccine!(
    ax::Axis, counterfactual::Vector{<:Real}, originaldata::Vector{<:Real}; 
    color=COLOURVECTOR[1], markersize=3,
)
    scatter!(
        ax, originaldata, 100 .* (counterfactual .- originaldata) ./ originaldata; 
        color, markersize,
    )
end

function plotcounterfactualvaccine!(
    ax::Axis, counterfactual::Dict{<:AbstractString, <:Any}, originaldata::NamedTuple; 
    kwargs...
)
    plotcounterfactualvaccine!(
        ax, counterfactual["mediantotaldiagnoses"], originaldata.mediantotaldiagnoses; 
        kwargs...
    )
end

function plotcounterfactualvaccine!(
    gl::GridLayout, cf_m2, cf_m1, cf_p1, cf_p2, originaldata; 
    kwargs...
)
    axs = [ Axis(gl[1, i]) for i ∈ 1:4 ]
    for (ax, cf) ∈ zip(axs, [ cf_m2, cf_m1, cf_p1, cf_p2 ])
        plotcounterfactualvaccine!(ax, cf, originaldata)
    end
    for (i, ax) ∈ enumerate(axs)
        hlines!(ax, 0; color=:black, linestyle=:dot)
        formataxis!(ax; hidey=(i != 1))
    end
    linkyaxes!(axs...)
    Label(gl[0, 1], "2 months earlier"; fontsize=11.84, tellwidth=false)
    Label(gl[0, 2], "1 month earlier"; fontsize=11.84, tellwidth=false)
    Label(gl[0, 3], "1 month later"; fontsize=11.84, tellwidth=false)
    Label(gl[0, 4], "2 months later"; fontsize=11.84, tellwidth=false)
    Label(
        gl[1, 0], "Proportional change in infections, %"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        gl[2, 1:4], "Mean infections per healthcare worker, after July 2021"; 
        fontsize=11.84, tellwidth=false
    )

    colgap!(gl, 1, 5)
    for r ∈ [ 1, 2 ] rowgap!(gl, r, 5) end
end


function calculatemodelleddiagnoses()


end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random sample of hospitals
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jseries = [
    82, 38, 23, 25, 74, 6, 9, 75, 57, 65, 7, 59, 
    81, 31, 3, 50, 10, 94, 15, 42, 19, 8, 47, 21, 27
]


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

unboosteddf_subset_omega180 = loadchainsperhospitaldf(
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556"; 
    jseries, omega=0.00556
)
filter!(:chain => x -> x ∈ [ 2, 4 ], unboosteddf_subset_omega180)

plotchains(
    unboosteddf_subset_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

# hospital-specific parameters, unboosted immunity lasts 180 days

unboostedoutputsperhospital_subset_omega180 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556", 
    vaccinated, 
    jseries; 
    dateid=:t, omega=0.00556, selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_subset_omega180_forcepsi0 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556", 
    vaccinated, 
    jseries; 
    forcepsi=0.0, dateid=:t, omega=0.00556, selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_subset_omega180_diagnosesafterjuly = predicttotaldiagnoses(
    unboostedoutputsperhospital_subset_omega180["chaindf"], 
    unboostedoutputsperhospital_subset_omega180["patients"], 
    vaccinated, 
    unboostedoutputsperhospital_subset_omega180["community"], 
    unboostedoutputsperhospital_subset_omega180["vpd"], 
    unboostedoutputsperhospital_subset_omega180["psb"], 
    unboostedoutputsperhospital_subset_omega180["stringency"], 
    unboostedoutputsperhospital_subset_omega180["ndates"], 
    jseries;
    daterange=470:831,
)

unboostedoutputsperhospital_subset_omega180_m2 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_minus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_subset_omega180_m1 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_minus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_subset_omega180_p1 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_plus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_subset_omega180_p2 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_plus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 4 ],
)


# hospital-specific parameters, unboosted immunity lasts 100 days

unboosteddf_subset_omega100 = loadchainsperhospitaldf(
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01"; 
    jseries, omega=0.01
)
filter!(:chain => x -> x ∈ [ 3, 4 ], unboosteddf_subset_omega100)

plotchains(
    unboosteddf_subset_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

unboostedoutputsperhospital_subset_omega100 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01", 
    vaccinated, 
    jseries; 
    dateid=:t, omega=0.01, selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_subset_omega100_forcepsi0 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01", 
    vaccinated, 
    jseries; 
    forcepsi=0.0, dateid=:t, omega=0.01, selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_subset_omega100_diagnosesafterjuly = predicttotaldiagnoses(
    unboostedoutputsperhospital_subset_omega100["chaindf"], 
    unboostedoutputsperhospital_subset_omega100["patients"], 
    vaccinated, 
    unboostedoutputsperhospital_subset_omega100["community"], 
    unboostedoutputsperhospital_subset_omega100["vpd"], 
    unboostedoutputsperhospital_subset_omega100["psb"], 
    unboostedoutputsperhospital_subset_omega100["stringency"], 
    unboostedoutputsperhospital_subset_omega100["ndates"], 
    jseries;
    daterange=470:831,
)

unboostedoutputsperhospital_subset_omega100_m2 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_minus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_subset_omega100_m1 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_minus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_subset_omega100_p1 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_plus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_subset_omega100_p2 = processoutputsperhospital(
    simulations["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_plus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 3, 4 ],
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# hospital-specific parameters, unboosted immunity lasts 180 days

boosteddf_subset_omega180 = loadchainsperhospitaldf(
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556"; 
    jseries, omega=0.00556
)
filter!(:chain => x -> x ∈ [ 2, 3, 4 ], boosteddf_subset_omega180)

plotchains(
    boosteddf_subset_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

boostedoutputsperhospital_subset_omega180 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556", 
    vaccinated, 
    jseries; 
    dateid=:t, omega=0.00556, selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_subset_omega180_forcepsi0 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556", 
    vaccinated, 
    jseries; 
    forcepsi=0.0, dateid=:t, omega=0.00556, selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_subset_omega180_diagnosesafterjuly = predicttotaldiagnoses(
    boostedoutputsperhospital_subset_omega180["chaindf"], 
    boostedoutputsperhospital_subset_omega180["patients"], 
    vaccinated, 
    boostedoutputsperhospital_subset_omega180["community"], 
    boostedoutputsperhospital_subset_omega180["vpd"], 
    boostedoutputsperhospital_subset_omega180["psb"], 
    boostedoutputsperhospital_subset_omega180["stringency"], 
    boostedoutputsperhospital_subset_omega180["ndates"], 
    jseries;
    daterange=470:831,
)

boostedoutputsperhospital_subset_omega180_m2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_minus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_subset_omega180_m1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_minus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_subset_omega180_p1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_plus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_subset_omega180_p2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.00556", 
    altvaccinated_plus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=[ 2, 3, 4 ],
)

fitfig_subset_omega180 = let 
    unboosted_estimatedeffect = estimateeffectofboosting(
        unboostedoutputsperhospital_subset_omega180["totaldiagnoses"], 
        unboostedoutputsperhospital_subset_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals=length(jseries)
    )
    boosted_estimatedeffect = estimateeffectofboosting(
        boostedoutputsperhospital_subset_omega180["totaldiagnoses"], 
        boostedoutputsperhospital_subset_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals=length(jseries)
    )
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:2 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:2 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ 
                unboostedoutputsperhospital_subset_omega180, 
                boostedoutputsperhospital_subset_omega180 
            ],
            [ unboosted_estimatedeffect, boosted_estimatedeffect ]
        )
            plotoutputs!(ax1, data)
            scatter!(
                ax2, 
                data["totalinfections"], 
                effect.medianeffectofboosting; 
                color=COLOURVECTOR[1], markersize=3
            )
            rangebars!(
                ax2, 
                data["totalinfections"], 
                effect.lceffectofboosting, 
                effect.uceffectofboosting; 
                color=( COLOURVECTOR[1], 0.1 ),
            )
        end

        linkaxes!(axs1...)
        linkaxes!(axs2...)

        fig    
    end
    fig
end

boostedtotalsperhospital_subset_omega180changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 550 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[2, 1])
        plotcounterfactualvaccine!(
            ga, 
            unboostedoutputsperhospital_subset_omega180_m2, 
            unboostedoutputsperhospital_subset_omega180_m1, 
            unboostedoutputsperhospital_subset_omega180_p1, 
            unboostedoutputsperhospital_subset_omega180_p2, 
            unboostedoutputsperhospital_subset_omega180_diagnosesafterjuly
        )
        plotcounterfactualvaccine!(
            gb, 
            boostedoutputsperhospital_subset_omega180_m2, 
            boostedoutputsperhospital_subset_omega180_m1, 
            boostedoutputsperhospital_subset_omega180_p1, 
            boostedoutputsperhospital_subset_omega180_p2, 
            boostedoutputsperhospital_subset_omega180_diagnosesafterjuly
        )
        fig    
    end
    fig
end


## hospital-specific parameters, unboosted immunity lasts 100 days

boosteddf_subset_omega100 = loadchainsperhospitaldf(
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
filter!(:chain => x -> x <= 2, boosteddf_subset_omega100)

plotchains(
    boosteddf_subset_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

boostedoutputsperhospital_subset_omega100 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01", 
    vaccinated, 
    jseries; 
    dateid=:t, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_subset_omega100_forcepsi0 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01", 
    vaccinated, 
    jseries; 
    forcepsi=0.0, dateid=:t, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_subset_omega100_diagnosesafterjuly = predicttotaldiagnoses(
    boostedoutputsperhospital_subset_omega100["chaindf"], 
    boostedoutputsperhospital_subset_omega100["patients"], 
    vaccinated, 
    boostedoutputsperhospital_subset_omega100["community"], 
    boostedoutputsperhospital_subset_omega100["vpd"], 
    boostedoutputsperhospital_subset_omega100["psb"], 
    boostedoutputsperhospital_subset_omega100["stringency"], 
    boostedoutputsperhospital_subset_omega100["ndates"], 
    jseries;
    daterange=470:831,
)

boostedoutputsperhospital_subset_omega100_m2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_minus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_subset_omega100_m1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_minus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_subset_omega100_p1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_plus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_subset_omega100_p2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_subset_omega_0.01", 
    altvaccinated_plus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)


fitfig_subset_omega100 = let 
    unboosted_estimatedeffect = estimateeffectofboosting(
        unboostedoutputsperhospital_subset_omega100["totaldiagnoses"], 
        unboostedoutputsperhospital_subset_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals=100
    )
    boosted_estimatedeffect = estimateeffectofboosting(
        boostedoutputsperhospital_subset_omega100["totaldiagnoses"], 
        boostedoutputsperhospital_subset_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals=100
    )
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:2 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:2 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ unboostedoutputsperhospital_subset_omega100, boostedoutputsperhospital_subset_omega100 ],
            [ unboosted_estimatedeffect, boosted_estimatedeffect ]
        )
            plotoutputs!(ax1, data)
            scatter!(
                ax2, 
                data["totalinfections"], 
                effect.medianeffectofboosting; 
                color=COLOURVECTOR[1], markersize=3
            )
            rangebars!(
                ax2, 
                data["totalinfections"], 
                effect.lceffectofboosting, 
                effect.uceffectofboosting; 
                color=( COLOURVECTOR[1], 0.1 ),
            )
        end

        linkaxes!(axs1...)
        linkaxes!(axs2...)

        fig    
    end
    fig
end

boostedtotalsperhospital_subset_omega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 550 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[2, 1])
        plotcounterfactualvaccine!(
            ga, 
            unboostedoutputsperhospital_subset_omega100_m2, 
            unboostedoutputsperhospital_subset_omega100_m1, 
            unboostedoutputsperhospital_subset_omega100_p1, 
            unboostedoutputsperhospital_subset_omega100_p2, 
            unboostedoutputsperhospital_subset_omega100_diagnosesafterjuly
        )
        plotcounterfactualvaccine!(
            gb, 
            boostedoutputsperhospital_subset_omega100_m2, 
            boostedoutputsperhospital_subset_omega100_m1, 
            boostedoutputsperhospital_subset_omega100_p1, 
            boostedoutputsperhospital_subset_omega100_p2, 
            boostedoutputsperhospital_subset_omega100_diagnosesafterjuly
        )
        fig    
    end
    fig
end


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

# hospital-specific parameters, unboosted immunity lasts 180 days

dataoutputsperhospital_omega180 = loadchainsperhospitaldf(
    "fittedvalues_coviddataperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x == 3, boosteddf_omega180)

plotchains(
    dataoutputsperhospital_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

boostedoutputsperhospital_omega180 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    jseries; 
    dateid=:t, omega=0.00556, selectchains=3,
)

boostedoutputsperhospital_omega180_forcepsi0 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    vaccinated, 
    jseries; 
    forcepsi=0.0, dateid=:t, omega=0.00556, selectchains=3,
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
    jseries;
    daterange=470:831,
)

boostedoutputsperhospital_omega180_m2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=3,
)

boostedoutputsperhospital_omega180_m1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=3,
)

boostedoutputsperhospital_omega180_p1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=3,
)

boostedoutputsperhospital_omega180_p2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.00556, selectchains=3,
)

fitfig_omega180 = let 
    unboosted_estimatedeffect = estimateeffectofboosting(
        unboostedoutputsperhospital_omega180["totaldiagnoses"], 
        unboostedoutputsperhospital_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals=100
    )
    boosted_estimatedeffect = estimateeffectofboosting(
        boostedoutputsperhospital_omega180["totaldiagnoses"], 
        boostedoutputsperhospital_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals=100
    )
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:2 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:2 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ unboostedoutputsperhospital_omega180, boostedoutputsperhospital_omega180 ],
            [ unboosted_estimatedeffect, boosted_estimatedeffect ]
        )
            plotoutputs!(ax1, data)
            scatter!(
                ax2, 
                data["totalinfections"], 
                effect.medianeffectofboosting; 
                color=COLOURVECTOR[1], markersize=3
            )
            rangebars!(
                ax2, 
                data["totalinfections"], 
                effect.lceffectofboosting, 
                effect.uceffectofboosting; 
                color=( COLOURVECTOR[1], 0.1 ),
            )
        end

        linkaxes!(axs1...)
        linkaxes!(axs2...)

        fig    
    end
    fig
end

boostedtotalsperhospital_omega180changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 550 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[2, 1])
        plotcounterfactualvaccine!(
            ga, 
            unboostedoutputsperhospital_omega180_m2, 
            unboostedoutputsperhospital_omega180_m1, 
            unboostedoutputsperhospital_omega180_p1, 
            unboostedoutputsperhospital_omega180_p2, 
            unboostedoutputsperhospital_omega180_diagnosesafterjuly
        )
        plotcounterfactualvaccine!(
            gb, 
            boostedoutputsperhospital_omega180_m2, 
            boostedoutputsperhospital_omega180_m1, 
            boostedoutputsperhospital_omega180_p1, 
            boostedoutputsperhospital_omega180_p2, 
            boostedoutputsperhospital_omega180_diagnosesafterjuly
        )
        fig    
    end
    fig
end


## hospital-specific parameters, unboosted immunity lasts 100 days

boosteddf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_boostedsimulationperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
filter!(:chain => x -> x <= 2, boosteddf_omega100)

plotchains(
    boosteddf_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

boostedoutputsperhospital_omega100 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    jseries; 
    dateid=:t, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_forcepsi0 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    vaccinated, 
    jseries; 
    forcepsi=0.0, dateid=:t, omega=0.01, selectchains=[ 1, 2 ],
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
    jseries;
    daterange=470:831,
)

boostedoutputsperhospital_omega100_m2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_m1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_p1 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus1month, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_p2 = processoutputsperhospital(
    simulations["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus2months, 
    jseries; 
    dateid=:t, daterange=470:831, omega=0.01, selectchains=[ 1, 2 ],
)


fitfig_omega100 = let 
    unboosted_estimatedeffect = estimateeffectofboosting(
        unboostedoutputsperhospital_omega100["totaldiagnoses"], 
        unboostedoutputsperhospital_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals=100
    )
    boosted_estimatedeffect = estimateeffectofboosting(
        boostedoutputsperhospital_omega100["totaldiagnoses"], 
        boostedoutputsperhospital_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals=100
    )
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:2 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:2 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ unboostedoutputsperhospital_omega100, boostedoutputsperhospital_omega100 ],
            [ unboosted_estimatedeffect, boosted_estimatedeffect ]
        )
            plotoutputs!(ax1, data)
            scatter!(
                ax2, 
                data["totalinfections"], 
                effect.medianeffectofboosting; 
                color=COLOURVECTOR[1], markersize=3
            )
            rangebars!(
                ax2, 
                data["totalinfections"], 
                effect.lceffectofboosting, 
                effect.uceffectofboosting; 
                color=( COLOURVECTOR[1], 0.1 ),
            )
        end

        linkaxes!(axs1...)
        linkaxes!(axs2...)

        fig    
    end
    fig
end

boostedtotalsperhospital_omega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 550 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[2, 1])
        plotcounterfactualvaccine!(
            ga, 
            unboostedoutputsperhospital_omega100_m2, 
            unboostedoutputsperhospital_omega100_m1, 
            unboostedoutputsperhospital_omega100_p1, 
            unboostedoutputsperhospital_omega100_p2, 
            unboostedoutputsperhospital_omega100_diagnosesafterjuly
        )
        plotcounterfactualvaccine!(
            gb, 
            boostedoutputsperhospital_omega100_m2, 
            boostedoutputsperhospital_omega100_m1, 
            boostedoutputsperhospital_omega100_p1, 
            boostedoutputsperhospital_omega100_p2, 
            boostedoutputsperhospital_omega100_diagnosesafterjuly
        )
        fig    
    end
    fig
end

# no hospital-specific parameters, unboosted immunity lasts 180 days
#=
dataoutputs_omega180 = let 
    datadf = loadchainsdf("fittedvalues_coviddata_omega_0.00556"; omega=0.00556)
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(finaldata, datadf, vaccinated)    
end

plotchains(dataoutputs_omega180["chaindf"])

dataoutputs_omega180fig = plotoutputs(dataoutputs_omega180)


# no hospital-specific parameters, unboosted immunity lasts 100 days

dataoutputs_omega100 = let 
    datadf = loadchainsdf("fittedvalues_coviddata_omega_0.01"; omega=0.01)
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputs(finaldata, datadf, vaccinated)    
end

plotchains(dataoutputs_omega100["chaindf"])

dataoutputs_omega100fig = plotoutputs(dataoutputs_omega100)

=#
# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 180 days

dataoutputsperhospital_omega180 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_omega_0.00556"; 
        jseries=1:nhospitals, omega=0.00556
    )
    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, vaccinated, 1:nhospitals)    
end

plotchains(dataoutputsperhospital_omega180["chaindf"]; size=( 400, 4800 ))
#=
dataoutputsperhospital_omega180_forcepsi0 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    for i ∈ axes(datadf, 1) 
        datadf.ψ[i] = 0.0 
    end

    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseries; )    
 
end

dataoutputsperhospital_omega180fig = let 
    raweffectofboosting = zeros(
        size(dataoutputsperhospital_omega180["totaldiagnoses"], 2), 
        length(jseries)
    )
    for i ∈ axes(raweffectofboosting, 1), j ∈ eachindex(jseries) 
        boosted = dataoutputsperhospital_omega180["totaldiagnoses"][j, i]
        unboosted = dataoutputsperhospital_omega180_forcepsi0["totaldiagnoses"][j, i]
        raweffectofboosting[i, j] = (boosted - unboosted) #/ boosted
    end

    medianeffectofboosting = [ 
        quantile(raweffectofboosting[:, j], 0.5) 
        for j ∈ eachindex(jseries) 
    ]
    lceffectofboosting = [ 
        quantile(raweffectofboosting[:, j], 0.05) 
        for j ∈ eachindex(jseries) 
    ]
    ucneffectofboosting = [ 
        quantile(raweffectofboosting[:, j], 0.95) 
        for j ∈ eachindex(jseries) 
    ]

    fig = Figure(; size=( 400, 500 ))
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[2, 1])
    
    plotoutputs!(ax1, dataoutputsperhospital_omega180)
    #plotoutputs!(ax2, dataoutputsperhospital_omega180_forcepsi0)

    scatter!(
        ax2, 
        dataoutputsperhospital_omega180["totalinfections"], 
        medianeffectofboosting; 
        color=:blue, markersize=3
    )
    rangebars!(
        ax2, 
        dataoutputsperhospital_omega180["totalinfections"], 
        lceffectofboosting, 
        ucneffectofboosting; 
        color=( :blue, 0.1 ),
    )

    fig
end


# modified vaccination schedule 
=#
#=
hospitaloutputs = plothospitaloutputs(dataoutputsperhospital_omega180; jseries, suppliedextra=false)

safesave(
    plotsdir("hospitaloutputs.svg"), hospitaloutputs
)
=#


#=
dataoutputsperhospital_omega180dr = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseries; daterange=470:831)    
end
=#



diagnosesafterjuly = predicttotaldiagnoses(
    dataoutputsperhospital_omega180["chaindf"], 
    dataoutputsperhospital_omega180["patients"], 
    vaccinated, 
    dataoutputsperhospital_omega180["community"], 
    dataoutputsperhospital_omega180["vpd"], 
    dataoutputsperhospital_omega180["psb"], 
    dataoutputsperhospital_omega180["stringency"], 
    dataoutputsperhospital_omega180["ndates"], 
    jseries;
    daterange=470:831,
)

dataoutputsperhospital_omega180_m2 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, altvaccinated_minus2months, jseries)    
end

dataoutputsperhospital_omega180_m1 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, altvaccinated_minus1month, jseries)    
end

dataoutputsperhospital_omega180_p1 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, altvaccinated_plus1month, jseries)    
end

dataoutputsperhospital_omega180_p2 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.00556"; 
        jseries, omega=0.00556
    )
    # remove chain that did not mix with others 
    filter!(:chain => x -> x in [ 2, 3 ], datadf)
    processoutputsperhospital(finaldata, datadf, altvaccinated_plus2months, jseries)    
end

changevaccinationdatefig = let 
    fig = Figure(; size=( 500, 300 ))
    #ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 1])
    #ax3 = Axis(fig[1, 2])    
    ax4 = Axis(fig[1, 2])
    #ax5 = Axis(fig[1, 3])
    ax6 = Axis(fig[1, 3])    
    #ax7 = Axis(fig[1, 4])
    ax8 = Axis(fig[1, 4])
    
    #plotoutputs!(ax1, dataoutputsperhospital_omega180_m2)
    scatter!(
        ax2, 
        diagnosesafterjuly.mediantotaldiagnoses, 
        (
            100 *
            (
                dataoutputsperhospital_omega180_m2["mediantotaldiagnoses"] .- 
                dataoutputsperhospital_omega180["mediantotaldiagnoses"]
            ) ./
            diagnosesafterjuly.mediantotaldiagnoses #diagnosesafterjuly.mediantotaldiagnoses
        ); 
        color=:blue, markersize=3
    )
    #plotoutputs!(ax3, dataoutputsperhospital_omega180_m1)
    scatter!(
        ax4, 
        diagnosesafterjuly.mediantotaldiagnoses, #dataoutputsperhospital_omega180["mediantotaldiagnoses"], 
        (
            100 *
            (
                dataoutputsperhospital_omega180_m1["mediantotaldiagnoses"] .- 
                dataoutputsperhospital_omega180["mediantotaldiagnoses"]
            ) ./
            diagnosesafterjuly.mediantotaldiagnoses #dataoutputsperhospital_omega180["mediantotaldiagnoses"]
        ); 
        color=:blue, markersize=3
    )
    #plotoutputs!(ax5, dataoutputsperhospital_omega180_p2)
    scatter!(
        ax6, 
        diagnosesafterjuly.mediantotaldiagnoses, #dataoutputsperhospital_omega180["mediantotaldiagnoses"], 
        (
            100 * 
            (
                dataoutputsperhospital_omega180_p1["mediantotaldiagnoses"] .- 
                dataoutputsperhospital_omega180["mediantotaldiagnoses"]
            ) ./
            diagnosesafterjuly.mediantotaldiagnoses #dataoutputsperhospital_omega180["mediantotaldiagnoses"]
        ); 
        color=:blue, markersize=3
    )
    #plotoutputs!(ax7, dataoutputsperhospital_omega180_p2)
    scatter!(
        ax8, 
        diagnosesafterjuly.mediantotaldiagnoses, #dataoutputsperhospital_omega180["mediantotaldiagnoses"], 
        (
            100 * 
            (
                dataoutputsperhospital_omega180_p2["mediantotaldiagnoses"] .- 
                dataoutputsperhospital_omega180["mediantotaldiagnoses"]
            ) ./
            dataoutputsperhospital_omega180["mediantotaldiagnoses"]
        ); 
        color=:blue, markersize=3
    )

    for ax ∈ [ ax2, ax4, ax6, ax8 ]
        hlines!(ax, 0; color=:black, linestyle=:dot)
    end

    linkyaxes!(ax2, ax4, ax6, ax8)
    #formataxis!(ax1; hidex=true, hidexticks=true, hidespines=( :b, :t, :r ) )
    formataxis!(ax2)

    #for ax ∈ [ ax3, ax5, ax7 ]
    #    formataxis!(
    #        ax; 
    #        hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
    #        hidespines=( :l, :b, :t, :r ) 
    #    )
    #end

    for ax ∈ [ ax4, ax6, ax8 ]
        formataxis!(
            ax; 
            hidey=true, hideyticks=true, 
            hidespines=( :l, :t, :r ) 
        )
    end

    Label(fig[0, 1], "2 months earlier"; fontsize=11.84, tellwidth=false)
    Label(fig[0, 2], "1 month earlier"; fontsize=11.84, tellwidth=false)
    Label(fig[0, 3], "1 month later"; fontsize=11.84, tellwidth=false)
    Label(fig[0, 4], "2 months later"; fontsize=11.84, tellwidth=false)
    Label(
        fig[1, 0], "Proportional change in infections, %"; 
        fontsize=11.84, rotation=π/2, tellheight=false
    )
    Label(
        fig[2, 1:4], "Mean infections per healthcare worker, after July 2021"; 
        fontsize=11.84, tellwidth=false
    )

    colgap!(fig.layout, 1, 5)
    for r ∈ [ 1, 2 ] rowgap!(fig.layout, r, 5) end

    fig 
end

safesave(
    plotsdir("changevaccinationdatefig.svg"), changevaccinationdatefig
)

# hospital-specific parameters for a subset of hospitals, unboosted immunity lasts 100 days
#=
dataoutputsperhospital_omega100 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.01"; 
        jseries, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseries; )    
end
#=
plotchains(dataoutputsperhospital_omega100["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospital_omega100fig = plotoutputs(dataoutputsperhospital_omega100)
=#
dataoutputsperhospital_omega100_forcepsi0 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_subset_omega_0.01"; 
        jseries, omega=0.01
    )
    for i ∈ axes(datadf, 1) 
        datadf.ψ[i] = 0.0 
    end

    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseries; )    
 
end

raweffectofboosting = zeros(10_000, length(jseries))
for i ∈ 1:10_000, j ∈ eachindex(jseries) 
    boosted = sample(dataoutputsperhospital_omega100["totaldiagnoses"][j, :])
    unboosted = sample(dataoutputsperhospital_omega100_forcepsi0["totaldiagnoses"][j, :])
    raweffectofboosting[i, j] = (boosted - unboosted) / boosted
end

medianeffectofboosting = [ quantile(raweffectofboosting[:, j], 0.5) for j ∈ eachindex(jseries) ]
lceffectofboosting = [ quantile(raweffectofboosting[:, j], 0.05) for j ∈ eachindex(jseries) ]
ucneffectofboosting = [ quantile(raweffectofboosting[:, j], 0.95) for j ∈ eachindex(jseries) ]

dataoutputsperhospital_omega100fig = plotoutputs(dataoutputsperhospital_omega100)

ax2 = Axis(dataoutputsperhospital_omega100fig[2, 1])
scatter!(
    ax2, 
    dataoutputsperhospital_omega100["totalinfections"], 
    medianeffectofboosting; 
    color=:blue, markersize=3
)
rangebars!(
    ax2, 
    dataoutputsperhospital_omega100["totalinfections"], 
    lceffectofboosting, 
    ucneffectofboosting; 
    color=( :blue, 0.1 ),
)

dataoutputsperhospital_omega100fig

safesave(
    plotsdir("dataoutputsperhospital_omega100fig.svg"), dataoutputsperhospital_omega100fig
)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 180 days

jseriesallhospitals = 1:nhospitals

dataoutputsperhospitalall_omega180 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_omega_0.00556"; 
        jseries=jseriesallhospitals, omega=0.00556
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseriesallhospitals)    
end
#=
plotchains(dataoutputsperhospitalall_omega180["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospitalall_omega180fig = plotoutputs(dataoutputsperhospitalall_omega180)
=#
dataoutputsperhospitalall_omega180_forcepsi0 = let 
    @unpack boostedsimulation = simulations
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_omega_0.00556"; 
        jseries=jseriesallhospitals, omega=0.00556
    )
    for i ∈ axes(datadf, 1) 
        datadf.ψ[i] = 0.0 
    end

    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseriesallhospitals)    
 
end

raweffectofboosting = zeros(10_000, nhospitals)
for i ∈ 1:10_000, j ∈ 1:nhospitals 
    boosted = sample(dataoutputsperhospitalall_omega180["totaldiagnoses"][j, :])
    unboosted = sample(dataoutputsperhospitalall_omega180_forcepsi0["totaldiagnoses"][j, :])
    raweffectofboosting[i, j] = boosted - unboosted 
end

medianeffectofboosting = [ quantile(raweffectofboosting[:, j], 0.5) for j ∈ 1:nhospitals ]
lceffectofboosting = [ quantile(raweffectofboosting[:, j], 0.05) for j ∈ 1:nhospitals ]
ucneffectofboosting = [ quantile(raweffectofboosting[:, j], 0.95) for j ∈ 1:nhospitals ]

dataoutputsperhospitalall_omega180fig = plotoutputs(dataoutputsperhospitalall_omega180)

ax2 = Axis(dataoutputsperhospitalall_omega180fig[2, 1])
scatter!(
    ax2, 
    dataoutputsperhospitalall_omega180["totalinfections"], 
    medianeffectofboosting; 
    color=:blue, markersize=3
)
rangebars!(
    ax2, 
    dataoutputsperhospitalall_omega180["totalinfections"], 
    lceffectofboosting, 
    ucneffectofboosting; 
    color=( :blue, 0.1 ),
)

dataoutputsperhospitalall_omega180fig

safesave(
    plotsdir("dataoutputsperhospitalall_omega180fig.svg"), 
    dataoutputsperhospitalall_omega180fig
)


# hospital-specific parameters for all hospitals, unboosted immunity lasts 100 days

dataoutputsperhospitalall_omega100 = let 
    datadf = loadchainsperhospitaldf(
        "fittedvalues_coviddataperhospital_omega_0.01"; 
        jseries=jseriesallhospitals, omega=0.01
    )
    # remove chain that did not mix with others 
    #filter!(:chain => x -> x in [ 2, 3 ], unboosteddf)
    processoutputsperhospital(finaldata, datadf, vaccinated, jseriesallhospitals)    
end

plotchains(dataoutputsperhospitalall_omega100["chaindf"]; size=( 400, 4800 ))

dataoutputsperhospitalall_omega100fig = plotoutputs(dataoutputsperhospitalall_omega100)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prior samples 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

priorsamples = sample(
    fitmodel(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    ), 
    Prior(), 
    1000
)

priorsamples_psi0 = sample(
    fitmodel(
        patients, staff, vaccinated, community, vpd, psb, stringency, ndates, nhospitals;
        psiprior=0
    ), 
    Prior(), 
    1000
)

priorsamplefig, priorsampleparameterfig = let 
    priorsamplesdf = DataFrame(priorsamples)
    priorsamplespredicteddiagnoses = predicttotaldiagnoses(
        priorsamplesdf, 
        patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
    priorsamplesbetas = calculatebetas(priorsamplesdf, vpd, psb, nhospitals)
    priorsampleslambdacs = calculatelambdacs(priorsamplesdf, stringency, community) 

    priorsamplesdf_psi0 = DataFrame(priorsamples_psi0)
    priorsamplespredicteddiagnoses_psi0 = predicttotaldiagnoses(
        priorsamplesdf_psi0, 0, 
        patients, vaccinated, community, vpd, psb, stringency, ndates, nhospitals
    )
    priorsamplesbetas_psi0 = calculatebetas(priorsamplesdf_psi0, vpd, psb, nhospitals)
    priorsampleslambdacs_psi0 = calculatelambdacs(priorsamplesdf_psi0, stringency, community) 

    priorsamplefig = Figure()
    axs = [ Axis(priorsamplefig[i, 1]) for i ∈ 1:3 ]
    scatter!(
        axs[1], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses.mediantotaldiagnoses; 
        color=:blue, markersize=3
    )
    rangebars!(
        axs[1], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses.lcitotaldiagnoses, 
        priorsamplespredicteddiagnoses.ucitotaldiagnoses; 
        color=( :blue, 0.1 ),
    )
    lines!(
        axs[1], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ]
    )

    scatter!(
        axs[2], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses_psi0.mediantotaldiagnoses; 
        color=:blue, markersize=3
    )
    rangebars!(
        axs[2], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        priorsamplespredicteddiagnoses_psi0.lcitotaldiagnoses, 
        priorsamplespredicteddiagnoses_psi0.ucitotaldiagnoses; 
        color=( :blue, 0.1 ),
    )
    lines!(
        axs[2], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ], 
        [ extrema([ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10)... ]
    )

    scatter!(
        axs[3], 
        [ sum(@view staff[:, i]) for i ∈ axes(staff, 2) ] ./ 10, 
        (
            priorsamplespredicteddiagnoses.mediantotaldiagnoses .- 
            priorsamplespredicteddiagnoses_psi0.mediantotaldiagnoses
        ) ./ priorsamplespredicteddiagnoses.mediantotaldiagnoses; 
        color=:blue, markersize=3,
    )

    ( priorsamplefig , 1 )
end

## Unboosted sims 

#unboostedsimulation = simulations["unboostedsimulation"]

@unpack unboostedsimulation = simulations

nhospitals = counthospitals(unboostedsimulation)
ndates = countdates(unboostedsimulation; dateid=:t)

# simulated numbers vaccinated  
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for t ∈ axes(vaccinated, 1), j ∈ axes(vaccinated, 2)
        vaccinated[300, j] = 0.8
    end
    vaccinated
end

@unpack newstaff, patients, staff = datamatrices(unboostedsimulation, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(unboostedsimulation)

stringency = finaldata.StringencyIndex_Average[1:ndates]
community = unboostedsimulation.CommunityCases[1:ndates] ./ 56_000_000

unboosteddf = loadchainsdf("fittedvalues_unboostedsimulation")
plotchains(unboosteddf)

@unpack βp, βh, βc = calculatebetas(unboosteddf, vpd, psb, stringency)

totalinfections = [ sum(@view newstaff[:, i]) for i ∈ 1:nhospitals ]

predictedinfections = summarizepredictedinfections(
    unboosteddf, patients, staff, community, vaccinated, vpd, psb, stringency
)

predictedinfectionswithoutboost = summarizepredictedinfections(
    unboosteddf, patients, staff, community, vaccinated, vpd, psb, stringency; ψ=0
)

medianbetas = let
    @unpack βp, βh, βc = calculatebetas(unboosteddf, vpd, psb, stringency)
    βclcris = zeros(ndates)
    βcmedians = zeros(ndates)
    βcucris = zeros(ndates)
    for i ∈ 1:ndates 
        vals = quantile([ b[i] for b ∈ βc ], [ 0.05, 0.5, 0.95 ])
        βclcris[i] = vals[1]
        βcmedians[i] = vals[2]
        βcucris[i] = vals[3]
    end
    βhlcris = zeros(nhospitals)
    βhmedians = zeros(nhospitals)
    βhucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βh ], [ 0.05, 0.5, 0.95 ])
        βhlcris[i] = vals[1]
        βhmedians[i] = vals[2]
        βhucris[i] = vals[3]
    end
    βplcris = zeros(nhospitals)
    βpmedians = zeros(nhospitals)
    βpucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βp ], [ 0.05, 0.5, 0.95 ])
        βplcris[i] = vals[1]
        βpmedians[i] = vals[2]
        βpucris[i] = vals[3]
    end
    @ntuple βclcris βcmedians βcucris βhlcris βhmedians βhucris βplcris βpmedians βpucris
end

fig = Figure()
axs1 = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
axs2 = [ Axis(fig[i, 2]) for i ∈ 1:3 ]
scatter!(axs1[1], totalinfections, predictedinfections.medians; color=:blue, markersize=3)
rangebars!(
    axs1[1], totalinfections, predictedinfections.lcris, predictedinfections.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[1], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[2], totalinfections, predictedinfectionswithoutboost.medians; 
    color=:blue, markersize=3,
)
rangebars!(
    axs1[2], 
    totalinfections, 
    predictedinfectionswithoutboost.lcris,
    predictedinfectionswithoutboost.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[2], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[3], 
    totalinfections, 
    (predictedinfections.medians .- predictedinfectionswithoutboost.medians) ./ 
        predictedinfections.medians; 
    color=:blue, markersize=3,
)
#scatter!(ax3, totalinfections, predictedinfections.means .- predictedinfectionswithoutboost.means; color=:blue, markersize=3)

scatter!(axs2[2], [ b[1] for b ∈ allbetas ], medianbetas[:βhmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[2], 
    [ b[1] for b ∈ allbetas ], 
    medianbetas[:βhlcris],
    medianbetas[:βhucris]; 
    color=( :blue, 0.1 ),
)
lines!(
    axs2[2], 
    [ extrema([ b[1] for b ∈ allbetas ])... ], 
    [ extrema([ b[1] for b ∈ allbetas ])... ])


scatter!(
    axs2[3], [ b[2] for b ∈ allbetas ], medianbetas[:βpmedians]; 
    color=:blue, markersize=3
)
rangebars!(
    axs2[3], 
    [ b[2] for b ∈ allbetas ], 
    medianbetas[:βplcris],
    medianbetas[:βpucris]; 
    color=( :blue, 0.1 ),
)
lines!(
    axs2[3], 
    [ extrema([ b[2] for b ∈ allbetas ])... ], 
    [ extrema([ b[2] for b ∈ allbetas ])... ]
)


linkxaxes!(axs1...)
linkyaxes!(axs1[1], axs1[2])

fig



## Boosted sims 

boostedsimulation = simulations["boostedsimulation"]

## Covid data 


nhospitals = counthospitals(finaldata)
ndates = countdates(finaldata)



# numbers vaccinated currently simulated 
vaccinated = let
    vaccinated = zeros(ndates, nhospitals)
    for t ∈ axes(vaccinated, 1), j ∈ axes(vaccinated, 2)
        if t ∈ [ 300, 450, 650, 800 ]
            vaccinated[t, j] = 0.8
        end
    end
    vaccinated
end

@unpack newstaff, patients, staff = datamatrices(finaldata, ndates, nhospitals)
@unpack vpd, psb = hospitalconditionmatrices(finaldata)

stringency = finaldata.StringencyIndex_Average[1:ndates]
community = finaldata.weeklycases[1:ndates] ./ 56_000_000

coviddf = loadchainsdf("fittedvalues_coviddata")
plotchains(coviddf)

@unpack βp, βh, βc = calculatebetas(coviddf, vpd, psb, stringency)

totalinfections = [ sum(@view newstaff[:, i]) for i ∈ 1:nhospitals ]

predictedinfections = summarizepredictedinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency
)

predictedinfectionswithoutboost = summarizepredictedinfections(
    coviddf, patients, staff, community, vaccinated, vpd, psb, stringency; 
    ψ=0,
)

medianbetas = let
    @unpack βp, βh, βc = calculatebetas(coviddf, vpd, psb, stringency)
    βclcris = zeros(ndates)
    βcmedians = zeros(ndates)
    βcucris = zeros(ndates)
    for i ∈ 1:ndates 
        vals = quantile([ b[i] for b ∈ βc ], [ 0.05, 0.5, 0.95 ])
        βclcris[i] = vals[1]
        βcmedians[i] = vals[2]
        βcucris[i] = vals[3]
    end
    βhlcris = zeros(nhospitals)
    βhmedians = zeros(nhospitals)
    βhucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βh ], [ 0.05, 0.5, 0.95 ])
        βhlcris[i] = vals[1]
        βhmedians[i] = vals[2]
        βhucris[i] = vals[3]
    end
    βplcris = zeros(nhospitals)
    βpmedians = zeros(nhospitals)
    βpucris = zeros(nhospitals)
    for i ∈ 1:nhospitals 
        vals = quantile([ b[i] for b ∈ βp ], [ 0.05, 0.5, 0.95 ])
        βplcris[i] = vals[1]
        βpmedians[i] = vals[2]
        βpucris[i] = vals[3]
    end
    @ntuple βclcris βcmedians βcucris βhlcris βhmedians βhucris βplcris βpmedians βpucris
end

fig = Figure()
axs1 = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
axs2 = [ Axis(fig[i, 2]) for i ∈ 1:3 ]
scatter!(axs1[1], totalinfections, predictedinfections.medians; color=:blue, markersize=3)
rangebars!(
    axs1[1], totalinfections, predictedinfections.lcris, predictedinfections.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[1], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[2], totalinfections, predictedinfectionswithoutboost.medians; 
    color=:blue, markersize=3,
)
rangebars!(
    axs1[2], 
    totalinfections, 
    predictedinfectionswithoutboost.lcris,
    predictedinfectionswithoutboost.ucris; 
    color=( :blue, 0.1 ),
)
lines!(axs1[2], [ extrema(totalinfections)... ], [ extrema(totalinfections)... ])

scatter!(
    axs1[3], 
    totalinfections, 
    (predictedinfections.medians .- predictedinfectionswithoutboost.medians) ./ 
        predictedinfections.medians; 
    color=:blue, markersize=3,
)
#scatter!(ax3, totalinfections, predictedinfections.means .- predictedinfectionswithoutboost.means; color=:blue, markersize=3)

scatter!(axs2[1], finaldata.Date[1:832], medianbetas[:βcmedians]; color=:blue, markersize=3)

#scatter!(axs2[2], ordinalrank(medianbetas[:βhmedians]), medianbetas[:βhmedians]; color=:blue, markersize=3)
scatter!(axs2[2], totalinfections, medianbetas[:βhmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[2], 
    totalinfections,
    #ordinalrank(medianbetas[:βhmedians]),
    medianbetas[:βhlcris],
    medianbetas[:βhucris]; 
    color=( :blue, 0.1 ),
)

#scatter!(axs2[3], ordinalrank(medianbetas[:βpmedians]), medianbetas[:βpmedians]; color=:blue, markersize=3)
scatter!(axs2[3], totalinfections, medianbetas[:βpmedians]; color=:blue, markersize=3)
rangebars!(
    axs2[3], 
    totalinfections,
    #ordinalrank(medianbetas[:βpmedians]),
    medianbetas[:βplcris],
    medianbetas[:βpucris]; 
    color=( :blue, 0.1 ),
)


linkxaxes!(axs1...)
linkyaxes!(axs1[1], axs1[2])

fig



##############################################################################################

# hospital RCD #25

using StaticArrays

idnumber = 80
id = unique(finaldata.StringCodes)[idnumber]

h20data = filter(:StringCodes => x -> x == id, finaldata)

@unpack βp, βh, βc = calculatebetas(
    coviddf, 
    h20data.VolumePerBed[1], 
    h20data.ProportionSingleBeds[1], 
    h20data.StringencyIndex_Average    
)

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, h20data.CovidAbsences; color=:black, markersize=3)
for i ∈ axes(coviddf, 1)
    predictedinfections = zeros(ndates) 
    immunevector = SizedVector{3}(zeros(Float64, 3)) 
    for t ∈ 2:ndates 
        foi = βc[i][t-1] * community[t-1] + βp[i][1] * h20data.PatientsProportion[t-1] + βh[i][1] * predictedinfections[t-1]
        v = vaccinated[(t - 1), idnumber]
        immune10 = predictedinfections[t-1] + v * (1 - sum(immunevector) - predictedinfections[t-1])
        # probability of boosting from natural immune boosting plus vaccination
        pb = (1 - exp(-coviddf.ψ[i] * foi)) * (1 - v) + v  
        for x ∈ 1:2
            immune10 += pb * immunevector[x]
            immunevector[x] += -(pb + 2 * coviddf.ω[i]) * immunevector[x] + 
                2 * coviddf.ω[i] * (1 - pb) * immunevector[x+1]
        end
        immunevector[2] += -(2 * coviddf.ω[i] * (1 - pb)) * 
            immunevector[2] + 
            immune10
        predictedinfections[t] = (1 - sum(immunevector)) * (1 - exp(-foi))
    end
    lines!(ax, predictedinfections; color=( :blue, 0.1 ))
end



fig




####

fig, ax = lines([ 1 - exp(-x * 0) for x ∈ 0:0.01:2 ])
lines!(ax, [ (1 - (1 - exp(-x)))^0 for x ∈ 0:0.01:2 ])


fig

=#

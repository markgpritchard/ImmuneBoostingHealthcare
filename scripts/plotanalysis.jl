
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
   
    medianeffectofboosting = [ quantile(raweffectofboosting[j, :], 0.5) for j ∈ 1:nhospitals ]
    lceffectofboosting = [ quantile(raweffectofboosting[j, :], CrI[1]) for j ∈ 1:nhospitals ]
    uceffectofboosting = [ quantile(raweffectofboosting[j, :], CrI[2]) for j ∈ 1:nhospitals ]
    
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
# Vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

altvaccinated_minus2months = [ vaccinatestaff(t; boostdateoffset=-62) for t ∈ 1:ndates ] 
altvaccinated_minus1month = [ vaccinatestaff(t; boostdateoffset=-31) for t ∈ 1:ndates ] 
altvaccinated_plus1month = [ vaccinatestaff(t; boostdateoffset=30) for t ∈ 1:ndates ] 
altvaccinated_plus2months = [ vaccinatestaff(t; boostdateoffset=61) for t ∈ 1:ndates ] 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations without natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unboosteddf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x ∈ [ 2, 4 ], unboosteddf_omega180)

plotchains(
    unboosteddf_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

# hospital-specific parameters, unboosted immunity lasts 180 days

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

unboostedoutputsperhospital_omega180_m2 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556,# selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_omega180_m1 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556,# selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_omega180_p1 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 4 ],
)

unboostedoutputsperhospital_omega180_p2 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 4 ],
)


# hospital-specific parameters, unboosted immunity lasts 100 days

unboosteddf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_unboostedsimulationperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x ∈ [ 3, 4 ], unboosteddf_omega100)

plotchains(
    unboosteddf_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

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

unboostedoutputsperhospital_omega100_m2 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_omega100_m1 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_omega100_p1 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 3, 4 ],
)

unboostedoutputsperhospital_omega100_p2 = processoutputsperhospital(
    unboostedsimulation["unboostedsimulation"], 
    finaldata, 
    "fittedvalues_unboostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 3, 4 ],
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with natural immune boosting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# hospital-specific parameters, unboosted immunity lasts 180 days

boosteddf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_boostedsimulationperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x ∈ [ 2, 3, 4 ], boosteddf_omega180)

plotchains(
    boosteddf_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

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

boostedoutputsperhospital_omega180_m2 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_omega180_m1 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_omega180_p1 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

boostedoutputsperhospital_omega180_p2 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

fitfig_omega180 = let 
    unboosted_estimatedeffect = estimateeffectofboosting(
        unboostedoutputsperhospital_omega180["totaldiagnoses"], 
        unboostedoutputsperhospital_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals
    )
    boosted_estimatedeffect = estimateeffectofboosting(
        boostedoutputsperhospital_omega180["totaldiagnoses"], 
        boostedoutputsperhospital_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals
    )
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:2 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:2 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ 
                unboostedoutputsperhospital_omega180, 
                boostedoutputsperhospital_omega180 
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
#filter!(:chain => x -> x <= 2, boosteddf_omega100)

plotchains(
    boosteddf_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

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

boostedoutputsperhospital_omega100_m2 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_m1 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01,# selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_p1 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)

boostedoutputsperhospital_omega100_p2 = processoutputsperhospital(
    boostedsimulation["boostedsimulation"], 
    finaldata, 
    "fittedvalues_boostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)


fitfig_omega100 = let 
    unboosted_estimatedeffect = estimateeffectofboosting(
        unboostedoutputsperhospital_omega100["totaldiagnoses"], 
        unboostedoutputsperhospital_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals
    )
    boosted_estimatedeffect = estimateeffectofboosting(
        boostedoutputsperhospital_omega100["totaldiagnoses"], 
        boostedoutputsperhospital_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulations with "mid'level" natural immune boosting (ψ = 0.5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# hospital-specific parameters, unboosted immunity lasts 180 days

midboosteddf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x ∈ [ 2, 3, 4 ], midboosteddf_omega180)

plotchains(
    midboosteddf_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

midboostedoutputsperhospital_omega180 = processoutputsperhospital(
    midboostedsimulation["boostedsimulation"], 
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

midboostedoutputsperhospital_omega180_m2 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

midboostedoutputsperhospital_omega180_m1 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

midboostedoutputsperhospital_omega180_p1 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)

midboostedoutputsperhospital_omega180_p2 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.00556", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=[ 2, 3, 4 ],
)


## hospital-specific parameters, unboosted immunity lasts 100 days

midboosteddf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_midboostedsimulationperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x <= 2, midboosteddf_omega100)

plotchains(
    midboosteddf_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

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

midboostedoutputsperhospital_omega100_m2 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)

midboostedoutputsperhospital_omega100_m1 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01,# selectchains=[ 1, 2 ],
)

midboostedoutputsperhospital_omega100_p1 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=[ 1, 2 ],
)

midboostedoutputsperhospital_omega100_p2 = processoutputsperhospital(
    midboostedsimulation["midboostedsimulation"], 
    finaldata, 
    "fittedvalues_midboostedsimulationperhospital_omega_0.01", 
    altvaccinated_plus2months, 
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

# hospital-specific parameters, unboosted immunity lasts 180 days

datadf_omega180 = loadchainsperhospitaldf(
    "fittedvalues_coviddataperhospital_omega_0.00556"; 
    jseries=1:nhospitals, omega=0.00556
)
#filter!(:chain => x -> x == 3, boosteddf_omega180)

plotchains(
    datadf_omega180; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

dataoutputsperhospital_omega180 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    vaccinated, 
    1:nhospitals; 
    dateid=:t, omega=0.01, #selectchains=[ 1, 2 ],
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

dataoutputsperhospital_omega180_m2 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=3,
)

dataoutputsperhospital_omega180_m1 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=3,
)

dataoutputsperhospital_omega180_p1 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556,# selectchains=3,
)

dataoutputsperhospital_omega180_p2 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.00556", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.00556, #selectchains=3,
)

fitfig_dataomega180 = let 
    data_estimatedeffect = estimateeffectofboosting(
        dataoutputsperhospital_omega180["totaldiagnoses"], 
        dataoutputsperhospital_omega180_forcepsi0["totaldiagnoses"]; 
        nhospitals
    )

    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:1 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:1 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ dataoutputsperhospital_omega180 ],
            [ data_estimatedeffect ]
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
        #gb = GridLayout(fig[2, 1])
        plotcounterfactualvaccine!(
            ga, 
            dataoutputsperhospital_omega180_m2, 
            dataoutputsperhospital_omega180_m1, 
            dataoutputsperhospital_omega180_p1, 
            dataoutputsperhospital_omega180_p2, 
            dataoutputsperhospital_omega180_diagnosesafterjuly
        )

        fig    
    end
    fig
end


## hospital-specific parameters, unboosted immunity lasts 100 days

datadf_omega100 = loadchainsperhospitaldf(
    "fittedvalues_coviddataperhospital_omega_0.01"; 
    jseries=1:nhospitals, omega=0.01
)
#filter!(:chain => x -> x == 3, boosteddf_omega100)

plotchains(
    datadf_omega100; 
    columns=[ 
        "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
        "sigma2", "hsigma2", "psigma2", "log_density" 
    ]
)

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

dataoutputsperhospital_omega100_m2 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    altvaccinated_minus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=3,
)

dataoutputsperhospital_omega100_m1 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    altvaccinated_minus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=3,
)

dataoutputsperhospital_omega100_p1 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    altvaccinated_plus1month, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01,# selectchains=3,
)

dataoutputsperhospital_omega100_p2 = processoutputsperhospital(
    finaldata, 
    "fittedvalues_coviddataperhospital_omega_0.01", 
    altvaccinated_plus2months, 
    1:nhospitals; 
    dateid=:t, daterange=470:831, omega=0.01, #selectchains=3,
)

fitfig_dataomega100 = let 
    data_estimatedeffect = estimateeffectofboosting(
        dataoutputsperhospital_omega100["totaldiagnoses"], 
        dataoutputsperhospital_omega100_forcepsi0["totaldiagnoses"]; 
        nhospitals
    )

    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 300 ))
        axs1 = [ Axis(fig[1, i]) for i ∈ 1:1 ]
        axs2 = [ Axis(fig[2, i]) for i ∈ 1:1 ]

        for (ax1, ax2, data, effect) ∈ zip(
            axs1, 
            axs2,
            [ dataoutputsperhospital_omega100 ],
            [ data_estimatedeffect ]
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
        #gb = GridLayout(fig[2, 1])
        plotcounterfactualvaccine!(
            ga, 
            dataoutputsperhospital_omega100_m2, 
            dataoutputsperhospital_omega100_m1, 
            dataoutputsperhospital_omega100_p1, 
            dataoutputsperhospital_omega100_p2, 
            dataoutputsperhospital_omega100_diagnosesafterjuly
        )

        fig    
    end
    fig
end

omega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 500, 750 ))
        ga = GridLayout(fig[1, 1])
        gb = GridLayout(fig[2, 1])
        gc = GridLayout(fig[3, 1])
        gd = GridLayout(fig[4, 1])
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
            midboostedoutputsperhospital_omega100_m2, 
            midboostedoutputsperhospital_omega100_m1, 
            midboostedoutputsperhospital_omega100_p1, 
            midboostedoutputsperhospital_omega100_p2, 
            midboostedoutputsperhospital_omega100_diagnosesafterjuly
        )
        plotcounterfactualvaccine!(
            gc, 
            boostedoutputsperhospital_omega100_m2, 
            boostedoutputsperhospital_omega100_m1, 
            boostedoutputsperhospital_omega100_p1, 
            boostedoutputsperhospital_omega100_p2, 
            boostedoutputsperhospital_omega100_diagnosesafterjuly
        )
        plotcounterfactualvaccine!(
            gd, 
            dataoutputsperhospital_omega100_m2, 
            dataoutputsperhospital_omega100_m1, 
            dataoutputsperhospital_omega100_p1, 
            dataoutputsperhospital_omega100_p2, 
            dataoutputsperhospital_omega100_diagnosesafterjuly
        )
        fig    
    end
    fig
end

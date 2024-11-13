
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions
import ImmuneBoostingHealthcare: Automatic, automatic

#include("processanalysis.jl")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chains
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const COLSFORCHAINPLOTS = [ 
    "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
    "sigma2", "hsigma2", "psigma2", "log_density" 
]

## Model without boosting

plotchains(unboosteddf_omega180; columns=COLSFORCHAINPLOTS)
plotchains(unboosteddf_omega100; columns=COLSFORCHAINPLOTS)

## Model with mid-level boosting (ψ = 0.5)

plotchains(midboosteddf_omega180; columns=COLSFORCHAINPLOTS)
plotchains(midboosteddf_omega100; columns=COLSFORCHAINPLOTS)

## Model with double-strength boosting (ψ = 2)

plotchains(boosteddf_omega180; columns=COLSFORCHAINPLOTS)
plotchains(boosteddf_omega100; columns=COLSFORCHAINPLOTS)

## Applied to covid-19 data

plotchains(datadf_omega180; columns=COLSFORCHAINPLOTS)
plotchains(datadf_omega100; columns=COLSFORCHAINPLOTS)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Changes in numbers of cases with changes in vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unboostedoutputs180 = load(datadir("sims", "unboostedoutputs180.jld2"))

omega180changevaccinationdatefig = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 587, 500 ))
    plotcaseswithchangedvaccinationdates!(
        fig, 
        470:831,  # dates
        [unboostedoutputs180, midboostedoutputs180, boostedoutputs180, dataoutputs180 ]
    )
    fig
end
safesave(plotsdir("omega180changevaccinationdatefig.pdf"), omega180changevaccinationdatefig)

omega100changevaccinationdatefig = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 587, 500 ))
    plotcaseswithchangedvaccinationdates!(
        fig, 
        470:831,  # dates
        [unboostedoutputs100, midboostedoutputs100, boostedoutputs100, dataoutputs100 ]
    )
    fig
end
safesave(plotsdir("omega100changevaccinationdatefig.pdf"), omega100changevaccinationdatefig)


omega180changevaccinationdatefig = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 587, 500 ))
    plotcaseswithchangedvaccinationdates!(
        fig, 
        470:831,  # dates
        [  # observedcasesvector
            unboostedobserveddiagnosesafterjuly,
            midboostedobserveddiagnosesafterjuly,
            boostedobserveddiagnosesafterjuly,
            dataobserveddiagnosesafterjuly
        ],
        [  # modelledcasesvector
            unboostedoutputsperhospital_omega180_diagnosesafterjuly.predicteddiagnoses,
            midboostedoutputsperhospital_omega180_diagnosesafterjuly.predicteddiagnoses,
            boostedoutputsperhospital_omega180_diagnosesafterjuly.predicteddiagnoses,
            dataoutputsperhospital_omega180_diagnosesafterjuly.predicteddiagnoses
        ],
        [  # counterfactualsvector
            unboostedoutputsperhospital_omega180_counterfactuals,
            midboostedoutputsperhospital_omega180_counterfactuals,
            boostedoutputsperhospital_omega180_counterfactuals,
            dataoutputsperhospital_omega180_counterfactuals 
        ]
    )
    fig
end

safesave(plotsdir("omega180changevaccinationdatefig.pdf"), omega180changevaccinationdatefig)



plotcumulativecounterfactualvaccine!(
    gl::GridLayout, 
    dates, 
    modelledoriginal::Array{<:Number, 3}, 
    cf_m2::Array{<:Number, 3}, 
    cf_m1::Array{<:Number, 3}, 
    cf_p1::Array{<:Number, 3}, 
    cf_p2::Array{<:Number, 3};
    observeddiagnoses=automatic,
    xticks=Makie.automatic,
    vline=nothing,
    vspan=nothing,
)

function plotcounterfactualvaccine!(
    ax::Axis, 
    counterfactual::Vector{<:Real}, 
    modelledoriginal::Vector{<:Real}, 
    originaldata_x::Vector{<:Real}; 
    color=COLOURVECTOR[1], markersize=3,
)


unboostedomega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 587, 800 ))
        ga = GridLayout(fig[1, 1])
        plotcumulativecounterfactualvaccine!(
            ga, 
            470:831, 
            unboostedoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[470:831, highinds, :] , 
            unboostedoutputsperhospital_omega100_m2["predictdiagnoses"][470:831, highinds, :], 
            unboostedoutputsperhospital_omega100_m1["predictdiagnoses"][470:831, highinds, :], 
            unboostedoutputsperhospital_omega100_p1["predictdiagnoses"][470:831, highinds, :], 
            unboostedoutputsperhospital_omega100_p2["predictdiagnoses"][470:831, highinds, :];
            observeddiagnoses=unboostedobserveddiagnosesafterjuly[highinds],
            xticks=(
                [ 469, 561, 653, 743, 834 ],
                [ "July", "Oct.", "Jan.", "April", "July" ]
            ),
            vline=531,
            vspan = [ (531 + x):(621 + x) for x ∈ [ -62, -31, 30, 61 ] ],
        )
        fig    
    end
    fig
end



fig



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
            unboostedoutputsperhospital_omega100_diagnosesafterjuly,
            unboostedobserveddiagnosesafterjuly 
        )
        plotcounterfactualvaccine!(
            gb, 
            midboostedoutputsperhospital_omega100_m2, 
            midboostedoutputsperhospital_omega100_m1, 
            midboostedoutputsperhospital_omega100_p1, 
            midboostedoutputsperhospital_omega100_p2, 
            midboostedoutputsperhospital_omega100_diagnosesafterjuly,
            midboostedobserveddiagnosesafterjuly 
        )
        plotcounterfactualvaccine!(
            gc, 
            boostedoutputsperhospital_omega100_m2, 
            boostedoutputsperhospital_omega100_m1, 
            boostedoutputsperhospital_omega100_p1, 
            boostedoutputsperhospital_omega100_p2, 
            boostedoutputsperhospital_omega100_diagnosesafterjuly,
            boostedobserveddiagnosesafterjuly 
        )
        plotcounterfactualvaccine!(
            gd, 
            dataoutputsperhospital_omega100_m2, 
            dataoutputsperhospital_omega100_m1, 
            dataoutputsperhospital_omega100_p1, 
            dataoutputsperhospital_omega100_p2, 
            dataoutputsperhospital_omega100_diagnosesafterjuly,
            dataobserveddiagnosesafterjuly 
        )
        fig    
    end
    fig
end

# identify hospital with greatest number of cases 

fig = Figure(; size=( 800, 1600))
axs1 = [ Axis(fig[i, 1]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs2 = [ Axis(fig[i, 2]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs3 = [ Axis(fig[i, 3]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs4 = [ Axis(fig[i, 4]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs5 = [ Axis(fig[i, 5]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs6 = [ Axis(fig[i, 6]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
for (i, ax) ∈ enumerate(axs1)
    lines!(ax, 470:831, staff[470:831, i])
    formataxis!(ax; hidex=(i!=23))
end
for (i, ax) ∈ enumerate(axs2)
    med = [ quantile(dataoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, med)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
for (i, ax) ∈ enumerate(axs3)
    med = [ quantile(dataoutputsperhospital_omega100_m2["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, med)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
for (i, ax) ∈ enumerate(axs4)
    med = [ quantile(dataoutputsperhospital_omega100_m1["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, med)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
for (i, ax) ∈ enumerate(axs5)
    med = [ quantile(dataoutputsperhospital_omega100_p1["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, med)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
for (i, ax) ∈ enumerate(axs6)
    med = [ quantile(dataoutputsperhospital_omega100_p2["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, med)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
linkaxes!(axs1..., axs2..., axs3..., axs4..., axs5..., axs6...)

fig



fig = Figure(; size=( 800, 1600))
axs1 = [ Axis(fig[i, 1]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs2 = [ Axis(fig[i, 2]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs3 = [ Axis(fig[i, 3]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
axs4 = [ Axis(fig[i, 4]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
#axs5 = [ Axis(fig[i, 5]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
#axs6 = [ Axis(fig[i, 6]; xticks=[ 531, 592, 653, 712 ]) for i ∈ 1:23 ]
#for (i, ax) ∈ enumerate(axs1)
#    lines!(ax, 470:831, staff[470:831, i])
#    formataxis!(ax; hidex=(i!=23))
#end
firstmed = [ quantile(dataoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[j, i, :], 0.5) for i ∈ 1:23, j ∈ 470:831 ]

#for (i, ax) ∈ enumerate(axs2)
#    lines!(ax, 470:831, firstmed)
#    formataxis!(ax; hidex=(i!=23), hidey=true)
#end
for (i, ax) ∈ enumerate(axs1)
    med = [ quantile(dataoutputsperhospital_omega100_m2["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, cumsum(med .- firstmed[i, :]) ./ 10)
    hlines!(ax, 0; color=:black, linestyle=:dot)
    formataxis!(ax; hidex=(i!=23))
end
for (i, ax) ∈ enumerate(axs2)
    med = [ quantile(dataoutputsperhospital_omega100_m1["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, cumsum(med .- firstmed[i, :]) ./ 10)
    hlines!(ax, 0; color=:black, linestyle=:dot)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
for (i, ax) ∈ enumerate(axs3)
    med = [ quantile(dataoutputsperhospital_omega100_p1["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, cumsum(med .- firstmed[i, :]) ./ 10)
    hlines!(ax, 0; color=:black, linestyle=:dot)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
for (i, ax) ∈ enumerate(axs4)
    med = [ quantile(dataoutputsperhospital_omega100_p2["predictdiagnoses"][j, i, :], 0.5) for j ∈ 470:831 ]
    lines!(ax, 470:831, cumsum(med .- firstmed[i, :]) ./ 10)
    hlines!(ax, 0; color=:black, linestyle=:dot)
    formataxis!(ax; hidex=(i!=23), hidey=true)
end
linkaxes!(axs1..., axs2..., axs3..., axs4...)

fig




highinds = findall(x -> x > median(unboostedobserveddiagnosesafterjuly), unboostedobserveddiagnosesafterjuly)

unboostedomega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 587, 800 ))
        ga = GridLayout(fig[1, 1])
        plotcumulativecounterfactualvaccine!(
            ga, 
            470:831, 
            unboostedoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[470:831, highinds, :] , 
            unboostedoutputsperhospital_omega100_m2["predictdiagnoses"][470:831, highinds, :], 
            unboostedoutputsperhospital_omega100_m1["predictdiagnoses"][470:831, highinds, :], 
            unboostedoutputsperhospital_omega100_p1["predictdiagnoses"][470:831, highinds, :], 
            unboostedoutputsperhospital_omega100_p2["predictdiagnoses"][470:831, highinds, :];
            observeddiagnoses=unboostedobserveddiagnosesafterjuly[highinds],
            xticks=(
                [ 469, 561, 653, 743, 834 ],
                [ "July", "Oct.", "Jan.", "April", "July" ]
            ),
            vline=531,
            vspan = [ (531 + x):(621 + x) for x ∈ [ -62, -31, 30, 61 ] ],
        )
        fig    
    end
    fig
end

midboostedomega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 800, 1600 ))
        ga = GridLayout(fig[1, 1])
        plotcumulativecounterfactualvaccine!(
            ga, 
            470:831, 
            midboostedoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[470:831, :, :] , 
            midboostedoutputsperhospital_omega100_m2["predictdiagnoses"][470:831, :, :], 
            midboostedoutputsperhospital_omega100_m1["predictdiagnoses"][470:831, :, :], 
            midboostedoutputsperhospital_omega100_p1["predictdiagnoses"][470:831, :, :], 
            midboostedoutputsperhospital_omega100_p2["predictdiagnoses"][470:831, :, :];
            observeddiagnoses=midboostedobserveddiagnosesafterjuly,
            xticks=[ 531, 592, 653, 712 ],
        )
        fig    
    end
    fig
end

boostedomega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 800, 1600 ))
        ga = GridLayout(fig[1, 1])
        plotcumulativecounterfactualvaccine!(
            ga, 
            470:831, 
            boostedoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[470:831, :, :] , 
            boostedoutputsperhospital_omega100_m2["predictdiagnoses"][470:831, :, :], 
            boostedoutputsperhospital_omega100_m1["predictdiagnoses"][470:831, :, :], 
            boostedoutputsperhospital_omega100_p1["predictdiagnoses"][470:831, :, :], 
            boostedoutputsperhospital_omega100_p2["predictdiagnoses"][470:831, :, :];
            observeddiagnoses=boostedobserveddiagnosesafterjuly,
            xticks=[ 531, 592, 653, 712 ],
        )
        fig    
    end
    fig
end

dataomega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 800, 1600 ))
        ga = GridLayout(fig[1, 1])
        plotcumulativecounterfactualvaccine!(
            ga, 
            470:831, 
            dataoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[470:831, :, :] , 
            dataoutputsperhospital_omega100_m2["predictdiagnoses"][470:831, :, :], 
            dataoutputsperhospital_omega100_m1["predictdiagnoses"][470:831, :, :], 
            dataoutputsperhospital_omega100_p1["predictdiagnoses"][470:831, :, :], 
            dataoutputsperhospital_omega100_p2["predictdiagnoses"][470:831, :, :];
            observeddiagnoses=dataobserveddiagnosesafterjuly,
            xticks=[ 531, 592, 653, 712 ],
        )
        fig    
    end
    fig
end


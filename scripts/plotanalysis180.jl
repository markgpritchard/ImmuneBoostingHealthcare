
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions
import ImmuneBoostingHealthcare: Automatic, automatic


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chains
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unboostedoutputs180 = load(datadir("sims", "unboostedoutputs180.jld2"))
midboostedoutputs180 = load(datadir("sims", "midboostedoutputs180.jld2"))
boostedoutputs180 = load(datadir("sims", "boostedoutputs180.jld2"))
dataoutputs180 = load(datadir("sims", "dataoutputs180.jld2"))

const COLSFORCHAINPLOTS = [ 
    "α1", "α2", "α3", "α4", "α5", "α6", "α7", "α8", "ω", "θ", "ψ", 
    "sigma2", "hsigma2", "psigma2", "log_density" 
]
CHAINPLOTYLABELS = [
    L"$\alpha_1$", 
    L"$\alpha_2$", 
    L"$\alpha_3$", 
    L"$\alpha_4$", 
    L"$\alpha_5$", 
    L"$\alpha_6$", 
    L"$\alpha_7$", 
    L"$\alpha_8$", 
    L"$\omega$",
    L"$\theta$",
    L"$\psi$",
    L"$\sigma^2$",
    L"$h_{\sigma^2}$",
    L"$p_{\sigma^2}$",
    "Log density"
]

## Model without boosting

plotchains(unboostedoutputs180["df"]; columns=COLSFORCHAINPLOTS)

## Model with mid-level boosting (ψ = 0.5)

plotchains(midboostedoutputs180["df"]; columns=COLSFORCHAINPLOTS)

## Model with double-strength boosting (ψ = 2)

plotchains(boostedoutputs180["df"]; columns=COLSFORCHAINPLOTS)

## Applied to covid-19 data

datachainfig180 = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 700 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    plotchains!(
        ga, dataoutputs180["df"]; 
        columns=COLSFORCHAINPLOTS[[ 1:6; [ 15 ] ]], 
        ylabels=CHAINPLOTYLABELS[[ 1:6; [ 15 ] ]], 
        yticks=WilkinsonTicks(3)
    )
    plotchains!(
        gb, dataoutputs180["df"]; 
        columns=COLSFORCHAINPLOTS[7:14], 
        ylabels=CHAINPLOTYLABELS[7:14], 
        yticks=WilkinsonTicks(3)
    )

    fig 
end
safesave(plotsdir("datachainfig180.pdf"), datachainfig180)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Changes in numbers of cases with changes in vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

omega180changevaccinationdatefig = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 587, 500 ))
    plotcaseswithchangedvaccinationdates!(
        fig, 
        470:831,  # dates
        [ unboostedoutputs180, midboostedoutputs180, boostedoutputs180, dataoutputs180 ]
    )
    fig
end
safesave(plotsdir("omega180changevaccinationdatefig.pdf"), omega180changevaccinationdatefig)

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

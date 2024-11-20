
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions
import ImmuneBoostingHealthcare: Automatic, automatic


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chains
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unboostedoutputs100 = load(datadir("sims", "unboostedoutputs100.jld2"))
midboostedoutputs100 = load(datadir("sims", "midboostedoutputs100.jld2"))
boostedoutputs100 = load(datadir("sims", "boostedoutputs100.jld2"))
dataoutputs100 = load(datadir("sims", "dataoutputs100.jld2"))

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

plotchains(unboostedoutputs100["df"]; columns=COLSFORCHAINPLOTS)

## Model with mid-level boosting (ψ = 0.5)

plotchains(midboostedoutputs100["df"]; columns=COLSFORCHAINPLOTS)

## Model with double-strength boosting (ψ = 2)

plotchains(boostedoutputs100["df"]; columns=COLSFORCHAINPLOTS)

## Applied to covid-19 data

datachainfig100 = with_theme(theme_latexfonts()) do 
    fig = Figure(; size=( 500, 700 ))
    ga = GridLayout(fig[1, 1])
    gb = GridLayout(fig[1, 2])

    plotchains!(
        ga, dataoutputs100["df"]; 
        columns=COLSFORCHAINPLOTS[[ 1:6; [ 15 ] ]], 
        ylabels=CHAINPLOTYLABELS[[ 1:6; [ 15 ] ]], 
        yticks=WilkinsonTicks(3)
    )
    plotchains!(
        gb, dataoutputs100["df"]; 
        columns=COLSFORCHAINPLOTS[7:14], 
        ylabels=CHAINPLOTYLABELS[7:14], 
        yticks=WilkinsonTicks(3)
    )

    fig 
end
safesave(plotsdir("datachainfig100.pdf"), datachainfig100)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Changes in numbers of cases with changes in vaccination times
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

omega100changevaccinationdatefig = with_theme(theme_latexfonts()) do
    fig = Figure(; size=( 587, 500 ))
    plotcaseswithchangedvaccinationdates!(
        fig, 
        470:831,  # dates
        [ unboostedoutputs100, midboostedoutputs100, boostedoutputs100, dataoutputs100 ]
    )
    fig
end
safesave(plotsdir("omega100changevaccinationdatefig.pdf"), omega100changevaccinationdatefig)

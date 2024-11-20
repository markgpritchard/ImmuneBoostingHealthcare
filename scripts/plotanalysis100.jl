
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions
import ImmuneBoostingHealthcare: Automatic, automatic

include("analysedatasetup.jl")

using PlotFormatting


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
# Cases per hospital
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


function modhcwseiirrr(  # version where vaccination is given as a parameter
    u, p::HCWSEIIRRRp, t::Integer, λc::Number, patients::Number, vaccinated::Number
) 
    S, E, I, I′1, I′2, I′3, I′4, I′5, I′6, I′7, I′8, I′9, I′10, R1, R2, R3, = u 

    λ = 1 - exp(-(λc + p.βp * patients + p.βh * I))
    λψ = vaccinated + (1 - exp(-p.ψ * (λc + p.βp * patients + p.βh * I))) * (1 - vaccinated)

    new_u = [
        S * (1 - λ - vaccinated * (1 - λ)) + R3 * 3 * p.ω * (1 - λψ),  # S 
        E * (1 - p.η) + λ * S,  # E 
        I * (1 - p.θ - p.γ * (1 - p.θ)) + p.η * E,  # I 
        p.θ * I,  # I′1
        I′1,  # I′2
        I′2,  # I′3
        I′3,  # I′4
        I′4,  # I′5
        I′5,  # I′6
        I′6,  # I′7
        I′7,  # I′8
        I′8,  # I′9
        I′9,  # I′10
        # R1:
        R1 * (1 - 3 * p.ω * (1 - λψ)) +  # previous R1 that has not waned
            S * (vaccinated * (1 - λ)) +   # vaccinated from S 
            p.γ * (1 - p.θ) * I +  # recovered 
            I′10 +  # ended isolation 
            λψ * (R2 + R3),  # boosted by exposure or vaccination from R2 and R3 
        R2 * (1 - 3 * p.ω * (1 - λψ) - λψ) + R1 * 3 * p.ω * (1 - λψ),  # R3
        R3 * (1 - 3 * p.ω * (1 - λψ) - λψ) + R2 * 3 * p.ω * (1 - λψ),  # R3
        λ,
        λψ,
    ]
    return new_u
end

totalstaff = [  
    (
        d = filter(:StringCodes => x -> x == h, finaldata);
        d.StaffTotal[1]
    )
    for h ∈ unique(finaldata.StringCodes)
]

lambda = zeros(832, 23, 16_384)
forceofboosting = zeros(832, 23, 16_384)
susceptible = zeros(832, 23, 16_384)
predictedinfections = zeros(832, 23, 16_384)

for k ∈ axes(lambda, 3)
    λc = [ 
        max(zero(Float64), df.α7[k] + df.α8[k] * (100 - s)) 
        for s ∈ stringency 
    ] .* community
    for j ∈ axes(lambda, 2)
        u = zeros(18)
        u[1] = staff[j]
        p = HCWSEIIRRRp(
            getproperty(df, "betahs[$j]")[k],  # βh 
            getproperty(df, "betaps[$j]")[k],  # βp 
            0.5,  # η 
            0.2,  # γ 
            df.ψ[k],  
            df.ω[k],  
            df.θ[k]
        )
        for t ∈ axes(lambda, 1)
            u = modhcwseiirrr(u, p, t, λc[t], patients[t, j], vaccinated[t])
            lambda[t, j, k] = u[17]
            forceofboosting[t, j, k] = u[18]
            susceptible[t, j, k] = u[1]
            predictedinfections[t, j, k] = sum(@view u[4:13])
        end
    end
end





mediandiagnoses = zeros(832, 23)
lcidiagnoses = zeros(832, 23) 
ucidiagnoses = zeros(832, 23)

bandcolour = ( COLOURVECTOR[1], 0.33)

for t ∈ axes(mediandiagnoses, 1), j ∈ axes(mediandiagnoses, 2)
    l, m, u = quantile(dataoutputs100["modelledoutput"]["predictdiagnoses"][t, j, :], [ 0.05, 0.5, 0.95 ])
    lcidiagnoses[t, j] = l
    mediandiagnoses[t, j] = m
    ucidiagnoses[t, j] = u
end

fig = Figure(; size=( 500, 700 ))

axs = [ Axis(fig[i, j]) for i ∈ 1:6, j ∈ 1:4 ]

for (i, hosp) ∈ enumerate(unique(finaldata.StringCodes))
    d = filter(:StringCodes => x -> x == hosp, finaldata)
    scatter!(axs[i], d.t, d.StaffProportion; color=:black, markersize=2)
    lines!(axs[i], d.t, mediandiagnoses[:, i]; color=COLOURVECTOR[1])
    band!(axs[i], d.t, lcidiagnoses[:, i], ucidiagnoses[:, i]; color=bandcolour)
end

linkaxes!(axs...)

fig

lambda = zeros(832, 23, 16_384)
medianlambda = zeros(832, 23)
lcilambda = zeros(832, 23) 
ucilambda = zeros(832, 23)

df = dataoutputs100["df"]
for k ∈ axes(lambda, 3)
    λc = [ 
        max(zero(Float64), df.α7[k] + df.α8[k] * (100 - s)) 
        for s ∈ stringency 
    ] .* community
    for j ∈ axes(medianlambda, 2)

    end
end 

for t ∈ axes(medianlambda, 1), j ∈ axes(medianlambda, 2)

    #end

end

λc = [ 
            max(zero(Float64), df.α7[i] + df.α8[i] * (100 - s)) 
            for s ∈ stringency 
        ] .* community

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

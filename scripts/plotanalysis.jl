
using DrWatson

@quickactivate :ImmuneBoostingHealthcare

include(srcdir("PlottingFunctions.jl"))
using CairoMakie, CategoricalArrays, CSV, DataFrames, Dates, Pigeons, StatsBase
using .PlottingFunctions
import ImmuneBoostingHealthcare: Automatic, automatic

include("analysesimssetup.jl")



function subsetdiagnoses(dataset, vaccinated, jseries=1:dataset["nhospitals"]; daterange)
    return predicttotaldiagnoses(
        dataset["chaindf"], 
        dataset["patients"], dataobserveddiagnosesafterjuly 
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
    ax::Axis, 
    counterfactual::Vector{<:Real}, 
    modelledoriginal::Vector{<:Real}, 
    originaldata_x::Vector{<:Real}; 
    color=COLOURVECTOR[1], markersize=3,
)
    scatter!(
        ax, originaldata_x, 100 .* (counterfactual .- modelledoriginal) ./ modelledoriginal; 
        color, markersize,
    )
end

function plotcounterfactualvaccine!(
    ax::Axis, 
    counterfactual::Dict{<:AbstractString, <:Any}, 
    modelledoriginal, 
    originaldata_x; 
    kwargs...
)
    plotcounterfactualvaccine!(
        ax, counterfactual["mediantotaldiagnoses"], modelledoriginal, originaldata_x; 
        kwargs...
    )
end

function plotcounterfactualvaccine!(
    ax::Axis, 
    counterfactual::Vector{<:Real}, 
    modelledoriginal::NamedTuple, 
    originaldata_x; 
    kwargs...
)
    plotcounterfactualvaccine!(
        ax, counterfactual, modelledoriginal.mediantotaldiagnoses, originaldata_x; 
        kwargs...
    )
end

function plotcounterfactualvaccine!(ax::Axis, counterfactual, originaldata; kwargs...)
    plotcounterfactualvaccine!(ax, counterfactual, originaldata, originaldata; kwargs...)
end

function plotcounterfactualvaccine!(
    gl::GridLayout, cf_m2, cf_m1, cf_p1, cf_p2, modelledoriginal, originaldata_x; 
    kwargs...
)
    axs = [ Axis(gl[1, i]) for i ∈ 1:4 ]
    for (ax, cf) ∈ zip(axs, [ cf_m2, cf_m1, cf_p1, cf_p2 ])
        plotcounterfactualvaccine!(ax, cf, modelledoriginal, originaldata_x)
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

function plotcounterfactualvaccine!(
    gl::GridLayout, cf_m2, cf_m1, cf_p1, cf_p2, modelledoriginal; 
    kwargs...
)
    plotcounterfactualvaccine!(
        gl, cf_m2, cf_m1, cf_p1, cf_p2, modelledoriginal, modelledoriginal; 
        kwargs...
    )
end

function _cumulativecounterfactualvaccine(
    counterfactual::Matrix{<:Real}, originaldata::Matrix{<:Real}
)
    cumulativedifference = zeros(size(counterfactual))
    for j ∈ axes(counterfactual, 2)
        cumulativedifference[:, j] = cumsum(counterfactual[:, j] .- originaldata[:, j])
    end
    return cumulativedifference 
end

function _quantilecumulativecounterfactualvaccine(counterfactual, originaldata; kwargs...)
    cumulativedifference = _cumulativecounterfactualvaccine(counterfactual, originaldata)
    return _quantilecumulativecounterfactualvaccine(cumulativedifference; kwargs...)
end

function _quantilecumulativecounterfactualvaccine(
    cumulativedifference; 
    quantiles=[ 0.05, 0.5, 0.95 ]
)
    return [ 
        quantile(cumulativedifference[i, :], quantiles) 
        for i ∈ axes(cumulativedifference, 1) 
    ]
end

function plotcumulativecounterfactualvaccine!(
    ax, dates, counterfactual::Matrix{<:Real}, originaldata::Matrix{<:Real};
    quantiles=[ 0.05, 0.5, 0.95 ], kwargs...
)
    quantilevalues = _quantilecumulativecounterfactualvaccine(
        counterfactual, originaldata; 
        quantiles
    )
    plotcumulativecounterfactualvaccine!(ax, dates, quantilevalues; kwargs...)
end

function plotcumulativecounterfactualvaccine!(
    ax, dates, quantiles::AbstractVector{<:AbstractVector{<:Real}};
    color=COLOURVECTOR[1]
)
    lines!(ax, dates, [ quantiles[i][2] for i ∈ eachindex(quantiles) ]; color)
    band!(
        ax, 
        dates, 
        [ quantiles[i][1] for i ∈ eachindex(quantiles) ], 
        [ quantiles[i][3] for i ∈ eachindex(quantiles) ];
        color=( color, 0.25),
    )
end

function plotcumulativecounterfactualvaccine!(
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
    # rank hospitals so those with most cases are plotted first
    rankindices = _rankindicesoftotaldiagnoses(modelledoriginal, observeddiagnoses)
    nhospitals = size(modelledoriginal, 2)

    axs1 = [ Axis(gl[i, 1]; xticks) for i ∈ 1:nhospitals ]
    axs2 = [ Axis(gl[i, 2]; xticks) for i ∈ 1:nhospitals ]
    axs3 = [ Axis(gl[i, 3]; xticks) for i ∈ 1:nhospitals ]
    axs4 = [ Axis(gl[i, 4]; xticks) for i ∈ 1:nhospitals ]
    for (i, ax) ∈ enumerate(axs1)
        k = rankindices[i]
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_m2[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
            # Values are divided by 10 as we are plotting a cumulative prevalence of isolating
            # from work. Those who isolate do so for 10 days, so the cumulative incidence of
            # isolating is approximately this value over 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 1)
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 1)
        formataxis!(ax; hidex=(i!=nhospitals))
    end
    for (i, ax) ∈ enumerate(axs2)
        k = rankindices[i]
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_m1[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 2)
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 2)
        formataxis!(ax; hidex=(i!=nhospitals), hidey=true)
    end
    for (i, ax) ∈ enumerate(axs3)
        k = rankindices[i]
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_p1[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 3)
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 3)
        formataxis!(ax; hidex=(i!=nhospitals), hidey=true)
    end
    for (i, ax) ∈ enumerate(axs4)
        k = rankindices[i]
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_p2[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 4)
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 4)
        formataxis!(ax; hidex=(i!=nhospitals), hidey=true)
    end
    linkaxes!(axs1..., axs2..., axs3..., axs4...)
end

function _plotcumulativecounterfactualvaccinevline!(
    ax, x::Number, ::Any; 
    color=:black, linestyle=:dot
)
    vlines!(ax, x; color, linestyle)
end

function _plotcumulativecounterfactualvaccinevline!(
    ax, x::AbstractVector{<:Number}, i::Integer; 
    kwargs...
)
    _plotcumulativecounterfactualvaccinevline!(ax, x[i], nothing; kwargs...)
end

_plotcumulativecounterfactualvaccinevline!(::Any, ::Nothing, ::Any) = nothing

function _plotcumulativecounterfactualvaccinevspan!(
    ax, xs::AbstractVector{<:Real}, ::Any; 
    color=( :grey, 10 ),
)
    vspan!(ax, xs; color)
end

function _plotcumulativecounterfactualvaccinevspan!(
    ax, xs::AbstractVector{<:AbstractVector}, i::Integer; 
    kwargs...
)
    _plotcumulativecounterfactualvaccinevspan!(ax, xs[i], nothing; kwargs...)
end

_plotcumulativecounterfactualvaccinevspan!(::Any, ::Nothing, ::Any) = nothing

function _rankindicesoftotaldiagnoses(totaldiagnosis, ::Automatic)
    rankvector = ordinalrank(totaldiagnosis; rev=true)
    return [ findfirst(x -> x == i, rankvector) for i ∈ eachindex(rankvector) ]
end

function _rankindicesoftotaldiagnoses(::Any, observeddiagnosis::AbstractVector)
    return _rankindicesoftotaldiagnoses(observeddiagnosis, automatic)
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






unboostedomega100changevaccinationdatefig = let 
    fig = with_theme(theme_latexfonts()) do
        fig = Figure(; size=( 587, 800 ))
        ga = GridLayout(fig[1, 1])
        plotcumulativecounterfactualvaccine!(
            ga, 
            470:831, 
            unboostedoutputsperhospital_omega100_diagnosesafterjuly.predicteddiagnoses[470:831, :, :] , 
            unboostedoutputsperhospital_omega100_m2["predictdiagnoses"][470:831, :, :], 
            unboostedoutputsperhospital_omega100_m1["predictdiagnoses"][470:831, :, :], 
            unboostedoutputsperhospital_omega100_p1["predictdiagnoses"][470:831, :, :], 
            unboostedoutputsperhospital_omega100_p2["predictdiagnoses"][470:831, :, :];
            observeddiagnoses=unboostedobserveddiagnosesafterjuly,
            xticks=(
                [ 469, 561, 653, 743, 834 ],
                [ "July", "Oct.", "Jan.", "April", "July" ]
            ),
            vline=531,
            vspan = [ [ 531 + x, 621 + x ] for x ∈ [ -62, -31, 30, 61 ] ],
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


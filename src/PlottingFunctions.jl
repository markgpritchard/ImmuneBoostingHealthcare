
# Plotting functions are kept in a separate module so that CairoMakie does not need to be
# compiled on the servers running MCMC

module PlottingFunctions 

using DrWatson
using CairoMakie, DataFrames, PlotFormatting, StatsBase
import ImmuneBoostingHealthcare: Automatic, automatic

export COLOUR_I, COLOUR_R, COLOUR_S, estimateeffectofboosting,  
    plotcaseswithchangedvaccinationdates!, plotchains, plotchains!, 
    plotcounterfactualvaccine!, plotcumulativecounterfactualvaccine!, plothospitaloutputs, 
    plotoutputs, plotoutputs!

# Consistent colour scheme across plots 

const COLOUR_S = COLOURVECTOR[1]
const COLOUR_I = COLOURVECTOR[2]
const COLOUR_R = COLOURVECTOR[3]

# outputs from MCMC that are not plotted 
const _NOPLOTNAMES = [ 
    "iteration", "chain", "lp", "n_steps", "is_accept", "acceptance_rate", "log_density", 
    "hamiltonian_energy", "hamiltonian_energy_error", "max_hamiltonian_energy_error", 
    "tree_depth", "numerical_error", "step_size", "nom_step_size"
]

function plotchains(data::DataFrame; size=( 400, 1200 ), kwargs...)
    fig = Figure(; size)
    plotchains!(fig, data; kwargs...)
    return fig
end

function plotchains!(fig::Figure, data::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plotchains!(gl, data; kwargs...)
end

function plotchains!(
    gl::GridLayout, data::DataFrame; 
    columns=automatic, ylabels=automatic, yticks=Makie.automatic, kwargs...
)
    @unpack colnames, plotnames_ind = _processplotchains(data, columns; kwargs...)
        
    ax = [ Axis(gl[i, 1]; yticks) for i ∈ eachindex(plotnames_ind) ]
    for (j, chainid) ∈ enumerate(unique(data.chain))
        inds = findall(x -> x == chainid, data.chain)
        for (i, k) ∈ enumerate(plotnames_ind) 
            lines!(ax[i], getproperty(data, colnames[k])[inds]; color=COLOURVECTOR[j])
            _plotchainsylabel!(gl, colnames, ylabels, i, k)
        end
    end
end

function _plotchainsylabel!(gl, colnames, ::Automatic, i, k)
    Label(gl[i, 0], "$(colnames[k])"; rotation=π/2, tellheight=false)
end

function _plotchainsylabel!(gl, ::Any, ylabels::AbstractVector{<:AbstractString}, i, k)
    Label(gl[i, 0], "$(ylabels[k])"; rotation=π/2, tellheight=false)
end

function _processplotchains(data, ::Automatic; kwargs...)
    colnames = names(data)
    return _processplotchains(data, colnames; kwargs...)
end

function _processplotchains(data, colnames; logdensity="log_density")
    # "log_density" is the label given by `Pigeons` output. Turing labels it "lp".
    lp_ind = findall(x -> x == logdensity, colnames)
    _plotnames_ind = findall(x -> x ∉ _NOPLOTNAMES, colnames)
    plotnames_ind = [ lp_ind; _plotnames_ind ]
    return @ntuple colnames plotnames_ind
end

function plotoutputs(outputs)
    fig = Figure()
    plotoutputs!(fig, outputs)
    return fig
end

function plotoutputs!(fig::Figure, outputs)
    ax = Axis(fig[1, 1])
    plotoutputs!(ax, outputs)
end

function plotoutputs!(ax::Axis, outputs)
    scatter!(
        ax, 
        outputs["totalinfections"], 
        outputs["mediantotaldiagnoses"]; 
        color=COLOURVECTOR[1], markersize=3
    )
    rangebars!(
        ax, 
        outputs["totalinfections"], 
        outputs["lcitotaldiagnoses"], 
        outputs["ucitotaldiagnoses"]; 
        color=( COLOURVECTOR[1], 0.1 ),
    )
    lines!(
        ax, 
        [ extrema(outputs["totalinfections"])... ], 
        [ extrema(outputs["totalinfections"])... ];
        color=:black,
    )
end

function plothospitaloutputs(
    outputs; 
    firstplot=1, jseries=firstplot:(firstplot + 25), 
    size=( 800, 800 ),
    xticks=Makie.automatic,
    suppliedextra=true,
)
    # default for 25 plots 
    fig = Figure(; size)
    axs = [ Axis(fig[i, j]; xticks) for i ∈ 1:5, j ∈ 1:5 ]
    plotind = 0 
    for (j, code) ∈ enumerate(unique(outputs["data"].Code))
        j ∉ jseries && continue
        plotind += 1
        inds = findall(x -> x == code, outputs["data"].Code)
        for k ∈ axes(outputs["predictdiagnoses"], 3)
            lines!(
                axs[plotind],
                outputs["data"].t[inds],
                suppliedextra ? 
                    outputs["predictdiagnoses"][:, j, k] : 
                    outputs["predictdiagnoses"][:, plotind, k];
                color=( :blue, 0.1 )
            )
        end
        if j <= 19 
            if (j - 1) / 5 == round(Int, (j - 1) / 5)
                formataxis!(
                    axs[plotind]; 
                    hidespines=( :r, :t, :b, ), 
                    hidex=true, hidexticks=true,  
                )
            else
                formataxis!(
                    axs[plotind]; 
                    hidespines=( :r, :t, :l, :b, ), 
                    hidex=true, hidexticks=true, hidey=true, hideyticks=true, 
                )
            end
        else
            if (j - 1) / 5 == round(Int, (j - 1) / 5)
                formataxis!(
                    axs[plotind]; 
                    hidespines=( :r, :t ), 
                    hidex=true,
                )
            else
                formataxis!(
                    axs[plotind]; 
                    hidespines=( :r, :t, :b, ), 
                    hidey=true, hideyticks=true, 
                    hidex=true,
                )
            end
        end
        
        scatter!(
            axs[plotind], 
            outputs["data"].t[inds], 
            outputs["data"].CovidAbsences[inds];
            color=:black, markersize=3,
        )
    end

    Label(fig[1:5, 0], "Prevalence"; rotation=π/2, tellheight=false)
    Label(fig[6, 1:5], "Date"; fontsize=11.84, tellwidth=false)

    colgap!(fig.layout, 1, 5)
    rowgap!(fig.layout, 5, 5)
    
    fig
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
    color=COLOURVECTOR[1], colorrange=Makie.automatic,
)
    lines!(ax, dates, [ quantiles[i][2] for i ∈ eachindex(quantiles) ]; color, colorrange)
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
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 1)
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_m2[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
            # Values are divided by 10 as we are plotting a cumulative prevalence of isolating
            # from work. Those who isolate do so for 10 days, so the cumulative incidence of
            # isolating is approximately this value over 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 1)
        formataxis!(ax; hidex=(i!=nhospitals))
    end
    for (i, ax) ∈ enumerate(axs2)
        k = rankindices[i]
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 2)
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_m1[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 2)
        formataxis!(ax; hidex=(i!=nhospitals), hidey=true)
    end
    for (i, ax) ∈ enumerate(axs3)
        k = rankindices[i]
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 3)
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_p1[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 3)
        formataxis!(ax; hidex=(i!=nhospitals), hidey=true)
    end
    for (i, ax) ∈ enumerate(axs4)
        k = rankindices[i]
        _plotcumulativecounterfactualvaccinevspan!(ax, vspan, 4)
        plotcumulativecounterfactualvaccine!(
            ax, dates, cf_p2[:, k, :] ./ 10, modelledoriginal[:, k, :] ./ 10
        )
        hlines!(ax, 0; color=:black, linestyle=:dot)
        _plotcumulativecounterfactualvaccinevline!(ax, vline, 4)
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
    vspan!(ax, xs[1], last(xs); color)
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

function plotcaseswithchangedvaccinationdates!(
    fig::Figure, dates, observedcasesvector, modelledcasesvector, counterfactualsvector;
    xticks=(
        [ 469, 561, 653, 743, 834 ],
        [ "July", "Oct.", "Jan.", "April", "July" ]
    ),
    kwargs...
)
    axs1 = [ Axis(fig[1, i]; xticks) for i ∈ 1:4 ]
    axs2 = [ Axis(fig[2, i]; xticks) for i ∈ 1:4 ]
    axs3 = [ Axis(fig[3, i]; xticks) for i ∈ 1:4 ]
    axs4 = [ Axis(fig[4, i]; xticks) for i ∈ 1:4 ]

    for (axs, observedcases, modelledcases, counterfactuals, hidex) ∈ zip(
        [ axs1, axs2, axs3, axs4 ],
        observedcasesvector,
        modelledcasesvector,
        counterfactualsvector,
        [ true, true, true, false ]
    )
        plotcaseswithchangedvaccinationdates!(
            axs, dates, observedcases, modelledcases, counterfactuals;
            hidex, kwargs...
        )
    end
end

function plotcaseswithchangedvaccinationdates!(
    fig::Figure, dates, outputsdicts::Vector{<:Dict};
    kwargs...
)
    observedcasesvector = Vector{Vector{Float64}}(undef, 4)
    modelledcasesvector = Vector{Array{Float64, 3}}(undef, 4)

    for i ∈ 1:4
        observedcasesvector[i] = outputsdicts[i]["observationssincejuly"]
        modelledcasesvector[i] = outputsdicts[i]["modelledoutput"]["predictdiagnoses"]
    end

    plotcaseswithchangedvaccinationdates!(
        fig, dates, observedcasesvector, modelledcasesvector, outputsdicts;
        kwargs...
    )
end

function plotcaseswithchangedvaccinationdates!(
    axs::Vector{<:Axis}, dates, observedcases, modelledcases, counterfactuals;
    hidex=false,
    maxcolour=COLOURVECTOR[1],
    mincolour=COLOURVECTOR[2],
    otherhospitalscolour=( :gray, 0.1 ),
    nhospitals=23,
    vaccinationtimes=( 531, 621 ),
    vspancolour=( :gray, 0.1),
)
    # greatest and least number of cases 
    minobs, = findmax(observedcases)
    maxobs, = findmin(observedcases)

    @unpack m2, m1, p1, p2 = counterfactuals
        
    for (i, v) ∈ enumerate([ m2, m1, p1, p2 ])
        vspan!(
            axs[i],
            vaccinationtimes[1] + [ -62, -31, 30, 61 ][i],
            vaccinationtimes[2] + [ -62, -31, 30, 61 ][i];
            color=vspancolour
        )
        #=plotcumulativecounterfactualvaccine!(
            axs[i], 
            dates,
            v["predictdiagnoses"][dates, maxind, :],
            modelledcases[dates, maxind, :];
            color=maxcolour 
        )
        plotcumulativecounterfactualvaccine!(
            axs[i], 
            dates,
            v["predictdiagnoses"][dates, minind, :],
            modelledcases[dates, minind, :];
            color=mincolour
        )=#
        for j ∈ 1:nhospitals 
            plotcumulativecounterfactualvaccine!(
                axs[i], 
                dates,
                v["predictdiagnoses"][dates, j, :],
                modelledcases[dates, j, :];
                color=observedcases[j], colorrange=[ minobs, maxobs ]
            )
        end
        hlines!(axs[i], 0; color=:black, linestyle=:dot)
        vlines!(axs[i], vaccinationtimes[1]; color=:black, linestyle=:dot)
        formataxis!(axs[i]; hidex, hidey=(i != 1))
    end
    linkaxes!(axs...)
end
    
end  # module PlottingFunctions 

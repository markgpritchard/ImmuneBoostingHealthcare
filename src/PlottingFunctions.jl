
# Plotting functions are kept in a separate module so that CairoMakie does not need to be
# compiled on the servers running MCMC

module PlottingFunctions 

using DrWatson
using CairoMakie, DataFrames
import ImmuneBoostingHealthcare: Automatic, automatic

export COLOUR_I, COLOUR_R, COLOUR_S, COLOURVECTOR, formataxis!, formataxishidespines!, 
    labelplots!, plotchains, plothospitaloutputs, plotoutputs, plotoutputs!, setorigin!, 
    setvalue!

# Consistent colour scheme across plots 

const COLOURVECTOR = [
    RGBf(33 / 255, 145 / 255, 140 / 255),
    RGBf(68 / 255, 57 / 255, 131 / 255),
    RGBf(253 / 255, 231 / 255, 37 / 255),
    RGBf(53 / 255, 183 / 255, 121 / 255),
    RGBf(49 / 255, 104 / 255, 142 / 255),
    RGBf(68 / 255, 1 / 255, 84 / 255),
    RGBf(144 / 255, 215 / 255, 67 / 255),
]

const COLOUR_S = COLOURVECTOR[1]
const COLOUR_I = COLOURVECTOR[2]
const COLOUR_R = COLOURVECTOR[3]

# outputs from MCMC that are not plotted 
const _NOPLOTNAMES = [ 
    "iteration", "chain", "lp", "n_steps", "is_accept", "acceptance_rate", "log_density", 
    "hamiltonian_energy", "hamiltonian_energy_error", "max_hamiltonian_energy_error", 
    "tree_depth", "numerical_error", "step_size", "nom_step_size"
]

function plotchains(data::DataFrame; size=( 400, 1200 ), columns=automatic, kwargs...)
    @unpack colnames, plotnames_ind = _processplotchains(data, columns; kwargs...)
        
    fig = Figure(; size)
    ax = [ Axis(fig[i, 1]) for i ∈ eachindex(plotnames_ind) ]
    for (j, chainid) ∈ enumerate(unique(data.chain))
        inds = findall(x -> x == chainid, data.chain)
        for (i, k) ∈ enumerate(plotnames_ind) 
            lines!(ax[i], getproperty(data, colnames[k])[inds]; color=COLOURVECTOR[j])
            Label(fig.layout[i, 0], "$(colnames[k])"; rotation=π/2, tellheight=false)
        end
    end
    
    return fig
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


""" 
    formataxis!(axis::Axis, width = 800; <keyword arguments>)
    formataxis!(axis::Axis3, width = 800; setorigin = false)
    formataxis!(cb::Colorbar, width = 800; horizontal = false)
    formataxis!(label::Label, width = 800)
    formataxis!(legend::Legend, width = 800; horizontal = true)
    formataxis!(axes::Array, args...; kwargs...)

Apply consistent formatting to components of figures.

## Arguments 
The first argument is the element to be formatted. 

`width` refers to the width of the whole figure. Text sizes and line widths are formatted 
    to be a constant proportion of this width.

## Keyword arguments
Most keyword arguments are only available when formatting an `Axis` 
* `hidespines = ( :r, :t )`: which sides of the plot are not to be displayed: accepts 
    a tuple of symbols, `:b` for bottom, `:l` for left, `:r` for right, `:t` for top
* `hidex = false`: whether to hide values on x-axis
* `hidexticks = false`: whether to hide tick marks on x-axis (can only be hidden if `hidex = true`)
* `hidey = false`: whether to hide values on y-axis
* `hideyticks = false`: whether to hide tick marks on y-axis (can only be hidden if `hidey = true`)
* `setorigin = false`: whether axes should be extended to a value of `0`
* `setpoint = nothing`: can take a number or tuple and calls `setvalue!`
* `trimspines = true`: whether to trim to axes at the minimum and maximum tick values
The additional keyword argument `horizontal` is available when formatting a `Colorbar`
    or `Legend`, and sets whether that item should be oriented horizontally.
""" 
function formataxis!(
    axis::Axis, width=800; 
    hidex=false, hidexticks=false, hidey=false, hideyticks=false, 
    setorigin=false, setpoint=nothing, hidespines=nothing, trimspines = false
)
    if !isnothing(hidespines)
        @warn "Current preference to avoid using hidespines"
        formataxishidespines!(axis, hidespines)
    end
    axis.spinewidth = width / 800
    axis.xtrimspine = trimspines; axis.ytrimspine = trimspines
    axis.xgridvisible = false; axis.ygridvisible = false
    axis.xtickwidth = width / 800; axis.ytickwidth = width / 800
    axis.xlabelsize = width / 67; axis.ylabelsize = width / 67
    axis.xticklabelsize = width / 80; axis.yticklabelsize = width / 80
    axis.titlealign = :left; axis.titlesize = width / 65
    if setorigin setorigin!(axis) end 
    setvalue!(axis, setpoint)
    if hidex 
        hidexdecorations!(axis; ticks = hidexticks) 
    else
        if hidexticks 
            @info "Function `formataxis!` cannot hide ticks on x axis unless `hidex` is true" 
        end 
    end 
    if hidey 
        hideydecorations!(axis; ticks = hideyticks) 
    else 
        if hideyticks 
            @info "Function `formataxis!` cannot hide ticks on y axis unless `hidey` is true" 
        end 
    end 
end 

function formataxis!(axis::Axis3, width = 800; setorigin = false)
    axis.xspinewidth = width / 800; axis.yspinewidth = width / 800; axis.zspinewidth = width / 800; 
    axis.xgridvisible = false; axis.ygridvisible = false; axis.zgridvisible = false
    axis.xtickwidth = width / 800; axis.ytickwidth = width / 800; axis.ztickwidth = width / 800
    axis.xlabelsize = width / 67; axis.ylabelsize = width / 67; axis.zlabelsize = width / 67
    axis.xticklabelsize = width / 80; axis.yticklabelsize = width / 80; axis.zticklabelsize = width / 80
    axis.titlealign = :left; axis.titlesize = width / 65
    if setorigin setorigin!(axis) end
end 

function formataxis!(legend::Legend, width = 800; horizontal = true)
    legend.framevisible = false
    legend.labelsize = width / 80; legend.titlesize = width / 80
    legend.patchsize = (width / 40, width / 40)
    if horizontal
        legend.orientation = :horizontal
        legend.titleposition = :left
    else 
        legend.margin = (10, 10, 10, 10)
    end 
end 

function formataxis!(cb::Colorbar, width = 800; horizontal = false)
    cb.ticklabelsize = width / 80; cb.labelsize = width / 67
    if horizontal 
        cb.height = width / 80 
    else 
        cb.width = width / 80 
    end
end 

function formataxis!(label::Label, width = 800)
    label.fontsize  = width / 67
end 

function formataxis!(axes::Array, width = 800; kwargs...)
    for ax ∈ axes formataxis!(ax, width; kwargs...) end 
end 

"""
    setvalue!(axis, <additional arguments>)

Extends axes to include the provided value. 

## Additional arguments 
Separate arguments can be given for the `x`, `y` and `z` (for `Axis3`) axes. They 
    may also be provided as a `Tuple`.  

For an `Axis3`, the `z` value may be provided alone, which assumes `x = y = 0`. 
    For an `Axis`, the `y` value may be provided alone, which assumes `x = 0`. If 
    no additional arguments are given, the origin, is added. If `x` is supplied as 
    `nothing`, no extension to the axes is performed.

"""
setvalue!(axis::Axis, y::Real = 0) = setvalue!(axis, 0, y)
setvalue!(axis::Axis, x, y) = scatter!(axis, [x], [y], markersize = 0)
setvalue!(axis::Axis3, z::Real = 0) = setvalue!(axis, 0, 0, z)
setvalue!(axis::Axis3, x, y, z) = scatter!(axis, [x], [y], [z], markersize = 0)
setvalue!(axis, x::Nothing) = nothing 
setvalue!(axis, xy::Tuple) = setvalue!(axis, xy...)
setorigin!(axis) = setvalue!(axis)

# Function to hide spines. Not exported.

formataxishidespines!(axis, hidespines::Nothing) = nothing
formataxishidespines!(axis, hidespines::Symbol) = hidespines!(axis, hidespines)

function formataxishidespines!(axis, hidespines)
    for d ∈ hidespines hidespines!(axis, d) end
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Label subplots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

""" 
    labelplots!(labels, layouts; <keyword arguments>)

Applies labels to plots within a larger figure.

`labels` and `layouts` can each refer to an individual plot or can be a vector referring 
    to several plots.

# Keyword arguments 
* `cols = 0`: which column within the plots should have the label applied. Can be 
    provided as one integer for all plots or a vector of integers for each plot.
* `rows = 0`: which row within the plots should have the label applied. Can be provided 
    as one integer for all plots or a vector of integers for each plot.
* `font = "TeX Gyre Heros Bold"`: font for the labels
* `fontsize = 14`: label fontsize
* `halign = :left`: label alignment
* `padding = ( 0, 5, 5, 0 )`: label padding, provided as a `Tuple` in the order left, 
    right, bottom, top
"""
function labelplots!(labels, layouts; cols = 0, rows = 0, kwargs...)
    return _labelplots!(labels, layouts, rows, cols; kwargs...)
end 

function _labelplots!(labels::Vector{String}, layouts, rows::Int, cols; kwargs...)
    rowvector = zeros(Int, length(labels)) .+ rows 
    return _labelplots!(labels, layouts, rowvector, cols; kwargs...)
end 

function _labelplots!(labels::Vector{String}, layouts, rows::Vector{<:Int}, cols::Int; kwargs...)
    colvector = zeros(Int, length(labels)) .+ cols 
    return _labelplots!(labels, layouts, rows, colvector; kwargs...)
end 

function _labelplots!(labels::Vector{String}, layouts::Vector, rows::Vector{<:Int}, cols::Vector{<:Int};
        kwargs...
    )
    @assert length(labels) == length(layouts)
    for (row, col, label, layout) ∈ zip(rows, cols, labels, layouts)
        _labelplots!(label, layout, row, col; kwargs...)
    end 
end

function _labelplots!(labels::Vector{String}, layout, rows::Vector{<:Int}, cols::Vector{<:Int};
        kwargs...
    )
    for (row, col, label) ∈ zip(rows, cols, labels)
        _labelplots!(label, layout, row, col; kwargs...)
    end 
end

function _labelplots!(label::String, layout, row::Int, col::Int;
        font = "TeX Gyre Heros Bold", fontsize = 14, halign = :left, padding = ( 0, 5, 5, 0 )
    )
    Label(layout[row, col, TopLeft()], label; font, fontsize, halign, padding)
end 
    
end  # module PlottingFunctions 


# Plotting functions are kept in a separate module so that CairoMakie does not need to be
# compiled on the servers running MCMC

module PlottingFunctions 

using DrWatson
using CairoMakie, DataFrames

export COLOUR_I, COLOUR_R, COLOUR_S, COLOURVECTOR, plotchains, plotoutputs

# Consistent colour scheme across plots 

const COLOUR_S = :blue
const COLOUR_I = :darkgoldenrod1
const COLOUR_R = :seagreen4
const COLOURVECTOR = [ 
    COLOUR_S, COLOUR_I, COLOUR_R, :plum, :brown2, :dodgerblue3, :skyblue2, :lightgray 
]

# outputs from MCMC that are not plotted 
const _NOPLOTNAMES = [ 
    "iteration", "chain", "lp", "n_steps", "is_accept", "acceptance_rate", "log_density", 
    "hamiltonian_energy", "hamiltonian_energy_error", "max_hamiltonian_energy_error", 
    "tree_depth", "numerical_error", "step_size", "nom_step_size"
]

function plotchains(data::DataFrame; size=( 400, 1200 ), kwargs...)
    @unpack colnames, plotnames_ind = _processplotchains(data; kwargs...)
        
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

function _processplotchains(data; logdensity="log_density")
    # "log_density" is the label given by `Pigeons` output. Turing labels it "lp".
    colnames = names(data)
    lp_ind = findall(x -> x == logdensity, colnames)
    _plotnames_ind = findall(x -> x ∉ _NOPLOTNAMES, colnames)
    plotnames_ind = [ lp_ind; _plotnames_ind ]
    return @ntuple colnames plotnames_ind
end

function plotoutputs(outputs)
    fig = Figure()
    axs = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
    scatter!(
        axs[1], 
        outputs[:totalinfections], 
        outputs[:mediantotaldiagnoses]; 
        color=:blue, markersize=3
    )
    rangebars!(
        axs[1], 
        outputs[:totalinfections], 
        outputs[:lcitotaldiagnoses], 
        outputs[:ucitotaldiagnoses]; 
        color=( :blue, 0.1 ),
    )
    lines!(
        axs[1], 
        [ extrema(outputs[:totalinfections])... ], 
        [ extrema(outputs[:totalinfections])... ]
    )
    
    scatter!(
        axs[2], 
        outputs[:totalinfections], 
        outputs[:mediantotaldiagnosesnoboost]; 
        color=:blue, markersize=3
    )
    rangebars!(
        axs[2], 
        outputs[:totalinfections], 
        outputs[:lcitotaldiagnosesnoboost],
        outputs[:ucitotaldiagnosesnoboost]; 
        color=( :blue, 0.1 ),
    )
    lines!(
        axs[2], 
        [ extrema(outputs[:totalinfections])... ], 
        [ extrema(outputs[:totalinfections])... ]
    )
    
    scatter!(
        axs[3], 
        outputs[:totalinfections],
        (outputs[:mediantotaldiagnoses] .- outputs[:mediantotaldiagnosesnoboost]) ./ 
            outputs[:mediantotaldiagnoses]; 
        color=:blue, markersize=3,
    )
        
    return fig
end
    
end  # module PlottingFunctions 

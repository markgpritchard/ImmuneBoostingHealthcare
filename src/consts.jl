
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

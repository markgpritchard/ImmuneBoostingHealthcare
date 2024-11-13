#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=4
module load Julia/1.9.3-linux-x86_64

n_rounds=12
omega=0.00556

julia scripts/analysesimsperhospital.jl "4" "$n_rounds" "boostedsimulation" "0.00556" &
julia scripts/analysesimsperhospital.jl "4" "$n_rounds" "unboostedsimulation" "0.00556" &
julia scripts/analysesimsperhospital.jl "4" "$n_rounds" "midboostedsimulation" "0.00556" &
julia scripts/analysedataperhospital.jl "4" "$n_rounds" "0.00556" &
julia scripts/analysesimsperhospital.jl "4" "$n_rounds" "boostedsimulation" "0.01" &
julia scripts/analysesimsperhospital.jl "4" "$n_rounds" "unboostedsimulation" "0.01" &
julia scripts/analysesimsperhospital.jl "4" "$n_rounds" "midboostedsimulation" "0.01" &
julia scripts/analysedataperhospital.jl "4" "$n_rounds" "0.01" &

wait




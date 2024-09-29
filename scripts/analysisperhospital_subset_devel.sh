#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=4
module load Julia/1.9.3-linux-x86_64

n_rounds=8
omega=0.00556

for n in {1..4}
do
	julia scripts/analysesimsperhospital_subset.jl "$n" "$n_rounds" "boostedsimulation" "$omega" &
	julia scripts/analysesimsperhospital_subset.jl "$n" "$n_rounds" "unboostedsimulation" "$omega" &
        julia scripts/analysedataperhospital_subset.jl "$n" "$n_rounds" "$omega" &
done

wait




#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=4
module load Julia/1.9.3-linux-x86_64

n_rounds=12
omega=0.00556

for n in {1..4}
do
        julia scripts/analysedatafixomega.jl "$n" "$n_rounds" "$omega" &
	julia scripts/analysesimsfixomega.jl "$n" "$n_rounds" "unboostedsimulation" "$omega" &
	julia scripts/analysesimsfixomega.jl "$n" "$n_rounds" "boostedsimulation" "$omega" &
done

wait




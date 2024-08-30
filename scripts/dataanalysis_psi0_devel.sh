#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.pritchard@ndm.ox.ac.uk

export JULIA_NUM_THREADS=10
module load Julia/1.8.5-linux-x86_64

n_rounds=8

for n in {1..4}
do
	julia scripts/analysedata_psi0.jl "$n" "$n_rounds" &
done

wait




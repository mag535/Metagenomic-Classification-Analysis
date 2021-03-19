#!/bin/bash

#SBATCH --mail-user=mag535@drexel.edu
#SBATCH -A rosenGrp
### 1 core = 1 slot = 1 task
#SBATCH -N 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem=40G
#SBATCH -p def

. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc
module load slurm

DBNAME=/beegfs/scratch/mag535/KRAKEN/--_kraken_db
years=(1999 2004 2009 2014 2019 2020)

for year in ${years[@]}
do
	start=`date +%s`
	kraken2-build --build --db ${DBNAME/--/${year}}
	end=`date +%s`
	runtime=$((end-start))
	echo "The database is ${DBNAME/--/${year}}"
	echo "The runtime to build is $runtime"
done

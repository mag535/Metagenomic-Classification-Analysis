#!/bin/bash

#SBATCH --mail-user=mag535@drexel.edu
#SBATCH -A rosenGrp
### 1 core = 1 slot = 1 task
#SBATCH -N 1
### n-tasks-per-node only available in default partition
#SBATCH --ntasks-per-node 16
#SBATCH -t 02:00:00
#SBATCH --mem=100G
#SBATCH -p def

. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc
module load slurm

years=(1999 2004 2009 2014 2019 2020)
scratch_space=/beegfs/scratch/mag535
base_db=bac_kraken_db

#loop over years
for year in ${years[@]}
do
	# create year db (<year>_kraken_db)
	mkdir ${scratch_space}/KRAKEN/${year}_kraken_db
	# copy /beegfs/scratch/mag535/KRAKEN/bac_kraken_db/* into <year>_kraken_db [FIXME!!!]
	cp -r ${scratch_space}/KRAKEN/${base_db}/* ${scratch_space}/KRAKEN/${year}_kraken_db
	# rm *.k2d
	rm ${scratch_space}/KRAKEN/${year}_kraken_db/*.k2d
	# rm <year>_kraken_db/library/bacteria/libray.fna
	rm ${scratch_space}/KRAKEN/${year}_kraken_db/library/bacteria/library.fna
	# mv <year>/library.fna <year>_kraken_db/library/bacteria/
	cp ${scratch_space}/${year}/library.fna ${scratch_space}/KRAKEN/${year}_kraken_db/library/bacteria
	# signal sub-database has been built
	echo "${year}_kraken_db has been set up in ${scratch_space}/KRAKEN"
	echo "Don't forget to build it!"
done

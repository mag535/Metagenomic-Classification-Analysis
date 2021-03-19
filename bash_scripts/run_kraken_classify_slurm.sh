#!/bin/bash

#sbatch --mail-user=mag535@drexel.edu
#sbatch -A rosenGrp
### 1 core = 1 slot = 1 task
#sbatch -N 1
#sbatch --ntasks-per-node=16
### memory is per node
#sbatch --mem=210G
#sbatch -p bm


. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc
module load slurm

### load kraken2?

DBNAME=/beegfs/scratch/mag535/KRAKEN/--_kraken_db
sample=/beegfs/scratch/mag535/ZZ_real_sample.fa
years=(1999 2004 2009 2014 2019 2020)

for year in ${years[@]}
do
	start=`date +%s`
	kraken2 --db ${DBNAME/--/$year} $sample > /beegfs/scratch/mag535/kraken_${year}_real_classify_output.txt
	end=`date +%s`
	runtime=$((end-start))
	echo "The database is ${DBNAME/--/$year}"
	echo "The sample is $sample"
	echo "Classification time is $runtime"
done

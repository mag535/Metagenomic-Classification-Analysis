#!/bin/bash

#SBATCH --mail-user=mag535@drexel.edu
#SBATCH -A rosenGrp
### 1 core = 1 slot = 1 task
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

#SBATCH --mem=40G
#SBATCH -p def

. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc
module load slurm

# Clark directory with bacteria directory in it; this is in picotte's version of scratch
DBNAME=/beegfs/scratch/mag535/CLARK/yy_clark_db
script=/ifs/groups/rosenGrp/mag535/CLARKSCV1.2.6.1/set_targets.sh
type=custom
###years=(1999 2004 2009 2014 2019 2020)
year=2009

start=`date +%s`
$script ${DBNAME/yy/$year} $type
end=`date +%s`
runtime=$((end-start))
echo "The databse is ${DBNAME/yy/$year}"
echo "The runtime for building the CLARK $type database is $runtime seconds"

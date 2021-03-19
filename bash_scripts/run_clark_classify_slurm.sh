#!/bin/bash

#sbatch --mail-user=mag535@drexel.edu
#sbatch -A rosenGrp
### 1 core = 1 slot = 1 task
#sbatch -N 1
### n-tasks-per-node only available in default partition
#sbatch --mem=512G
#sbatch -p bm

. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc
module load slurm

year=2009
DBNAME=/beegfs/scratch/mag535/CLARK/yy_clark_db
script=/ifs/groups/rosenGrp/mag535/CLARKSCV1.2.6.1/classify_metagenome.sh
sample=/beegfs/scratch/mag535/ZZ_real_sample.fa
output=/beegfs/scratch/mag535/clark_yy_bac_real_classify_output

start=`date +%s`
$script -O $sample -R ${output/yy/$year} -n 24
end=`date +%s`

runtime=$((end-start))
echo "The database is ${DBNAME/yy/$year}"
echo "The sample is $sample"
echo "Classification time is $runtime"

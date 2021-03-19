#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -M mag535@drexel.edu
#$ -P rosenPrj
#$ -l h_rt=24:00:00
#$ -l h_vmem=50G
#$ -l m_mem_free=40G
### vendor
#$ -q all.q

. /etc/profile.d/modules.sh
. ~/.bashrc

### basics
module load shared
module load proteus
module load sge/univa
module load gcc

module load rosenGrp-anaconda3

step=1
py_file=/mnt/HA/groups/rosenGrp/mag535/python_files/kraken_db.py
parent_dir=/lustre/scratch/mag535
### 1999-20006 have already been done
years=(2011 2016 2019)
folds=(1 2 3 4 5)

for year in ${years[@]}
do
	for f in ${folds[@]}
	do
		python $py_file $parent_dir $year $f $step
	done
done

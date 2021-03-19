#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -M mag535@drexel.edu
#$ -P rosenPrj
#$ -l h_rt=24:00:00
#$ -l h_vmem=50G
#$ -l m_mem_free=40G
### -t 1999-2006:1
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

fold=1
###year=$SGE_TASK_ID
years=(2011 2016 2019)
py_file=/mnt/HA/groups/rosenGrp/mag535/python_files/known.py
classification_results=/lustre/scratch/mag535/yy_kraken_fold${fold}_classify_output.txt
library_headers=/lustre/scratch/mag535/Kraken/yy_kraken_db/library/bacteria/library_headers.fna
output_path=/lustre/scratch/mag535

for year in ${years[@]}
do
	python $py_file ${classification_results/yy/$year} ${library_headers/yy/$year} $output_path $year $fold
	echo "Ran $python for $year"
done

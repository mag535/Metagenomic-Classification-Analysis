#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -M mag535@drexel.edu
#$ -P rosenPrj
#$ -l h_rt=24:00:00
#$ -l h_vmem=60G
#$ -l m_mem_free=50G
### vendor
#$ -q all.q

. /etc/profile.d/modules.sh
. ~/.bashrc

### basics
module load shared
module load proteus
module load sge/univa
module load gcc

### others
module load rosenGrp-anaconda3

### to do
py_file=/mnt/HA/groups/rosenGrp/mag535/python_files/cross_validation.py
### 2020 not in ZZ's dir
years=(1999 2000 2001 2002 2003 2004 2005 2006 2011 2016 2019)
### should be only 5 folds
folds=(1 2 3 4 5)
core=/lustre/scratch/mag535/year_exp/{}/{}
headers=/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria/library_headers.fna
output_path=/lustre/scratch/mag535

for year in ${years[@]}
do
	for f in ${folds[@]}
	do
		start=`date +%s`
		python $py_file $year "fold$f" $core $headers $output_path
		end=`date +%s`
		runtime=$((end-start))
		echo "Creating ${year}_fold${f} training data took $runtime"
	done
done

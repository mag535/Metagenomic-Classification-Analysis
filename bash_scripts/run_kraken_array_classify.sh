#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$-j y
#$ -M mag535@drexel.edu
#$ -P rosenPrj
#$-l h_rt=05:00:00
#$-l h_vmem=210G
#$ -l m_mem_free=200G
### -t 1999-2006:1
### -l vendor=amd
#$ -q all.q

. /etc/profile.d/modules.sh
. ~/.bashrc

### basic
module load shared
module load proteus
module load sge/univa
module load gcc

### others
module load rosenGrp-kraken2

### to do
fold=5
###year=$SGE_TASK_ID
years=(2011 2016 2019)
parent_dir=/lustre/scratch/mag535
DBNAME=${parent_dir}/Kraken/yy_kraken_db
sample=/lustre/scratch/mag535/ZZ_testing_data/testing_data_fold_${fold}.fna

for year in ${years[@]}
do
	start=`date +%s`
	kraken2 --db ${DBNAME/yy/$year} $sample > ${parent_dir}/${year}_kraken_fold${fold}_classify_output.txt
	end=`date +%s`
	runtime=$((end-start))

	echo "The database is ${DBNAME/yy/$year}"
	echo "The sample is $sample"
	echo "The excluded fold is fold$fold"
	echo "Classification time is $runtime"
done

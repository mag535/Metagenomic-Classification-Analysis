#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -M mag535@drexel.edu
#$ -P rosenPrj
#$ -pe shm 15
#$ -l h_rt=48:00:00
#$ -l h_vmem=7G
#$ -l m_mem_free=6G
### -t 1999-2006:1
### vendor
#$ -q all.q

. /etc/profile.d/modules.sh
. ~/.bashrc

### These four modules must ALWAYS be loaded
module load shared
module load proteus
module load sge/univa
module load gcc

module load rosenGrp-kraken2
module load ncbi-blast/gcc/64/2.6.0

fold=1
###year=$SGE_TASK_ID
years=(2011 2016 2019)
DIR=/lustre/scratch/mag535
DBNAME=${DIR}/Kraken/yy_kraken_db

for year in ${years[@]}
do
	rm ${DBNAME/yy/$year}/*.k2d
	cp ${DIR}/${year}_kraken/combined_data_${year}_fold${fold}.fna ${DBNAME/yy/$year}/library/bacteria/library.fna
	grep ">" ${DBNAME/yy/$year}/library/bacteria/library.fna > ${DBNAME/yy/$year}/library/bacteria/library_headers.fna

	start=`date +%s`
	kraken2-build --build --db ${DBNAME/yy/$year} --threads 15
	end=`date +%s`
	runtime=$((end-start))

	echo "Database is ${DBNAME/yy/$year}"
	echo "Fold excluded: fold$fold"
	echo "Building is $runtime seconds"
done

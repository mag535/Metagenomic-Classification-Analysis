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

### to do
kraken_db=/lustre/scratch/mag535/Kraken/1999_kraken_db
destination=/lustre/scratch/mag535/Kraken/--_kraken_db
## no 2020 dir in ZZ's scratch
years=(2000 2001 2002 2003 2004 2005 2006 2011 2016 2019)

for year in ${years[@]}
do
	mkdir ${destination/--/$year}
	cp -r ${kraken_db}/* ${destination/--/$year}
done

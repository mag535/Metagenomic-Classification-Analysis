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

cp -r /lustre/scratch/mag535/Kraken/bac_kraken_db /lustre/scratch/mag535/Kraken/1999_kraken_db
rm /lustre/scratch/mag535/Kraken/1999_kraken_db/library/bacteria/library.fna
rm /lustre/scratch/mag535/Kraken/1999_kraken_db/*.k2d

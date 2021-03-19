#!/bin/bash

#sbatch --mail-user=mag535@drexel.edu
#sbatch -A rosenGrp
### 1 core = 1 slot = 1 task
#sbatch -N 1
#sbatch --n-tasks-per-node=16
#sbatch --mem=6G
#sbatch -p def

. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc

KRAKEN_DIR=/beegfs/scratch/KRAKEN

### install kraken2
kraken2/kraken2-master/install_kraken2.sh $KRAKEN_DIR
### copy important scripts
cp $KRAKEN_DIR/kraken2{,-build,-inspect} $HOME/bin

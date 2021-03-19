#!/bin/bash

#sbatch --mail-user=mag535@drexel.edu
#sbatch -A rosenGrp
### 1 core = 1 slot = 1 task
#sbatch -N 1
#sbatch --n-tasks-per-node=16
#sbatch --mem=40G
#sbatch -p bm


. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc

# Clark directory with bacteria directory in it; this is in picotte's version of scratch
DBNAME=/beegfs/scratch/KRAKEN/bac_kraken_db
type=bacteria

start=`date +%s`
###kraken2-build --standard --threads 24 --db $DBNAME
kraken2-build --download-taxonomy --db $DBNAME
kraken2-build --download-library $type --db $DBNAME
kraken2-build --build --db $DBNAME
end=`date +%s`

runtime=$((end-start))
echo "The runtime for building the Kraken 2 $type database is $runtime seconds"

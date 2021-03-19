#!/bin/bash

#sbatch --mail-user=mag535@drexel.edu
#sbatch -A rosenGrp
### 1 core = 1 slot = 1 task
#sbatch -N 1
#sbatch --ntasks-per-node=16
### memory is per node
#sbatch --mem=60G
#sbatch -p def


. /etc/profile.d/modules.sh
. ~/.bashrc

module load shared
module load gcc
module load slurm

## python?
module load python37
### did these, but don't know if needs to be done each time I log in
# /cm/local/apps/python37/bin/python3 -m pip install --upgrade pip
# pip install pandas

py_file=/ifs/groups/rosenGrp/mag535/python_files/clark_db.py
year=2019
end_date=12/31${year}
start_date=-1/-1/-1
clark_db=/beegfs/scratch/mag535/CLARK/${year}_clark_db
output_path=/beegfs/scratch/mag535/CLARK

python $py_file $end_date $start_date $clark_db $output_path
### run setup_clark_sub-databases_slurm.sh after to copy genomes into Custom of sub-database

#!/bin/bash

#sbatch --mail-user=mag535@drexel.edu
#sbatch -A rosenGrp
### 1 core = 1 slot = 1 task
#sbatch -N 1
#sbatch --ntasks-per-node=16
### memory is per node
#sbatch --mem=60G
#sbatch -p def

module load shared
module load gcc
module load slurm

## python?
module load python37
### did these, but don't know if needs to be done each time I log in
# /cm/local/apps/python37/bin/python3 -m pip install --upgrade pip
# pip install pandas

py_file=/ifs/groups/rosenGrp/mag535/python_files/clark_db.py
years=(1999 2004 2009 2014 2019 2020)
DIR=/beegfs/scratch/mag535/CLARK
clark_db=${DIR}/bac_clark_db
end=12/31/yy
start=-1/-1/-1

for year in ${years[@]}
do
	mkdir ${DIR}/${year}_clark_db
	mkdir ${DIR}/${year}_clark_db/Custom
	python $py_file ${end/yy/$year} $start $clark_db $DIR
	input=${DIR}/sub_db_scrap/${year}_genome_paths.txt
	n=0
	while read -r line
	do
		cp $line ${DIR}/${year}_clark_db/Custom
		n=$((n+1))
	done < $input
	echo "$n genomes copied"
	echo "$year sub-database has been set up"
done

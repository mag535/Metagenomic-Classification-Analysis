#!/bin/bash

#SBATCH --mail-user=mag535@drexel.edu
#SBATCH --account=rosenMRIPrj
### 1 core = 1 slot = 1 task
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 48:00:00
#SBATCH --mem=160G
#SBATCH -p def

. /etc/profile.d/modules.sh
. ~/.bashrc

# DEFAULTS
module load shared
module load gcc
module load slurm

# OTHERS
module load python37
# PIP, for the python modules I needs
###/cm/local/apps/python37/bin/python3 -m pip install --upgrade pip
###pip install pandas
###pip install numpy
###pip install matplotlib
###pip install openpyxl
###pip install xlsxwriter
###pip install scipy
###pip install seaborn

# variables
DIR=/ifs/group/rosenMRIGrp/mag535
py_file=${DIR}/python_files/clark_graphing.py
input_path=${DIR}/classification_results
output_path=${DIR}/clark_graphing_outputs
year=(1999 2004 2009 2014 2019 2020)
level=species
who=c
uncl=1
excel=0
plot=2

start=`date +%s`
python $py_file $input_path $output_path $level $who $uncl $excel $plot ${year[@]}
end=`date +%s`

runtime=$((end-start))

if [$who -eq "c"]
then
        echo "This took $runtime seconds for CLARK"
elif [$who -eq "k"]
then
        echo "This took $runtime seconds for Kraken 2"
else
        echo "Something went wrong in printing who: $who"
fi

echo "Years:  ${year[@]}"
echo "Rank: $level"

if [$plot -eq 1]
then
        echo "Tested function: plot_relative_abundance()"
elif [$plot -eq 2]
then
        echo "Tested function: plot_triangular_bray_curtis()"
else
        echo "Something went wrong in printing tested function: $plot"
fi

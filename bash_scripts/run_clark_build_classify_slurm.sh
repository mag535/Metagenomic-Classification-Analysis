### bash array job with dependency
### 'afterany' can be used instead of 'afterok', but will run second job even if first exits with error

build_clark_db=/ifs/groups/rosenGrp/mag535/scripts/run_clark_download_slurm.sh
classify_with_clark_db=/ifs/groups/rosenGrp/mag535/scripts/run_clark_classify_slurm.sh

build_ID=$(sbatch $build_clark_db | sed 's/Submitted batch job //')
echo "Build ID: $build_ID"
classify_ID=$(sbatch --depend=afterok:$build_ID $classify_with_clark_db | sed 's/Submitted batch job //')
echo "First Classify ID: $classify_ID"
second_classify_ID=$(sbatch --depend=afterok:$classify_ID $classify_with_clark_db | sed 's/Submitted batch job //')
echo "Second Classify ID: $second_classify_ID"

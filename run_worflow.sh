#!/bin/sh
### Job name
#SBATCH --job-name=annotation

### Requirements
#SBATCH --partition=long



### Output
#SBATCH --output=/shared/home/tdurand/annotation_fijiensism2/workflow_effectors/log_effector.out
#SBATCH --error=/shared/home/tdurand/annotation_fijiensism2/workflow_effectors/log_effector.err
module purge
module load r
module load python/3.7


snakemake --profile effector --use-envmodules

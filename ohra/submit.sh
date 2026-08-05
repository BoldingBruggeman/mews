#!bin/bash   

#SBATCH --ntasks-per-node=10
#SBATCH --nodes=1
#SBATCH --job-name=pygetm-run-1 
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/work/%u/logs/%x-%j.log

module load Conda/25.3.1
module load foss/2024a
conda activate pygetm-env
python3 run.py

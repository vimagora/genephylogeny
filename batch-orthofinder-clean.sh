#!/bin/bash

#SBATCH --job-name=orthofinder_clean              # Job name
#SBATCH --output=logs/%x_%j.out         # Output file for stdout
#SBATCH --error=logs/%x_%j.err          # Output file for stderr
#SBATCH --ntasks=1                     # Number of tasks (cores)
#SBATCH --cpus-per-task=40             # Number of cores per task
#SBATCH --time=3-00:00:00              # Time limit hrs:min:sec
#SBATCH --account=project_2002833      # Project number
#SBATCH --mem-per-cpu=2G               # Memory to reserve
#SBATCH --partition=small              # Job queue (partition)


# Define variables for paths
ORTHOFINDER_CMD="/scratch/project_2002833/VG/software/OrthoFinder/bin/orthofinder"     # Path to the OrthoFinder executable
PROTEOMES_DIR="local_data/proteomes/clean/"         # Path to the directory containing the proteomes

# Print paths for debugging purposes (optional)
echo "OrthoFinder command: $ORTHOFINDER_CMD"
echo "Proteomes directory: $PROTEOMES_DIR"
echo "Output directory: $OUTPUT_DIR"

# Load the biopython module and biokit
module load biopythontools
module load biokit

# Run OrthoFinder
$ORTHOFINDER_CMD -f $PROTEOMES_DIR -t 40 -a 10

echo "OrthoFinder run finished."
echo "Results can be found in: $OUTPUT_DIR"

#!/bin/bash -e
 
#SBATCH --job-name=2022_08_17.SVCalling_starlingwgs_delly_step1.sl
#SBATCH --account=uoa02613
#SBATCH --time=00-24:00:00
#SBATCH --mem=10GB
#SBATCH --output=%x_%j_a%.errout
#SBATCH --mail-user=katarina.stuart@auckland.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --profile task
#SBATCH --array=1-55

# load modules
module purge
module load Delly/1.1.3

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/mapped_reads/starling_wgs_mapped/sample_names_starling_55_wgs.txt)
echo "working with sample:" ${SAMPLE}

# set paths
DIR=/nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/mapped_reads/starling_wgs_mapped
BAM=${DIR}/${SAMPLE}.sorted.dup.bam
GENOME=/nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/resources/Svulgaris_vAU_1.0.fasta
OUT_DIR=/nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/variant_calling/SV_delly

#Run delly step 1
delly call -o ${OUT_DIR}/${SAMPLE}.bcf -g ${GENOME} ${BAM}

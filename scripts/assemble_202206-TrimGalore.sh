#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=8  # number of cpus requesting on node <=20
#SBATCH --time=23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=TrimGalore
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/TrimGalore-20220616-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/TrimGalore-20220616-%A-%a.err
#SBATCH --mem=24000
#
#SBATCH --array=1-38

# Job array size will be half of the number of lines in the file

PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"
cd $SCRATCH_DIR/trimmed_reads

SECOND=$((SLURM_ARRAY_TASK_ID*2))
FIRST="$(($SECOND - 1))"

R1=$(sed -n "$FIRST"p "$PROJECT_DIR"scripts/assemble_202206_TrimGalore_filepaths.txt) # 76 files = 38 array size
R2=$(sed -n "$SECOND"p "$PROJECT_DIR"scripts/assemble_202206_TrimGalore_filepaths.txt)

$LOCAL_SRC/TrimGalore-0.6.6/trim_galore --paired --phred33 --output_dir $SCRATCH_DIR/trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 -j 8 $R1 $R2

#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=FiltUncorrP
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/FiltUncorrP-20220616-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/FiltUncorrP-20220616-%A-%a.err
#SBATCH --mem=2000
#
#SBATCH --array=1-38

# Job array size will be half of the number of lines in the file

PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
cd $PROJECT_DIR/assemble/Rcorrector/

SECOND=$((SLURM_ARRAY_TASK_ID*2))
FIRST="$(($SECOND - 1))"

R1=$(sed -n "$FIRST"p "$PROJECT_DIR"scripts/assemble_202206_FiltUncorrPE_filepaths.txt) # 76 files = 38 array size
R2=$(sed -n "$SECOND"p "$PROJECT_DIR"scripts/assemble_202206_FiltUncorrPE_filepaths.txt)
SAMPLE_NAME=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$PROJECT_DIR"scripts/assemble_202206_FiltUncorrPE_SampleNames.txt)

python $LOCAL_SRC/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 $R1 -2 $R2 -s $SAMPLE_NAME

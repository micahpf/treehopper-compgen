#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=Rcorrector
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/Rcorrector-20220616-%k-%j.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/Rcorrector-20220616-%k-%j.err
#SBATCH --mem=24000
#
#SBATCH --array=1-38

# Job array size will be half of the number of lines in the file

PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
cd $PROJECT_DIR/assemble/Rcorrector/

SECOND=$((SLURM_ARRAY_TASK_ID*2))
FIRST="$(($SECOND - 1))"

LINE1=$(sed -n "$FIRST"p "$PROJECT_DIR"scripts/assemble_202206_rcorrector_filepaths.txt) # 76 files = 38 array size
LINE2=$(sed -n "$SECOND"p "$PROJECT_DIR"scripts/assemble_202206_rcorrector_filepaths.txt)

perl $LOCAL_SRC/rcorrector/run_rcorrector.pl -t 12 -1 $LINE1 -2 $LINE2
#perl $LOCAL_SRC/rcorrector/run_rcorrector.pl -t 12 -1 $PROJECT_DIR/reads/fastq/INHS/Agcon_ACTTGA_L001_R1_001.fastq.gz -2 $PROJECT_DIR/reads/fastq/INHS/Agcon_ACTTGA_L001_R1_001.fastq.gz

#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=6-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=rRNA_blacklist
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/rRNA_blacklist-20220623-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/rRNA_blacklist-20220623-%A-%a.err
#SBATCH --mem=24000
#
#SBATCH --array=22,4

# Job array size will be half of the number of lines in the file

PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"
cd $SCRATCH_DIR/rRNA_blacklist

SECOND=$((SLURM_ARRAY_TASK_ID*2))
FIRST="$(($SECOND - 1))"

R1=$(sed -n "$FIRST"p "$PROJECT_DIR"scripts/assemble_202206-rRNA_blacklist_filepaths.txt) # 76 files = 38 array size
R2=$(sed -n "$SECOND"p "$PROJECT_DIR"scripts/assemble_202206-rRNA_blacklist_filepaths.txt)
SAMPLE_NAME=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$PROJECT_DIR"scripts/assemble_202206-rRNA_blacklist_SampleNames.txt)

# $R1 = R1 fastq file
# $R2 = R2 fastq file
# $SAMPLE_NAME = sample_id

bowtie2 --quiet --very-sensitive-local --phred33 --threads 20 \
    -x $PROJECT_DIR/../../databases/silva/SILVA_BOTH_NR99 \
    -1 $R1 -2 $R2 --met-file ${SAMPLE_NAME}_bowtie2_metrics.txt \
    --al-conc-gz blacklist_paired_aligned_${SAMPLE_NAME}.fq.gz \
    --un-conc-gz blacklist_paired_unaligned_${SAMPLE_NAME}.fq.gz  \
    --al-gz blacklist_unpaired_aligned_${SAMPLE_NAME}.fq.gz \
    --un-gz blacklist_unpaired_unaligned_${SAMPLE_NAME}.fq.gz

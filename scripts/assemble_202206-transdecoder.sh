#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=1  # number of cpus requesting on node <=20
#SBATCH --time=6-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=transdecoder
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/transdecoder-20220624i-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/transdecoder-20220624i-%A-%a.err
#SBATCH --mem=5000
#
#SBATCH --array=14

export PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"

export SAMPLE_NAME=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$PROJECT_DIR"scripts/assemble_202206-assembly_names.txt)
cd $SCRATCH_DIR/trinity/${SAMPLE_NAME}/

#TransDecoder.LongOrfs -t ${SAMPLE_NAME}.Trinity.fasta 
#TransDecoder.Predict -t ${SAMPLE_NAME}.Trinity.fasta

TransDecoder.LongOrfs -t trinity_genes.fasta
TransDecoder.Predict -t trinity_genes.fasta

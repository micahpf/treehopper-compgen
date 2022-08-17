#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=10  # number of cpus requesting on node <=20
#SBATCH --time=6-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=super
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/trinity-super-20220623-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/trinity-super-20220623-%A-%a.err
#SBATCH --mem=100000
#
#SBATCH --array=14

export PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"

export SAMPLE_NAME=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$PROJECT_DIR"scripts/assemble_202206-assembly_names.txt)
cd $SCRATCH_DIR/trinity/${SAMPLE_NAME}/reads

#R1=$(ls -1p *.fq.1.gz | grep -v / | xargs echo | sed 's/ /,/g')
#R2=$(ls -1p *.fq.2.gz | grep -v / | xargs echo | sed 's/ /,/g')

R1=$(readlink -f *.fq.1.gz | xargs echo | sed 's/ /,/g')
R2=$(readlink -f *.fq.2.gz | xargs echo | sed 's/ /,/g')

cd $SCRATCH_DIR/trinity/${SAMPLE_NAME}

#Trinity \
#          --seqType fq \
#          --left $SCRATCH_DIR/trinity/${SAMPLE_NAME}/reads/$R1  \
#          --right $SCRATCH_DIR/trinity/${SAMPLE_NAME}/reads/$R2 \
#          --max_memory 100G --CPU 20 \
#          --output $SCRATCH_DIR/trinity/${SAMPLE_NAME}/trinity_out

export SINGULARITY_BINDPATH="/scratch/tmp/micahf"

#singularity exec -e /Genomics/kocherlab/micahf/singularity_images/trinity/trinityrnaseq.v2.14.0.simg Trinity \
#    --verbose --seqType fq \
#    --left $R1 \
#    --right $R2 \
#    --max_memory 100G --CPU 10 --output $SCRATCH_DIR/trinity/${SAMPLE_NAME}/trinity_out

singularity exec -e /Genomics/kocherlab/micahf/singularity_images/trinity/trinityrnaseq.v2.14.0.simg \
    /usr/local/bin/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py \
    --trinity_fasta ${SAMPLE_NAME}.Trinity.fasta

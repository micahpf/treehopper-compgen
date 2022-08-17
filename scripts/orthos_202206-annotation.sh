#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=23:00:00 --qos=1day # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=EnTAP-tr1
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/EnTAP-20220725-tree1-%A.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/EnTAP-20220725-tree1-%A.err
#SBATCH --mem=10000

## 2022/07/22

# Annotate orthogroups using EnTAP
# - Pick a single representative sequence from each orthogroup (based on % gaps in cds) using custom python script
# - Translate into pep using (same) custom script
# - Annotate using EnTAP

export PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"
export DB_DIR="/Genomics/kocherlab/micahf/databases/"
cd ${SCRATCH_DIR}

## Get representative sequence for each orthogroup
cd ${SCRATCH_DIR}/EnTAP
#cut -f1 ../selection_pipeline/aaml_compiled/selection_202206-tree1_t10_aamls.txt > selection_202206-tree1_OG_names.txt

# OG2gene.py [OG_dir] [OG_IDs_file] [OG_seqs_out]
# outputs a fasta with cds for each OG 'selection_202206-tree1_OG_seqs.cds'
#python ${PROJECT_DIR}/scripts/python_utils/OG2gene.py \
#    ${SCRATCH_DIR}/selection_pipeline/selection_202206-tree1_fsa_coding/ \
#    ${SCRATCH_DIR}/EnTAP/selection_202206-tree1_OG_names.txt \
#    ${SCRATCH_DIR}/EnTAP/selection_202206-tree1_OG_seqs.cds \
#    ${SCRATCH_DIR}/EnTAP/selection_202206-tree1_OG_seqs.pep

## EnTAP
# No need to run the first command, configged before but here for ref
#$LOCAL_SRC/EnTAP-v0.9.1-beta/EnTAP --config \
#    --paths ${SCRATCH_DIR}/EnTAP/entap_config.txt \
#    -d ${DB_DIR}/nr/nr.fasta \
#    -d ${DB_DIR}/uniprot_sprot/uniprot_sprot.fasta \
#    --out-dir ${SCRATCH_DIR}/EnTAP \
#    -t 20

# This is here in case I need to re-configure
#$LOCAL_SRC/EnTAP-v0.9.1-beta/EnTAP --runP \
#    -i ${SCRATCH_DIR}/EnTAP/selection_202206-tree1_OG_seqs.cds.transdecoder.pep \
#    --paths ${SCRATCH_DIR}/EnTAP/entap_config.txt \
#    -d ${SCRATCH_DIR}/EnTAP/bin/nr.dmnd \
#    -d ${SCRATCH_DIR}/EnTAP/bin/uniprot_sprot.dmnd \
#    --out-dir ${SCRATCH_DIR}/EnTAP/EnTAP_20220722 \
#    -t 20

# RUN THIS
module load hyphy

$LOCAL_SRC/EnTAP-v0.9.1-beta/EnTAP --runP \
    -i ${SCRATCH_DIR}/EnTAP/selection_202206-tree1_OG_seqs.pep \
    --paths /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/EnTAP/entap_config.txt \
    -d /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/EnTAP/bin/nr.dmnd \
    -d /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/EnTAP/bin/uniprot_sprot.dmnd \
    --out-dir ${SCRATCH_DIR}/EnTAP/EnTAP_20220722 \
    -t 20

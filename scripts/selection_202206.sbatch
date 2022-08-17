#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=2-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --mail-user=micahf@princeton.edu
#SBATCH --job-name=BUSCO
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/BUSCO_Llanq-20220616-%j.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/BUSCO_Llanq-20220616-%j.err
#SBATCH --mem=75000

## Add a leafhopper outgroup (Agallia) and Melizoderid to the dataset
# Adding the outgroup solves the problem of arbitrarily bifurcating the (AETA,everything_else) branch, plus gives a reference for dN/dS calculations
# Adding Melizoderid (if assembly is high quality) might not be necessary, as it would increase number of ant foreground without increasing mom foreground (6 vs 5), but I want to check BUSCO

# First steps:
# - Get BUSCO scores for both assemblies
# - Copy and rename assemblies into new directory for this project
# - OrthoFinder for full list of assemblies, check how many useable orthogroups
# - Run RAxML and check that species tree makes sense
# - If all goes well, move on to RERconverge and RELAX

# Setup
DATA_HOME="/Genomics/kocherlab/micahf/data"
cd /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206

# BUSCO

module load ncbi-blast
module load augustus
export PYTHONPATH="${PYTHONPATH}:/Genomics/kocherlab/micahf/local/src/BUSCO/src/"
BUSCO_CONFIG_FILE="/Genomics/kocherlab/micahf/local/src/BUSCO/config/config.ini"
cd /scratch/tmp/micahf/membracoid_transcriptomes/BUSCO
#python $LOCAL_SRC/BUSCO/scripts/run_BUSCO.py \
#       --in /scratch/tmp/micahf/membracoid_transcriptomes/treehopper_Transdecoder/Agallia_constricta.free.accepted.fa.transdecoder.pep \
#       --out Agallia_constrica_prot \
#       --lineage_path $LOCAL_SRC/BUSCO/insecta_odb9/ \
#       --mode proteins \
#       -f

#python $LOCAL_SRC/BUSCO/scripts/run_BUSCO.py \
#       --in /scratch/tmp/micahf/membracoid_transcriptomes/treehopper_Transdecoder/Llanquihuea_pilosa.free.accepted.fa.transdecoder.pep \
#       --out Llanquihuea_pilosa_prot \
#       --lineage_path $LOCAL_SRC/BUSCO/insecta_odb9/ \
#       --mode proteins \
#       -f

# Both have middling BUSCO scores (C:75%), but I will proceed with Orthofinder and see how it performs
# I wonder if Trinity could assemble these better than the original assemblies using SOAPdenovo ca 2016?
# As a dry-run to see if reassembling from raw reads generates better assemblies I'll try to do that for Agallia



#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=1  # number of cpus requesting on node <=20
#SBATCH --time=6-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=fix_ids
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/fix_ids-20220624-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/fix_ids-20220624-%A-%a.err
#SBATCH --mem=5000
#
#SBATCH --array=1-25

export PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"

THIRD=$((SLURM_ARRAY_TASK_ID*3))
SECOND=$(($THIRD - 1))
FIRST="$(($SECOND - 1))"

SP_ID=$(sed -n "$FIRST"p "$PROJECT_DIR"scripts/orthos_202206-assembly_names+paths.txt) # 75 lines = 25 array size
CDS=$(sed -n "$SECOND"p "$PROJECT_DIR"scripts/orthos_202206-assembly_names+paths.txt)
PEP=$(sed -n "$THIRD"p "$PROJECT_DIR"scripts/orthos_202206-assembly_names+paths.txt)
cd $SCRATCH_DIR/fixed_ids/genes

#python ${PROJECT_DIR}scripts/python_utils/rename_genes.py $SP_ID ${SCRATCH_DIR}trinity/${CDS} ${SCRATCH_DIR}trinity/${PEP}

if [[ $SP_ID != 'UMBO' ]]
then
    python ${PROJECT_DIR}scripts/python_utils/filter_braker.py 300 no_braker ${SP_ID}.fixed_ids.cds ${SP_ID}.fixed_ids.pep # array=1-24
elif [[ $SP_ID = 'UMBO' ]]
then
    python ${PROJECT_DIR}scripts/python_utils/filter_braker.py 300 braker ${SP_ID}.fixed_ids.cds ${SP_ID}.fixed_ids.pep # array=25
else
    echo "if/else didn't work..."
fi

#python ${PROJECT_DIR}scripts/python_utils/filter_braker.py 300 no_braker ${SP_ID}.fixed_ids.cds ${SP_ID}.fixed_ids.pep # array=1-24
#python ${PROJECT_DIR}scripts/python_utils/filter_braker.py 300 braker ${SP_ID}.fixed_ids.cds ${SP_ID}.fixed_ids.pep # array=25

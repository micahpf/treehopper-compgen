#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=10  # number of cpus requesting on node <=20
#SBATCH --time=223:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=busco
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/busco-20220624-%A-%a.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/busco-20220624-%A-%a.err
#SBATCH --mem=50000
#
#SBATCH --array=14

export PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"

export SAMPLE_NAME=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$PROJECT_DIR"scripts/assemble_202206-assembly_names.txt)
cd $SCRATCH_DIR/trinity/${SAMPLE_NAME}/

module load ncbi-blast
module load augustus
export PYTHONPATH="${PYTHONPATH}:/Genomics/kocherlab/micahf/local/src/BUSCO/src/"
BUSCO_CONFIG_FILE="/Genomics/kocherlab/micahf/local/src/BUSCO/config/config.ini"

python $LOCAL_SRC/BUSCO/scripts/run_BUSCO.py --in ${SAMPLE_NAME}.Trinity.fasta --out ${SAMPLE_NAME}.Trinity \
    --lineage_path $LOCAL_SRC/BUSCO/insecta_odb9/ --mode transcriptome -f -c 10

python $LOCAL_SRC/BUSCO/scripts/run_BUSCO.py --in ${SAMPLE_NAME}.Trinity.fasta.transdecoder.pep --out ${SAMPLE_NAME}.Trinity.fasta.transdecoder.pep \
    --lineage_path $LOCAL_SRC/BUSCO/insecta_odb9/ --mode proteins -f -c 10

python $LOCAL_SRC/BUSCO/scripts/run_BUSCO.py --in trinity_genes.fasta.transdecoder.pep --out ${SAMPLE_NAME}.SuperTranscripts \
    --lineage_path $LOCAL_SRC/BUSCO/insecta_odb9/ --mode proteins -f -c 10

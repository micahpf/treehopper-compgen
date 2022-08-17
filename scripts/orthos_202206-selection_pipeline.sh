#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=6:00:00 --qos=1day # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=rer-tr1
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/selection-20220721-rer-tree1-%A.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/selection-20220721-rer-tree1-%A.err
#SBATCH --mem=20000

### Run Ben's selection pipeline to get from orthofinder results to alignments and gene trees used in RELAX and RERConverge

export PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"
cd ${SCRATCH_DIR}
source activate python2.7

## 2022/07/15
# The first step to using Ben's pipeline is to convert the HOG output from the new version of OrthoFinder into the same format as in the older version of OrthoFinder that Ben's pipeline can read:
#cd ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jul13_1/
#cut -f1,4-23 Phylogenetic_Hierarchical_Orthogroups/N0.tsv | tr ", " " " | tr -s "\t" " " | sed "s/./&:/13" > HOG_N0.txt # if write_orthos can handle ID that start N1.HOG0000000 rather than OG0000000
#cut -f1,4-23 Phylogenetic_Hierarchical_Orthogroups/N0.tsv | tr ", " " " | tr -s "\t" " " | sed "s/./&:/13" | cut -c 5- | tail -n +2 > HOG_N0.txt # if need OG000000
#mkdir ${PROJECT_DIR}/orthofinder/tree1/
#cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv ${PROJECT_DIR}/orthofinder/tree1/
#cp HOG_N0.txt ${PROJECT_DIR}/orthofinder/tree1/
#
#cd ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jul13_2/
#cut -f1,4-23 Phylogenetic_Hierarchical_Orthogroups/N0.tsv | tr ", " " " | tr -s "\t" " " | sed "s/./&:/13" > HOG_N0.txt # if write_orthos can handle ID that start N1.HOG0000000 rather than OG0000000
#cut -f1,4-23 Phylogenetic_Hierarchical_Orthogroups/N0.tsv | tr ", " " " | tr -s "\t" " " | sed "s/./&:/13" | cut -c 5- | tail -n +2 > HOG_N0.txt # if need OG000000
#mkdir ${PROJECT_DIR}/orthofinder/tree2/
#cp Phylogenetic_Hierarchical_Orthogroups/N0.tsv ${PROJECT_DIR}/orthofinder/tree2/
#cp HOG_N0.txt ${PROJECT_DIR}/orthofinder/tree2/
#
#cp -r ${SCRATCH_DIR}orthofinder/Results_20220627_n20 ${PROJECT_DIR}/orthofinder/Results_20220627_n20 &

# Now we need a params file with the paths to the cds files for each of our assemblies. I think our pep files need to be in the same directory with the same file stems.
# To keep things backed-up in case the scratch dir gets erased, cp the assemblies to the project dir
#cd $PROJECT_DIR
#mkdir seqs
#cp /scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/fixed_ids/genes/* seqs/
# params file:
#cd seqs
#ls -d "$PWD"/*filtbrak.cds > selection_202206.params # manually edited so the species ID is the first field followed by the path after a space
#mv selection_202206.params ../scripts

## 2022/07/14 
#mkdir selection_pipeline
cd ${SCRATCH_DIR}/selection_pipeline/
#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree1 \
#    -r ${PROJECT_DIR}/orthofinder/tree1/HOG_N0.txt \
#    -t 2 \
#    -a write_orthos \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params

#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree1 \
#    -r ${PROJECT_DIR}/orthofinder/tree1/HOG_N0.txt \
#    -t 2 \
#    -p 20 \
#    -a align_coding \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params

#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree1 \
#    -r ${PROJECT_DIR}/orthofinder/tree1/HOG_N0.txt \
#    -t 10 \
#    -p 20 \
#    -a alignment_filter \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params \
#    --nogap_min_count 6

#/Genomics/kocherlab/berubin/local/src/raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAWAG \
#    -o AETA \
#    -T 20 \
#    -s /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/treehoppers_20191108.afa \
#    -n treehoppers_20191108

#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910 \
#    -o treehoppers_20191108 \
#    -r /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/orthofinder_seqs_20191108/Results_Nov08/Orthogroups.txt \
#    -p 20 -t 10 \
#    -a rer_converge \
#    -d treehoppers_20191108.params \
#    --outputfile treehoppers_20191108_t10_aamls.txt \
#    --taxa_inclusion /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/empty_file.txt \
#    -e /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/RAxML_bestTree.treehoppers_20191108 &> treehoppers_20191108_t10_aamls.log

## 2022/07/17 
#mkdir selection_pipeline
cd ${SCRATCH_DIR}/selection_pipeline/
#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree2 \
#    -r ${PROJECT_DIR}/orthofinder/tree2/HOG_N0.txt \
#    -t 2 \
#    -a write_orthos \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params

#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree2 \
#    -r ${PROJECT_DIR}/orthofinder/tree2/HOG_N0.txt \
#    -t 2 \
#    -p 20 \
#    -a align_coding \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params

#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree2 \
#    -r ${PROJECT_DIR}/orthofinder/tree2/HOG_N0.txt \
#    -t 10 \
#    -p 20 \
#    -a alignment_filter \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params \
#    --nogap_min_count 6

## 2022/07/20
cd ${SCRATCH_DIR}/selection_pipeline/
# Using the predefined tree topology to contrain gene trees

# tree2
# Manually copied tree used for orthofinder and changed node names into the 4-char ID, and added branch lengths = 1 to create tree2_20220713_shortnames_b1.txt
# (AGAL:1,((AETA:1,LOPH:1):1,(((LYCO:1,MICR:1):1,(NESS:1,CENT:1):1):1,(LLAN:1,(PROC:1,((NOTO:1,(UMBO:1,(MEMB:1,(TYLO:1,CAMP:1):1):1):1):1,((CYPH:1,STIC:1):1,(CHEL:1,(ENTY:1,AMAS:1):1):1):1):1):1):1):1):1);

#python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
#    -b ${SCRATCH_DIR}/selection_pipeline/ \
#    -o selection_202206-tree2 \
#    -r ${PROJECT_DIR}/orthofinder/tree2/HOG_N0.txt \
#    -t 10 -p 20 \
#    -a rer_converge \
#    -d ${PROJECT_DIR}/scripts/selection_202206.params \
#    --outputfile selection_202206-tree2_t10_aamls.txt \
#    --taxa_inclusion ${SCRATCH_DIR}/selection_pipeline/empty_file.txt \
#    -e ${PROJECT_DIR}/selection/tree2_20220713_shortnames_b1.txt

# tree1

python2.7 $LOCAL_SRC/ComparativeGenomics-devel/selection_pipeline.py \
    -b ${SCRATCH_DIR}/selection_pipeline/ \
    -o selection_202206-tree1 \
    -r ${PROJECT_DIR}/orthofinder/tree1/HOG_N0.txt \
    -t 10 -p 20 \
    -a rer_converge \
    -d ${PROJECT_DIR}/scripts/selection_202206.params \
    --outputfile selection_202206-tree1_t10_aamls.txt \
    --taxa_inclusion ${SCRATCH_DIR}/selection_pipeline/empty_file.txt \
    -e ${PROJECT_DIR}/selection/tree1_20220713_shortnames_b1.txt

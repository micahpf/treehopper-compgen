#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=2-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --mail-user=micahf@princeton.edu
#SBATCH --job-name=BUSCO
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/BUSCO_Llanq-20220616-%j.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/BUSCO_Llanq-20220616-%j.err
#SBATCH --mem=75000

## 2022/06/24

## Add a leafhopper outgroup (Agallia) and Melizoderid to the dataset
# Adding the outgroup solves the problem of arbitrarily bifurcating the (AETA,everything_else) branch, plus gives a reference for dN/dS calculations
# Adding Melizoderid (if assembly is high quality) might not be necessary, as it would increase number of ant foreground without increasing mom foreground (6 vs 5), but I want to check BUSCO

# First steps:
# - Get BUSCO scores for both assemblies
# - Copy and rename assemblies into new directory for this project
# - **OrthoFinder for full list of assemblies, check how many useable orthogroups
# - Run RAxML and check that species tree makes sense
# - If all goes well, move on to RERconverge and RELAX

# Setup
DATA_HOME="/Genomics/kocherlab/micahf/data"
PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"
cd /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206

#### Collect files, rename sequence IDs to be compatible with Ben's selection pipeline
### Note that Umbonia file is taken from a previous analysis (25 June 2018)

sbatch orthos_202206-fix_ids.sh

#### Run orthofinder

# OrthoFinder (2.2.7)
#cd /Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/orthofinder_seqs
#cp seqs/filtbrak_nr100/*.pep orthofinder_results/
sbatch orthos_202206-orthofinder.sh

# This is the latest version using the new HOG method 
sbatch orthos_202206-orthofinder.sh

## 2022/06/27
# The previous runs used 4 outgroup species and two species with poor assemblies (HETE and TOLA) which may have caused significant ortholog drop out
# HOLD also represents a species for which we have no behavioral data
# I may also need to apply an a priori species tree based on the AHE phlyo to fix some misplaced branches

# New run after removing the three outgroups besides AGAL, HETE, TOLA and HOLD
# The result is a tree with 18 ingroup species and 1 outgroup (AGAL)
# 1 care, no ants (UMBO)
# 5 care, ants (AETA, LOPH, ENTY, AMAS, CHEL) (If LOPH is good, I may remove CHEL later because of dubious evidence of care)
# 7 no care, ants (LYCO, CENT, NOTO, MEMB, TYLO, CAMP, CYPH)
# 5 no care, no ants (MICR, NESS, LLAN, PROC, STIC)
cd ${SCRATCH_DIR}orthofinder/
source activate scipy-numpy
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -f . -S diamond -t 20 -a 2

# Finally, I will re-run with a fixed species tree (and remove HOLD for lack of behavior data?)
# tree: (AGAL,((AETA,LOPH),(LLAN,((LYCO,(MICR,(NESS,CENT))),(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));
$LOCAL_SRC/OrthoFinder_source/orthofinder.py -f . -S diamond -t 20 -a 2 \
    -ft ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jun27/ \
    -s ${SCRATCH_DIR}orthofinder/Species_Tree_Fixed_20220627.txt

# The result is a tree with 18 ingroup species and 1 outgroup (AGAL)
# 1 care, no ants (UMBO)
# 5 care, ants (AETA, LOPH, ENTY, AMAS, CHEL) (If LOPH is good, I may remove CHEL later because of dubious evidence of care)
# 7 no care, ants (LYCO, CENT, NOTO, MEMB, TYLO, CAMP, CYPH)
# 5 no care, no ants (MICR, NESS, LLAN, PROC, STIC)

## 2022/07/13
# I wan't to account for two likely phylogenetic scenarios: (Aet,(Mel,Mem)) and (Aet,(Mel+Mem))
# The former is the expectation from Dietrich 2017 and from monophyly of families
# The latter is the (poorly supported) topology of Yanghui's draft tree, with Mel arising within Membracidae
# I'm not sure which will end up having more support in the final phylogeny dataset/analysis but I will inclube both scenarios so I don't have to redo everything later
# Plus it will be interesting/important to note the sensitivity to tree topology (it might not matter much at all in the end)
# Also I fixed the fact that Stegaspidinae wasn't monophyletic in the previous version

cd ${SCRATCH_DIR}orthofinder/
source activate scipy-numpy

# Option one: (AGAL,((AETA,LOPH),(LLAN,(((LYCO,MICR),(NESS,CENT)),(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));
mkdir tree1
echo '(AGAL,((AETA,LOPH),(LLAN,(((LYCO,MICR),(NESS,CENT)),(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));' > tree1/tree1_20220713.txt
# I manually added .fixed_ids_filtbrak to the end of each ID in the tree file before proceeding
$LOCAL_SRC/OrthoFinder_source/orthofinder.py -S diamond -t 20 -a 2 \
    -ft ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jun27 \
    -s ${SCRATCH_DIR}orthofinder/tree1/tree1_20220713.txt
mv Results_Jul13_1 ../tree1

# Option two: (AGAL,((AETA,LOPH),(((LYCO,MICR),(NESS,CENT)),(LLAN,(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));
mkdir tree2
echo '(AGAL,((AETA,LOPH),(((LYCO,MICR),(NESS,CENT)),(LLAN,(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));' > tree2/tree2_20220713.txt
# I manually added .fixed_ids_filtbrak to the end of each ID in the tree file before proceeding
$LOCAL_SRC/OrthoFinder_source/orthofinder.py -S diamond -t 20 -a 2 \
    -ft ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jun27 \
    -s ${SCRATCH_DIR}orthofinder/tree2/tree2_20220713.txt
mv Results_Jul13_2 ../tree2

# Option three: just let orthofinder estimate the phylogeny (maybe it will turn out to be one or the other scenario and will simplify methods/justification
mkdir tree3
$LOCAL_SRC/OrthoFinder_source/orthofinder.py -f . -S diamond -t 20 -a 2

## 2022/07/15
# tree3 creates a terrible tree, so I'm abandoning that approach
# For tree1 and tree2, we need to use the N1 HOG (i.e. excluding AGAL) for downstream analyses, although if we remove AGAL at this stage and run write_orthos, align_orthos and create aRAxML phylo, we again won't be using an outgroup for the phylogeny, which defeats the purpose of this reanalysis. 
# Since I alrady have a species tree I want to use, should I skip those steps and just jump to the selection_pipeline -a raxml command?
# That depends on whether I need all the fitering steps from align_coding and align_filter for relax. I think align_coding and align_filter are just for creating the concatenated alignment for RAxML, and relax uses CodeML to generate gene trees constrained by topology of the resulting species tree. I'm not sure how the relax function approaches alignment and filtering. 
# If the previous steps do the alignment and filtering for each orthogroup and subsequently concatenates them, then I can run those and just skip the RAxML step and use the pre-defined species trees: tree1 and tree2. I think this is how it works.

# So I should actually use N0 to ensure that the gene trees have an outgroup? But what if AGAL is missing for a particular orthogroup? Wait, are CodeML trees even rooted? Should I just include AGAL in all the analyses? Maybe I can remove AGAL later in the analysis? Using N0 for now.



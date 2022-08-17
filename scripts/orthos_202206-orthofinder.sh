#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=6-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --job-name=orthofinder
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/orthofinder-20220713_tree1_hogs-%A.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/orthofinder-20220713_tree1_hogs-%A.err
#SBATCH --mem=100000

export SCRATCH_DIR="/scratch/tmp/micahf/membracoid_transcriptomes/selection_202206/"
cd ${SCRATCH_DIR}orthofinder/

#$LOCAL_SRC/OrthoFinder-2.2.7/orthofinder -f . -S diamond -t 20 -a 2 # this is the older (deprecated) version using the less accurate ortho clustering method

source activate scipy-numpy
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -f . -S diamond -t 20 -a 2 # This is the latest version using the new HOG method

# New run after removing the three outgroups besides AGAL, HETE, TOLA and HOLD
# The result is a tree with 18 ingroup species and 1 outgroup (AGAL)
# 1 care, no ants (UMBO)
# 5 care, ants (AETA, LOPH, ENTY, AMAS, CHEL) (If LOPH is good, I may remove CHEL later because of dubious evidence of care)
# 7 no care, ants (LYCO, CENT, NOTO, MEMB, TYLO, CAMP, CYPH)
# 5 no care, no ants (MICR, NESS, LLAN, PROC, STIC)
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -f . -S diamond -t 20 -a 2 \
#    -o Results_20220627_n20

# Finally, I will re-run with a fixed species tree (and remove HOLD for lack of behavior data?)
# tree: (AGAL,((AETA,LOPH),(LLAN,((LYCO,(MICR,(NESS,CENT))),(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -S diamond -t 20 -a 2 \
#    -ft ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jun27/ \
#    -s ${SCRATCH_DIR}orthofinder/Species_Tree_Fixed_20220627.txt

## 2022/07/13
# I wan't to account for two likely phylogenetic scenarios: (Aet,(Mel,Mem)) and (Aet,(Mel+Mem))
# The former is the expectation from Dietrich 2017 and from monophyly of families
# The latter is the (poorly supported) topology of Yanghui's draft tree, with Mel arising within Membracidae
# I'm not sure which will end up having more support in the final phylogeny dataset/analysis but I will inclube both scenarios so I don't have to redo everything later
# Plus it will be interesting/important to note the sensitivity to tree topology (it might not matter much at all in the end)
# Also I fixed the fact that Stegaspidinae wasn't monophyletic in the previous version


# Option one: (AGAL,((AETA,LOPH),(LLAN,(((LYCO,MICR),(NESS,CENT)),(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));
#mkdir tree1
#echo '(AGAL,((AETA,LOPH),(LLAN,(((LYCO,MICR),(NESS,CENT)),(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));' > tree1/tree1_20220713.txt
# I manually added .fixed_ids_filtbrak to the end of each ID in the tree file before proceeding
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -S diamond -t 20 -a 2 \
#    -ft ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jun27 \
#    -s ${SCRATCH_DIR}orthofinder/tree1/tree1_20220713.txt

# Need to prep fasta files from hogs before moving on to Ben's pipeline
#python $LOCAL_SRC/OrthoFinder_source/tools/create_files_for_hogs.py Results_20220627_n20/Results_Jul13_1 Results_20220627_n20/Results_Jul13_1/write_hogs/ 'N1'

# Option two: (AGAL,((AETA,LOPH),(((LYCO,MICR),(NESS,CENT)),(LLAN,(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));
#mkdir tree2
#echo '(AGAL,((AETA,LOPH),(((LYCO,MICR),(NESS,CENT)),(LLAN,(PROC,(((NOTO,(UMBO,(MEMB,(TYLO,CAMP))))),((CYPH,STIC),(CHEL,(ENTY,AMAS)))))))));' > tree2/tree2_20220713.txt
# I manually added .fixed_ids_filtbrak to the end of each ID in the tree file before proceeding
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -S diamond -t 20 -a 2 \
#    -ft ${SCRATCH_DIR}orthofinder/Results_20220627_n20/Results_Jun27 \
#    -s ${SCRATCH_DIR}orthofinder/tree2/tree2_20220713.txt

python $LOCAL_SRC/OrthoFinder_source/tools/create_files_for_hogs.py Results_20220627_n20/Results_Jul13_2 Results_20220627_n20/Results_Jul13_2/write_hogs/ 'N1'

# Option three: just let orthofinder estimate the phylogeny (maybe it will turn out to be one or the other scenario and will simplify methods/justification)
#mkdir tree3
#$LOCAL_SRC/OrthoFinder_source/orthofinder.py -f . -S diamond -t 20 -a 2


#!/bin/bash
#SBATCH -n 1  # almost always one node
#SBATCH --cpus-per-task=20  # number of cpus requesting on node <=20
#SBATCH --time=2-23:00 --qos=1wk # time defaults to one hour; shorter jobs usually run sooner
#SBATCH --mail-user=micahf@princeton.edu
#SBATCH --job-name=BUSCO
#SBATCH --output=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_out/BUSCO_Llanq-20220616-%j.out
#SBATCH --error=/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/sbatch_err/BUSCO_Llanq-20220616-%j.err
#SBATCH --mem=75000

## 2022/06/16

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
PROJECT_DIR="/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_202206/"
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

##### Assemble transcriptomes
# Pipeline based on: https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
cd $PROJECT_DIR

#### Read quality
### - Run fastqc on raw reads
# fastqc ran previously, see results summarized using multiqc
ls -R reads/fastqc
source activate multiqc
cd reads/fastqc/multiqc
multiqc ../INHS ../NCBI ../Princeton

### - Remove erroneous k-mers from Illumina paired-end reads
cd $PROJECT_DIR

mkdir assemble assemble/Rcorrector 
cd assemble/Rcorrector

# Used a separate script to run an array job: assemble_202206-rcorrector.sh, which uses assemble_202206_rcorrector_filepaths.txt to get paths to each file
# The core command looks like this:
# perl $LOCAL_SRC/rcorrector/run_rcorrector.pl -t 12 -1 $PROJECT_DIR/reads/fastq/INHS/Agcon_ACTTGA_L001_R1_001.fastq.gz -2 $PROJECT_DIR/reads/fastq/INHS/Agcon_ACTTGA_L001_R1_001.fastq.gz

module load gcc/4.9.4
sbatch $PROJECT_DIR/scripts/assemble_202206-rcorrector.sh

### - Discard read pairs for which one of the reads is deemed unfixable
# This is from https://github.com/harvardinformatics/TranscriptomeAssemblyTools.git
# Again calling another script to run this in parallel
# python /PATH/TO/FilterUncorrectabledPEfastq.py -1 $1 -2 $2 -s $3

cd $PROJECT_DIR/assemble/Rcorrector/
source activate python2.7
sbatch $PROJECT_DIR/scripts/assemble_202206-FilterUncorrectabledPEfastq.sh

### - Trim adapter and low quality bases from fastq files
# From https://github.com/FelixKrueger/TrimGalore/archive/0.6.0.zip

cd $PROJECT_DIR/assemble/Rcorrector/
mkdir ../trimmed_reads
cd trimmed_reads
source activate cutadaptenv
sbatch $PROJECT_DIR/scripts/assemble_202206-TrimGalore.sh

### - Map trimmed reads to a blacklist to remove unwanted (rRNA reads) 
# Download rRNA database, concatenate long strand and short strand and replace U with T to get into DNA space
cd $PROJECT_DIR/../../databases/silva
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSRef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/README.txt

gunzip *.gz
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVA_138.1_BOTH_NR99_tax_silva.fasta
sed '/^[^>]/s/U/T/g' SILVA_138.1_BOTH_NR99_tax_silva.fasta > SILVA_138.1_BOTH_NR99_tax_silva_DNA.fasta

# There was a problem with the base environment re tbb and bowtie2
conda create -n bt2 bioconda bowtie2
source activate bt2
bowtie2-build --threads 20 SILVA_138.1_BOTH_NR99_tax_silva_DNA.fasta SILVA_BOTH_NR99

mkdir $SCRATCH_DIR/rRNA_blacklist
cd $SCRATCH_DIR/rRNA_blacklist

sbatch $PROJECT_DIR/scripts/assemble_202206-rRNA_blacklist.sh


### - Run fastqc on your processed reads that pass qc and filtering from the above steps
cd $SCRATCH_DIR
mkdir rRNA_removed_fastqc
cd rRNA_removed_fastqc
fastqc --threads 20 $SCRATCH_DIR/rRNA_blacklist/blacklist_paired_unaligned_* -o ./
source activate multiqc
multiqc .

# 2022/06/22
#### Run trinity
# Make a separate dir for each transcriptome assembly, with all reads linked in a subdir and another subdir for assembly
cd $SCRATCH_DIR
mkdir trinity; cd trinity
cd $SCRATCH_DIR/trinity/
for LINE in $(cat $PROJECT_DIR/scripts/assemble_202206-assembly_names.txt)
do mkdir $LINE; mkdir $LINE/reads; mkdir $LINE/trinity_out;
done

# I then manually soft-linked all read files for each assembly into trinity/reads
# run trinity from script
cd $PROJECT_DIR/scripts
module purge
module load singularity
sbatch assemble_202206-Trinity.sh

# 2022/06/23
# LOPH run errored out complaining about an unexpected end of file in rRNA_blackist file. There is no problem with trimmed file so the corruption happened while blacklisting
# Trying rRNA blacklisting again:

source activate bt2
sbatch $PROJECT_DIR/scripts/assemble_202206-rRNA_blacklist.sh # array=22
# I accidentally removed all the output files in this directory D:
# Re-run all of them so I have them for future reference
sbatch $PROJECT_DIR/scripts/assemble_202206-rRNA_blacklist.sh # array=1-38

### Run BUSCO on de novo trinity assemblies
# First I manually changed the names of the trinity assemblies to correspond to the $SAMPLE_NAMES I've been using for directory names, like so:
#mv UMBO/trinity_out.Trinity.fasta UMBO/UMBO.Trinity.fasta
#mv UMBO/trinity_out.Trinity.fasta.gene_trans_map UMBO/UMBO.Trinity.fasta.gene_trans_map

# Then set up busco and run in parallel
sbatch assemble_202206-BUSCO.sh

# This was for the nucleotide output of Trinity. To be comparable to previous runs I need to first use TransDecoder to find ORFs and predict prot sequences
source activate transdecoder
sbatch assemble_202206-transdecoder.sh

sbatch assemble_202206-BUSCO.sh

### Get SuperTranscripts
# Collapse transcripts into SuperTranscripts representing single genes. This may help reduce duplicates and spurious paralog calls in BUSCO or OrthoFinder
# Editted the Trinity script to do this
module purge
module load singularity
sbatch assemble_202206-Trinity.sh

# And now run BUSCO on this output (after transdecoder)
source activate transdecoder
sbatch assemble_202206-transdecoder.sh

sbatch assemble_202206-BUSCO.sh

## UMBO assembly was created by using transcriptome assembly to annotate the geneome assembly using BRAKER and taking the predicted pep file
# This was performed previously (25 June 2018) using the following code.

#cd /scratch/tmp/micahf/membracoid_transcriptomes/Umbo+Tylo/annotation/braker
#
#module add BRAKER1
#module add exonerate
#module add ncbi-blast
#module rm bioperl

# Build index for Umbonia genome (already made; move to working dir)
#cp /scratch/tmp/micahf/membracoid_transcriptomes/Umbo+Tylo/Umbo_trans_genome-guided/UCRA_genome_v1.0.FINAL.hisat2-idx* \
#   /scratch/tmp/micahf/membracoid_transcriptomes/Umbo+Tylo/annotation/

# Map reads to Umbonia genome and convert to sorted bam
#$HISAT2_HOME/hisat2 --max-intronlen 100000 -p 10 \
#       -x /scratch/tmp/micahf/membracoid_transcriptomes/Umbo+Tylo/annotation/UCRA_genome_v1.0.FINAL.hisat2-idx --phred33 \
#       -1 $DATA_HOME/transcriptomes/Membracoidea/reads_Kocher/Umbonia/Umbo_180327_R1_fixednames.fastq,\
#       $DATA_HOME/transcriptomes/Membracoidea/reads_Dietrich/aetalionidae+membracidae/Umcra_ACTTGA_L008_R1_001_fixednames.fastq \
#       -2 $DATA_HOME/transcriptomes/Membracoidea/reads_Kocher/Umbonia/Umbo_180327_R2_fixednames.fastq,\
#       $DATA_HOME/transcriptomes/Membracoidea/reads_Dietrich/aetalionidae+membracidae/Umcra_ACTTGA_L008_R2_001_fixednames.fastq \
#       | samtools view -bS - | samtools sort -o /scratch/tmp/micahf/membracoid_transcriptomes/Umbo+Tylo/annotation/braker/Umbonia_old+new_sorted.bam

#/Genomics/kocherlab/micahf/local/src/BRAKER_v2.0.4/braker.pl \
#        --genome=/Genomics/kocherlab/micahf/data/genomes/umbonia_genome/UCRA_genome_v1.0.FINAL.fasta \
#        --bam=/scratch/tmp/micahf/membracoid_transcriptomes/Umbo+Tylo/annotation/braker/Umbonia_old+new_sorted.bam \
#        --species=Umbonia_crassicornis_BRAKER_v2.0.4  --gff3 --cores=16 \
#        --AUGUSTUS_CONFIG_PATH=/Genomics/kocherlab/micahf/local/src/augustus-3.3.1/config/ \
#        --GENEMARK_PATH=/usr/local/GeneMark/4.10 \
#        --overwrite

## Use augustus.codingseq (=SuperTranscripts, roughly) and augustus.aa (=SuperTranscripts, roughly)


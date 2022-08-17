#!/usr/bin/env python3
# Input list of ortho groups and return fasta of one representative seq for each orthogroup (based on % missingness)
import sys
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


#OG_dir = '/scratch/tmp/micahf/membracoid_transcriptomes/selection/treehoppers_MSPMout_fsa_coding/'
#OG_IDs_file = 'treehopper_OG_names_201802.csv'
#OG_dir = '/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/treehoppers_20191108_fsa_coding/'
#OG_IDs_file = '/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/OG_names_20191108.txt'

OG_dir = sys.argv[1]
OG_IDs_file = sys.argv[2]
OG_seqs_cds = sys.argv[3]
OG_seqs_pep = sys.argv[4]

os.chdir(OG_dir)

OG_ID_prefix_len = len('OG_')

with open(OG_IDs_file, 'r') as f:
        reader = csv.reader(f)
        OG_IDs_list = list(reader)

OG_IDs_list = OG_IDs_list[1:]
ogids = [x[0][OG_ID_prefix_len:] for x in OG_IDs_list]

reprGenes_cds = list()
#reprGenes_pep = list()

for ogid in ogids:
    completeness = list()
    afa_name = "og_cds_%s.afa" % ogid
    print(ogid)
    # save the most complete ortholog as the representative gene
    for rec_idx, seq_record in enumerate(SeqIO.parse(afa_name, "fasta")):
        completeness.append(seq_record.seq.count("-")/len(seq_record.seq))
    
    max_completeness = max(completeness)
    max_index = completeness.index(max_completeness)
    #print(max_index)
    most_complete_seq = list(SeqIO.parse(afa_name, "fasta"))[max_index].seq
    
    # omit the first run of missing values (blast doesn't like seqs starting w missings)
    # 	if str(most_complete_seq).find == 0:
    # 		first_nonmiss = [i for i, ltr in enumerate(str(most_complete_seq)) if ltr.isalpha()][0]
    # 		most_complete_seq = most_complete_seq[first_nonmiss:]
    
    reprGene_cds = SeqRecord(Seq(str(list(SeqIO.parse(afa_name, "fasta"))[max_index].seq)),
            id="OG_%s|c%s_g1_i1" % (ogid, ogid),
            name="OG_%s|c%s_g1_i1" % (ogid, ogid),
            description="OG_%s|c%s_g1_i1" % (ogid, ogid))
    
    # ungap entire sequence (does this work better for transdecoder and annotating?) [why did I do this? wouldn't it introduce frameshifts?]
    reprGene_cds.seq = reprGene_cds.seq.ungap('-')
    reprGenes_cds.append(reprGene_cds)
    
    #reprGene_pep = SeqRecord(seq=reprGene_cds.seq.translate(to_stop=True),
    #                         id="trans_" + reprGenes_cds.id,
    #                         description="translation of CDS, using default table",
    #                        )
    #reprGenes_pep.append(reprGenes_pep)


#SeqIO.write(reprGenes, "/Genomics/kocherlab/micahf/membracoid_transcriptomes/selection_201910/trinotate_20191108/OG_seqs_20191108.cds", "fasta")
SeqIO.write(reprGenes_cds, OG_seqs_cds, "fasta")

# Translate coding seqs
def make_protein_record(cds_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(
        seq=cds_record.seq.translate(to_stop=True),
        id="trans_" + cds_record.id,
        description="translation of CDS, using default table",
    )

reprGenes_pep = (
    make_protein_record(cds_rec)
    for cds_rec in SeqIO.parse(OG_seqs_cds, "fasta")
)
SeqIO.write(reprGenes_pep, OG_seqs_pep, "fasta")

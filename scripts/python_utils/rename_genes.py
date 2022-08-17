#!/usr/bin/env python3
# Rename fasta sequence IDs to be compatible with Ben's ComparativeGenomics pipeline
# Write renamed files to pwd
# Run: python rename_genes.py [sp_id] [cds_file] [pep_file]

import sys
from Bio import SeqIO
from timeit import default_timer
from time import strftime, localtime, asctime

#sp_ids = list(('ARET', 'AMAS', 'CLAT', 'CSP1', 'CCLA', 'HSPM', 'LSPM', 'LBUR', 'MTEC', 'MCAR',
#	'NGIB', 'NSP1', 'PSP1', 'SBIS'))

#input_files = list(('Aetalion_reticulatum_accepted_assembly.fsa.transdecoder.pep',
#	'Amastris_sp_accepted_assembly.fsa.transdecoder.pep',
#	'Campylencha_lat.free.accepted.fa.transdecoder.pep',
#	'Chelyoidea_sp_accepted_assembly.fsa.transdecoder.pep',
#	'Cyphonia_clavata_accepted_assembly.fsa.transdecoder.pep',
#	'Heteronotus_spMe_032416.fsa.transdecoder.pep',
#	'Lophyraspis_spM.free.accepted.fa.transdecoder.pep',
#	'Lycoderes_burmeisteri_accepted_assembly.fsa.transdecoder.pep',
#	'Membracis_tectigera_accepted_assembly.fsa.transdecoder.pep',
#	'Microcentrus_caryae_032416.fsa.transdecoder.pep',
#	'Nessorhinus_gib.free.accepted.fa.transdecoder.pep',
#	'Notocera_sp_accepted_assembly.fsa.transdecoder.pep',
#	'Procyrta_sp_accepted_assembly.fsa.transdecoder.pep',
#	'Stictocephala_bis.free.accepted.fa.transdecoder.pep'))

#sp_ids = list(('EFAB', 'LPIL', 'MSPM', 'XSPC'))

#input_files = list(('Empoasca_fab.free.accepted.fa.transdecoder.cds',
#	'Llanquihuea_pilosa.free.accepted.fa.transdecoder.cds',
#	'Mapuchea_spM8.free.accepted.fa.transdecoder.cds',
#	'Xestocephalus_spCC_032416.fsa.transdecoder.cds'))

# sp_ids = list(('ACUT', 'AMAS', 'CAMP', 'CHEL', 'CYPH', 'HETE', 'MEMB', 'MICR', 'NOTO', 'PROC', 'STIC', 'TOLA', 'TYLO', 'UMBO'))

# input_files = list(('ACUT.pep','AMAS.pep', 'CAMP.pep', 'CHEL.pep', 'CYPH.pep', 'HETE.pep', 'MEMB.pep', 'MICR.pep', 'NOTO.pep', 'PROC.pep', 'STIC.pep', 'TOLA.pep', 'TYLO.pep', 'UMBO.pep'))

#sp_ids = list(('ENTY', 'AETA', 'LYCO', 'NESS', 'MICR', 'LYCO', 'LOPH'))

#input_files = list(('ENTY.pep', 'AETA.pep', 'LYCO.pep', 'NESS.pep', 'MICR.pep', 'LYCO.pep', 'LOPH.pep'))

sp_id = sys.argv[1]
cds_file = sys.argv[2]
pep_file = sys.argv[3]

#file_ext = '.pep'
file_type = 'fasta'	
input_files = [cds_file, pep_file]
file_exts = [".cds", ".pep"]

for count1, input_file in enumerate(input_files):
    start_time = default_timer()
    file_ext = file_exts[count1]
    
    renamed_file = str(sp_id + '.fixed_ids' + file_ext)
    
    with open(input_file) as original, open(renamed_file, 'w') as renamed:
	    records = list(SeqIO.parse(input_file, file_type))
	    for count2, record in enumerate(records):
    #		print(record.id)
		    new_record = str(sp_id + '_') + str(count2+1).zfill(5)
		    record.id = new_record
		    record.description = new_record
    #		print(record.id)
	    SeqIO.write(records, renamed, file_type)
    
    print('Renaming took ' + str(default_timer() - start_time) + 's')

print('Done')

#!/usr/bin/env python3
# Filter out seqs below a 'min_length' threshold, and fix braker cds fuckery
### Execute by: python filter_braker.py min_length=[integer] from_braker=[braker/no_braker] fastafile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import numpy as np
import os

# input = open('Amastris_test.fsa', 'r')
#input = open(sys.argv[1], 'r')

min_length = int(sys.argv[1])
from_braker = sys.argv[2]
cds_file = sys.argv[3]
pep_file = sys.argv[4]

# get number of seqs in cds file
#num_seqs = 0
#for record in SeqIO.parse(cds_file):
#    num_seqs += 1



if from_braker == 'braker':
    print("Fixing braker's fuckery")
    # find seqs where stop codon is in cds but no corresponding stop symbol in pep
    stop_codons = []
    for count,record in enumerate(SeqIO.parse(cds_file, 'fasta')):
        in_seq = record.seq.upper()
        if in_seq[-3:] == 'TAA' or in_seq[-3:] == 'TAG' or in_seq[-3:] == 'TGA':
            stop_codons.append(True)
        elif in_seq[-3:] != 'TAA' and in_seq[-3:] != 'TAG' and in_seq[-3:] != 'TGA':
            stop_codons.append(False)
            #print(record.id + ': no stop codon detected')
        else:
            print('Something went wrong while checking '+ record.id +' for STOP codons!')

    stops_fixed = []
    # check that * is indeed missing in pep and add if necessary
    for count,record in enumerate(SeqIO.parse(pep_file, 'fasta')):
        if record.seq[-1:] == '*':
            stop_fixed = record
            stops_fixed.append(stop_fixed)
        elif record.seq[-1:] != '*':
            stop_fixed = record
            # only add * if not already a * in pep
            if stop_codons[count]:
                stop_fixed.seq = record.seq + '*'
                #print(record.id + ': stop codon added')
            elif stop_codons[count] != True:
                pass
            else:
                print('Failed to properly add stop to '+ record.id +'!')
            stops_fixed.append(stop_fixed)
        else:
            print('Failed to properly add stop to '+ record.id +'!')

    ## Write the filtered seqs
    pep_in_filename, pep_ext = os.path.splitext(pep_file)

    # Create a file in the same dir where you ran this script
    with open(pep_in_filename + "_stopsfixed" + pep_ext, "w+") as stopsfixed_file:
        for record in stops_fixed:
            stopsfixed_file.write(record.format("fasta"))
    print('Stops fixed')
    
    # use the fixed pep file in the future
    original_pep_file = pep_file
    pep_file = pep_in_filename + "_stopsfixed" + pep_ext

else:
    print("Skipping fixing braker's fuckery")
    original_pep_file = pep_file

## Check that the braker cds seqs are in the right reading frame
cds_records = list(SeqIO.parse(cds_file, 'fasta'))
pep_records = list(SeqIO.parse(pep_file, 'fasta'))

funky_reading_frame = []
for count,cds_record in enumerate(cds_records):
    #print("ID %s" % cds_record.id)
    #print("Sequence\n%s" % cds_record.seq)
    if cds_record.seq.translate() != pep_records[count].seq:
        #print('Translation of ' + cds_record.id + ' cds != pep!')
        #print('cds_record.seq.translate(): \n' + cds_record.seq.translate())
        #print('pep_record: \n' + pep_records[count])
        funky_reading_frame.append(True)
    elif cds_record.seq.translate() == pep_records[count].seq:
        print('Translation of ' + cds_record.id + ' cds matches pep.')
        funky_reading_frame.append(False)
    else:
        print('Something went wrong while translating cds!')
        break

print('Number of seqs with funky reading frames: ' + str(np.count_nonzero(funky_reading_frame)))

# Check cds_file for partial codons and seqs below min_length
print('Checking ' + cds_file + ' for partial codons and seqs below min_length')

partial_codons = []
length_below_thresh = []
for cds in SeqIO.parse(cds_file, 'fasta'):
    # check for partial codons
    if len(cds.seq)%3 != 0:
       #print(record.id + ': partial codon detected! Seq length not a multiple of 3')
       partial_codons.append(True)
    elif len(cds.seq)%3 == 0:
       partial_codons.append(False)
    else:
        print('Failed to check for partial codon in '+ record.id +'!')
        break

    # check for seqs below the length threshold
    if len(cds.seq) >= min_length:
        length_below_thresh.append(False)
    elif len(cds.seq) < min_length:
        length_below_thresh.append(True)
    else:
        print('Failed to check seq length for '+ record.id +'!')
        break


# simplify removal process by combining the funky_reading_frame and partial_codons vectors
partial_codons = partial_codons + funky_reading_frame

# Filter cds seqs
print('Filtering ' + cds_file)
filt_cds = []
# remove cds and pep seqs with partial codons and with lengths below thresh
for count,cds in enumerate(SeqIO.parse(cds_file, 'fasta')):
    if partial_codons[count] == True or length_below_thresh[count] == True:
        pass
    elif partial_codons[count] == False and length_below_thresh[count] == False:
        filt_cds.append(cds)
    else:
        print('Failed to filter out ' + record.id + ' in ' + cds_file + '!')

# Save filtered cds file
cds_in_filename, cds_ext = os.path.splitext(cds_file)
cds_out_filename = cds_in_filename + "_filtbrak" + cds_ext
print('Saving filtered seqs in ' + cds_out_filename)
with open(cds_out_filename, "w+") as output_file:
  for record in filt_cds:
      output_file.write(record.format("fasta"))


# Filter pep seqs
print('Filtering ' + pep_file)
filt_peps = []
# remove cds and pep seqs with partial codons and with lengths below thresh
for count,pep in enumerate(SeqIO.parse(pep_file, 'fasta')):
    if partial_codons[count] == True or length_below_thresh[count] == True:
        pass
    elif partial_codons[count] == False and length_below_thresh[count] == False:
        filt_peps.append(pep)
    else:
        print('Failed to filter out '+ record.id + ' in ' + pep_file + '!')

# Save filtered pep file
pep_in_filename, pep_ext = os.path.splitext(original_pep_file)
pep_out_filename = pep_in_filename + "_filtbrak" + pep_ext
print('Saving filtered seqs in ' + pep_out_filename)
with open(pep_out_filename, "w+") as output_file:
  for record in filt_peps:
      output_file.write(record.format("fasta"))

# Done
print('Done')

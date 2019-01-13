"""
@title: GFF3, Unpacks transcriptome files from GFF3 file.
@author: alexander lucaci
@description: Turns GFF3 file to its individual file, in this case I am looking for the transcript files. 

GFF3 File Format: https://useast.ensembl.org/info/website/upload/gff3.html

#Requires: gffutils

<usage>
python GFF3_unpack.py <input file> > <output filename>


https://docs.python.org/3/library/itertools.html#itertools-recipes

testing on spider genome files:
test1: https://www.ncbi.nlm.nih.gov/genome/54232
test2: https://www.ncbi.nlm.nih.gov/genome/12925
control (has transcript file): https://www.ncbi.nlm.nih.gov/genome/13270 
"""

# =============================================================================
# Imports
# =============================================================================
import sys
from datetime import datetime
from optparse import OptionParser
from Bio import Entrez, SeqIO
from multiprocessing import Pool
import itertools
from itertools import zip_longest
# =============================================================================
# Declares
# =============================================================================
startTime = datetime.now()
GB_IDList = []
Entrez.email="alexander.lucaci@temple.edu"
in_file, out_file = "GCA_000611955.2_Stegodyphus_mimosarum_v1_genomic.gff", "GFF3_UNPACK_GCA_000611955.2_Stegodyphus_mimosarum_v1_genomic.gff"

# =============================================================================
# Helper Functions
# =============================================================================  
def singleEntry(ID_List):   #the singleID is the accession number
    #print(ID_List)
    #print(type(','.join(ID_List)))
    
    #handle = Entrez.efetch(db='nucleotide',id=singleID, rettype = 'fasta', retmode= 'text')
    handle = Entrez.efetch(db='nucleotide', id=','.join(ID_List), rettype = 'fasta', retmode= 'text')
    #handle = Entrez.efetch(db='nucleotide', rettype="fasta", id=','.join(accs))
    
    #append a file.
    #f = open('%s.fasta' % singleID, 'a+')
    #f.write(handle.read())
    #handle.close()
    #f.close()
    
    print()
    print("handle", handle.read())
    
    #f = open(out_file, 'a+')
    #f.write(handle.read())
    #handle.close()
    #f.close()

def grouper(iterable, n, fillvalue=""):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
# =============================================================================
# Main loop
# =============================================================================
    
def main_loop(gff_file, output_file):
    global GB_IDList, startTime
    
    with open(gff_file, "r") as f:
        #for n in range(250):
        while True:
            line = f.readline()
            if line == "":
                break
            try:
                transcript = line.split("\t")
                if transcript[2] == "mRNA":
                    print(transcript)
                    GB_IDList += [transcript[0]]
            except:
                pass
            
    print("Number of mRNA transcripts:", len(GB_IDList))
    
    print("Grabbing Genbank sequence")
    
    #for line in grouper(set(GB_IDList), 100):
    #    singleEntry(line)
        
    #Initialize file.
    print("Creating output file")
    output_file = open(output_file, 'w')

# =============================================================================
# Main Program Starts Here
# =============================================================================   
print("# --- Starting, GFF3 file unpacker (transcripts) --- #")
print("# --- Version 0.01                               --- #")
print("# --- Processing:", in_file)
print() 

main_loop(in_file, out_file)
#singleEntry("KK111985.1")
#singleEntry("KK111986.1")

#print(len(set(GB_IDList)))
print("Unpacking complete. Have a nice day!")
print("Runtime:", datetime.now() - startTime) 


# =============================================================================
# End of file. 
# =============================================================================
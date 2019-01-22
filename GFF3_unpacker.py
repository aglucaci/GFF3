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
#in_file= "GCA_000611955.2_Stegodyphus_mimosarum_v1_genomic.gff"
#in_file = sys.argv[1]
in_file = "Stegodyphus_mimosarum.Stegodyphus_mimosarum_v1.42.gff3"
out_file = "GFF3_UNPACK_" + in_file
#"GFF3_UNPACK_GCA_000611955.2_Stegodyphus_mimosarum_v1_genomic.gff"

# =============================================================================
# Helper Functions
# =============================================================================  
def singleEntry(singleID, From, To):   #the singleID is the accession number
    #print(ID_List)
    #print(type(','.join(ID_List)))
    
    handle = Entrez.efetch(db='nucleotide',id=singleID, rettype = 'fasta', retmode= 'text')
    #handle = Entrez.efetch(db='nucleotide', id=','.join(ID), rettype = 'fasta', retmode= 'text')
    #handle = Entrez.efetch(db='nucleotide', rettype="fasta", id=','.join(accs))
    
    #append a file.
    #f = open('%s.fasta' % singleID, 'a+')
    #f.write(handle.read())
    #handle.close()
    #f.close()
    
    #print()
    a = handle.read().split("\n")
    print(a[0])

    try:
        #print(handle.read().split("\n")[1:])
        #print(type(a[1:]), From, To, type(From), type(To))
        print(''.join(a[1:])[int(From)+1:int(To)+1]) #works!
        #print(''.join(a[1:])[0])
    except Exception as e:
        print(e)
    #print("seq", "".join(handle.read().split("\n")[1:])[From+1:To+1])
    #print(handle.read().split("\n")[1])


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
                    #print(transcript)
                    #GB_IDList[transcript[0]] = transcript[3] + "," + transcript[4]
                    GB_IDList.append(transcript[0] + "," + transcript[3] + "," + transcript[4])
                    #print("Grabbing", transcript[0], "From", transcript[3], "to", transcript[4])
                    singleEntry(transcript[0], transcript[3], transcript[4])

            except:
                pass
    #print(GB_IDList[0])
    #print(GB_IDList[GB_IDList.keys()[0]])
    #for k in GB_IDList: print(k, GB_IDList[k], len(GB_IDList))
    print("Number of mRNA transcripts:", len(GB_IDList))
    #Sanity check, len GB_IDlist and biopythons, number of sequences/entries/ whatever they call it. on the output file.
    #print("Grabbing Genbank sequence")
    #for line in grouper(GB_IDList, 100):
        #singleEntry(line)
    #    print(str(line).split(","), "\n")
        
    #Initialize file.
    #print("Creating output file")
    #output_file = open(output_file, 'w')

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

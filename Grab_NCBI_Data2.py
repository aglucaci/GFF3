#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 15:47:18 2018

@author: iamtokr

https://github.com/mttmartin/transcript2gene/blob/master/transcript2gene.py
https://insidedna.io/tutorials/view/biopython-retrieving-fasta-sequences-ncbi-databases
https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists/7961390

"""
#!/usr/bin/env python3
# =============================================================================
# Imports
# =============================================================================
import sys, gffutils, time
from optparse import OptionParser
from Bio import Entrez, SeqIO

# =============================================================================
# Declares
# =============================================================================
GB_IDList = []
Entrez.email="alexander.lucaci@temple.edu"
in_file = "GCA_002102615.1_NepCla1.0_genomic.gff"
out_file = "gff3_to_transcript1.gff"

# =============================================================================
# Helper Functions
# =============================================================================
def singleEntry(singleID):   #the singleID is the accession number
    handle = Entrez.efetch(db='nucleotide',id=singleID, rettype = 'fasta', retmode= 'text')
    
    #append a file.
    #f = open('%s.fasta' % singleID, 'a+')
    #f.write(handle.read())
    #handle.close()
    #f.close()
    
    f = open(out_file, 'a+')
    f.write(handle.read())
    handle.close()
    f.close()
    
# =============================================================================
# Main Program Starts Here
# =============================================================================
def main_loop(gff_file, output_file):
    print("looping")
    global GB_IDList
    
    #Computationally heavy.
    db = gffutils.create_db(gff_file, dbfn="gff3.db", merge_strategy="merge", force=True)

    #initialize file.
    output_file = open(output_file, 'w')
    output_file.write('transcript,gene\n')
    
    #print(len(db))
    print("entering for")
    for i in db.features_of_type('mRNA'):
        #print(type(i))
        a = str(i)
        output_file.write(a + "\n")
        b = a.split("\t")
        
        #Grab ID's for later. We will download fasta's
        GB_IDList += [b[0]]
      
    print("Outputting file, GFF3")    
    #output_file.close()
    #print(GB_IDList)   
    
    print("Grabbing Genbank sequence")
    for c in range(len(set(GB_IDList))):
        singleEntry(GB_IDList[c])
    
def Main():   
    main_loop(in_file, out_file)
    
    print("Done")
    print("Runtime:", datetime.now() - startTime)
    
startTime = datetime.now()
print("Here")
Main()

# =============================================================================
# End of file. 
# =============================================================================

"""
I have all the fasta files from NepCla but they are in separate files. (~3400 files)
-script to combine them all into one big file. 
--This is the transcript file.

-take the gff file output from above. It has the protein product names and sequence sites
-Grab the actual sequence for each thing.

Run this against the protein file.
However, they are not lined up like the common house spider
Have to search for mRNA product name to match protein file.


Tonight, set up the 3rd species to download fasta files. (untitled7_Github.py)


All vs all blast on the transcript files.
https://www.biostars.org/p/169706/

"""


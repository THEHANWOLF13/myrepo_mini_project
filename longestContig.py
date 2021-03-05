import os
from Bio import SeqIO

def longestContig():
    #grabbing our large contigs file
    inFile = open('LargeContigs.txt','r')

    data = SeqIO.parse('LargeContigs.txt','fasta')

    
    longest = None
    longest_rec = None
    #determining the longest contig via updating 'longest' as the file is parsed
    for record in data:
        x = len(record.seq)

        if (longest ==None):
            longest = x
            longest_rec = record
        if x > longest:
            longest = x
            longest_rec = record
    inFile.close()
    #writing into a new file which will be used in the BLAST
    SeqIO.write(longest_rec, 'LongestContig.fasta','fasta')

longestContig()
        
            

    

    
            
        

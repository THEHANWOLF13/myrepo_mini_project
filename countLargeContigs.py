import os
from Bio import SeqIO

log = open('miniproject.log' , 'a')

def countLargeContigs():
    #this function sums the total base pairs in the assembly

    contFile = open('LargeContigs.txt','r')

    #parsing through the file with our large contigs in fasta format

    data = SeqIO.parse('LargeContigs.txt','fasta')

    contigLength = []
    
    for record in data:
        x = len(record.seq)
        contigLength.append(int(x))


    total = sum(contigLength)
    #adding our total to our log file
    log.write('\n'+'There are '+ str(total)+' bp in the assembly.')
    contFile.close()
countLargeContigs()



    

    

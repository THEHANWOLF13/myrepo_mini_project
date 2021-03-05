import os
from Bio import SeqIO

log = open('miniproject.log' , 'a')
def numLargeContigs():

    #we will be counting contings that are larger than 1000 base pairs


    contigFile = open('LargeContigs.txt','w')
    count = 0

    #initializing count and updating as each large contig is found
    #must supply this path because Spades stores the contigs in a subfolder
    data = SeqIO.parse('./SpadesAssembly/contigs.fasta','fasta')

    for record in data:
        x = len(record.seq)

        if x > 1000:
            count = count+1

            contigFile.write('> '+str(record.id) +'\n' +str(record.seq) + '\n')

    contigFile.close()
    log.write('There are '+str(count)+' contigs > 1000 bp in the assembly.')

numLargeContigs()

    

import os

seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
log = open('miniproject.log','a')
def numReads(seq):
    header =''
    # keeping track of which donor we will write into the log file
    if seq == 'SRR5660030':
        header == 'Donor 1(2dpi)'
    elif seq = 'SRR5660033':
        header == 'Donor 1(6dpi)'
    elif seq = 'SRR5660044':
        header == 'Donor 3 (2dpi)'
    else:
        header == 'Donor 3 (6dpi)'

    # grabbing our originial fastq's
    srrFile1 = open(seq +'.1_1.fastq')
    srrFile2 = open(seq +'.1_2.fastq')

    count = 0
    counter = 0
    # counting reads before bowtie
    for line in srrFile1:
        count = count + 1
    for line in srrFile2:
        counter = counter +1

    beforeCount = (count+counter)/8
    #opening the files that the bowtie2 created
    afterFile1 = open(seq+'_mapped_1.fastq')
    afterFile2 = open(seq+'_mapped_2.fastq')

    sec_count = 0
    third_count = 0
    #counting reads after bowtie
    for line in afterFile1:
        sec_count = sec_count + 1

    for line in afterFile2:
        third_count = third_count + 1

    afterCount = (sec_count + third_count)/8
    #updating the log file with the results
    log.write(str(header)+' had ' +str(beforeCount) + ' read pairs before Bowtie2 filtering and '+ str(afterCount)+ ' pairs after.')

for i in seqs:
    numReads(i)
        
        

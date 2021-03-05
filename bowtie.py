import os
seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

def bowtie2(seq):
    #builing Bowtie index for the HCMV
    #followed the provided commands in class

    build_bowtie = 'bowtie2-build ./EF999921.fasta EF999921_index'
    os.system(build_bowtie)

    
    bowtie_index = 'bowtie2 -x EF999921_index -1 '+seq+ '.1_1.fastq -2 '+ seq +'.1_2.fastq -S tmp.sam --al-conc' +seq+ '_mapped_%.fastq'
    os.system(bowtie_index)



'''
def SAM_FASTQ_conversion(seq):
    #converting SAM files to fastq format
    conversion = 'bash samtofasq.sh ' + seq
    os.system(conversion)

'''

for i in seqs:
    bowtie2(i)



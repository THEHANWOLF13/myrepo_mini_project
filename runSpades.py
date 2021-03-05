import os

seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
log = open('miniproject.log' , 'a')

def runSpades(sequences):
    #running spades on the generate paired-end fastqs for each SRR

    seq1 = sequences[0]
    seq2 = sequences[1]
    seq3 = sequences[2]
    seq4 = sequences[3]

    spades_cmd ='spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 '+seq1+'_mapped_1.fastq --pe1-2 '+seq1+'_mapped_2.fastq --pe2-1 '+seq2+'_mapped_1.fastq --pe2-2 '+seq2+'_mapped_2.fastq --pe3-1 '+seq3+'_mapped_1.fastq --pe3-2 '+seq3+'_mapped_2.fastq --pe4-1 '+seq4+'_mapped_1.fastq --pe4-2 '+seq4+'_mapped_2.fastq -o SpadesAssembly/'
    os.system(spades_cmd)
    log.write(str(spades_cmd))


runSpades(seqs)

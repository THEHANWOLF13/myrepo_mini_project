import os

seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
def kallisto(srr):

    # following the kallisto commands in the powerpoint, we have the following
    kallisto_index = "time kallisto index -i index.idx EF999921_CDS.fasta"
    os.system(kallisto_index)

    kallisto_quant = "time kallisto quant -i index.idx -o ./" + srr + " -b 30 -t 4 " + srr + ".1_1.fastq " + srr + ".1_2.fastq"
    os.system(kallisto_quant)

for i in seqs:
    kallisto(i)

import os
import logging
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq

current_dir = os.getcwd()
os.chdir(current_dir)

log = open('miniproject.log' , 'a')

seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045'] # SRR numbers of desired sequences

def getSeq(seq): # grabbing sequences using the wget command
    grabSeq = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/' + str(seq) + '/' + str(seq) + '.1'

    rename = 'mv' + str(seq) + '.1' + str(seq)
    # converting them into fastq paired end files
    uncomp = 'fastq-dump -I --split-files' + str(seq)

    os.system(grabSeq)
    os.system(rename)
    os.system(uncomp)

'''
for i in seqs:
    getSeq(seqs)
'''
def makeTestData(files):
    c = os.getcwd()

    for file in files:
        inFile = open(c+'/'+file,'r')
        outFile = open(c+'/test_'+file,'w')

        inFileStrip = inFile.read().strip().split('\n')

        for i in range(40000):
            outFile.write(inFileStrip[i]+'\n')

testData = ['SRR5660030.1_2.fastq','SRR5660030.1_1.fastq', 'SRR5660033.1_1.fastq', 'SRR5660033.1_2.fastq', 'SRR5660044.1_1.fastq', 'SRR5660044.1_2.fastq','SRR5660045.1_1.fastq','SRR5660045.1_2.fastq']

makeTestData(testData)
def transcriptomeIndex():
    Entrez.email = 'rrajagopal@luc.edu'

    outFasta = open("EF999921.fasta",'w')
    outCDS = open("EF999921_CDS.fasta",'w')

    
    #grabbing the records
    handle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype = 'fasta')
    recs = list(SeqIO.parse(handle, 'fasta'))

    #writing in the desired FASTA format
    outFasta.write('>' + str(recs[0].description) + '\n' + str(recs[0].seq))
    outFasta.close()

    #retrieving GenBank records
    GBhandle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype = 'gb', retmode = 'text')

    count = 0
    # we are looking specifically for CDS features
    # count is updated as more CDS are found
    for record in SeqIO.parse(GBhandle, 'genbank'):
        for item in record.features:
            if item.type =='CDS':
                count = count + 1 
                # writing in the desired format
                outCDS.write('>' + str(item.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'", "") + '\n' + str(item.location.extract(record).seq) + '\n')
    outCDS.close()
    print(count)

result = transcriptomeIndex()

log.write('THE HCMV genome (EF999921) has ' + str(result)+ ' CDS. ' + '\n')
#print('THE HCMV genome (EF999921) has ' + str(result)+ ' CDS. ')


def kallisto(srr):

    #performing kallisto indexing and running using the given commands    
    kallisto_index = 'time kallisto index -i index.idx EF999921_CDS.fasta'
    os.system(kallisto_index)

    kallisto_quant = 'time kallisto quant -i index.idx -o ./' +'test_'+ srr + ' -b 30 -t 4 '+ 'test_' + srr + '.1_1.fastq ' + 'test_' + srr + '.1_2.fastq'
    os.system(kallisto_quant)

for i in seqs:
    kallisto(i)

#seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
    
def sleuth_input(srrs):
    infile = open('sleuth_input.txt', 'w')

    path = os.getcwd() # this path is not hardcoded
    

    first_cond = '2dpi'
    sec_cond = '6dpi'
    # first line of the file
    infile.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')

    for item in srrs: # assigning condition based on SRR number and writing to file accordingly
        if int(item[3:])%2==0:
            infile.write(item + '\t' + first_cond + '\t' + str(path)+'/'+ item + '\n')

        else:
            infile.write(item + '\t' + sec_cond + '\t' + str(path) + '/' + item + '\n')


    infile.close()


sleuth_input(seqs)

def run_sleuth():

    #will run the sleuth from R and will add to the log file

    sleuthRun = 'Rscript sleuth.r'

    os.system(sleuthRun)
    output = open('topsleuth.txt', 'r')

    line_list = output.readlines()

    for line in line_list:
        
        log.write(line)

run_sleuth()

def bowtie2(seq):
    #builing Bowtie index for the HCMV
    #followed the provided commands in class

    build_bowtie = 'bowtie2-build ./EF999921.fasta EF999921_index'
    os.system(build_bowtie)

    
    bowtie_index = 'bowtie2 -x EF999921_index -1 '+'test_' +seq+ '.1_1.fastq -2 '+ 'test_'+ seq +'.1_2.fastq -S tmp.sam --al-conc'+ ' test_'+seq+ '_mapped_%.fastq'
    os.system(bowtie_index)


for i in seqs:
    bowtie2(i)
    
def numReads(seq):
    header =''
    # keeping track of which donor we will write into the log file
    if seq == 'SRR5660030':
        header = 'Donor 1 (2dpi)'
    elif seq == 'SRR5660033':
        header = 'Donor 1 (6dpi)'
    elif seq == 'SRR5660044':
        header = 'Donor 3 (2dpi)'
    else:
        header = 'Donor 3 (6dpi)'

    # grabbing our originial fastq's
    srrFile1 = open('test_'+seq +'.1_1.fastq')
    srrFile2 = open('test_'+seq +'.1_2.fastq')

    count = 0
    counter = 0
    # counting reads before bowtie
    for line in srrFile1:
        count = count + 1
    for line in srrFile2:
        counter = counter +1

    beforeCount = (count+counter)/8
    #opening the files that the bowtie2 created
    afterFile1 = open('test_'+seq+'_mapped_1.fastq')
    afterFile2 = open('test_'+seq+'_mapped_2.fastq')

    sec_count = 0
    third_count = 0
    #counting reads after bowtie
    for line in afterFile1:
        sec_count = sec_count + 1

    for line in afterFile2:
        third_count = third_count + 1

    afterCount = (sec_count + third_count)/8
    #updating the log file with the results
    log.write(str(header)+' had ' +str(beforeCount) + ' read pairs before Bowtie2 filtering and '+ str(afterCount)+ ' pairs after.'+'\n')

for i in seqs:
    numReads(i)

def runSpades(sequences):
    #running spades on the generate paired-end fastqs for each SRR

    seq1 = sequences[0]
    seq2 = sequences[1]
    seq3 = sequences[2]
    seq4 = sequences[3]

    spades_cmd ='spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 '+'test_'+seq1+'_mapped_1.fastq --pe1-2 '+'test_'+seq1+'_mapped_2.fastq --pe2-1 '+'test_'+seq2+'_mapped_1.fastq --pe2-2 '+'test_'+seq2+'_mapped_2.fastq --pe3-1 '+'test_'+seq3+'_mapped_1.fastq --pe3-2 '+'test_'+seq3+'_mapped_2.fastq --pe4-1 '+'test_'+seq4+'_mapped_1.fastq --pe4-2 '+'test_'+seq4+'_mapped_2.fastq -o SpadesAssembly/'
    os.system(spades_cmd)
    log.write(str(spades_cmd))


runSpades(seqs)

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

def makeBlastDB():
    #following the commands on the powerpoints for blastn and grabbing the top 10 hits
    construct_db = 'makeblastdb -in sequence.fasta -out family_db -title Betaherpesvirinae -dbtype nucl'
    os.system(construct_db)

    blast_cmd = "blastn -query LongestContig.fasta -db family_db -out blast_results.csv -outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    os.system(blast_cmd)


makeBlastDB()

def logBlast():
    #tab delimited format
    log.write('\n'+'sacc'+'\t'+'pident'+'\t'+'length'+'\t'+'qstart'+'\t'+'qend'+'\t'+'sstart'+'\t'+'send'+'\t'+'bitscore'+'\t'+'evalue'+'\t'+'stitle' +'\n')
    numHit = 10
    handle = open('blast_results.csv','r')

    data = handle.readlines()[0:10]
    # we want the top 10 hits to be displayed in tab delimited format
    for line in data:
        log.write(line.replace(',','\t'))  
logBlast()
log.close()

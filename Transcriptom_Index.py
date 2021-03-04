import os
from Bio import SeqIO
from Bio import Entrez

log = open('miniproject.log' , 'a')

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
log.close()


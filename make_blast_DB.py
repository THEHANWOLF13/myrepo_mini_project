import os

def makeBlastDB():
    #following the commands on the powerpoints for blastn and grabbing the top 10 hits
    construct_db = 'makeblastdb -in sequence.fasta -out family_db -title Betaherpesvirinae -dbtype nucl'
    os.system(construct_db)

    blast_cmd = "blastn -query LongestContig.fasta -db family_db -out blast_results.csv -outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
    os.system(blast_cmd)


makeBlastDB()

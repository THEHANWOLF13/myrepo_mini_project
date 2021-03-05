import os

log = open('miniproject.log','a')


def logBlast():
    #tab delimited format
    log.write('\n'+'sacc'+'\t'+'pident'+'\t'+'length'+'\t'+'qstart'+'\t'+'qend'+'\t'+'sstart'+'\t'+'send'+'\t'+'bitscore'+'\t'+'evalue'+'\t'+'stitle' +'\n')
    numHit = 10
    handle = open('blast_results.csv','r')

    data = handle.readlines()[0:10]
    # we want the top 10 hits to be displayed in tab delimited format
    for line in data:
        log.write(line.replace(',','\t'))

    log.close()
logBlast()

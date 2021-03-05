import os

seqs = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
def sleuth_input(srrs):
    infile = open('sleuth_input.txt', 'w')

    path = os.getcwd()
    

    first_cond = '2dpi'
    sec_cond = '6dpi'
    # first line of the file
    infile.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')

    for item in srrs: # assigning condition based on SRR number and writing to file accordingly, writing in tab delimited format
        if int(item[3:])%2==0:
            infile.write(item + '\t' + first_cond + '\t' + str(path)+'/'+ item + '\n')

        else:
            infile.write(item + '\t' + sec_cond + '\t' + str(path) + '/' + item + '\n')


    infile.close()


sleuth_input(seqs)

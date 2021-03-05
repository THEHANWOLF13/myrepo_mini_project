import os


log = open('miniproject.log' ,'a')
def run_sleuth():

    #will run the sleuth from R and will add to the log file

    sleuthRun = 'Rscript sleuth.r'

    os.system(sleuthRun)
    output = open('topsleuth.txt', 'r')

    line_list = output.readlines()

    for line in line_list:
        
        log.write(line)

run_sleuth()

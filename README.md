# myrepo_mini_project COMP 383 Loyola University Chicago

We are interested in analyzing transcriptomes of Human Herpes Virus 5

Softwares Required:
  Linux/Unix,
  Python3,
  Biopython,
  Kallisto,
  SPAdes,
  Bowtie2
  
 Run miniproject_main.py for the desired output:
  this pipeline contains test data for the first 40,000 lines of transcriptome data for HCMV (first 10000 reads)
  
  
To run the script, first clone the repo via the following command
  'git clone https://github.com/THEHANWOLF13/myrepo_mini_project.git'
  
This will create a local repository.
  
Within the local repo, run the following command:

'python3 miniproject_main.py'
  
 Notable outputs from miniproject_main.py:
 
 miniproject.log -> contains results of the full pipeline
 EF999921.fasta -> contains HCMV genome
 EF999921_CDS.fasta -> contains coding sequences for the HCMV genome
 LargeContig.txt -> shows all contigs exceeding 1000 bp
 
 
 Notable files included in this repo:
 
 Transcriptom_index.py -> used to generate input for Kaliisto,
 Sleuth_input.py -> converts Kallisto output into suitable input for Sleuth,
 sleuth.r -> performs sleuth,
 runSleuth.py -> runs sleuth in python
 
 
 
 
 
 
 
  
  
  
  




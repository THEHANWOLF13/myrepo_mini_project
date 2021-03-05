#library(devtools)
library(sleuth)
library(data.table)

#source("http://bioconductor.org/biocLite.R")
#biocLite("devtools")    # only if devtools not yet installed
#biocLite("pachterlab/sleuth")


stab <- read.table("sleuth_input.txt", header = TRUE, stringsAsFactors = FALSE, sep = '\t')

#sleuth object
so <- sleuth_prep(stab)

#fitting model to compare conditions

so <- sleuth_fit(so, ~condition, 'full')

#reduced model
so <- sleuth_fit(so, ~1, 'reduced')


#likelihood ratio test for differential expression between conditions
so<- sleuth_lrt(so, 'reduced','full')

#extracting results
library(dplyr)
sleuth_table <- sleuth_results(so, 'reduced:full','lrt', show_all=FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <=0.05) %>% dplyr::arrange(pval)

signif_sleuth <- sleuth_significant %>% select(target_id, test_stat, pval, qval)

#writing target id, test stat, pval, and qval for each significant transcript

# this will include the header and will be in tab-delimited format
write.table(signif_sleuth, file = "topsleuth.txt",quote = FALSE,row.names = FALSE)





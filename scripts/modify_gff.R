library(dplyr)

#setwd("~/Documents/mason/browser/pnag/static/data/2022_08_18_05_51_45/Ga0066243_109818/reference_sequences") #

options(echo = F)
args <- commandArgs(trailingOnly = T)

#args <- c("./full_transcripts_2606217223.gff", "./genelengths.tsv")
print(args)

# import gff and fasta lengths:
gff <-  read.delim(args[1], header = FALSE) 

fastas <- read.delim(args[2], header = FALSE, row.names = 1)


# fet length of fasta for each gff entry:
gff$RL <- fastas[gff$V1,]

errorgenes <- gff[(gff$RL < gff$V5 | gff$V4 < 1),]$V3

if (length(errorgenes > 0)) {
  texterror <- paste0("Following genes were not screened for off-targets, as their 5'UTR is too short (<30 nt): \n", 
                      paste0(errorgenes, collapse = "; "))
  write(texterror, file=args[3], append = T) # write(texterror, file="../../warnings.txt", append = T)
}

#change gff if RL is bigger than end:
#gff$V5 <- ifelse(gff$RL < gff$V5, gff$RL, gff$V5)
#gff$V4 <- ifelse(gff$V4 < 1, 1, gff$V4)

gff <- gff[!(duplicated(gff$V4) & duplicated(gff$V5) & duplicated(gff$V3)),]
# write gff that is modified:
write.table(gff[!(gff$RL < gff$V5 | gff$V4 < 1),-10],
          file = args[1],  
          sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



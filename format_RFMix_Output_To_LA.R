library(methods)
library(data.table)
library(dplyr)
library(ggplot2)
library(DelayedArray)


#for(chr in 1:22)# loop through chromosomes
#{

argv <- commandArgs(trailingOnly = TRUE)
chr <- argv[1]
print(paste0('Chromosome ',chr))
SNP_Identifiers<-fread(paste0(### path and file name to SNP annotations ### 
                              ),stringsAsFactors=F, header=F,data.table = F)$V3
count_2_allelePairs<-data.frame()
  for(un in 1:9) # loop through the chuncks
  {
    print(paste0('Chunck ',un))
    # Read RFMIx viterbi output as well as SNP identifiers from RFMIx input map files
    RFMix_Viterbi<-fread(paste0(### path and file name for the Viterbi outcome of RFMix ###
                                ),stringsAsFactors=F, header=F,data.table = F)   
    
    alleleACols<-seq(1,ncol(RFMix_Viterbi),by=2) # sequence of odd numbered columns
    #loop through pairs of columns and count the number of occurrence of "2" in each pair of columns
    for(y in alleleACols)
    {
      if(un==1 && y==1) # initialize allele count data frame at first y of first chunck
      {  
        count_2_allelePairs<- as.data.frame(rowCounts(data.matrix(RFMix_Viterbi)[,(y:(y+1))],value=2))
      } else # concatenate the allele count data frame with the count of every two column in the Viterbi matrix for all chuncks
      {  
        count_2_allelePairs<- cbind(count_2_allelePairs, as.data.frame(rowCounts(data.matrix(RFMix_Viterbi)[,(y:(y+1))],value=2)))
      }
    }
    
  } 
colnames(count_2_allelePairs)<-1:length(count_2_allelePairs)
RFMixToLA_output<-cbind(SNP_Identifiers, count_2_allelePairs) #add SNP identifiers
colnames(RFMixToLA_output)<-c('rsid',1:length(count_2_allelePairs))
filename_RFtoLA<-paste0(###Path and File name to save the LA information####)
fwrite(RFMixToLA_output, file = filename_RFtoLA, append = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
#}

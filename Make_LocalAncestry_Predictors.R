library(methods)
library(data.table)
library(dplyr)
source("/projects/b1047/LNahlawi/CNV_AA_RFMix/check_flipped_alleles.R")



##reading command line arguments
#uncomment when running code on command line
argv <- commandArgs(trailingOnly = TRUE)
LA_file<- argv[1] #full path and file name for local ancestry
Geno_file<-argv[2] # full path and file name for genotypes
check_flip_alleles_Flag<-argv[3] # 0 or 1
chr<-argv[4]
## uncomment when testing the code in Rstudio
# LA_file<-### path and file name for LA data ### 
# Geno_file<-### path and file name for genotype dosage ###
# check_flip_alleles_Flag<- # 1 or 0

Geno_chrom_flip_SNPs<-data.frame()

local_ancestry <-fread(LA_file, stringsAsFactors = F, header=T, data.table = F) ##local ancestry file
print("Done reading RFmix output matrix")


Geno_chrom <- fread(Geno_file, stringsAsFactors=F, data.table = F)
print("Done reading the genotype of the chromosome")
colnames(Geno_chrom)[2]<-'rsid'
## check if there are alleles flipped between the genotype data and model predictors
if(check_flip_alleles_Flag==1){
  
  #uncomment when running in command line
  Model_file<-argv[5] # full path and filename for the model summary
  Model_weight_file<-argv[6] # full path and filename for the predictors weights
  predictors_annotation<-argv[7] # full path and filename for the predictors annotations
  
  #uncomment to test code in RStudio
  # Model_file<- ### path  and file name for the LA-GEM model summary ###
  # Model_weight_file<- ### path and filename for the weights of Predictors in LA-GEM model ###
  # predictors_annotation<- ### path and filename for predictors (SNP/LA terms/LA-dosage interaction terms) annotations ###
  
  #call the function that chech allele flips
  Geno_chrom_flip_SNPs<-check_flipped_Alleles_func(Geno_file,Model_file,Model_weight_file, predictors_annotation)
  #fix the dosage for the flipped alleles 
  Geno_chrom[which(Geno_chrom$rsid %in% Geno_chrom_flip_SNPs$rsid_dbSNP150),7:length(Geno_chrom)] <- 2- Geno_chrom[which(Geno_chrom$rsid %in% Geno_chrom_flip_SNPs$rsid_dbSNP150),7:length(Geno_chrom)]
  #flip the ref and alt alleles according to the data frame read above to match the heaptocytes AA
  Geno_chrom[which(Geno_chrom$rsid %in% Geno_chrom_flip_SNPs$rsid_dbSNP150),]<-Geno_chrom[which(Geno_chrom$rsid %in% Geno_chrom_flip_SNPs$rsid_dbSNP150),c('V1','rsid','V3','V5','V4', x<-paste('V',6:length(Geno_chrom),sep = ''))]
  # save the Geno files after flipping alleles
  fwrite(Geno_chrom, paste0(gsub("\\..*","",Geno_file),"wfixed_flipAlleles.txt"), append = F, sep = " ", row.names = F, col.names = F)
  Prefix_columns_geno<-Geno_chrom[,2:6]
}

Geno_chrom<-Geno_chrom[,c('rsid',x<-paste('V',7:length(Geno_chrom),sep = ''))]
## get all local ancestry sequences for current chromosome 
non_unique_LA_chrom<-subset(local_ancestry, rsid %in% Geno_chrom$rsid)
print("Done getting local ancestry for the chromosome")
Prefix_columns_NoUniq_LA<-subset(Prefix_columns_geno, rsid %in% non_unique_LA_chrom$rsid)


##get the UNIQUE local ancestry sequences for current chromosome to make the additional candidate predictors from local ancestry information
unique_LA_chrom<-non_unique_LA_chrom[!duplicated(non_unique_LA_chrom[,2:length(non_unique_LA_chrom)]), ]
Prefix_columns_UNIQ_LA<-subset(Prefix_columns_NoUniq_LA, rsid %in% unique_LA_chrom$rsid)

Prefix_columns_NoUniq_LA[,1]<-paste0(Prefix_columns_NoUniq_LA[,1],'_INT')
Prefix_columns_UNIQ_LA[,1]<-paste0(Prefix_columns_UNIQ_LA[,1],'_LA')
print("Done getting unique local ancestry for the chromosome")


##use non-unique local ancestry sequences to calculcate the interaction term between genotype and local ancestry
Inter_Geno_LA<-subset(Geno_chrom, rsid %in% non_unique_LA_chrom$rsid)[,2:length(non_unique_LA_chrom)]*non_unique_LA_chrom[,2:length(non_unique_LA_chrom)]
Inter_Geno_LA<-cbind(non_unique_LA_chrom$rsid, Inter_Geno_LA)
print("Done multiplying genotype dosage with local ancestry sequences, Interaction term")


## Make the new dataFrame of candidate predictors by concatenating the genotype, unique local ancestry sequences and
#       the interaction term (multiplication of all ancestry sequences for the specific chromosome  with genotype of the chromosome)
Geno_LA_Inter <- cbind(transpose(Geno_chrom[-1,]), transpose(unique_LA_chrom[-1,]), transpose(Inter_Geno_LA[-1,])) ## for hepatocyte data
LA_Inter <- cbind(transpose(unique_LA_chrom[-1,]), transpose(Inter_Geno_LA[-1,]))

## Assign the correct colnames to the new dataFrame of candidate predictors
colnames(Geno_LA_Inter)<-Geno_LA_Inter[1, ]
Geno_LA_Inter<-Geno_LA_Inter[-1,]

## Add _LA suffix to the rsid colnames to differentiate between genotype columns and local ancestry columns
LA_Boundary_1<-dim(Geno_chrom)[1]
LA_Boundary_2<-(dim(Geno_chrom)[1]-1+dim(unique_LA_chrom)[1]-1)
colnames(Geno_LA_Inter)[LA_Boundary_1:LA_Boundary_2]<-paste(colnames(Geno_LA_Inter)[LA_Boundary_1:LA_Boundary_2],"_LA",sep="")

## Add _INT suffix to the rsid colnames to identify them as columns of the interaction term (the multiplication between genotype and local ancestry terms)
INTER_Boundary_1<-(dim(Geno_chrom)[1]+dim(unique_LA_chrom)[1]-1)
INTER_Boundary_2<-dim(Geno_chrom)[1]-1+dim(unique_LA_chrom)[1]-1+dim(Inter_Geno_LA)[1]-1
colnames(Geno_LA_Inter)[INTER_Boundary_1:INTER_Boundary_2]<-paste(colnames(Geno_LA_Inter)[INTER_Boundary_1:INTER_Boundary_2],"_INT",sep="")

## Assign the sample names as the rownames of the dataFrame of candidate predictors
rownames(Geno_LA_Inter)<-paste(colnames(Geno_chrom)[-1])
print("Done making the predictors data frame")

# save the dataFrame of candidate predictors
file_Geno_LA_Inter <- paste0(gsub("\\..*","",Geno_file),"Geno_LA_Inter_Predictors.txt")
file_LA_Inter <- paste0(gsub("\\..*","",Geno_file),"LA_Inter_Predictors.txt")
fwrite(LA_Inter, file = file_LA_Inter, append = FALSE, quote = FALSE, sep = " ", eol = "\n", na = "NA", dec = ".", col.names = TRUE)
fwrite(Geno_LA_Inter, file = file_Geno_LA_Inter, append = FALSE, quote = FALSE, sep = " ", eol = "\n", na = "NA", dec = ".", col.names = TRUE)


## Transpose the predictors data-frame (genotype, local ancestry and interaction terms) testing/using LA-GEM
# -> samples are columns and predictors are rows
geLaInter_file_for_testing<-paste0(gsub("\\..*","",Geno_file),"Geno_LA_Inter_Predictors_TESTING.txt")
colnamesBeforeTrans<-colnames(Geno_LA_Inter)
geLaInter_chr<-cbind(rbind(Prefix_columns_geno[-1,],Prefix_columns_UNIQ_LA[-1,], Prefix_columns_NoUniq_LA[-1,]), transpose(Geno_LA_Inter))
# uncomment if dataframe has rownames
#colnames(geLaInter_chr)<-geLaInter_chr[1,] 
#geLaInter_chr<- geLaInter_chr[-1,]
fwrite(geLaInter_chr, file = geLaInter_file_for_testing, append = FALSE, quote = FALSE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = TRUE)




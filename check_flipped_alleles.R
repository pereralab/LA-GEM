## Check how many LA and INT predictors of LA-GEM model have flipped ref/alt alleles in the tested dataset
# Arguments: tested_genotype->filename & path of genotype to be tested
#            model_file-> filename & path of model summary
#            model_weight_file->filename & path of predictors weights 
#            model_Predictors_annoantions->filename & path of predictors annotations

check_flipped_Alleles_func<-function(tested_genotype, model_file, model_weight_file,model_Predictors_annotations){
  
  p_value_thresh<-0.05
  rho_thresh<-0.1

  LA_INT_terms<-data.frame()
  LA_INT_terms_filtered<-data.frame()
  mergeLAIN_Flip<-data.frame()
  mergeLAIN_Flip_filtered<-data.frame()

  # Read model weights and model summary
  weights_chr<-fread(model_weight_file,stringsAsFactors=F, header=T,sep=",", data.table=F)
  model_summary_chr<-subset(fread(model_file, stringsAsFactors = F, header = T, data.table = F), gene_id %in% weights_chr$gene_id)
  
  # Fetch predictors of well-predicted genes
  model_summary_chr_filtered<-model_summary_chr[which(model_summary_chr$rho_avg>=rho_thresh & model_summary_chr$zscore_pval<=p_value_thresh),]
  weights_chr<-subset(weights_chr, gene_id %in% model_summary_chr$gene_id)
  weights_chr_filtered<-subset(weights_chr, gene_id %in% model_summary_chr_filtered$gene_id)
  
  # Read predictors annotations and tested genotype annotations
  annotations_chr<-fread(model_Predictors_annotations, stringsAsFactors=F,header=T, data.table=F)
  annotations_chr_tested<-fread(tested_genotype,stringsAsFactors = F, header = F, data.table=F)[,1:6]
  
  # Fetch alleles for LA and INT terms from annotations and filtered predictors weights
  # using unfiltered weights
  LA_INT_terms<- rbind(LA_INT_terms,annotations_chr[which(grepl(paste0(annotations_chr$dbSNP150,'_LA'), weights_chr$rsid, ignore.case = T)==T),4:6])
  LA_INT_terms<- rbind(LA_INT_terms,annotations_chr[which(grepl(paste0(annotations_chr$dbSNP150,'_INT'), weights_chr$rsid, ignore.case = T)==T),4:6])
  #using filtered weights
  LA_INT_terms_filtered<- rbind(LA_INT_terms_filtered,annotations_chr[which(grepl(paste0(annotations_chr$dbSNP150,'_LA'), weights_chr_filtered$rsid, ignore.case = T)==T),4:6])
  LA_INT_terms_filtered<- rbind(LA_INT_terms_filtered,annotations_chr[which(grepl(paste0(annotations_chr$dbSNP150,'_INT'), weights_chr_filtered$rsid, ignore.case = T)==T),4:6])
  print(paste0('number of LA & INT terms = ',nrow(LA_INT_terms)))
  
  #Fetch the annotations of the LA and INT terms in the tested genotypes
  LA_INT_tested<-subset(annotations_chr_tested, V2 %in% LA_INT_terms$rsid_dbSNP150)[,c("V2","V4","V5")]
  LA_INT_tested_filtered<-subset(annotations_chr_tested, V2 %in% LA_INT_terms_filtered$rsid_dbSNP150)[,c("V2","V4","V5")]
  print(paste0('number of overlap with LA & INT terms = ',nrow(LA_INT_tested)))
  
  #Merge data about flipped alleles to save in a file
  mergeLAIN_Flip<- rbind(mergeLAIN_Flip, merge(LA_INT_terms, LA_INT_tested, by.x=c('rsid_dbSNP150', 'ref_vcf'), by.y = c('V2', 'V5')))
  print(paste0('size of the merge matrix = ',nrow(merge(LA_INT_terms, LA_INT_tested, by.x=c('rsid_dbSNP150', 'ref_vcf'), by.y = c('V2', 'V5')))))
  mergeLAIN_Flip_filtered<- rbind(mergeLAIN_Flip_filtered, merge(LA_INT_terms_filtered, LA_INT_tested_filtered, by.x=c('rsid_dbSNP150', 'ref_vcf'), by.y = c('V2', 'V5')))
  
  return(mergeLAIN_Flip_filtered)

}



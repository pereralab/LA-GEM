# LA-GEM
Local Ancestry based Gene Expression Models (LA-GEM)

This is an extension of PredixCan gene expression prediction framework o incorporate loci-specific inferred Local Ancestry (LA) into the prediction
model. We train a linear model per gene to map genotype to GE levels.
For each model, we generate 3 sets of predictors: dosage data for cis-SNPs, LA data for the respective loci, and
interaction terms consisting of the product of dosage and LA data for each locus.

To Run Train LA-GEM on another data from another tissue, you need to use the PrediXcan model training code (https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7) and replace their "gtex_v7_nested_cv_elnet.R" file
with the modified version available in this LA-GEM repository. You'll also need to run RFMix (https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip) on your genotype data to infer Local Ancestry (LA) for all of your variants.

Steps To run:

1- Make local ancestry Predictors
  - Run "Make_LocalAncestry_Predictors.R" and make sure to set the flag to run the function that check flipped alleles "check_flipped_alleles.R".

2- Format LA into LA-GEM acceptable format
  - Run "format_RFMix_Output_To_LA.R" to change the format of the RFMix output. RFMix infers LA for both alleles of each variants. This formatting function counts the number of alleles that has African Ancestry. 

3- Train the LA-GEM Model
  - Run PrediXcan training pipeline but using "gtex_v7_nested_cv_elnet.R" provided in La-GEM Repository.


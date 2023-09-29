# LA-GEM
LA-GEM: Local Ancestry Guided Expression Modeling

Description
LA-GEM is a specialized computational model developed for the purpose of gene expression prediction. It innovatively integrates SNP-based Local Ancestry (LA) information with cis-regional genetic variants. This approach creates a powerful and comprehensive prediction tool suitable for applications in genomics and bioinformatics.

Key Feature
SNP-Based Local Ancestry Integration: Employs single-nucleotide polymorphism (SNP)-based local ancestry data as a distinct set of predictors for gene expression prediction.

Prerequisites
R
RFMix v 1.5.4

Getting Started
To utilize LA-GEM effectively, you must have calculated SNP-based local ancestry using RFMix v 1.5.4. Below are the steps to run LA-GEM:

Steps To Run LA-GEM

1- Make Local Ancestry Predictors
 - Execute Make_LocalAncestry_Predictors.R and set the flag to run the function that checks flipped alleles, check_flipped_alleles.R.

2- Format LA into LA-GEM Acceptable Format
 - Run format_RFMix_Output_To_LA.R to modify the output format from RFMix. RFMix infers LA for both alleles of each variant. This formatting function tallies the number of alleles that have African Ancestry.

3- Train the LA-GEM Model
 - Execute the PrediXcan training pipeline but use the modified gtex_v7_nested_cv_elnet.R file provided in the LA-GEM repository.

## Additional Resources

- [PrediXcan model training code](https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7)
- [RFMix v 1.5.4](https://www.dropbox.com/s/cmq4saduh9gozi9/RFMix_v1.5.4.zip)

###Bidirectional two-sample MR analyses: Common infections and platelet function###

#Install relevant packages 
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("ggplot2")

#Load packages
library(TwoSampleMR)

##Two-sample MR: Common infections on platelet function##

#Read in infection exposure data
#EBV EA-D
data <- read.table("/projects/MRC-IEU/research/projects/ieu1/wp1/028/working/data/results/UKBB/Kachuri_infection_GWAS/ukb_seroreact_EBV_EAD_5E08.txt", header = T)
exposure <- format_data(data, type ="exposure", header = TRUE, snp_col = "ID", chr_col = "CHR", pos_col = "POS", effect_allele_col = "EFFECT_ALLELE", other_allele_col = "OTHER_ALLELE", samplesize_col = "N", beta_col = "BETA", se_col = "SE", pval_col = "P", eaf_col = "EAF") 

#EBV-EBNA 
data <- read.table("/projects/MRC-IEU/research/projects/ieu1/wp1/028/working/data/results/UKBB/Kachuri_infection_GWAS/ukb_seroreact_EBV_EBNA_5E08.txt", header = T)
exposure <- format_data(data, type ="exposure", header = TRUE, snp_col = "ID", chr_col = "CHR", pos_col = "POS", effect_allele_col = "EFFECT_ALLELE", other_allele_col = "OTHER_ALLELE", samplesize_col = "N", beta_col = "BETA", se_col = "SE", pval_col = "P", eaf_col = "EAF") 

#BKV VP1
data <- read.table("/projects/MRC-IEU/research/projects/ieu1/wp1/028/working/data/results/UKBB/Kachuri_infection_GWAS/ukb_seroreact_BKV_5E08.txt", header = T)
exposure <- format_data(data, type ="exposure", header = TRUE, snp_col = "ID", chr_col = "CHR", pos_col = "POS", effect_allele_col = "EFFECT_ALLELE", other_allele_col = "OTHER_ALLELE", samplesize_col = "N", beta_col = "BETA", se_col = "SE", pval_col = "P", eaf_col = "EAF") 

#VSV gEgI
data <- read.table("/projects/MRC-IEU/research/projects/ieu1/wp1/028/working/data/results/UKBB/Kachuri_infection_GWAS/ukb_seroreact_VZV_5E08.txt", header = T) 

#LD clumping 
try(exposure <- clump_data(exposure)) 

#Extract outcome data
ao<-available_outcomes()
id.cvd.out <- c("ieu-a-1008", "ieu-a-1006")
outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = id.cvd.out)

#Harmonise exposure and outcome data 
dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome_dat)

#MR analysis 
Res <- mr(dat) #IVW (default), Egger, WM, MODE, Wald ratio 
Res_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test 
Res_hetero <- mr_heterogeneity(dat) #
Res_single <- mr_singlesnp(dat) #single SNP analysis
or_results_res <- generate_odds_ratios(Res) #Odds ratio 
Res_MR_PRESSO <- run_mr_presso(dat) #Run if there is evidence of heterogeneity in Res_hetero findings 
Res_mr_leaveoneout <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)

#Write results 
write.table(Res,"MR_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_pleio,"MR_pleio_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_hetero,"MR_hetero_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_single,"MR_single_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(or_results_res,"MR_OR_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_MR_PRESSO, "MR_PRESSO_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F) 
write.table(Res_mr_leaveoneout, "Res_mr_leaveoneout_exposure_outcome", sep="\t",col.names=T,row.names=F,quote=F) 

##Two-sample MR: Platelet function on common infections##

#Extract platelet function exposures
ao<-available_outcomes()
id.platelet.exp <- c("ieu-a-1008", "ieu-a-1006")
exposure <- extract_instruments(id.platelet.exp)

#LD clumping 
try(exposure <- clump_data(exposure)) 

#Read in outcome data 
outcome_dat <- read_outcome_data(snps = exposure$SNP, filename = "outcome_data", sep = "\t", snp_col = "ID", beta_col = "BETA", se_col = "SE", pval_col = "P", effect_allele_col = "ALT", other_allele_col = "REF")
                                 
#Harmonise exposure and outcome data 
dat <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome_dat)

#MR analysis 
Res <- mr(dat) #IVW (default), Egger, WM, MODE, Wald ratio 
Res_pleio <- mr_pleiotropy_test(dat) # MR-Egger intercept test 
Res_hetero <- mr_heterogeneity(dat) #
Res_single <- mr_singlesnp(dat) #single SNP analysis
or_results_res <- generate_odds_ratios(Res) #Odds ratio 
Res_MR_PRESSO <- run_mr_presso(dat) #Run if there is evidence of heterogeneity in Res_hetero findings 
Res_mr_leaveoneout <- mr_leaveoneout(dat, parameters = default_parameters(), method = mr_ivw)

#Write results 
write.table(Res,"MR_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_pleio,"MR_pleio_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_hetero,"MR_hetero_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_single,"MR_single_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(or_results_res,"MR_OR_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F)
write.table(Res_MR_PRESSO, "MR_PRESSO_exposure_outcome",sep="\t",col.names=T,row.names=F,quote=F) 
write.table(Res_mr_leaveoneout, "Res_mr_leaveoneout_exposure_outcome", sep="\t",col.names=T,row.names=F,quote=F) 

##Check instrument strength## 

#Rename required columns
exposure$BetaXG<- exposure$beta.exposure
exposure$seBetaXG<- exposure$se.exposure
BetaXG   = exposure$BetaXG
seBetaXG = exposure$seBetaXG 
seBetaYG<-dat$se.outcome

BXG             = abs(BetaXG)         # gene-exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted = Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance. This can then be used to indicate whether SIMEX correction is needed for MR-Egger. If I2 is <90%, then this is typically recommended. The code for running SIMEX is below: 




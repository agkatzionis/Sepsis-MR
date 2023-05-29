
############################################################
##########   MENDELIAN RANDOMIZATION FOR SEPSIS   ##########
############################################################

## We repeat the main analysis conducted in the "Sepsis MR.R" 
## file without palindromic SNPs, in response to an issue 
## raised in GitHub. The analysis reported here is NOT included
## in the published paper.

## We work with the "MendelianRandomization" R package.
## install.packages("MendelianRandomization")
library(MendelianRandomization)

## Specify the working directory.
setwd("H:/.../Sepsis MR")

## Load the GWAS data for the various traits.
bmi <- read.delim("bmi_i.csv", header = TRUE, sep = ",")
ldl <- read.delim("ldl_i.csv", header = TRUE, sep = ",")
sbp <- read.delim("sbp_i.csv", header = TRUE, sep = ",")
smoking <- read.delim("smoking_i.csv", header = TRUE, sep = ",")
t2dm <- read.delim("t2dm_i.csv", header = TRUE, sep = ",")

## Load data for sepsis risk.
ukbb <- read.delim("Sepsis Data from UKBB.txt", header = TRUE, sep = "\t")
hunt <- read.delim("Sepsis Data from HUNT.txt", header = TRUE, sep = " ")

## Detect and remove palindromic SNPs in UKBB.
ukbb.palindromic <- which(ukbb$ALLELE0 == "A" & ukbb$ALLELE1 == "T" |
                          ukbb$ALLELE0 == "T" & ukbb$ALLELE1 == "A" |
                          ukbb$ALLELE0 == "G" & ukbb$ALLELE1 == "C" |
                          ukbb$ALLELE0 == "C" & ukbb$ALLELE1 == "G" )
length(ukbb.palindromic)
ukbb <- ukbb[- ukbb.palindromic, ]
dim(ukbb)

## Do the same for the HUNT dataset.
hunt.palindromic <- which(hunt$OA == "A" & hunt$EA == "T" |
                          hunt$OA == "T" & hunt$EA == "A" |
                          hunt$OA == "G" & hunt$EA == "C" |
                          hunt$OA == "C" & hunt$EA == "G" )
length(hunt.palindromic)
hunt <- hunt[- hunt.palindromic, ]
dim(hunt)

############################################################

## First, we run the analysis using UK Biobank data.

##########   UKBB - PRE-PROCESS THE DATA   ##########

## Subset UKBB data, get G-Y associations per trait.
sepsis.ukbb.bmi <- ukbb[as.character(ukbb$SNP) %in% as.character(bmi$SNP), ]
sepsis.ukbb.ldl <- ukbb[as.character(ukbb$SNP) %in% as.character(ldl$SNP), ]
sepsis.ukbb.sbp <- ukbb[as.character(ukbb$SNP) %in% as.character(sbp$SNP), ]
sepsis.ukbb.smk <- ukbb[as.character(ukbb$SNP) %in% as.character(smoking$SNP), ]
sepsis.ukbb.t2d <- ukbb[as.character(ukbb$SNP) %in% as.character(t2dm$SNP), ]

## Subset G-X associations to exclude any SNPs not included in UK Biobank.
## Since we have already discarded palindromic SNPs from UKBB, this will
## get rid of the palindromic SNPs in the other datasets.
bmi1 <- bmi[as.character(bmi$SNP) %in% as.character(sepsis.ukbb.bmi$SNP), ]
ldl1 <- ldl[as.character(ldl$SNP) %in% as.character(sepsis.ukbb.ldl$SNP), ]
sbp1 <- sbp[as.character(sbp$SNP) %in% as.character(sepsis.ukbb.sbp$SNP), ]
smoking1 <- smoking[as.character(smoking$SNP) %in% as.character(sepsis.ukbb.smk$SNP), ]
t2dm1 <- t2dm[as.character(t2dm$SNP) %in% as.character(sepsis.ukbb.t2d$SNP), ]

## Order all datasets by increasing rsID.
bmi1 <- bmi1[order(as.character(bmi1$SNP)), ]; sepsis.ukbb.bmi <- sepsis.ukbb.bmi[order(as.character(sepsis.ukbb.bmi$SNP)), ]
ldl1 <- ldl1[order(as.character(ldl1$SNP)), ]; sepsis.ukbb.ldl <- sepsis.ukbb.ldl[order(as.character(sepsis.ukbb.ldl$SNP)), ]
sbp1 <- sbp1[order(as.character(sbp1$SNP)), ]; sepsis.ukbb.sbp <- sepsis.ukbb.sbp[order(as.character(sepsis.ukbb.sbp$SNP)), ]
smoking1 <- smoking1[order(as.character(smoking1$SNP)), ]; sepsis.ukbb.smk <- sepsis.ukbb.smk[order(as.character(sepsis.ukbb.smk$SNP)), ]
t2dm1 <- t2dm1[order(as.character(t2dm1$SNP)), ]; sepsis.ukbb.t2d <- sepsis.ukbb.t2d[order(as.character(sepsis.ukbb.t2d$SNP)), ]



## Sanity check. Different numbers of SNPs for BMI and T2DM?
nrow(bmi1) == nrow(sepsis.ukbb.bmi)
nrow(ldl1) == nrow(sepsis.ukbb.ldl)
nrow(sbp1) == nrow(sepsis.ukbb.sbp)
nrow(smoking1) == nrow(sepsis.ukbb.smk)
nrow(t2dm1) == nrow(sepsis.ukbb.t2d)

## The variant rs78896587 is triallelic in the Biobank data for T2DM
## We manually discard the allele that does not appear in the G-X dataset.
sepsis.ukbb.t2d[155:156, ]
t2dm1[155, ]
sepsis.ukbb.t2d <- sepsis.ukbb.t2d[-155, ]



## Align the reference alleles across datasets.

## For BMI.
sepsis.ukbb.bmi2 <- sepsis.ukbb.bmi
for (i in 1:nrow(bmi1)) {
  if (as.character(sepsis.ukbb.bmi$ALLELE1[i]) != as.character(bmi1$bmi_ea[i])) {
    sepsis.ukbb.bmi2$BETA[i] <- - sepsis.ukbb.bmi$BETA[i]
    sepsis.ukbb.bmi2$A1FREQ[i] <- 1 - sepsis.ukbb.bmi$A1FREQ[i]
    sepsis.ukbb.bmi2$ALLELE1[i] <- sepsis.ukbb.bmi$ALLELE0[i]
    sepsis.ukbb.bmi2$ALLELE0[i] <- sepsis.ukbb.bmi$ALLELE1[i]
  }
}
sepsis.ukbb.bmi <- sepsis.ukbb.bmi2
rm(sepsis.ukbb.bmi2)

## For LDL-C.
sepsis.ukbb.ldl2 <- sepsis.ukbb.ldl
for (i in 1:nrow(ldl1)) {
  if (as.character(sepsis.ukbb.ldl$ALLELE1[i]) != as.character(ldl1$ldl_ea[i])) {
    sepsis.ukbb.ldl2$BETA[i] <- - sepsis.ukbb.ldl$BETA[i]
    sepsis.ukbb.ldl2$A1FREQ[i] <- 1 - sepsis.ukbb.ldl$A1FREQ[i]
    sepsis.ukbb.ldl2$ALLELE1[i] <- sepsis.ukbb.ldl$ALLELE0[i]
    sepsis.ukbb.ldl2$ALLELE0[i] <- sepsis.ukbb.ldl$ALLELE1[i]
  }
}
sepsis.ukbb.ldl <- sepsis.ukbb.ldl2
rm(sepsis.ukbb.ldl2)

## For SBP.
sepsis.ukbb.sbp2 <- sepsis.ukbb.sbp
for (i in 1:nrow(sbp1)) {
  if (as.character(sepsis.ukbb.sbp$ALLELE1[i]) != as.character(sbp1$sbp_ea[i])) {
    sepsis.ukbb.sbp2$BETA[i] <- - sepsis.ukbb.sbp$BETA[i]
    sepsis.ukbb.sbp2$A1FREQ[i] <- 1 - sepsis.ukbb.sbp$A1FREQ[i]
    sepsis.ukbb.sbp2$ALLELE1[i] <- sepsis.ukbb.sbp$ALLELE0[i]
    sepsis.ukbb.sbp2$ALLELE0[i] <- sepsis.ukbb.sbp$ALLELE1[i]
  }
}
sepsis.ukbb.sbp <- sepsis.ukbb.sbp2
rm(sepsis.ukbb.sbp2)

## For Smoking.
sepsis.ukbb.smk2 <- sepsis.ukbb.smk
for (i in 1:nrow(smoking1)) {
  if (as.character(sepsis.ukbb.smk$ALLELE1[i]) != as.character(smoking1$smoking_ea[i])) {
    sepsis.ukbb.smk2$BETA[i] <- - sepsis.ukbb.smk$BETA[i]
    sepsis.ukbb.smk2$A1FREQ[i] <- 1 - sepsis.ukbb.smk$A1FREQ[i]
    sepsis.ukbb.smk2$ALLELE1[i] <- sepsis.ukbb.smk$ALLELE0[i]
    sepsis.ukbb.smk2$ALLELE0[i] <- sepsis.ukbb.smk$ALLELE1[i]
  }
}
sepsis.ukbb.smk <- sepsis.ukbb.smk2
rm(sepsis.ukbb.smk2)

## For T2DM.
sepsis.ukbb.t2d2 <- sepsis.ukbb.t2d
for (i in 1:nrow(t2dm1)) {
  if (as.character(sepsis.ukbb.t2d$ALLELE1[i]) != as.character(t2dm1$t2dm_ea[i])) {
    sepsis.ukbb.t2d2$BETA[i] <- - sepsis.ukbb.t2d$BETA[i]
    sepsis.ukbb.t2d2$A1FREQ[i] <- 1 - sepsis.ukbb.t2d$A1FREQ[i]
    sepsis.ukbb.t2d2$ALLELE1[i] <- sepsis.ukbb.t2d$ALLELE0[i]
    sepsis.ukbb.t2d2$ALLELE0[i] <- sepsis.ukbb.t2d$ALLELE1[i]
  }
}
sepsis.ukbb.t2d <- sepsis.ukbb.t2d2
rm(sepsis.ukbb.t2d2)


## Sanity check. Confirms that data are aligned.
all(as.character(sepsis.ukbb.bmi$SNP) == as.character(bmi1$SNP)); all(as.character(sepsis.ukbb.bmi$ALLELE1) == as.character(bmi1$bmi_ea))
all(as.character(sepsis.ukbb.ldl$SNP) == as.character(ldl1$SNP)); all(as.character(sepsis.ukbb.ldl$ALLELE1) == as.character(ldl1$ldl_ea))
all(as.character(sepsis.ukbb.sbp$SNP) == as.character(sbp1$SNP)); all(as.character(sepsis.ukbb.sbp$ALLELE1) == as.character(sbp1$sbp_ea))
all(as.character(sepsis.ukbb.smk$SNP) == as.character(smoking1$SNP)); all(as.character(sepsis.ukbb.smk$ALLELE1) == as.character(smoking1$smoking_ea))
all(as.character(sepsis.ukbb.t2d$SNP) == as.character(t2dm1$SNP)); all(as.character(sepsis.ukbb.t2d$ALLELE1) == as.character(t2dm1$t2dm_ea))

###############################################

##########   UKBB - MR ANALYSIS   ##########

## We use functions from the R package MendelianRandomization.
## We implement the following MR methods:
## IVW, MR-Egger, mode, median and the contamination mixture.

## Note: mr_allmethods performs IVW using a random-effects model
## by default when using more than 3 SNPs, which is what we want.

## BMI.
bmi1.main <- mr_allmethods( mr_input(bx = bmi1$bmi_beta, bxse = bmi1$bmi_se, by = sepsis.ukbb.bmi$BETA, byse = sepsis.ukbb.bmi$SE), method = "main" )
bmi1.mode1 <- mr_mbe( mr_input(bx = bmi1$bmi_beta, bxse = bmi1$bmi_se, by = sepsis.ukbb.bmi$BETA, byse = sepsis.ukbb.bmi$SE), weighting = "unweighted" )
bmi1.mode2 <- mr_mbe( mr_input(bx = bmi1$bmi_beta, bxse = bmi1$bmi_se, by = sepsis.ukbb.bmi$BETA, byse = sepsis.ukbb.bmi$SE), weighting = "weighted" )
bmi1.conmix <- mr_conmix( mr_input(bx = bmi1$bmi_beta, bxse = bmi1$bmi_se, by = sepsis.ukbb.bmi$BETA, byse = sepsis.ukbb.bmi$SE), CIStep = 0.001 )

## LDL-C.
ldl1.main <- mr_allmethods( mr_input(bx = ldl1$ldl_beta, bxse = ldl1$ldl_se, by = sepsis.ukbb.ldl$BETA, byse = sepsis.ukbb.ldl$SE), method = "main" )
ldl1.mode1 <- mr_mbe( mr_input(bx = ldl1$ldl_beta, bxse = ldl1$ldl_se, by = sepsis.ukbb.ldl$BETA, byse = sepsis.ukbb.ldl$SE), weighting = "unweighted" )
ldl1.mode2 <- mr_mbe( mr_input(bx = ldl1$ldl_beta, bxse = ldl1$ldl_se, by = sepsis.ukbb.ldl$BETA, byse = sepsis.ukbb.ldl$SE), weighting = "weighted" )
ldl1.conmix <- mr_conmix( mr_input(bx = ldl1$ldl_beta, bxse = ldl1$ldl_se, by = sepsis.ukbb.ldl$BETA, byse = sepsis.ukbb.ldl$SE), CIStep = 0.001 )

## SBP.
sbp1.main <- mr_allmethods( mr_input(bx = sbp1$sbp_beta, bxse = sbp1$sbp_se, by = sepsis.ukbb.sbp$BETA, byse = sepsis.ukbb.sbp$SE), method = "main" )
sbp1.mode1 <- mr_mbe( mr_input(bx = sbp1$sbp_beta, bxse = sbp1$sbp_se, by = sepsis.ukbb.sbp$BETA, byse = sepsis.ukbb.sbp$SE), weighting = "unweighted" )
sbp1.mode2 <- mr_mbe( mr_input(bx = sbp1$sbp_beta, bxse = sbp1$sbp_se, by = sepsis.ukbb.sbp$BETA, byse = sepsis.ukbb.sbp$SE), weighting = "weighted" )
sbp1.conmix <- mr_conmix( mr_input(bx = sbp1$sbp_beta, bxse = sbp1$sbp_se, by = sepsis.ukbb.sbp$BETA, byse = sepsis.ukbb.sbp$SE), CIStep = 0.001 )

## Smoking.
smk1.main <- mr_allmethods( mr_input(bx = smoking1$smoking_beta, bxse = smoking1$smoking_se, by = sepsis.ukbb.smk$BETA, byse = sepsis.ukbb.smk$SE), method = "main" )
smk1.mode1 <- mr_mbe( mr_input(bx = smoking1$smoking_beta, bxse = smoking1$smoking_se, by = sepsis.ukbb.smk$BETA, byse = sepsis.ukbb.smk$SE), weighting = "unweighted" )
smk1.mode2 <- mr_mbe( mr_input(bx = smoking1$smoking_beta, bxse = smoking1$smoking_se, by = sepsis.ukbb.smk$BETA, byse = sepsis.ukbb.smk$SE), weighting = "weighted" )
smk1.conmix <- mr_conmix( mr_input(bx = smoking1$smoking_beta, bxse = smoking1$smoking_se, by = sepsis.ukbb.smk$BETA, byse = sepsis.ukbb.smk$SE), CIStep = 0.001 )

## T2DM.
t2d1.main <- mr_allmethods( mr_input(bx = t2dm1$t2dm_beta, bxse = t2dm1$t2dm_se, by = sepsis.ukbb.t2d$BETA, byse = sepsis.ukbb.t2d$SE), method = "main" )
t2d1.mode1 <- mr_mbe( mr_input(bx = t2dm1$t2dm_beta, bxse = t2dm1$t2dm_se, by = sepsis.ukbb.t2d$BETA, byse = sepsis.ukbb.t2d$SE), weighting = "unweighted" )
t2d1.mode2 <- mr_mbe( mr_input(bx = t2dm1$t2dm_beta, bxse = t2dm1$t2dm_se, by = sepsis.ukbb.t2d$BETA, byse = sepsis.ukbb.t2d$SE), weighting = "weighted" )
t2d1.conmix <- mr_conmix( mr_input(bx = t2dm1$t2dm_beta, bxse = t2dm1$t2dm_se, by = sepsis.ukbb.t2d$BETA, byse = sepsis.ukbb.t2d$SE), CIStep = 0.001 )

###############################################

##########   UKBB - RESULTS TABLES   ##########

## This creates results tables for all methods per trait.
## We report both default MR estimates (log-odds ratios)
## and results in an odds-ratio scale (exponentiated).


## Label the methods.
method.names <- c(bmi1.main$Values[, 1], "Simple mode", "Weighted mode", "Con mix")

## BMI.
bmi1.est <- round(c(bmi1.main$Values[, 2], bmi1.mode1$Estimate, bmi1.mode2$Estimate, bmi1.conmix$Estimate), digits = 3)
bmi1.se <- round(c(bmi1.main$Values[, 3], bmi1.mode1$StdError, bmi1.mode2$StdError, NA), digits = 3)
bmi1.ci1 <- round(c(bmi1.main$Values[, 4], bmi1.mode1$CILower, bmi1.mode2$CILower, bmi1.conmix$CILower), digits = 3)
bmi1.ci2 <- round(c(bmi1.main$Values[, 5], bmi1.mode1$CIUpper, bmi1.mode2$CIUpper, bmi1.conmix$CIUpper), digits = 3)
bmi1.p <- round(c(bmi1.main$Values[, 6], bmi1.mode1$Pvalue, bmi1.mode2$Pvalue, bmi1.conmix$Pvalue), digits = 3)
bmi1.results <- data.frame(method.names, bmi1.est, bmi1.se, bmi1.ci1, bmi1.ci2, bmi1.p)
bmi1.results <- cbind(bmi1.results, round(exp(bmi1.results[, c(2, 4, 5)]), digits = 3))
bmi1.results[5, 7:9] <- NA
names(bmi1.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
bmi1.results

## LDL-C.
ldl1.est <- round(c(ldl1.main$Values[, 2], ldl1.mode1$Estimate, ldl1.mode2$Estimate, ldl1.conmix$Estimate), digits = 3)
ldl1.se <- round(c(ldl1.main$Values[, 3], ldl1.mode1$StdError, ldl1.mode2$StdError, NA), digits = 3)
ldl1.ci1 <- round(c(ldl1.main$Values[, 4], ldl1.mode1$CILower, ldl1.mode2$CILower, ldl1.conmix$CILower), digits = 3)
ldl1.ci2 <- round(c(ldl1.main$Values[, 5], ldl1.mode1$CIUpper, ldl1.mode2$CIUpper, ldl1.conmix$CIUpper), digits = 3)
ldl1.p <- round(c(ldl1.main$Values[, 6], ldl1.mode1$Pvalue, ldl1.mode2$Pvalue, ldl1.conmix$Pvalue), digits = 3)
ldl1.results <- data.frame(method.names, ldl1.est, ldl1.se, ldl1.ci1, ldl1.ci2, ldl1.p)
ldl1.results <- cbind(ldl1.results, round(exp(ldl1.results[, c(2, 4, 5)]), digits = 3))
ldl1.results[5, 7:9] <- NA
names(ldl1.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
ldl1.results

## SBP.
sbp1.est <- round(c(sbp1.main$Values[, 2], sbp1.mode1$Estimate, sbp1.mode2$Estimate, sbp1.conmix$Estimate), digits = 3)
sbp1.se <- round(c(sbp1.main$Values[, 3], sbp1.mode1$StdError, sbp1.mode2$StdError, NA), digits = 3)
sbp1.ci1 <- round(c(sbp1.main$Values[, 4], sbp1.mode1$CILower, sbp1.mode2$CILower, sbp1.conmix$CILower), digits = 3)
sbp1.ci2 <- round(c(sbp1.main$Values[, 5], sbp1.mode1$CIUpper, sbp1.mode2$CIUpper, sbp1.conmix$CIUpper), digits = 3)
sbp1.p <- round(c(sbp1.main$Values[, 6], sbp1.mode1$Pvalue, sbp1.mode2$Pvalue, sbp1.conmix$Pvalue), digits = 3)
sbp1.results <- data.frame(method.names, sbp1.est, sbp1.se, sbp1.ci1, sbp1.ci2, sbp1.p)
sbp1.results <- cbind(sbp1.results, round(exp(sbp1.results[, c(2, 4, 5)]), digits = 3))
sbp1.results[5, 7:9] <- NA
names(sbp1.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
sbp1.results

## Smoking.
smk1.est <- round(c(smk1.main$Values[, 2], smk1.mode1$Estimate, smk1.mode2$Estimate, smk1.conmix$Estimate), digits = 3)
smk1.se <- round(c(smk1.main$Values[, 3], smk1.mode1$StdError, smk1.mode2$StdError, NA), digits = 3)
smk1.ci1 <- round(c(smk1.main$Values[, 4], smk1.mode1$CILower, smk1.mode2$CILower, smk1.conmix$CILower), digits = 3)
smk1.ci2 <- round(c(smk1.main$Values[, 5], smk1.mode1$CIUpper, smk1.mode2$CIUpper, smk1.conmix$CIUpper), digits = 3)
smk1.p <- round(c(smk1.main$Values[, 6], smk1.mode1$Pvalue, smk1.mode2$Pvalue, smk1.conmix$Pvalue), digits = 3)
smk1.results <- data.frame(method.names, smk1.est, smk1.se, smk1.ci1, smk1.ci2, smk1.p)
smk1.results <- cbind(smk1.results, round(exp(smk1.results[, c(2, 4, 5)]), digits = 3))
smk1.results[5, 7:9] <- NA
names(smk1.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
smk1.results

## T2DM.
t2d1.est <- round(c(t2d1.main$Values[, 2], t2d1.mode1$Estimate, t2d1.mode2$Estimate, t2d1.conmix$Estimate), digits = 3)
t2d1.se <- round(c(t2d1.main$Values[, 3], t2d1.mode1$StdError, t2d1.mode2$StdError, NA), digits = 3)
t2d1.ci1 <- round(c(t2d1.main$Values[, 4], t2d1.mode1$CILower, t2d1.mode2$CILower, t2d1.conmix$CILower), digits = 3)
t2d1.ci2 <- round(c(t2d1.main$Values[, 5], t2d1.mode1$CIUpper, t2d1.mode2$CIUpper, t2d1.conmix$CIUpper), digits = 3)
t2d1.p <- round(c(t2d1.main$Values[, 6], t2d1.mode1$Pvalue, t2d1.mode2$Pvalue, t2d1.conmix$Pvalue), digits = 3)
t2d1.results <- data.frame(method.names, t2d1.est, t2d1.se, t2d1.ci1, t2d1.ci2, t2d1.p)
t2d1.results <- cbind(t2d1.results, round(exp(t2d1.results[, c(2, 4, 5)]), digits = 3))
t2d1.results[5, 7:9] <- NA
names(t2d1.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
t2d1.results

###############################################

## Now we run the analysis using HUNT data.

##########   HUNT - PRE-PROCESS THE DATA   ##########

## Subset HUNT data, get G-Y associations per trait.
sepsis.hunt.bmi <- hunt[as.character(hunt$SNP) %in% as.character(bmi$SNP), ]
sepsis.hunt.ldl <- hunt[as.character(hunt$SNP) %in% as.character(ldl$SNP), ]
sepsis.hunt.sbp <- hunt[as.character(hunt$SNP) %in% as.character(sbp$SNP), ]
sepsis.hunt.smk <- hunt[as.character(hunt$SNP) %in% as.character(smoking$SNP), ]
sepsis.hunt.t2d <- hunt[as.character(hunt$SNP) %in% as.character(t2dm$SNP), ]

## Subset G-X associations to exclude any SNPs not included in HUNT.
## Since we have already discarded palindromic SNPs from HUNT, this will
## get rid of the palindromic SNPs in the other datasets.
bmi2 <- bmi[as.character(bmi$SNP) %in% as.character(sepsis.hunt.bmi$SNP), ]
ldl2 <- ldl[as.character(ldl$SNP) %in% as.character(sepsis.hunt.ldl$SNP), ]
sbp2 <- sbp[as.character(sbp$SNP) %in% as.character(sepsis.hunt.sbp$SNP), ]
smoking2 <- smoking[as.character(smoking$SNP) %in% as.character(sepsis.hunt.smk$SNP), ]
t2dm2 <- t2dm[as.character(t2dm$SNP) %in% as.character(sepsis.hunt.t2d$SNP), ]

## Order all datasets by increasing rsID.
bmi2 <- bmi2[order(as.character(bmi2$SNP)), ]; sepsis.hunt.bmi <- sepsis.hunt.bmi[order(as.character(sepsis.hunt.bmi$SNP)), ]
ldl2 <- ldl2[order(as.character(ldl2$SNP)), ]; sepsis.hunt.ldl <- sepsis.hunt.ldl[order(as.character(sepsis.hunt.ldl$SNP)), ]
sbp2 <- sbp2[order(as.character(sbp2$SNP)), ]; sepsis.hunt.sbp <- sepsis.hunt.sbp[order(as.character(sepsis.hunt.sbp$SNP)), ]
smoking2 <- smoking2[order(as.character(smoking2$SNP)), ]; sepsis.hunt.smk <- sepsis.hunt.smk[order(as.character(sepsis.hunt.smk$SNP)), ]
t2dm2 <- t2dm2[order(as.character(t2dm2$SNP)), ]; sepsis.hunt.t2d <- sepsis.hunt.t2d[order(as.character(sepsis.hunt.t2d$SNP)), ]



## Sanity check. Same number of SNPs per trait - no manual adjustments needed.
nrow(bmi2) == nrow(sepsis.hunt.bmi)
nrow(ldl2) == nrow(sepsis.hunt.ldl)
nrow(sbp2) == nrow(sepsis.hunt.sbp)
nrow(smoking2) == nrow(sepsis.hunt.smk)
nrow(t2dm2) == nrow(sepsis.hunt.t2d)



## Align the reference alleles across datasets.

## Note: we are aligning according to the effect allele in the exposure GWAS 
## so data for UKBB and HUNT should be alinged in the same way.

## For BMI.
sepsis.hunt.bmi2 <- sepsis.hunt.bmi
for (i in 1:nrow(bmi2)) {
  if (as.character(sepsis.hunt.bmi$EA[i]) != as.character(bmi2$bmi_ea[i])) {
    sepsis.hunt.bmi2$Beta[i] <- - sepsis.hunt.bmi$Beta[i]
    sepsis.hunt.bmi2$EA[i] <- sepsis.hunt.bmi$OA[i]
    sepsis.hunt.bmi2$OA[i] <- sepsis.hunt.bmi$EA[i]
    sepsis.hunt.bmi2$EAF[i] <- 1 - sepsis.hunt.bmi$EAF[i]
  }
}
sepsis.hunt.bmi <- sepsis.hunt.bmi2
rm(sepsis.hunt.bmi2)

## For LDL-C.
sepsis.hunt.ldl2 <- sepsis.hunt.ldl
for (i in 1:nrow(ldl2)) {
  if (as.character(sepsis.hunt.ldl$EA[i]) != as.character(ldl2$ldl_ea[i])) {
    sepsis.hunt.ldl2$Beta[i] <- - sepsis.hunt.ldl$Beta[i]
    sepsis.hunt.ldl2$EA[i] <- sepsis.hunt.ldl$OA[i]
    sepsis.hunt.ldl2$OA[i] <- sepsis.hunt.ldl$EA[i]
    sepsis.hunt.ldl2$EAF[i] <- 1 - sepsis.hunt.ldl$EAF[i]
  }
}
sepsis.hunt.ldl <- sepsis.hunt.ldl2
rm(sepsis.hunt.ldl2)

## For SBP.
sepsis.hunt.sbp2 <- sepsis.hunt.sbp
for (i in 1:nrow(sbp2)) {
  if (as.character(sepsis.hunt.sbp$EA[i]) != as.character(sbp2$sbp_ea[i])) {
    sepsis.hunt.sbp2$Beta[i] <- - sepsis.hunt.sbp$Beta[i]
    sepsis.hunt.sbp2$EA[i] <- sepsis.hunt.sbp$OA[i]
    sepsis.hunt.sbp2$OA[i] <- sepsis.hunt.sbp$EA[i]
    sepsis.hunt.sbp2$EAF[i] <- 1 - sepsis.hunt.sbp$EAF[i]
  }
}
sepsis.hunt.sbp <- sepsis.hunt.sbp2
rm(sepsis.hunt.sbp2)

## For Smoking.
sepsis.hunt.smk2 <- sepsis.hunt.smk
for (i in 1:nrow(smoking2)) {
  if (as.character(sepsis.hunt.smk$EA[i]) != as.character(smoking2$smoking_ea[i])) {
    sepsis.hunt.smk2$Beta[i] <- - sepsis.hunt.smk$Beta[i]
    sepsis.hunt.smk2$EA[i] <- sepsis.hunt.smk$OA[i]
    sepsis.hunt.smk2$OA[i] <- sepsis.hunt.smk$EA[i]
    sepsis.hunt.smk2$EAF[i] <- 1 - sepsis.hunt.smk$EAF[i]
  }
}
sepsis.hunt.smk <- sepsis.hunt.smk2
rm(sepsis.hunt.smk2)

## For T2DM.
sepsis.hunt.t2d2 <- sepsis.hunt.t2d
for (i in 1:nrow(t2dm2)) {
  if (as.character(sepsis.hunt.t2d$EA[i]) != as.character(t2dm2$t2dm_ea[i])) {
    sepsis.hunt.t2d2$Beta[i] <- - sepsis.hunt.t2d$Beta[i]
    sepsis.hunt.t2d2$EA[i] <- sepsis.hunt.t2d$OA[i]
    sepsis.hunt.t2d2$OA[i] <- sepsis.hunt.t2d$EA[i]
    sepsis.hunt.t2d2$EAF[i] <- 1 - sepsis.hunt.t2d$EAF[i]
  }
}
sepsis.hunt.t2d <- sepsis.hunt.t2d2
rm(sepsis.hunt.t2d2)

## Sanity check. Confirm that data are aligned.
all(as.character(sepsis.hunt.bmi$SNP) == as.character(bmi2$SNP)); all(as.character(sepsis.hunt.bmi$EA) == as.character(bmi2$bmi_ea))
all(as.character(sepsis.hunt.ldl$SNP) == as.character(ldl2$SNP)); all(as.character(sepsis.hunt.ldl$EA) == as.character(ldl2$ldl_ea))
all(as.character(sepsis.hunt.sbp$SNP) == as.character(sbp2$SNP)); all(as.character(sepsis.hunt.sbp$EA) == as.character(sbp2$sbp_ea))
all(as.character(sepsis.hunt.smk$SNP) == as.character(smoking2$SNP)); all(as.character(sepsis.hunt.smk$EA) == as.character(smoking2$smoking_ea))
all(as.character(sepsis.hunt.t2d$SNP) == as.character(t2dm2$SNP)); all(as.character(sepsis.hunt.t2d$EA) == as.character(t2dm2$t2dm_ea))

###############################################

##########   HUNT - MR ANALYSIS   ##########

## Same as for the UKBB data.

## BMI.
bmi2.main <- mr_allmethods( mr_input(bx = bmi2$bmi_beta, bxse = bmi2$bmi_se, by = sepsis.hunt.bmi$Beta, byse = sepsis.hunt.bmi$SE), method = "main" )
bmi2.mode1 <- mr_mbe( mr_input(bx = bmi2$bmi_beta, bxse = bmi2$bmi_se, by = sepsis.hunt.bmi$Beta, byse = sepsis.hunt.bmi$SE), weighting = "unweighted" )
bmi2.mode2 <- mr_mbe( mr_input(bx = bmi2$bmi_beta, bxse = bmi2$bmi_se, by = sepsis.hunt.bmi$Beta, byse = sepsis.hunt.bmi$SE), weighting = "weighted" )
bmi2.conmix <- mr_conmix( mr_input(bx = bmi2$bmi_beta, bxse = bmi2$bmi_se, by = sepsis.hunt.bmi$Beta, byse = sepsis.hunt.bmi$SE), CIStep = 0.001 )

## LDL-C.
ldl2.main <- mr_allmethods( mr_input(bx = ldl2$ldl_beta, bxse = ldl2$ldl_se, by = sepsis.hunt.ldl$Beta, byse = sepsis.hunt.ldl$SE), method = "main" )
ldl2.mode1 <- mr_mbe( mr_input(bx = ldl2$ldl_beta, bxse = ldl2$ldl_se, by = sepsis.hunt.ldl$Beta, byse = sepsis.hunt.ldl$SE), weighting = "unweighted" )
ldl2.mode2 <- mr_mbe( mr_input(bx = ldl2$ldl_beta, bxse = ldl2$ldl_se, by = sepsis.hunt.ldl$Beta, byse = sepsis.hunt.ldl$SE), weighting = "weighted" )
ldl2.conmix <- mr_conmix( mr_input(bx = ldl2$ldl_beta, bxse = ldl2$ldl_se, by = sepsis.hunt.ldl$Beta, byse = sepsis.hunt.ldl$SE), CIStep = 0.001 )

## SBP.
sbp2.main <- mr_allmethods( mr_input(bx = sbp2$sbp_beta, bxse = sbp2$sbp_se, by = sepsis.hunt.sbp$Beta, byse = sepsis.hunt.sbp$SE), method = "main" )
sbp2.mode1 <- mr_mbe( mr_input(bx = sbp2$sbp_beta, bxse = sbp2$sbp_se, by = sepsis.hunt.sbp$Beta, byse = sepsis.hunt.sbp$SE), weighting = "unweighted" )
sbp2.mode2 <- mr_mbe( mr_input(bx = sbp2$sbp_beta, bxse = sbp2$sbp_se, by = sepsis.hunt.sbp$Beta, byse = sepsis.hunt.sbp$SE), weighting = "weighted" )
sbp2.conmix <- mr_conmix( mr_input(bx = sbp2$sbp_beta, bxse = sbp2$sbp_se, by = sepsis.hunt.sbp$Beta, byse = sepsis.hunt.sbp$SE), CIStep = 0.001 )

## Smoking.
smk2.main <- mr_allmethods( mr_input(bx = smoking2$smoking_beta, bxse = smoking2$smoking_se, by = sepsis.hunt.smk$Beta, byse = sepsis.hunt.smk$SE), method = "main" )
smk2.mode1 <- mr_mbe( mr_input(bx = smoking2$smoking_beta, bxse = smoking2$smoking_se, by = sepsis.hunt.smk$Beta, byse = sepsis.hunt.smk$SE), weighting = "unweighted" )
smk2.mode2 <- mr_mbe( mr_input(bx = smoking2$smoking_beta, bxse = smoking2$smoking_se, by = sepsis.hunt.smk$Beta, byse = sepsis.hunt.smk$SE), weighting = "weighted" )
smk2.conmix <- mr_conmix( mr_input(bx = smoking2$smoking_beta, bxse = smoking2$smoking_se, by = sepsis.hunt.smk$Beta, byse = sepsis.hunt.smk$SE), CIStep = 0.001 )

## T2DM.
t2d2.main <- mr_allmethods( mr_input(bx = t2dm2$t2dm_beta, bxse = t2dm2$t2dm_se, by = sepsis.hunt.t2d$Beta, byse = sepsis.hunt.t2d$SE), method = "main" )
t2d2.mode1 <- mr_mbe( mr_input(bx = t2dm2$t2dm_beta, bxse = t2dm2$t2dm_se, by = sepsis.hunt.t2d$Beta, byse = sepsis.hunt.t2d$SE), weighting = "unweighted" )
t2d2.mode2 <- mr_mbe( mr_input(bx = t2dm2$t2dm_beta, bxse = t2dm2$t2dm_se, by = sepsis.hunt.t2d$Beta, byse = sepsis.hunt.t2d$SE), weighting = "weighted" )
t2d2.conmix <- mr_conmix( mr_input(bx = t2dm2$t2dm_beta, bxse = t2dm2$t2dm_se, by = sepsis.hunt.t2d$Beta, byse = sepsis.hunt.t2d$SE), CIStep = 0.001 )

###############################################

##########   HUNT - RESULTS TABLES   ##########

## Again, we report both default MR estimates and odds ratios.


## Name the methods.
method.names <- c(bmi2.main$Values[, 1], "Simple mode", "Weighted mode", "Con mix")

## BMI.
bmi2.est <- round(c(bmi2.main$Values[, 2], bmi2.mode1$Estimate, bmi2.mode2$Estimate, bmi2.conmix$Estimate), digits = 3)
bmi2.se <- round(c(bmi2.main$Values[, 3], bmi2.mode1$StdError, bmi2.mode2$StdError, NA), digits = 3)
bmi2.ci1 <- round(c(bmi2.main$Values[, 4], bmi2.mode1$CILower, bmi2.mode2$CILower, bmi2.conmix$CILower), digits = 3)
bmi2.ci2 <- round(c(bmi2.main$Values[, 5], bmi2.mode1$CIUpper, bmi2.mode2$CIUpper, bmi2.conmix$CIUpper), digits = 3)
bmi2.p <- round(c(bmi2.main$Values[, 6], bmi2.mode1$Pvalue, bmi2.mode2$Pvalue, bmi2.conmix$Pvalue), digits = 3)
bmi2.results <- data.frame(method.names, bmi2.est, bmi2.se, bmi2.ci1, bmi2.ci2, bmi2.p)
bmi2.results <- cbind(bmi2.results, round(exp(bmi2.results[, c(2, 4, 5)]), digits = 3))
bmi2.results[5, 7:9] <- NA
names(bmi2.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
bmi2.results

## LDL-C.
ldl2.est <- round(c(ldl2.main$Values[, 2], ldl2.mode1$Estimate, ldl2.mode2$Estimate, ldl2.conmix$Estimate), digits = 3)
ldl2.se <- round(c(ldl2.main$Values[, 3], ldl2.mode1$StdError, ldl2.mode2$StdError, NA), digits = 3)
ldl2.ci1 <- round(c(ldl2.main$Values[, 4], ldl2.mode1$CILower, ldl2.mode2$CILower, ldl2.conmix$CILower), digits = 3)
ldl2.ci2 <- round(c(ldl2.main$Values[, 5], ldl2.mode1$CIUpper, ldl2.mode2$CIUpper, ldl2.conmix$CIUpper), digits = 3)
ldl2.p <- round(c(ldl2.main$Values[, 6], ldl2.mode1$Pvalue, ldl2.mode2$Pvalue, ldl2.conmix$Pvalue), digits = 3)
ldl2.results <- data.frame(method.names, ldl2.est, ldl2.se, ldl2.ci1, ldl2.ci2, ldl2.p)
ldl2.results <- cbind(ldl2.results, round(exp(ldl2.results[, c(2, 4, 5)]), digits = 3))
ldl2.results[5, 7:9] <- NA
names(ldl2.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
ldl2.results

## SBP.
sbp2.est <- round(c(sbp2.main$Values[, 2], sbp2.mode1$Estimate, sbp2.mode2$Estimate, sbp2.conmix$Estimate), digits = 3)
sbp2.se <- round(c(sbp2.main$Values[, 3], sbp2.mode1$StdError, sbp2.mode2$StdError, NA), digits = 3)
sbp2.ci1 <- round(c(sbp2.main$Values[, 4], sbp2.mode1$CILower, sbp2.mode2$CILower, sbp2.conmix$CILower), digits = 3)
sbp2.ci2 <- round(c(sbp2.main$Values[, 5], sbp2.mode1$CIUpper, sbp2.mode2$CIUpper, sbp2.conmix$CIUpper), digits = 3)
sbp2.p <- round(c(sbp2.main$Values[, 6], sbp2.mode1$Pvalue, sbp2.mode2$Pvalue, sbp2.conmix$Pvalue), digits = 3)
sbp2.results <- data.frame(method.names, sbp2.est, sbp2.se, sbp2.ci1, sbp2.ci2, sbp2.p)
sbp2.results <- cbind(sbp2.results, round(exp(sbp2.results[, c(2, 4, 5)]), digits = 3))
sbp2.results[5, 7:9] <- NA
names(sbp2.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
sbp2.results

## Smoking.
smk2.est <- round(c(smk2.main$Values[, 2], smk2.mode1$Estimate, smk2.mode2$Estimate, smk2.conmix$Estimate), digits = 3)
smk2.se <- round(c(smk2.main$Values[, 3], smk2.mode1$StdError, smk2.mode2$StdError, NA), digits = 3)
smk2.ci1 <- round(c(smk2.main$Values[, 4], smk2.mode1$CILower, smk2.mode2$CILower, smk2.conmix$CILower), digits = 3)
smk2.ci2 <- round(c(smk2.main$Values[, 5], smk2.mode1$CIUpper, smk2.mode2$CIUpper, smk2.conmix$CIUpper), digits = 3)
smk2.p <- round(c(smk2.main$Values[, 6], smk2.mode1$Pvalue, smk2.mode2$Pvalue, smk2.conmix$Pvalue), digits = 3)
smk2.results <- data.frame(method.names, smk2.est, smk2.se, smk2.ci1, smk2.ci2, smk2.p)
smk2.results <- cbind(smk2.results, round(exp(smk2.results[, c(2, 4, 5)]), digits = 3))
smk2.results[5, 7:9] <- NA
names(smk2.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
smk2.results

## T2DM.
t2d2.est <- round(c(t2d2.main$Values[, 2], t2d2.mode1$Estimate, t2d2.mode2$Estimate, t2d2.conmix$Estimate), digits = 3)
t2d2.se <- round(c(t2d2.main$Values[, 3], t2d2.mode1$StdError, t2d2.mode2$StdError, NA), digits = 3)
t2d2.ci1 <- round(c(t2d2.main$Values[, 4], t2d2.mode1$CILower, t2d2.mode2$CILower, t2d2.conmix$CILower), digits = 3)
t2d2.ci2 <- round(c(t2d2.main$Values[, 5], t2d2.mode1$CIUpper, t2d2.mode2$CIUpper, t2d2.conmix$CIUpper), digits = 3)
t2d2.p <- round(c(t2d2.main$Values[, 6], t2d2.mode1$Pvalue, t2d2.mode2$Pvalue, t2d2.conmix$Pvalue), digits = 3)
t2d2.results <- data.frame(method.names, t2d2.est, t2d2.se, t2d2.ci1, t2d2.ci2, t2d2.p)
t2d2.results <- cbind(t2d2.results, round(exp(t2d2.results[, c(2, 4, 5)]), digits = 3))
t2d2.results[5, 7:9] <- NA
names(t2d2.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
t2d2.results

###############################################

## Now create figures for the manuscript.

## We use the ".tiff" format for the figures (and "lzw" compression to save space).
## We plot point estimates and 95% CI and report numerical values next to each plot.

###############################################

##########   FIGURE 1   ##########

## Plot random-effects IVW estimates for UKBB and HUNT.

## UK Biobank data.
ivw.means1 <- exp(rev(c(bmi1.main$Values[3, 2], ldl1.main$Values[3, 2], sbp1.main$Values[3, 2], smk1.main$Values[3, 2], t2d1.main$Values[3, 2])))
ivw.ci11 <- exp(rev(c(bmi1.main$Values[3, 4], ldl1.main$Values[3, 4], sbp1.main$Values[3, 4], smk1.main$Values[3, 4], t2d1.main$Values[3, 4])))
ivw.ci21 <- exp(rev(c(bmi1.main$Values[3, 5], ldl1.main$Values[3, 5], sbp1.main$Values[3, 5], smk1.main$Values[3, 5], t2d1.main$Values[3, 5])))

## HUNT data.
ivw.means2 <- exp(rev(c(bmi2.main$Values[3, 2], ldl2.main$Values[3, 2], sbp2.main$Values[3, 2], smk2.main$Values[3, 2], t2d2.main$Values[3, 2])))
ivw.ci12 <- exp(rev(c(bmi2.main$Values[3, 4], ldl2.main$Values[3, 4], sbp2.main$Values[3, 4], smk2.main$Values[3, 4], t2d2.main$Values[3, 4])))
ivw.ci22 <- exp(rev(c(bmi2.main$Values[3, 5], ldl2.main$Values[3, 5], sbp2.main$Values[3, 5], smk2.main$Values[3, 5], t2d2.main$Values[3, 5])))

## Combine them.
ivw.means <- c(ivw.means2[1], ivw.means1[1], NA, ivw.means2[2], ivw.means1[2], NA, ivw.means2[3], ivw.means1[3], NA, ivw.means2[4], ivw.means1[4], NA, ivw.means2[5], ivw.means1[5], NA)
ivw.ci1 <- c(ivw.ci12[1], ivw.ci11[1], NA, ivw.ci12[2], ivw.ci11[2], NA, ivw.ci12[3], ivw.ci11[3], NA, ivw.ci12[4], ivw.ci11[4], NA, ivw.ci12[5], ivw.ci11[5], NA)
ivw.ci2 <- c(ivw.ci22[1], ivw.ci21[1], NA, ivw.ci22[2], ivw.ci21[2], NA, ivw.ci22[3], ivw.ci21[3], NA, ivw.ci22[4], ivw.ci21[4], NA, ivw.ci22[5], ivw.ci21[5], NA)
ivw.traits <- rev(c("", "BMI (Biobank)", "BMI (HUNT)", "", "LDL-C (Biobank)", "LDL-C (HUNT)", "", "SBP (Biobank)", "SBP (HUNT)", "", "Smoking (Biobank)", "Smoking (HUNT)", "", "T2DM (Biobank)", "T2DM (HUNT)"))
ivw.colors <- rep(c("darkblue", "darkorange", NA), times = 5)


## Do the plot.
tiff(file = "Sepsis_IVW_NoPal.tiff", width = 1800, height = 1200, pointsize = 11, res = 300, compression = 'lzw')
par(mar = c(4, 7, 2, 10))

## Plot point estimates and confidence intervals.
plot(x = ivw.means, y = 1:15, type = "p", axes = FALSE, main = "", xlab = "Odds Ratio", ylab = "", cex.lab = 1, xlim = c(0.5, 4), ylim = c(0.5, 14), pch = 19, col = ivw.colors)
for (i in 1:15) {
  lines(c(ivw.ci1[i], ivw.ci2[i]), c(i, i), col = ivw.colors[i])
  lines(c(ivw.ci1[i], ivw.ci1[i]), c(i - 0.1, i + 0.1), col = ivw.colors[i])
  lines(c(ivw.ci2[i], ivw.ci2[i]), c(i - 0.1, i + 0.1), col = ivw.colors[i])
}

## Add axes and vertical/horizontal lines.
abline(h = c(3, 6, 9, 12), lty = 2, col = "black")
abline(v = c(1, 2, 3, 4), lty = 2, col = "grey")
abline(v = 1, lty = 1, col = "brown")
axis(side = 2, at = c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14), labels = ivw.traits[c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], las = 1, cex.axis = 0.75)
axis(side = 1, at = c(0.5, 1, 2, 3, 4), labels = c("", "1", "2", "3", "4"), las = 1, cex.axis = 1)

## Report point estimates and CIs (2 decimal places).
ivw.means.print <- sprintf("%.2f", round(ivw.means[c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], digits = 2))
ivw.ci1.print <- sprintf("%.2f", round(ivw.ci1[c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], digits = 2))
ivw.ci2.print <- sprintf("%.2f", round(ivw.ci2[c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], digits = 2))
ivw.ci.print <- paste("(", ivw.ci1.print, ", ", ivw.ci2.print, ")", sep = "")
est.print <- paste(ivw.means.print, ivw.ci.print, sep = "  ")
mtext(text = "Estimate", side = 2, line = -22, at = 15.5, las = 1, font = 1, cex = 0.9, col = "black")
mtext(text = est.print, side = 2, line = -24, at = c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14), las = 1, font = 1, cex = 0.9)

dev.off()

###############################################

##########   FIGURE 2   ##########

## Sensitivity analysis plotting results from 
## various pleiotropy-robust methods.

## We use a logarithmic scale for the x-axis.


## UK Biobank data.
all.means1 <- exp(rev(c(NA, bmi1.main$Values[c(3, 4, 2), 2], bmi1.mode2$Estimate, bmi1.conmix$Estimate, 
                        NA, ldl1.main$Values[c(3, 4, 2), 2], ldl1.mode2$Estimate, ldl1.conmix$Estimate, 
                        NA, sbp1.main$Values[c(3, 4, 2), 2], sbp1.mode2$Estimate, sbp1.conmix$Estimate, 
                        NA, smk1.main$Values[c(3, 4, 2), 2], smk1.mode2$Estimate, smk1.conmix$Estimate, 
                        NA, t2d1.main$Values[c(3, 4, 2), 2], t2d1.mode2$Estimate, t2d1.conmix$Estimate)))
all.ci11 <- exp(rev(c(NA, bmi1.main$Values[c(3, 4, 2), 4], bmi1.mode2$CILower, bmi1.conmix$CILower, 
                      NA, ldl1.main$Values[c(3, 4, 2), 4], ldl1.mode2$CILower, ldl1.conmix$CILower, 
                      NA, sbp1.main$Values[c(3, 4, 2), 4], sbp1.mode2$CILower, sbp1.conmix$CILower, 
                      NA, smk1.main$Values[c(3, 4, 2), 4], smk1.mode2$CILower, smk1.conmix$CILower, 
                      NA, t2d1.main$Values[c(3, 4, 2), 4], t2d1.mode2$CILower, t2d1.conmix$CILower)))
all.ci21 <- exp(rev(c(NA, bmi1.main$Values[c(3, 4, 2), 5], bmi1.mode2$CIUpper, bmi1.conmix$CIUpper, 
                      NA, ldl1.main$Values[c(3, 4, 2), 5], ldl1.mode2$CIUpper, ldl1.conmix$CIUpper, 
                      NA, sbp1.main$Values[c(3, 4, 2), 5], sbp1.mode2$CIUpper, sbp1.conmix$CIUpper, 
                      NA, smk1.main$Values[c(3, 4, 2), 5], smk1.mode2$CIUpper, smk1.conmix$CIUpper, 
                      NA, t2d1.main$Values[c(3, 4, 2), 5], t2d1.mode2$CIUpper, t2d1.conmix$CIUpper)))
all.traits <- rev(rep(c("", "IVW", "MR-Egger", "Median", "Mode", "ConMix"), 5))

## HUNT data.
all.means2 <- exp(rev(c(NA, bmi2.main$Values[c(3, 4, 2), 2], bmi2.mode2$Estimate, bmi2.conmix$Estimate, 
                        NA, ldl2.main$Values[c(3, 4, 2), 2], ldl2.mode2$Estimate, ldl2.conmix$Estimate, 
                        NA, sbp2.main$Values[c(3, 4, 2), 2], sbp2.mode2$Estimate, sbp2.conmix$Estimate, 
                        NA, smk2.main$Values[c(3, 4, 2), 2], smk2.mode2$Estimate, smk2.conmix$Estimate, 
                        NA, t2d2.main$Values[c(3, 4, 2), 2], t2d2.mode2$Estimate, t2d2.conmix$Estimate)))
all.ci12 <- exp(rev(c(NA, bmi2.main$Values[c(3, 4, 2), 4], bmi2.mode2$CILower, bmi2.conmix$CILower, 
                      NA, ldl2.main$Values[c(3, 4, 2), 4], ldl2.mode2$CILower, ldl2.conmix$CILower, 
                      NA, sbp2.main$Values[c(3, 4, 2), 4], sbp2.mode2$CILower, sbp2.conmix$CILower, 
                      NA, smk2.main$Values[c(3, 4, 2), 4], smk2.mode2$CILower, smk2.conmix$CILower, 
                      NA, t2d2.main$Values[c(3, 4, 2), 4], t2d2.mode2$CILower, t2d2.conmix$CILower)))
all.ci22 <- exp(rev(c(NA, bmi2.main$Values[c(3, 4, 2), 5], bmi2.mode2$CIUpper, bmi2.conmix$CIUpper, 
                      NA, ldl2.main$Values[c(3, 4, 2), 5], ldl2.mode2$CIUpper, ldl2.conmix$CIUpper, 
                      NA, sbp2.main$Values[c(3, 4, 2), 5], sbp2.mode2$CIUpper, sbp2.conmix$CIUpper, 
                      NA, smk2.main$Values[c(3, 4, 2), 5], smk2.mode2$CIUpper, smk2.conmix$CIUpper, 
                      NA, t2d2.main$Values[c(3, 4, 2), 5], t2d2.mode2$CIUpper, t2d2.conmix$CIUpper)))


## Start plotting.
tiff(file = "Sepsis_All_Methods_NoPal.tiff", width = 2800, height = 1600, pointsize = 10, res = 300, compression = 'lzw')

## Set overall plot parameters.
par(oma = c(0, 5, 0, 0))
par(mfrow = c(1, 2))

## Set inner plot margins for the UKBB plot.
par(mar = c(4, 2, 1, 8))

## Plot point estimates and confidence intervals.
plot(x = log2(all.means1), y = 1:30, type = "p", axes = FALSE, main = "", xlab = "UK Biobank", ylab = "", xlim = c(-1, 5), cex.lab = 1, pch = 19)
for (i in 1:30) {
  lines(c(log2(all.ci11[i]), log2(all.ci21[i])), c(i, i))
  lines(c(log2(all.ci11[i]), log2(all.ci11[i])), c(i - 0.2, i + 0.2))
  lines(c(log2(all.ci21[i]), log2(all.ci21[i])), c(i - 0.2, i + 0.2))
}

## Add axes and vertical/horizontal lines.
abline(h = c(6, 12, 18, 24, 30), lty = 2, col = "black")
abline(v = c(-1, 0, 1, 2, 3, 4, 5), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")
axis(side = 2, at = 1:30, labels = all.traits, las = 1, cex.axis = 0.75)
axis(side = 1, at = c(-1, 0, 1, 2, 3, 4, 5), labels = c("0.5", "1", "2", "4", "8", "16", "32"), las = 1, cex.axis = 1)

## Add trait labels.
mtext(text = c("T2DM", "Smoking", "SBP", "LDL-C", "BMI"), side = 2, line = 3, at = c(1:5 * 6), las = 1, font = 2, cex = 1, col = "brown")

## Report point estimates and CIs (2 decimal places).
all.means1.print <- sprintf("%.2f", round(all.means1[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci11.print <- sprintf("%.2f", round(all.ci11[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci21.print <- sprintf("%.2f", round(all.ci21[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci1.print <- paste("(", all.ci11.print, ", ", all.ci21.print, ")", sep = "")
est1.print <- paste(all.means1.print, all.ci1.print, sep = "  ")
mtext(text = "Estimate", side = 2, line = -21, at = 31, las = 1, font = 1, cex = 1, col = "black")
mtext(text = est1.print, side = 2, line = -23, at = c(1:5, 7:11, 13:17, 19:23, 25:29), las = 1, font = 1, cex = 0.9)

## Set inner plot margins for the HUNT plot.
par(mar = c(4, 2, 1, 8))

## Plot point estimates and confidence intervals.
plot(x = log(all.means2, base = 4), y = 1:30, type = "p", axes = FALSE, main = "", xlab = "HUNT", ylab = "", xlim = c(-1, 4), cex.lab = 1, pch = 19)
for (i in 1:30) {
  lines(c(log(all.ci12[i], base = 4), log(all.ci22[i], base = 4)), c(i, i))
  lines(c(log(all.ci12[i], base = 4), log(all.ci12[i], base = 4)), c(i - 0.2, i + 0.2))
  lines(c(log(all.ci22[i], base = 4), log(all.ci22[i], base = 4)), c(i - 0.2, i + 0.2))
}

## Add axes and vertical/horizontal lines.
abline(h = c(6, 12, 18, 24, 30), lty = 2, col = "black")
abline(v = c(-1, 0, 1, 2, 3, 4), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")
#axis(side = 2, at = 1:30, labels = all.traits, las = 1, cex.axis = 0.75)
axis(side = 1, at = c(-1, 0, 1, 2, 3, 4), labels = c("0.25", "1", "4", "16", "64", "256"), las = 1, cex.axis = 1)

## Report point estimates and CIs (2 decimal places).
all.means2.print <- sprintf("%.2f", round(all.means2[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci12.print <- sprintf("%.2f", round(all.ci12[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci22.print <- sprintf("%.2f", round(all.ci22[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci2.print <- paste("(", all.ci12.print, ", ", all.ci22.print, ")", sep = "")
est2.print <- paste(all.means2.print, all.ci2.print, sep = "  ")
mtext(text = "Estimate", side = 2, line = -21, at = 31, las = 1, font = 1, cex = 1, col = "black")
mtext(text = est2.print, side = 2, line = -23, at = c(1:5, 7:11, 13:17, 19:23, 25:29), las = 1, font = 1, cex = 0.9)

dev.off()

###############################################

##########   FIGURE 3   ##########

## Multivariate MR analysis of lipid traits 
## and sepsis risk. This analysis was run by 
## Dipender Gill. Here we only plot the results.

## Here are the results from Dipender's analysis.
mvmr.means <- c(-0.039880, -0.089341, 0.000348, -0.03319, 0.10649, -0.06646)
mvmr.se <- c(0.058392, 0.050243, 0.036376, 0.11959, 0.10212, 0.0746)
mvmr.ci1 <- mvmr.means - 1.96 * mvmr.se
mvmr.ci2 <- mvmr.means + 1.96 * mvmr.se
mvmr.traits <- c("TG", "HDL-C", "LDL-C")

## Create the Figure.
tiff(file = "Sepsis_MVMR_Lipids_NoPal.tiff", width = 2400, height = 800, pointsize = 11, res = 300, compression = 'lzw')

## Set plotting parameters.
par(oma = c(0, 2, 2, 1))
par(mfrow = c(1, 2))

## Plot point estimates and confidence intervals for UKBB.
par(mar = c(4, 2, 1, 6))
plot(x = exp(mvmr.means[1:3]), y = 1:3, type = "p", axes = FALSE, main = "", xlab = "UK Biobank", ylab = "", ylim = c(0.5, 3), xlim = c(0.8, 1.2), pch = 19)
for (i in 1:3) {
  lines(c(exp(mvmr.ci1[i]), exp(mvmr.ci2[i])), c(i, i))
  lines(c(exp(mvmr.ci1[i]), exp(mvmr.ci1[i])), c(i - 0.1, i + 0.1))
  lines(c(exp(mvmr.ci2[i]), exp(mvmr.ci2[i])), c(i - 0.1, i + 0.1))
}

## Add axes and vertical/horizontal lines.
axis(side = 1, at = c(0.8, 0.9, 1, 1.1, 1.2), cex.axis = 0.9, las = 1)
abline(v = c(0.8, 0.9, 1, 1.1, 1.2), lty = 2, col = "grey")
abline(v = 1, lty = 1, col = "brown")
axis(side = 2, at = 1:3, labels = mvmr.traits, las = 1)

## Report point estimates and CIs (2 decimal places).
mvmr.means.print <- sprintf("%.2f", round(exp(mvmr.means), digits = 2))
mvmr.ci1.print <- sprintf("%.2f", round(exp(mvmr.ci1), digits = 2))
mvmr.ci2.print <- sprintf("%.2f", round(exp(mvmr.ci2), digits = 2))
mvmr.print1 <- paste(mvmr.means.print[1:3], "  (", mvmr.ci1.print[1:3], ", ", mvmr.ci2.print[1:3], ")", sep = "")
mvmr.print2 <- paste(mvmr.means.print[4:6], "  (", mvmr.ci1.print[4:6], ", ", mvmr.ci2.print[4:6], ")", sep = "")
mtext(text = "Estimate", side = 2, line = -17.5, at = 3.75, las = 1, font = 1, cex = 1, col = "black")
mtext(text = mvmr.print1, side = 2, line = -19, at = 1:3, las = 1, font = 1, cex = 0.9)


## Plot point estimates and confidence intervals for HUNT.
par(mar = c(4, 2, 1, 6))
plot(x = exp(mvmr.means[4:6]), y = 1:3, type = "p", axes = FALSE, main = "", xlab = "HUNT", ylab = "", ylim = c(0.5, 3), xlim = c(0.75, 1.5), pch = 19)
for (i in 1:3) {
  lines(c(exp(mvmr.ci1[i + 3]), exp(mvmr.ci2[i + 3])), c(i, i))
  lines(c(exp(mvmr.ci1[i + 3]), exp(mvmr.ci1[i + 3])), c(i - 0.1, i + 0.1))
  lines(c(exp(mvmr.ci2[i + 3]), exp(mvmr.ci2[i + 3])), c(i - 0.1, i + 0.1))
}

## Add axes and vertical/horizontal lines.
axis(side = 1, at = c(0.75, 1, 1.25, 1.5), cex.axis = 0.9, las = 1)
abline(v =  c(0.75, 1, 1.25, 1.5), lty = 2, col = "grey")
abline(v = 1, lty = 1, col = "brown")

## Report point estimates and CIs (2 decimal places).
mtext(text = "Estimate", side = 2, line = -17.5, at = 3.75, las = 1, font = 1, cex = 1, col = "black")
mtext(text = mvmr.print2, side = 2, line = -19, at = 1:3, las = 1, font = 1, cex = 0.9)

dev.off()


################################################

## Those were the days, my friend... ##

################################################
################################################
################################################


###########################################################
##########   MENDELIAN RANDOMIZATION FOR COVID   ##########
###########################################################

## This code runs an MR analysis between a range of 
## cardiometabolic traits and severe covid-19 infection 
## risk, using NEJM data for covid-19 infection.

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

## Load the NEJM data for severe covid-19 infection.
covid <- read.delim("Covid Data.txt", header = TRUE, sep = "\t")

## Store alleles as characters - helps with pre-processing.
covid[, 5] <- as.character(covid[, 5])
covid[, 6] <- as.character(covid[, 6])

############################################################

##########   COVID - PRE-PROCESS THE DATA   ##########

## Subset UKBB data, get G-Y associations per trait.
covid.bmi <- covid[as.character(covid$rs_id_dbSNP151_GRCh38p7) %in% as.character(bmi$SNP), ]
covid.ldl <- covid[as.character(covid$rs_id_dbSNP151_GRCh38p7) %in% as.character(ldl$SNP), ]
covid.sbp <- covid[as.character(covid$rs_id_dbSNP151_GRCh38p7) %in% as.character(sbp$SNP), ]
covid.smk <- covid[as.character(covid$rs_id_dbSNP151_GRCh38p7) %in% as.character(smoking$SNP), ]
covid.t2d <- covid[as.character(covid$rs_id_dbSNP151_GRCh38p7) %in% as.character(t2dm$SNP), ]

## Subset G-X associations to exclude any SNPs not included in UK Biobank.
bmi3 <- bmi[as.character(bmi$SNP) %in% as.character(covid.bmi$rs_id_dbSNP151_GRCh38p7), ]
ldl3 <- ldl[as.character(ldl$SNP) %in% as.character(covid.ldl$rs_id_dbSNP151_GRCh38p7), ]
sbp3 <- sbp[as.character(sbp$SNP) %in% as.character(covid.sbp$rs_id_dbSNP151_GRCh38p7), ]
smoking3 <- smoking[as.character(smoking$SNP) %in% as.character(covid.smk$rs_id_dbSNP151_GRCh38p7), ]
t2dm3 <- t2dm[as.character(t2dm$SNP) %in% as.character(covid.t2d$rs_id_dbSNP151_GRCh38p7), ]

## Order all datasets by increasing rsID.
bmi3 <- bmi3[order(as.character(bmi3$SNP)), ]; covid.bmi <- covid.bmi[order(as.character(covid.bmi$rs_id_dbSNP151_GRCh38p7)), ]
ldl3 <- ldl3[order(as.character(ldl3$SNP)), ]; covid.ldl <- covid.ldl[order(as.character(covid.ldl$rs_id_dbSNP151_GRCh38p7)), ]
sbp3 <- sbp3[order(as.character(sbp3$SNP)), ]; covid.sbp <- covid.sbp[order(as.character(covid.sbp$rs_id_dbSNP151_GRCh38p7)), ]
smoking3 <- smoking3[order(as.character(smoking3$SNP)), ]; covid.smk <- covid.smk[order(as.character(covid.smk$rs_id_dbSNP151_GRCh38p7)), ]
t2dm3 <- t2dm3[order(as.character(t2dm3$SNP)), ]; covid.t2d <- covid.t2d[order(as.character(covid.t2d$rs_id_dbSNP151_GRCh38p7)), ]



## Sanity check. Number of SNPs in each dataset is the same.
nrow(bmi3) == nrow(covid.bmi)
nrow(ldl3) == nrow(covid.ldl)
nrow(sbp3) == nrow(covid.sbp)
nrow(smoking3) == nrow(covid.smk)
nrow(t2dm3) == nrow(covid.t2d)



## Align the reference alleles across datasets.

## For BMI.
covid.bmi2 <- covid.bmi
for (i in 1:nrow(bmi3)) {
  if (as.character(covid.bmi$effect_allele[i]) != as.character(bmi3$bmi_ea[i])) {
    covid.bmi2$beta[i] <- - covid.bmi$beta[i]
    covid.bmi2$effect_allele[i] <- covid.bmi$other_allele[i]
    covid.bmi2$other_allele[i] <- covid.bmi$effect_allele[i]
  }
}
covid.bmi <- covid.bmi2
rm(covid.bmi2)

## For LDL-C.
covid.ldl2 <- covid.ldl
for (i in 1:nrow(ldl3)) {
  if (as.character(covid.ldl$effect_allele[i]) != as.character(ldl3$ldl_ea[i])) {
    covid.ldl2$beta[i] <- - covid.ldl$beta[i]
    covid.ldl2$effect_allele[i] <- covid.ldl$other_allele[i]
    covid.ldl2$other_allele[i] <- covid.ldl$effect_allele[i]
  }
}
covid.ldl <- covid.ldl2
rm(covid.ldl2)

## For SBP.
covid.sbp2 <- covid.sbp
for (i in 1:nrow(sbp3)) {
  if (as.character(covid.sbp$effect_allele[i]) != as.character(sbp3$sbp_ea[i])) {
    covid.sbp2$beta[i] <- - covid.sbp$beta[i]
    covid.sbp2$effect_allele[i] <- covid.sbp$other_allele[i]
    covid.sbp2$other_allele[i] <- covid.sbp$effect_allele[i]
  }
}
covid.sbp <- covid.sbp2
rm(covid.sbp2)

## For Smoking.
covid.smk2 <- covid.smk
for (i in 1:nrow(smoking3)) {
  if (as.character(covid.smk$effect_allele[i]) != as.character(smoking3$smoking_ea[i])) {
    covid.smk2$beta[i] <- - covid.smk$beta[i]
    covid.smk2$effect_allele[i] <- covid.smk$other_allele[i]
    covid.smk2$other_allele[i] <- covid.smk$effect_allele[i]
  }
}
covid.smk <- covid.smk2
rm(covid.smk2)

## For T2DM.
covid.t2d2 <- covid.t2d
for (i in 1:nrow(t2dm3)) {
  if (as.character(covid.t2d$effect_allele[i]) != as.character(t2dm3$t2dm_ea[i])) {
    covid.t2d2$beta[i] <- - covid.t2d$beta[i]
    covid.t2d2$effect_allele[i] <- covid.t2d$other_allele[i]
    covid.t2d2$other_allele[i] <- covid.t2d$effect_allele[i]
  }
}
covid.t2d <- covid.t2d2
rm(covid.t2d2)


## Sanity check. Confirms that  that data are aligned.
all(as.character(covid.bmi$SNP) == as.character(bmi3$SNP)); all(as.character(covid.bmi$ALLELE1) == as.character(bmi3$bmi_ea))
all(as.character(covid.ldl$SNP) == as.character(ldl3$SNP)); all(as.character(covid.ldl$ALLELE1) == as.character(ldl3$ldl_ea))
all(as.character(covid.sbp$SNP) == as.character(sbp3$SNP)); all(as.character(covid.sbp$ALLELE1) == as.character(sbp3$sbp_ea))
all(as.character(covid.smk$SNP) == as.character(smoking3$SNP)); all(as.character(covid.smk$ALLELE1) == as.character(smoking3$smoking_ea))
all(as.character(covid.t2d$SNP) == as.character(t2dm3$SNP)); all(as.character(covid.t2d$ALLELE1) == as.character(t2dm3$t2dm_ea))

###############################################

##########   COVID - MR ANALYSIS   ##########

## We use functions from the R package MendelianRandomization.
## We use IVW, MR-Egger, mode, median and the contamination mixture.

## Note: mr_allmethods performs IVW using a random-effects model
## by default when using more than 3 SNPs, which is what we want.

## BMI.
bmi3.main <- mr_allmethods( mr_input(bx = bmi3$bmi_beta, bxse = bmi3$bmi_se, by = covid.bmi$beta, byse = covid.bmi$standard_error), method = "main" )
bmi3.mode1 <- mr_mbe( mr_input(bx = bmi3$bmi_beta, bxse = bmi3$bmi_se, by = covid.bmi$beta, byse = covid.bmi$standard_error), weighting = "unweighted" )
bmi3.mode2 <- mr_mbe( mr_input(bx = bmi3$bmi_beta, bxse = bmi3$bmi_se, by = covid.bmi$beta, byse = covid.bmi$standard_error), weighting = "weighted" )
bmi3.conmix <- mr_conmix( mr_input(bx = bmi3$bmi_beta, bxse = bmi3$bmi_se, by = covid.bmi$beta, byse = covid.bmi$standard_error), CIStep = 0.001 )

## LDL-C.
ldl3.main <- mr_allmethods( mr_input(bx = ldl3$ldl_beta, bxse = ldl3$ldl_se, by = covid.ldl$beta, byse = covid.ldl$standard_error), method = "main" )
ldl3.mode1 <- mr_mbe( mr_input(bx = ldl3$ldl_beta, bxse = ldl3$ldl_se, by = covid.ldl$beta, byse = covid.ldl$standard_error), weighting = "unweighted" )
ldl3.mode2 <- mr_mbe( mr_input(bx = ldl3$ldl_beta, bxse = ldl3$ldl_se, by = covid.ldl$beta, byse = covid.ldl$standard_error), weighting = "weighted" )
ldl3.conmix <- mr_conmix( mr_input(bx = ldl3$ldl_beta, bxse = ldl3$ldl_se, by = covid.ldl$beta, byse = covid.ldl$standard_error), CIStep = 0.001 )

## SBP.
sbp3.main <- mr_allmethods( mr_input(bx = sbp3$sbp_beta, bxse = sbp3$sbp_se, by = covid.sbp$beta, byse = covid.sbp$standard_error), method = "main" )
sbp3.mode1 <- mr_mbe( mr_input(bx = sbp3$sbp_beta, bxse = sbp3$sbp_se, by = covid.sbp$beta, byse = covid.sbp$standard_error), weighting = "unweighted" )
sbp3.mode2 <- mr_mbe( mr_input(bx = sbp3$sbp_beta, bxse = sbp3$sbp_se, by = covid.sbp$beta, byse = covid.sbp$standard_error), weighting = "weighted" )
sbp3.conmix <- mr_conmix( mr_input(bx = sbp3$sbp_beta, bxse = sbp3$sbp_se, by = covid.sbp$beta, byse = covid.sbp$standard_error), CIStep = 0.001 )

## Smoking.
smk3.main <- mr_allmethods( mr_input(bx = smoking3$smoking_beta, bxse = smoking3$smoking_se, by = covid.smk$beta, byse = covid.smk$standard_error), method = "main" )
smk3.mode1 <- mr_mbe( mr_input(bx = smoking3$smoking_beta, bxse = smoking3$smoking_se, by = covid.smk$beta, byse = covid.smk$standard_error), weighting = "unweighted" )
smk3.mode2 <- mr_mbe( mr_input(bx = smoking3$smoking_beta, bxse = smoking3$smoking_se, by = covid.smk$beta, byse = covid.smk$standard_error), weighting = "weighted" )
smk3.conmix <- mr_conmix( mr_input(bx = smoking3$smoking_beta, bxse = smoking3$smoking_se, by = covid.smk$beta, byse = covid.smk$standard_error), CIStep = 0.001 )

## T2DM.
t2d3.main <- mr_allmethods( mr_input(bx = t2dm3$t2dm_beta, bxse = t2dm3$t2dm_se, by = covid.t2d$beta, byse = covid.t2d$standard_error), method = "main" )
t2d3.mode1 <- mr_mbe( mr_input(bx = t2dm3$t2dm_beta, bxse = t2dm3$t2dm_se, by = covid.t2d$beta, byse = covid.t2d$standard_error), weighting = "unweighted" )
t2d3.mode2 <- mr_mbe( mr_input(bx = t2dm3$t2dm_beta, bxse = t2dm3$t2dm_se, by = covid.t2d$beta, byse = covid.t2d$standard_error), weighting = "weighted" )
t2d3.conmix <- mr_conmix( mr_input(bx = t2dm3$t2dm_beta, bxse = t2dm3$t2dm_se, by = covid.t2d$beta, byse = covid.t2d$standard_error), CIStep = 0.001 )

###############################################

##########   COVID - RESULTS TABLES   ##########

## This explicitly reports results per method per trait.
## We report both MR estimates and (exponentiated) 
## results in an odds-ratio scale.


## Label the methods.
method.names <- c(bmi3.main$Values[, 1], "Simple mode", "Weighted mode", "Con mix")

## BMI.
bmi3.est <- round(c(bmi3.main$Values[, 2], bmi3.mode1$Estimate, bmi3.mode2$Estimate, bmi3.conmix$Estimate), digits = 3)
bmi3.se <- round(c(bmi3.main$Values[, 3], bmi3.mode1$StdError, bmi3.mode2$StdError, NA), digits = 3)
bmi3.ci1 <- round(c(bmi3.main$Values[, 4], bmi3.mode1$CILower, bmi3.mode2$CILower, bmi3.conmix$CILower), digits = 3)
bmi3.ci2 <- round(c(bmi3.main$Values[, 5], bmi3.mode1$CIUpper, bmi3.mode2$CIUpper, bmi3.conmix$CIUpper), digits = 3)
bmi3.p <- round(c(bmi3.main$Values[, 6], bmi3.mode1$Pvalue, bmi3.mode2$Pvalue, bmi3.conmix$Pvalue), digits = 3)
bmi3.results <- data.frame(method.names, bmi3.est, bmi3.se, bmi3.ci1, bmi3.ci2, bmi3.p)
bmi3.results <- cbind(bmi3.results, round(exp(bmi3.results[, c(2, 4, 5)]), digits = 3))
bmi3.results[5, 7:9] <- NA
names(bmi3.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
bmi3.results

## LDL-C.
ldl3.est <- round(c(ldl3.main$Values[, 2], ldl3.mode1$Estimate, ldl3.mode2$Estimate, ldl3.conmix$Estimate), digits = 3)
ldl3.se <- round(c(ldl3.main$Values[, 3], ldl3.mode1$StdError, ldl3.mode2$StdError, NA), digits = 3)
ldl3.ci1 <- round(c(ldl3.main$Values[, 4], ldl3.mode1$CILower, ldl3.mode2$CILower, ldl3.conmix$CILower), digits = 3)
ldl3.ci2 <- round(c(ldl3.main$Values[, 5], ldl3.mode1$CIUpper, ldl3.mode2$CIUpper, ldl3.conmix$CIUpper), digits = 3)
ldl3.p <- round(c(ldl3.main$Values[, 6], ldl3.mode1$Pvalue, ldl3.mode2$Pvalue, ldl3.conmix$Pvalue), digits = 3)
ldl3.results <- data.frame(method.names, ldl3.est, ldl3.se, ldl3.ci1, ldl3.ci2, ldl3.p)
ldl3.results <- cbind(ldl3.results, round(exp(ldl3.results[, c(2, 4, 5)]), digits = 3))
ldl3.results[5, 7:9] <- NA
names(ldl3.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
ldl3.results

## SBP.
sbp3.est <- round(c(sbp3.main$Values[, 2], sbp3.mode1$Estimate, sbp3.mode2$Estimate, sbp3.conmix$Estimate), digits = 3)
sbp3.se <- round(c(sbp3.main$Values[, 3], sbp3.mode1$StdError, sbp3.mode2$StdError, NA), digits = 3)
sbp3.ci1 <- round(c(sbp3.main$Values[, 4], sbp3.mode1$CILower, sbp3.mode2$CILower, sbp3.conmix$CILower), digits = 3)
sbp3.ci2 <- round(c(sbp3.main$Values[, 5], sbp3.mode1$CIUpper, sbp3.mode2$CIUpper, sbp3.conmix$CIUpper), digits = 3)
sbp3.p <- round(c(sbp3.main$Values[, 6], sbp3.mode1$Pvalue, sbp3.mode2$Pvalue, sbp3.conmix$Pvalue), digits = 3)
sbp3.results <- data.frame(method.names, sbp3.est, sbp3.se, sbp3.ci1, sbp3.ci2, sbp3.p)
sbp3.results <- cbind(sbp3.results, round(exp(sbp3.results[, c(2, 4, 5)]), digits = 3))
sbp3.results[5, 7:9] <- NA
names(sbp3.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
sbp3.results

## Smoking.
smk3.est <- round(c(smk3.main$Values[, 2], smk3.mode1$Estimate, smk3.mode2$Estimate, smk3.conmix$Estimate), digits = 3)
smk3.se <- round(c(smk3.main$Values[, 3], smk3.mode1$StdError, smk3.mode2$StdError, NA), digits = 3)
smk3.ci1 <- round(c(smk3.main$Values[, 4], smk3.mode1$CILower, smk3.mode2$CILower, smk3.conmix$CILower), digits = 3)
smk3.ci2 <- round(c(smk3.main$Values[, 5], smk3.mode1$CIUpper, smk3.mode2$CIUpper, smk3.conmix$CIUpper), digits = 3)
smk3.p <- round(c(smk3.main$Values[, 6], smk3.mode1$Pvalue, smk3.mode2$Pvalue, smk3.conmix$Pvalue), digits = 3)
smk3.results <- data.frame(method.names, smk3.est, smk3.se, smk3.ci1, smk3.ci2, smk3.p)
smk3.results <- cbind(smk3.results, round(exp(smk3.results[, c(2, 4, 5)]), digits = 3))
smk3.results[5, 7:9] <- NA
names(smk3.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
smk3.results

## T2DM.
t2d3.est <- round(c(t2d3.main$Values[, 2], t2d3.mode1$Estimate, t2d3.mode2$Estimate, t2d3.conmix$Estimate), digits = 3)
t2d3.se <- round(c(t2d3.main$Values[, 3], t2d3.mode1$StdError, t2d3.mode2$StdError, NA), digits = 3)
t2d3.ci1 <- round(c(t2d3.main$Values[, 4], t2d3.mode1$CILower, t2d3.mode2$CILower, t2d3.conmix$CILower), digits = 3)
t2d3.ci2 <- round(c(t2d3.main$Values[, 5], t2d3.mode1$CIUpper, t2d3.mode2$CIUpper, t2d3.conmix$CIUpper), digits = 3)
t2d3.p <- round(c(t2d3.main$Values[, 6], t2d3.mode1$Pvalue, t2d3.mode2$Pvalue, t2d3.conmix$Pvalue), digits = 3)
t2d3.results <- data.frame(method.names, t2d3.est, t2d3.se, t2d3.ci1, t2d3.ci2, t2d3.p)
t2d3.results <- cbind(t2d3.results, round(exp(t2d3.results[, c(2, 4, 5)]), digits = 3))
t2d3.results[5, 7:9] <- NA
names(t2d3.results) <- c("Method", "Estimate", "Std Error", "95% CI", "", "P-value", "Odds Ratio", "95% CI", "")
t2d3.results

###############################################

## Now create figures.

## We use the ".tiff" format as required by the 
## journal (and "lzw" compression to save space).
## We plot point estimates and 95\% CI per method 
## and report numerical values next to each plot.

###############################################

##########   FIGURE 1   ##########

## Plot random-effects IVW estimates for UKBB and HUNT.

## We use a logarithmic scale for the x-axis.

## UK Biobank data.
ivw.means3 <- exp(rev(c(bmi3.main$Values[3, 2], ldl3.main$Values[3, 2], sbp3.main$Values[3, 2], smk3.main$Values[3, 2], t2d3.main$Values[3, 2])))
ivw.ci13 <- exp(rev(c(bmi3.main$Values[3, 4], ldl3.main$Values[3, 4], sbp3.main$Values[3, 4], smk3.main$Values[3, 4], t2d3.main$Values[3, 4])))
ivw.ci23 <- exp(rev(c(bmi3.main$Values[3, 5], ldl3.main$Values[3, 5], sbp3.main$Values[3, 5], smk3.main$Values[3, 5], t2d3.main$Values[3, 5])))
ivw.traits3 <- rev(c("BMI", "LDL-C", "SBP", "Smoking", "T2DM"))
ivw.colors3 <- rep("darkblue", 5)

## Do the plot.
tiff(file = "Covid_IVW.tiff", width = 1500, height = 1000, pointsize = 11, res = 300, compression = 'lzw')
par(mar = c(4, 3.5, 3, 7))

## Plot point estimates and confidence intervals.
plot(x = log2(ivw.means3), y = 1:5, type = "p", axes = FALSE, main = "", xlab = "Odds Ratio", ylab = "", cex.lab = 1, xlim = c(-1, 4), ylim = c(0.5, 5), pch = 19, col = ivw.colors3)
for (i in 1:5) {
  lines(c(log2(ivw.ci13[i]), log2(ivw.ci23[i])), c(i, i), col = ivw.colors3[i])
  lines(c(log2(ivw.ci13[i]), log2(ivw.ci13[i])), c(i - 0.1, i + 0.1), col = ivw.colors3[i])
  lines(c(log2(ivw.ci23[i]), log2(ivw.ci23[i])), c(i - 0.1, i + 0.1), col = ivw.colors3[i])
}

## Add axes and vertical/horizontal lines.
abline(v = c(-1, 0, 1, 2, 3, 4), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")
axis(side = 2, at = c(1:5), labels = ivw.traits3, las = 1, cex.axis = 0.75)
axis(side = 1, at = c(-1, 0, 1, 2, 3, 4), labels = c("0.5", "1", "2", "4", "8", "16"), las = 1, cex.axis = 1)

## Report point estimates and CIs (2 decimal places).
ivw.means3.print <- sprintf("%.2f", round(ivw.means3, digits = 2))
ivw.ci13.print <- sprintf("%.2f", round(ivw.ci13, digits = 2))
ivw.ci23.print <- sprintf("%.2f", round(ivw.ci23, digits = 2))
ivw.ci3.print <- paste("(", ivw.ci13.print, ", ", ivw.ci23.print, ")", sep = "")
est3.print <- paste(ivw.means3.print, ivw.ci3.print, sep = "  ")
mtext(text = "Estimate", side = 2, line = -22, at = 5.75, las = 1, font = 1, cex = 0.9, col = "black")
mtext(text = est3.print, side = 2, line = -24, at = 1:5, las = 1, font = 1, cex = 0.9)

dev.off()

###############################################

##########   FIGURE 2   ##########

## Sensitivity analysis plotting results from 
## various pleiotropy-robust methods.

## We use a logarithmic scale for the x-axis.


## Collect the data.
all.means3 <- exp(rev(c(NA, bmi3.main$Values[c(3, 4, 2), 2], bmi3.mode2$Estimate, bmi3.conmix$Estimate, 
                        NA, ldl3.main$Values[c(3, 4, 2), 2], ldl3.mode2$Estimate, ldl3.conmix$Estimate, 
                        NA, sbp3.main$Values[c(3, 4, 2), 2], sbp3.mode2$Estimate, sbp3.conmix$Estimate, 
                        NA, smk3.main$Values[c(3, 4, 2), 2], smk3.mode2$Estimate, smk3.conmix$Estimate, 
                        NA, t2d3.main$Values[c(3, 4, 2), 2], t2d3.mode2$Estimate, t2d3.conmix$Estimate)))
all.ci13 <- exp(rev(c(NA, bmi3.main$Values[c(3, 4, 2), 4], bmi3.mode2$CILower, bmi3.conmix$CILower, 
                      NA, ldl3.main$Values[c(3, 4, 2), 4], ldl3.mode2$CILower, ldl3.conmix$CILower, 
                      NA, sbp3.main$Values[c(3, 4, 2), 4], sbp3.mode2$CILower, sbp3.conmix$CILower, 
                      NA, smk3.main$Values[c(3, 4, 2), 4], smk3.mode2$CILower, smk3.conmix$CILower, 
                      NA, t2d3.main$Values[c(3, 4, 2), 4], t2d3.mode2$CILower, t2d3.conmix$CILower)))
all.ci23 <- exp(rev(c(NA, bmi3.main$Values[c(3, 4, 2), 5], bmi3.mode2$CIUpper, bmi3.conmix$CIUpper, 
                      NA, ldl3.main$Values[c(3, 4, 2), 5], ldl3.mode2$CIUpper, ldl3.conmix$CIUpper, 
                      NA, sbp3.main$Values[c(3, 4, 2), 5], sbp3.mode2$CIUpper, sbp3.conmix$CIUpper, 
                      NA, smk3.main$Values[c(3, 4, 2), 5], smk3.mode2$CIUpper, smk3.conmix$CIUpper, 
                      NA, t2d3.main$Values[c(3, 4, 2), 5], t2d3.mode2$CIUpper, t2d3.conmix$CIUpper)))
all.traits3 <- rev(rep(c("", "IVW", "MR-Egger", "Median", "Mode", "ConMix"), 5))

## Start plotting.
tiff(file = "Covid_All_Methods.tiff", width = 2000, height = 1600, pointsize = 10, res = 300, compression = 'lzw')

## Set overall plot parameters.
par(oma = c(0, 5, 0, 0))
par(mar = c(4, 2, 1, 8))

## Plot point estimates and confidence intervals.
plot(x = log(all.means3, base = 4), y = 1:30, type = "p", axes = FALSE, main = "", xlab = "Odds Ratio", ylab = "", xlim = c(-4, 4), cex.lab = 1, pch = 19)
for (i in 1:30) {
  lines(c(log(all.ci13[i], base = 4), log(all.ci23[i], base = 4)), c(i, i))
  lines(c(log(all.ci13[i], base = 4), log(all.ci13[i], base = 4)), c(i - 0.2, i + 0.2))
  lines(c(log(all.ci23[i], base = 4), log(all.ci23[i], base = 4)), c(i - 0.2, i + 0.2))
}

## Add axes and vertical/horizontal lines.
abline(h = c(6, 12, 18, 24, 30), lty = 2, col = "black")
abline(v = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), lty = 2, col = "grey")
abline(v = 0, lty = 1, col = "brown")
axis(side = 2, at = 1:30, labels = all.traits3, las = 1, cex.axis = 0.75)
axis(side = 1, at = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), labels = c("1/256", "1/64", "1/16", "1/4", "1", "4", "16", "64", "256"), las = 1, cex.axis = 1)

## Add trait labels.
mtext(text = c("T2DM", "Smoking", "SBP", "LDL-C", "BMI"), side = 2, line = 3, at = c(1:5 * 6), las = 1, font = 2, cex = 1, col = "brown")

## Report point estimates and CIs (2 decimal places).
all.means3.print <- sprintf("%.2f", round(all.means3[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci13.print <- sprintf("%.2f", round(all.ci13[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci23.print <- sprintf("%.2f", round(all.ci23[- c(6, 12, 18, 24, 30)], digits = 2))
all.ci3.print <- paste("(", all.ci13.print, ", ", all.ci23.print, ")", sep = "")
est3.print <- paste(all.means3.print, all.ci3.print, sep = "  ")
mtext(text = "Estimate", side = 2, line = -31, at = 31, las = 1, font = 1, cex = 1, col = "black")
mtext(text = est3.print, side = 2, line = -33, at = c(1:5, 7:11, 13:17, 19:23, 25:29), las = 1, font = 1, cex = 0.9)

dev.off()

###############################################

##########   FIGURE 3   ##########

## Multivariate MR analysis of lipid traits 
## and covid-19 risk. This analysis was run by 
## Dipender Gill. Here we only plot the results.

## Here are the results from Dipender's analysis.
mvmr.means3 <- c(0.1126, -0.1188, -0.2293)
mvmr.se3 <- c(0.2587, 0.2213, 0.1689)
mvmr.ci13 <- mvmr.means3 - 1.96 * mvmr.se3
mvmr.ci23 <- mvmr.means3 + 1.96 * mvmr.se3
mvmr.traits3 <- c("TG", "HDL-C", "LDL-C")

## Create the Figure.
tiff(file = "Covid_MVMR_Lipids.tiff", width = 1500, height = 800, pointsize = 11, res = 300, compression = 'lzw')

## Set plotting parameters.
par(oma = c(0, 2, 2, 1))
par(mar = c(4, 2, 1, 7))

## Plot point estimates and confidence intervals for UKBB.
plot(x = exp(mvmr.means3), y = 1:3, type = "p", axes = FALSE, main = "", xlab = "Odds Ratio", ylab = "", ylim = c(0.5, 3), xlim = c(0.5, 2), pch = 19)
for (i in 1:3) {
  lines(c(exp(mvmr.ci13[i]), exp(mvmr.ci23[i])), c(i, i))
  lines(c(exp(mvmr.ci13[i]), exp(mvmr.ci13[i])), c(i - 0.1, i + 0.1))
  lines(c(exp(mvmr.ci23[i]), exp(mvmr.ci23[i])), c(i - 0.1, i + 0.1))
}

## Add axes and vertical/horizontal lines.
axis(side = 1, at = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2), cex.axis = 0.85, las = 1)
abline(v = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2), lty = 2, col = "grey")
abline(v = 1, lty = 1, col = "brown")
axis(side = 2, at = 1:3, labels = mvmr.traits3, las = 1)

## Report point estimates and CIs (2 decimal places).
mvmr.means3.print <- sprintf("%.2f", round(exp(mvmr.means3), digits = 2))
mvmr.ci13.print <- sprintf("%.2f", round(exp(mvmr.ci13), digits = 2))
mvmr.ci23.print <- sprintf("%.2f", round(exp(mvmr.ci23), digits = 2))
mvmr.print3 <- paste(mvmr.means3.print, "  (", mvmr.ci13.print, ", ", mvmr.ci23.print, ")", sep = "")
mtext(text = "Estimate", side = 2, line = -20, at = 3.75, las = 1, font = 1, cex = 1, col = "black")
mtext(text = mvmr.print3, side = 2, line = -22, at = 1:3, las = 1, font = 1, cex = 0.9)

dev.off()


################################################


################################################
################################################
################################################

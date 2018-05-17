###################################################
### THESIS. CHAPTER 4. ############################
### Schizophrenia Example Part 1. #################
### Age and sex NOT included in risk estimation ###
###################################################
### Description: ##################################
### In part 1, included variables for risk 
# estimation are: 10 CNVs, cannabis usage, paternal
# age at conception, migrant status, urbanicity, 
# childhood trauma, a PRS and the disease status of
# a parent (family history of disease).
### We show how to parameterise the LTMM.
### Risk estimation LTMMM methods A, B, and C 
# (see main thesis chapter for details of these) 
# are each used to calculate risk to show that the 
# approximate methods (B and C) provide a good 
# approximation to method A.
### The LL method is also presented. 
###################################################
### Libraries needed and source other scripts #####
###################################################
### INSTALL REQUIRED LIBRARIES:
library(BB)
library(mvtnorm)
### IMPORTANT: TO USE THIS SCRIPT, DOWNLOAD 
# R-SCRIPTS: riskLTMM_methodA.R, riskLTMM_methodB.R,
# riskLTMM_methodC.R and riskLL.R
# available from: 
# https://github.com/alexgillett/THESIS/tree/master/Chapter3
# THEN USE SOURCE FUNCTION (WITH FILE PATH)
source(file="/Users/alexg/Documents/Online_appendices_work/Chapter3/Rcode/riskLTMM_methodA.R")
source(file="/Users/alexg/Documents/Online_appendices_work/Chapter3/Rcode/riskLTMM_methodB.R")
source(file="/Users/alexg/Documents/Online_appendices_work/Chapter3/Rcode/riskLTMM_methodC.R")
source(file="/Users/alexg/Documents/Online_appendices_work/Chapter3/Rcode/riskLL.R")
###################################################
### LTMM risk estimation ##########################
###################################################
### STEP 1: LTMM parameterisation #################
###################################################
### STEP 1A = effect size estimates on L scale
### STEP 1B = calculate E[L]
### STEP 1C = calculate disease threshold
### STEP 1D = calculate variance contribution for:
# a. measured major genetic risk loci and b. the 
# measured, discrete environmental risk factors
### STEP 1E = define/ calculate other variables
# required for using LTMM risk estimation functions
###################################################
### STEP 2: Compare methods A, B and C ############
###################################################
### STEP 3: Risk estimation for a set of risk #####
### variables #####################################
###################################################
### STEP 1A: Effect size estimates on L scale #####
###################################################
### General known parameters:
K <- 0.01 ### prevalence
H2 <- 0.81 ### Broad-sense heritability
h2 <- 0.81 ### Narrow-sense heritability
VML <- 0.07 ### Var[PRS] on L-scale
VD <- H2-h2 ### (Quasi) Dominant variance component
###################################################
### CNVs:
###################################################
### G1: $22q11.21$ del*
### p(carrier)/ p(non-carrier) known
pG1_0 <- 0.9999207151
pG1_1 <- 1 - 0.9999207151
### use p(carrier)/ p(non-carrier) to calculate 
# `risk allele frequency': f1
f1 <- 1 - sqrt(pG1_0)
### Odds ratio (OR) known:
OR_G1nc <- 1
OR_G1c <- 67.7
### OR used, with K, to calculate penetrance function:
# pD_G1nc/ pD_G1c (only results presented, but see 
# thesis for details)
pD_G1nc <- 0.009968652
pD_G1c <- 0.405354107
### Re-capture beta0 and beta1 on the logit scale 
# using penetrance function
b0logitG1 <- log(pD_G1nc/(1-pD_G1nc))
b1logitG1 <- log(pD_G1c/(1-pD_G1c)) - b0logitG1
### Use beta0 and beta1 on the logit scale to calculate
# beta0 and beta1 on the probit scale
b0probitG1 <- qnorm(plogis(b0logitG1))
b1probitG1 <- qnorm(plogis(b0logitG1 + b1logitG1)) - b0probitG1
### tau == effect size estimates on the liability scale
### set tau0 equal to the reference category
tau0G1 <- 0
### Calculate tau1 (CNV effect size estimate)
tau1G1 <- b1probitG1/sqrt(1 + ((b1probitG1^2)*pG1_1*(1-pG1_1)))
### Check using raf instead
###b1probitG1/sqrt(1 + ((b1probitG1^2)*2*f1*(1-f1)))
###################################################
### G2: $16p11.2$, proximal, dup
pG2_0 <- 0.9996242038
pG2_1 <- 1 - 0.9996242038
f2 <- 1 - sqrt(pG2_0)
OR_G2nc <- 1
OR_G2c <- 9.4
pD_G2nc <- 0.009971246
pD_G2c <- 0.086485796
b0logitG2 <- log(pD_G2nc/(1-pD_G2nc))
b1logitG2 <- log(pD_G2c/(1-pD_G2c)) - b0logitG2
b0probitG2 <- qnorm(plogis(b0logitG2))
b1probitG2 <- qnorm(plogis(b0logitG2 + b1logitG2)) - b0probitG2
tau0G2 <- 0
tau1G2 <- b1probitG2/sqrt(1 + ((b1probitG2^2)*pG2_1*(1-pG2_1)))
###################################################
### G3: $2p16.3$ (NRXN1) del*
pG3_0 <- 0.9998365742
pG3_1 <- 1 - 0.9998365742
f3 <- 1 - sqrt(pG3_0)
OR_G3nc <- 1
OR_G3c <- 14.4
pD_G3nc <- 0.009980914
pD_G3c <- 0.126770351
b0logitG3 <- log(pD_G3nc/(1-pD_G3nc))
b1logitG3 <- log(pD_G3c/(1-pD_G3c)) - b0logitG3
b0probitG3 <- qnorm(plogis(b0logitG3))
b1probitG3 <- qnorm(plogis(b0logitG3 + b1logitG3)) - b0probitG3
tau0G3 <- 0
tau1G3 <- b1probitG3/sqrt(1 + ((b1probitG3^2)*pG3_1*(1-pG3_1)))
###################################################
### G4: $15q13.3$ del 
pG4_0 <- 0.9998888371
pG4_1 <- 1 - 0.9998888371
f4 <- 1 - sqrt(pG4_0)
OR_G4nc <- 1
OR_G4c <-15.6
pD_G4nc <- 0.009985996
pD_G4c <- 0.135959286
b0logitG4 <- log(pD_G4nc/(1-pD_G4nc))
b1logitG4 <- log(pD_G4c/(1-pD_G4c)) - b0logitG4
b0probitG4 <- qnorm(plogis(b0logitG4))
b1probitG4 <- qnorm(plogis(b0logitG4 + b1logitG4)) - b0probitG4
tau0G4 <- 0
tau1G4 <- b1probitG4/sqrt(1 + ((b1probitG4^2)*pG4_1*(1-pG4_1)))
###################################################
### G5: $1q21.1$ del$+$dup
pG5_0 <- 0.9992863332
pG5_1 <- 1 - 0.9992863332
f5 <- 1 - sqrt(pG5_0)
OR_G5nc <- 1
OR_G5c <- 3.8
pD_G5nc <- 0.009980792
pD_G5c <- 0.036895907
b0logitG5 <- log(pD_G5nc/(1-pD_G5nc))
b1logitG5 <- log(pD_G5c/(1-pD_G5c)) - b0logitG5
b0probitG5 <- qnorm(plogis(b0logitG5))
b1probitG5 <- qnorm(plogis(b0logitG5 + b1logitG5)) - b0probitG5
tau0G5 <- 0
tau1G5 <- b1probitG5/sqrt(1 + ((b1probitG5^2)*pG5_1*(1-pG5_1)))
###################################################
### G6: $16p11.2$, distal, del
pG6_0 <- 9.999458e-01
pG6_1 <- 1 - 9.999458e-01
f6 <- 1 - sqrt(pG6_0)
OR_G6nc <- 1
OR_G6c <- 20.6
pD_G6nc <- 0.00999122
pD_G6c <- 0.17211433
b0logitG6 <- log(pD_G6nc/(1-pD_G6nc))
b1logitG6 <- log(pD_G6c/(1-pD_G6c)) - b0logitG6
b0probitG6 <- qnorm(plogis(b0logitG6))
b1probitG6 <- qnorm(plogis(b0logitG6 + b1logitG6)) - b0probitG6
tau0G6 <- 0
tau1G6 <- b1probitG6/sqrt(1 + ((b1probitG6^2)*pG6_1*(1-pG6_1)))
###################################################
### G7: $7q11.23$ dup
pG7_0 <- 9.999435e-01
pG7_1 <- 1 - 9.999435e-01
f7 <- 1 - sqrt(pG7_0)
OR_G7nc <- 1
OR_G7c <- 16.1
pD_G7nc <- 0.009992663
pD_G7c <- 0.139789189
b0logitG7 <- log(pD_G7nc/(1-pD_G7nc))
b1logitG7 <- log(pD_G7c/(1-pD_G7c)) - b0logitG7
b0probitG7 <- qnorm(plogis(b0logitG7))
b1probitG7 <- qnorm(plogis(b0logitG7 + b1logitG7)) - b0probitG7
tau0G7 <- 0
tau1G7 <- b1probitG7/sqrt(1 + ((b1probitG7^2)*pG7_1*(1-pG7_1)))
###################################################
### G8: $22q11.21$ dup
pG8_0 <- 0.9992154661
pG8_1 <- 1 - 0.9992154661
f8 <- 1 - sqrt(pG8_0)
OR_G8nc <- 1
OR_G8c <- 0.15
pD_G8nc <- 0.010006663
pD_G8c <- 0.001513876
b0logitG8 <- log(pD_G8nc/(1-pD_G8nc))
b1logitG8 <- log(pD_G8c/(1-pD_G8c)) - b0logitG8
b0probitG8 <- qnorm(plogis(b0logitG8))
b1probitG8 <- qnorm(plogis(b0logitG8 + b1logitG8)) - b0probitG8
tau0G8 <- 0
tau1G8 <- b1probitG8/sqrt(1 + ((b1probitG8^2)*pG8_1*(1-pG8_1)))
###################################################
### G9: $15q11.2$ del*
pG9_0 <- 0.997506317
pG9_1 <- 1 - 0.997506317
f9 <- 1 - sqrt(pG9_0)
OR_G9nc <- 1
OR_G9c <- 1.8
pD_G9nc <- 0.009980444
pD_G9c <- 0.017822499
b0logitG9 <- log(pD_G9nc/(1-pD_G9nc))
b1logitG9 <- log(pD_G9c/(1-pD_G9c)) - b0logitG9
b0probitG9 <- qnorm(plogis(b0logitG9))
b1probitG9 <- qnorm(plogis(b0logitG9 + b1logitG9)) - b0probitG9
tau0G9 <- 0
tau1G9 <- b1probitG9/sqrt(1 + ((b1probitG9^2)*pG9_1*(1-pG9_1)))
###################################################
### G10: $7p36.3$ (VIPR2; WDR60), del$+$dup  
pG10_0 <- 0.9996968517
pG10_1 <- 1 - 0.9996968517
f10 <- 1 - sqrt(pG10_0)
OR_G10nc <- 1
OR_G10c <- 3.5
pD_G10nc <- 0.009992685
pD_G10c <- 0.034121973
b0logitG10 <- log(pD_G10nc/(1-pD_G10nc))
b1logitG10 <- log(pD_G10c/(1-pD_G10c)) - b0logitG10
b0probitG10 <- qnorm(plogis(b0logitG10))
b1probitG10 <- qnorm(plogis(b0logitG10 + b1logitG10)) - b0probitG10
tau0G10 <- 0
tau1G10 <- b1probitG10/sqrt(1 + ((b1probitG10^2)*pG10_1*(1-pG10_1)))
###################################################
### Environmental variables:
###################################################
### E1 = cannabis usage (variable 11).
### Number of categories contained within E1
levels_E1 <- 3
### probability distribution function (known)
pvecE1 <- c(0.70, 0.15, 0.15)
### OR estimates (known)
OR_E1_0 <- 1
OR_E1_1 <- 1.41
OR_E1_2 <- 2.78
### Penetrance function calculated using K and OR
# (just results presented)
pD_E1_0 <- 0.007562538
pD_E1_1 <- 0.010630218
pD_E1_2 <- 0.020744605
### Re-capture effect size estimates on the logit
# scale using the penetrance function
b0logitE1 <- log(pD_E1_0/(1-pD_E1_0))
b1logitE1 <- log(pD_E1_1/(1-pD_E1_1)) - b0logitE1
b2logitE1 <- log(pD_E1_2/(1-pD_E1_2)) - b0logitE1
### Calculate the effect size estimates on the probit
# scale using the logit scale betas
b0probitE1 <- qnorm(plogis(b0logitE1))
b1probitE1 <- qnorm(plogis(b0logitE1 + b1logitE1)) - b0probitE1
b2probitE1 <- qnorm(plogis(b0logitE1 + b2logitE1)) - b0probitE1
### Set reference category on L-scale
tau0E1 <- 0
### Write function, and use BB library to estimate
# effect size of E1 on L scale
### See main thesis for detail
V1E1 <- pvecE1[2]*(1-pvecE1[2])
V2E1 <- pvecE1[3]*(1-pvecE1[3])
functionE1 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (b1probitE1^2) - (b1probitE1^2)*V1E1*(p[1]^2) - (b1probitE1^2)*V2E1*(p[2]^2) + 2*(b1probitE1^2)*pvecE1[2]*pvecE1[3]*p[1]*p[2] - (p[1]^2)
	r[2] <- (b2probitE1^2) - (b2probitE1^2)*V1E1*(p[1]^2) - (b2probitE1^2)*V2E1*(p[2]^2) + 2*(b2probitE1^2)*pvecE1[2]*pvecE1[3]*p[1]*p[2] - (p[2]^2)
	r
}
p0 <- rep(0.001, 2)
outE1 <- dfsane(par=p0, fn=functionE1, control=list(trace=F))
outE1
tau1E1 <- abs(outE1$par[1])
tau2E1 <- abs(outE1$par[2])
###################################################
### E2 = Migration (variable 12).
levels_E2 <- 2
pvecE2 <- c(0.924, 1-0.924)
OR_E2_0 <- 1
OR_E2_1 <- 2.3
pD_E2_0 <- 0.009117829
pD_E2_1 <- 0.020725345
b0logitE2 <- log(pD_E2_0/(1-pD_E2_0))
b1logitE2 <- log(pD_E2_1/(1-pD_E2_1)) - b0logitE2
b0probitE2 <- qnorm(plogis(b0logitE2))
b1probitE2 <- qnorm(plogis(b0logitE2 + b1logitE2)) - b0probitE2
tau0E2 <- 0
### Calculate effect size estimate on the L-scale. Note 
# different method used compared to E1. Both give 
# the same answer
tau1E2 <- b1probitE2/sqrt(1 + ((b1probitE2^2)*pvecE2[2]*(1-pvecE2[2])))
###################################################
### E3 = Urbanicity (variable 13).
levels_E3 <- 3
pvecE3 <- c(0.25, 0.25, 0.5)
OR_E3_0 <- 1
OR_E3_1 <- 1.5
OR_E3_2 <- 2
pD_E3_0 <- 0.006181612
pD_E3_1 <- 0.009243848
pD_E3_2 <- 0.012287270
b0logitE3 <- log(pD_E3_0/(1-pD_E3_0))
b1logitE3 <- log(pD_E3_1/(1-pD_E3_1)) - b0logitE3
b2logitE3 <- log(pD_E3_2/(1-pD_E3_2)) - b0logitE3
b0probitE3 <- qnorm(plogis(b0logitE3))
b1probitE3 <- qnorm(plogis(b0logitE3 + b1logitE3)) - b0probitE3
b2probitE3 <- qnorm(plogis(b0logitE3 + b2logitE3)) - b0probitE3
tau0E3 <- 0
V1E3 <- pvecE3[2]*(1-pvecE3[2])
V2E3 <- pvecE3[3]*(1-pvecE3[3])
functionE3 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (b1probitE3^2) - (b1probitE3^2)*V1E3*(p[1]^2) - (b1probitE3^2)*V2E3*(p[2]^2) + 2*(b1probitE3^2)*pvecE3[2]*pvecE3[3]*p[1]*p[2] - (p[1]^2)
	r[2] <- (b2probitE3^2) - (b2probitE3^2)*V1E3*(p[1]^2) - (b2probitE3^2)*V2E3*(p[2]^2) + 2*(b2probitE3^2)*pvecE3[2]*pvecE3[3]*p[1]*p[2] - (p[2]^2)
	r
}
p0 <- rep(0.001, 2)
outE3 <- dfsane(par=p0, fn=functionE3, control=list(trace=F))
outE3
tau1E3 <- abs(outE3$par[1])
tau2E3 <- abs(outE3$par[2])
###################################################
### E4 = Paternal age (variable 14).
levels_E4 <- 7
pvecE4 <- c(0.342, 0.204, 0.252, 0.123, 0.052, 0.019, 0.008)
OR_E4_0 <- 1
OR_E4_1 <- 1.06
OR_E4_2 <- 1.06
OR_E4_3 <- 1.13
OR_E4_4 <- 1.22
OR_E4_5 <- 1.21
OR_E4_6 <- 1.66
pD_E4_0 <- 0.009404157
pD_E4_1 <- 0.009962785
pD_E4_2 <- 0.009962785
pD_E4_3 <- 0.010613722
pD_E4_4 <- 0.011449384
pD_E4_5 <- 0.011356603
pD_E4_6 <- 0.015514606
b0logitE4 <- log(pD_E4_0/(1-pD_E4_0))
b1logitE4 <- log(pD_E4_1/(1-pD_E4_1)) - b0logitE4
b2logitE4 <- log(pD_E4_2/(1-pD_E4_2)) - b0logitE4
b3logitE4 <- log(pD_E4_3/(1-pD_E4_3)) - b0logitE4
b4logitE4 <- log(pD_E4_4/(1-pD_E4_4)) - b0logitE4
b5logitE4 <- log(pD_E4_5/(1-pD_E4_5)) - b0logitE4
b6logitE4 <- log(pD_E4_6/(1-pD_E4_6)) - b0logitE4
b0probitE4 <- qnorm(plogis(b0logitE4))
b1probitE4 <- qnorm(plogis(b0logitE4 + b1logitE4)) - b0probitE4
b2probitE4 <- qnorm(plogis(b0logitE4 + b2logitE4)) - b0probitE4
b3probitE4 <- qnorm(plogis(b0logitE4 + b3logitE4)) - b0probitE4
b4probitE4 <- qnorm(plogis(b0logitE4 + b4logitE4)) - b0probitE4
b5probitE4 <- qnorm(plogis(b0logitE4 + b5logitE4)) - b0probitE4
b6probitE4 <- qnorm(plogis(b0logitE4 + b6logitE4)) - b0probitE4
tau0E4 <- 0
V1E4 <- pvecE4[2]*(1-pvecE4[2])
V2E4 <- pvecE4[3]*(1-pvecE4[3])
V3E4 <- pvecE4[4]*(1-pvecE4[4])
V4E4 <- pvecE4[5]*(1-pvecE4[5])
V5E4 <- pvecE4[6]*(1-pvecE4[6])
V6E4 <- pvecE4[7]*(1-pvecE4[7])
functionE4 <- function(p){
	r <- rep(NA, length(p))
	r[1] <- (b1probitE4^2) - (b1probitE4^2)*(V1E4*(p[1]^2) + V2E4*(p[2]^2) + V3E4*(p[3]^2) + V4E4*(p[4]^2) + V5E4*(p[5]^2) + V6E4*(p[6]^2) - 2*p[1]*p[2]*pvecE4[2]*pvecE4[3] - 2*p[1]*p[3]*pvecE4[2]*pvecE4[4] - 2*p[1]*p[4]*pvecE4[2]*pvecE4[5] - 2*p[1]*p[5]*pvecE4[2]*pvecE4[6] - 2*p[1]*p[6]*pvecE4[2]*pvecE4[7] - 2*p[2]*p[3]*pvecE4[3]*pvecE4[4] - 2*p[2]*p[4]*pvecE4[3]*pvecE4[5] - 2*p[2]*p[5]*pvecE4[3]*pvecE4[6] - 2*p[2]*p[6]*pvecE4[3]*pvecE4[7] - 2*p[3]*p[4]*pvecE4[4]*pvecE4[5] - 2*p[3]*p[5]*pvecE4[4]*pvecE4[6] - 2*p[3]*p[6]*pvecE4[4]*pvecE4[7] - 2*p[4]*p[5]*pvecE4[5]*pvecE4[6] - 2*p[4]*p[6]*pvecE4[5]*pvecE4[7] - 2*p[5]*p[6]*pvecE4[6]*pvecE4[7]) - (p[1]^2)
	r[2] <- (b2probitE4^2) - (b2probitE4^2)*(V1E4*(p[1]^2) + V2E4*(p[2]^2) + V3E4*(p[3]^2) + V4E4*(p[4]^2) + V5E4*(p[5]^2) + V6E4*(p[6]^2) - 2*p[1]*p[2]*pvecE4[2]*pvecE4[3] - 2*p[1]*p[3]*pvecE4[2]*pvecE4[4] - 2*p[1]*p[4]*pvecE4[2]*pvecE4[5] - 2*p[1]*p[5]*pvecE4[2]*pvecE4[6] - 2*p[1]*p[6]*pvecE4[2]*pvecE4[7] - 2*p[2]*p[3]*pvecE4[3]*pvecE4[4] - 2*p[2]*p[4]*pvecE4[3]*pvecE4[5] - 2*p[2]*p[5]*pvecE4[3]*pvecE4[6] - 2*p[2]*p[6]*pvecE4[3]*pvecE4[7] - 2*p[3]*p[4]*pvecE4[4]*pvecE4[5] - 2*p[3]*p[5]*pvecE4[4]*pvecE4[6] - 2*p[3]*p[6]*pvecE4[4]*pvecE4[7] - 2*p[4]*p[5]*pvecE4[5]*pvecE4[6] - 2*p[4]*p[6]*pvecE4[5]*pvecE4[7] - 2*p[5]*p[6]*pvecE4[6]*pvecE4[7]) - (p[2]^2)
	r[3] <- (b3probitE4^2) - (b3probitE4^2)*(V1E4*(p[1]^2) + V2E4*(p[2]^2) + V3E4*(p[3]^2) + V4E4*(p[4]^2) + V5E4*(p[5]^2) + V6E4*(p[6]^2) - 2*p[1]*p[2]*pvecE4[2]*pvecE4[3] - 2*p[1]*p[3]*pvecE4[2]*pvecE4[4] - 2*p[1]*p[4]*pvecE4[2]*pvecE4[5] - 2*p[1]*p[5]*pvecE4[2]*pvecE4[6] - 2*p[1]*p[6]*pvecE4[2]*pvecE4[7] - 2*p[2]*p[3]*pvecE4[3]*pvecE4[4] - 2*p[2]*p[4]*pvecE4[3]*pvecE4[5] - 2*p[2]*p[5]*pvecE4[3]*pvecE4[6] - 2*p[2]*p[6]*pvecE4[3]*pvecE4[7] - 2*p[3]*p[4]*pvecE4[4]*pvecE4[5] - 2*p[3]*p[5]*pvecE4[4]*pvecE4[6] - 2*p[3]*p[6]*pvecE4[4]*pvecE4[7] - 2*p[4]*p[5]*pvecE4[5]*pvecE4[6] - 2*p[4]*p[6]*pvecE4[5]*pvecE4[7] - 2*p[5]*p[6]*pvecE4[6]*pvecE4[7]) - (p[3]^2)
	r[4] <- (b4probitE4^2) - (b4probitE4^2)*(V1E4*(p[1]^2) + V2E4*(p[2]^2) + V3E4*(p[3]^2) + V4E4*(p[4]^2) + V5E4*(p[5]^2) + V6E4*(p[6]^2) - 2*p[1]*p[2]*pvecE4[2]*pvecE4[3] - 2*p[1]*p[3]*pvecE4[2]*pvecE4[4] - 2*p[1]*p[4]*pvecE4[2]*pvecE4[5] - 2*p[1]*p[5]*pvecE4[2]*pvecE4[6] - 2*p[1]*p[6]*pvecE4[2]*pvecE4[7] - 2*p[2]*p[3]*pvecE4[3]*pvecE4[4] - 2*p[2]*p[4]*pvecE4[3]*pvecE4[5] - 2*p[2]*p[5]*pvecE4[3]*pvecE4[6] - 2*p[2]*p[6]*pvecE4[3]*pvecE4[7] - 2*p[3]*p[4]*pvecE4[4]*pvecE4[5] - 2*p[3]*p[5]*pvecE4[4]*pvecE4[6] - 2*p[3]*p[6]*pvecE4[4]*pvecE4[7] - 2*p[4]*p[5]*pvecE4[5]*pvecE4[6] - 2*p[4]*p[6]*pvecE4[5]*pvecE4[7] - 2*p[5]*p[6]*pvecE4[6]*pvecE4[7]) - (p[4]^2)
	r[5] <- (b5probitE4^2) - (b5probitE4^2)*(V1E4*(p[1]^2) + V2E4*(p[2]^2) + V3E4*(p[3]^2) + V4E4*(p[4]^2) + V5E4*(p[5]^2) + V6E4*(p[6]^2) - 2*p[1]*p[2]*pvecE4[2]*pvecE4[3] - 2*p[1]*p[3]*pvecE4[2]*pvecE4[4] - 2*p[1]*p[4]*pvecE4[2]*pvecE4[5] - 2*p[1]*p[5]*pvecE4[2]*pvecE4[6] - 2*p[1]*p[6]*pvecE4[2]*pvecE4[7] - 2*p[2]*p[3]*pvecE4[3]*pvecE4[4] - 2*p[2]*p[4]*pvecE4[3]*pvecE4[5] - 2*p[2]*p[5]*pvecE4[3]*pvecE4[6] - 2*p[2]*p[6]*pvecE4[3]*pvecE4[7] - 2*p[3]*p[4]*pvecE4[4]*pvecE4[5] - 2*p[3]*p[5]*pvecE4[4]*pvecE4[6] - 2*p[3]*p[6]*pvecE4[4]*pvecE4[7] - 2*p[4]*p[5]*pvecE4[5]*pvecE4[6] - 2*p[4]*p[6]*pvecE4[5]*pvecE4[7] - 2*p[5]*p[6]*pvecE4[6]*pvecE4[7]) - (p[5]^2)
	r[6] <- (b6probitE4^2) - (b6probitE4^2)*(V1E4*(p[1]^2) + V2E4*(p[2]^2) + V3E4*(p[3]^2) + V4E4*(p[4]^2) + V5E4*(p[5]^2) + V6E4*(p[6]^2) - 2*p[1]*p[2]*pvecE4[2]*pvecE4[3] - 2*p[1]*p[3]*pvecE4[2]*pvecE4[4] - 2*p[1]*p[4]*pvecE4[2]*pvecE4[5] - 2*p[1]*p[5]*pvecE4[2]*pvecE4[6] - 2*p[1]*p[6]*pvecE4[2]*pvecE4[7] - 2*p[2]*p[3]*pvecE4[3]*pvecE4[4] - 2*p[2]*p[4]*pvecE4[3]*pvecE4[5] - 2*p[2]*p[5]*pvecE4[3]*pvecE4[6] - 2*p[2]*p[6]*pvecE4[3]*pvecE4[7] - 2*p[3]*p[4]*pvecE4[4]*pvecE4[5] - 2*p[3]*p[5]*pvecE4[4]*pvecE4[6] - 2*p[3]*p[6]*pvecE4[4]*pvecE4[7] - 2*p[4]*p[5]*pvecE4[5]*pvecE4[6] - 2*p[4]*p[6]*pvecE4[5]*pvecE4[7] - 2*p[5]*p[6]*pvecE4[6]*pvecE4[7]) - (p[6]^2)
	r
}
p0 <- rep(0, 6)
outE4 <- dfsane(par=p0, fn=functionE4, control=list(trace=F))
outE4
tau1E4 <- abs(outE4$par[1])
tau2E4 <- abs(outE4$par[2])
tau3E4 <- abs(outE4$par[3])
tau4E4 <- abs(outE4$par[4])
tau5E4 <- abs(outE4$par[5])
tau6E4 <- abs(outE4$par[6])
###################################################
### E5 = Childhood adversity (variable 15).
levels_E5 <- 2
pvecE5 <- c(0.73, 1-0.73)
OR_E5_0 <- 1
OR_E5_1 <- 2.78
pD_E5_0 <- 0.006795188
pD_E5_1 <- 0.018776172
b0logitE5 <- log(pD_E5_0/(1-pD_E5_0))
b1logitE5 <- log(pD_E5_1/(1-pD_E5_1)) - b0logitE5
b0probitE5 <- qnorm(plogis(b0logitE5))
b1probitE5 <- qnorm(plogis(b0logitE5 + b1logitE5)) - b0probitE5
tau0E5 <- 0
tau1E5 <- b1probitE5/sqrt(1 + ((b1probitE5^2)*pvecE5[2]*(1-pvecE5[2])))
###################################################
### STEP 1B: Calculate E[L] #######################
###################################################
### Done as a loop- but definitely more efficient
# ways to calculate this. Loop commented out and
# final result given to save run time!
pvecG1 <- c(pG1_0, pG1_1)
pvecG2 <- c(pG2_0, pG2_1)
pvecG3 <- c(pG3_0, pG3_1)
pvecG4 <- c(pG4_0, pG4_1)
pvecG5 <- c(pG5_0, pG5_1)
pvecG6 <- c(pG6_0, pG6_1)
pvecG7 <- c(pG7_0, pG7_1)
pvecG8 <- c(pG8_0, pG8_1)
pvecG9 <- c(pG9_0, pG9_1)
pvecG10 <- c(pG10_0, pG10_1)
# EHGE <- 0
# for(ga in 1:2){
	# for(gb in 1:2){
		# for(gc in 1:2){
			# for(gd in 1:2){
				# for(ge in 1:2){
					# for(gf in 1:2){
						# for(gg in 1:2){
							# for(gh in 1:2){
								# for(gi in 1:2){
									# for(gj in 1:2){
										# for(ea in 1:levels_E1){
											# for(eb in 1:levels_E2){
												# for(ec in 1:levels_E3){
													# for(ed in 1:levels_E4){
														# for(ee in 1:levels_E5){
															# x.iter <- rep(0, (10 + (levels_E1 - 1) + (levels_E2 - 1) + (levels_E3 - 1) + (levels_E4 - 1) + (levels_E5 - 1)))
				# x.iter[1] <- as.numeric(ga == 2)
				# x.iter[2] <- as.numeric(gb == 2)
				# x.iter[3] <- as.numeric(gc == 2)
				# x.iter[4] <- as.numeric(gd == 2)
				# x.iter[5] <- as.numeric(ge == 2)
				# x.iter[6] <- as.numeric(gf == 2)
				# x.iter[7] <- as.numeric(gg == 2)
				# x.iter[8] <- as.numeric(gh == 2)
				# x.iter[9] <- as.numeric(gi == 2)
				# x.iter[10] <- as.numeric(gj == 2)
				# x.iter[11] <- as.numeric(ea == 2)	
				# x.iter[12] <- as.numeric(ea == 3)										
				# x.iter[13] <- as.numeric(eb == 2)	
				# x.iter[14] <- as.numeric(ec == 2)
				# x.iter[15] <- as.numeric(ec == 3)
				# x.iter[16] <- as.numeric(ed == 2)
				# x.iter[17] <- as.numeric(ed == 3)
				# x.iter[18] <- as.numeric(ed == 4)
				# x.iter[19] <- as.numeric(ed == 5)
				# x.iter[20] <- as.numeric(ed == 6)	
				# x.iter[21] <- as.numeric(ed == 7)
				# x.iter[22] <- as.numeric(ee == 2)
				# hx.iter <- hGEfunction(x=x.iter, tau=tau)
				# prob.iter <- prod(c(pvecG1[ga], pvecG2[gb], pvecG3[gc], pvecG4[gd], pvecG5[ge], pvecG6[gf], pvecG7[gg], pvecG8[gh], pvecG9[gi], pvecG10[gj], pvecE1[ea], pvecE2[eb], pvecE3[ec], pvecE4[ed], pvecE5[ee]))
				# EHGE.iter <- hx.iter*prob.iter
				# EHGE <- EHGE + EHGE.iter

														# }
													# }
												# }
											# }
										# }
									# }
								# }
							# }
						# }
					# }
				# }
			# }
		# }
	# }
# }
# EHGE
EHGE <- 0.3856964
###################################################
### STEP 1C: Calculate disease threshold ##########
###################################################
### This is an approximation
D_T_joint <- EHGE - qnorm(K)
###################################################
### STEP 1D: Calculate the variance components ####
### for the above measured risk variables #########
###################################################
VXi_function <- function(taui, pi){
	Vp1 <- sum((taui^2)*(pi)*(1-pi))
#	Vp3 <- sum((taui^2)*(pi^2))
	Vp2 <- matrix(0, nrow=length(pi), ncol=length(pi))
	for(j in 1:length(pi)){
		for(k in 1:length(pi)){
			Vp2[j,k] <- taui[j]*taui[k]*pi[j]*pi[k]
		}
	}
	diag(Vp2) <- 0
	Vp2 <- sum(Vp2)
	Vp1 - Vp2
}
### Genetic major risk loci variance comp.
VhG <- c(VXi_function(taui=c(tau0G1, tau1G1), pi=pvecG1),
VXi_function(taui=c(tau0G2, tau1G2), pi=pvecG2),
VXi_function(taui=c(tau0G3, tau1G3), pi=pvecG3),
VXi_function(taui=c(tau0G4, tau1G4), pi=pvecG4),
VXi_function(taui=c(tau0G5, tau1G5), pi=pvecG5),
VXi_function(taui=c(tau0G6, tau1G6), pi=pvecG6),
VXi_function(taui=c(tau0G7, tau1G7), pi=pvecG7),
VXi_function(taui=c(tau0G8, tau1G8), pi=pvecG8),
VXi_function(taui=c(tau0G9, tau1G9), pi=pvecG9),
VXi_function(taui=c(tau0G10, tau1G10), pi=pvecG10))
VhG <- sum(VhG)
### Discrete env. variance comp.
VhE <- c(VXi_function(taui=c(tau0E1, tau1E1, tau2E1), pi=pvecE1),
VXi_function(taui=c(tau0E2, tau1E2), pi=pvecE2),
VXi_function(taui=c(tau0E3, tau1E3, tau2E3), pi=pvecE3),
VXi_function(taui=c(tau0E4, tau1E4, tau2E4, tau3E4, tau4E4, tau5E4, tau6E4), pi=pvecE4),
VXi_function(taui=c(tau0E5, tau1E5), pi=pvecE5))
VhE <- sum(VhE)
### Assume that all of genetic contribution is additive:
VAhG <- VhG
VDhG <- 0
###################################################
### STEP 1E: Define other variables to be used ####
###################################################
### Parent = family history of disease and so...
R <- "P"
# Coefficient of relatedness is:
r <- 0.5
# Coefficient of coancestry is:
theta = 0.25
### General input definitions:
D_threshold <- D_T_joint
VM <- VML
maf.vec <- c(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10)
tau <- c(tau1G1, 2*tau1G1, tau1G2,2*tau1G2, tau1G3, 2*tau1G3 , tau1G4, 2*tau1G4, tau1G5, 2*tau1G5, tau1G6, 2*tau1G6, tau1G7, tau1G7*2, tau1G8, tau1G8*2, tau1G9, 2*tau1G9, tau1G10, tau1G10*2, tau1E1, tau2E1, tau1E2, tau1E3, tau2E3, tau1E4, tau2E4, tau3E4, tau4E4, tau5E4, tau6E4, tau1E5)
levelsE <- c(levels_E1, levels_E2, levels_E3, levels_E4, levels_E5)
pEmat <- matrix(0, nrow=length(levelsE), ncol=max(levelsE))
pEmat[1, 1: levels_E1] <- pvecE1
pEmat[2, 1: levels_E2] <- pvecE2
pEmat[3, 1: levels_E3] <- pvecE3
pEmat[4, 1: levels_E4] <- pvecE4
pEmat[5, 1: levels_E5] <- pvecE5
n_E <- length(levelsE)
n_G <- length(maf.vec)
combosR_E <- expand.grid(1:(levelsE[1]), 1:(levelsE[2]), 1:(levelsE[3]), 1:(levelsE[4]), 1:(levelsE[5]))
###################################################
### STEP 2: Compare Methods A, B and C ############
###################################################
### For a restricted set of PRS values, and major 
# and environmental risk variable combinations we 
# compare the outputted results for methods A, B and 
# C. 
### This restricted risk variable set is defined by
# combinations of: 
### PRS vector
mIvec <- c(-0.18, 0, 0.18)
### x2 extreme environmental risk factor profiles
eImin <- c(1,1,1,1,1)
eImax <- levelsE
### x2 Schizophrenia CNV risk profiles
# 1. carry no CNVs
gImin <- rep(0,10)
# 2. carry risk variant at G1
g1I <- c(1, rep(0,9))
###################################################
### Generate output for method A...
###################################################
### NOTE: This code will take a long time to run. 
# We have therefore commented it out so it is not 
# accidentally run.
### If going to run please update pathway for writing 
# results to file...
# ### 1. Min risk major loci, min risk environmental:
# outA_work <- riskLTMM_methodA(mIvec=mIvec, g_Ivec=gImin, e_Ivec=eImin, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, D_threshold=D_threshold, tau=tau, VM=VM, VhG=VhG, VhE=VhE, H2=H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E=combosR_E, R="P")
### Create output matrix
# out_methodA <- cbind(rep(0, 6), rep(0,6), c(mIvec, mIvec), c(rep(1,3), rep(0,3)), c(outA_work[,1], outA_work[,2]))
# colnames(out_methodA) <- c("GI", "EI", "mI", "YR", "risk")
# ### 2. Min risk major loci, max risk environmental:
# outA_work <- riskLTMM_methodA(mIvec=mIvec, g_Ivec=gImin, e_Ivec=eImax, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, D_threshold=D_threshold, tau=tau, VM=VM, VhG=VhG, VhE=VhE, H2=H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E=combosR_E, R="P")
# ### Update output matrix:
# out_methodA <- rbind(cbind(rep(0, 6), rep(1,6), c(mIvec, mIvec), c(rep(1,3), rep(0,3)), c(outA_work[,1], outA_work[,2])), out_methodA)
# ### 3. G1 == 1, min risk environmental:
# outA_work <- riskLTMM_methodA(mIvec=mIvec, g_Ivec=g1I, e_Ivec=eImin, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, D_threshold=D_threshold, tau=tau, VM=VM, VhG=VhG, VhE=VhE, H2=H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E=combosR_E, R="P")
# ### Update output matrix:
# out_methodA <- rbind(out_methodA,cbind(rep(1, 6), rep(0,6), c(mIvec, mIvec), c(rep(1,3), rep(0,3)), c(outA_work[,1], outA_work[,2])))
# ### 4. G1 == 1, max risk environmental:
# outA_work <- riskLTMM_methodA(mIvec=mIvec, g_Ivec=g1I, e_Ivec=eImax, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, D_threshold=D_threshold, tau=tau, VM=VM, VhG=VhG, VhE=VhE, H2=H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E=combosR_E, R="P")
# out_methodA <- rbind(out_methodA,cbind(rep(1, 6), rep(1,6), c(mIvec, mIvec), c(rep(1,3), rep(0,3)), c(outA_work[,1], outA_work[,2])))
# colnames(out_methodA) <- c("GI", "EI", "mI", "YR", "risk")
# ### Write your output to file for future reference:
# write.csv(out_methodA, file="XXX/methodA_scz_example.csv")
###################################################
### Generate output for method B...
###################################################
### 1. Min risk major loci, min risk environmental:
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=gImin, e_Ivec=eImin, maf.vec=maf.vec, levelsE=levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_gImin_eImin_YR1MB <- riskB[,1]
risk_gImin_eImin_YR0MB <- riskB[,2]
### 2. Min risk major loci, max risk environmental:
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=gImin, e_Ivec=eImax, maf.vec=maf.vec, levelsE=levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_gImin_eImax_YR1MB <- riskB[,1]
risk_gImin_eImax_YR0MB <- riskB[,2]
### 3. G1 == 1, min risk environmental:
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g1I, e_Ivec=eImin, maf.vec=maf.vec, levelsE=levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g1I_eImin_YR1MB <- riskB[,1]
risk_g1I_eImin_YR0MB <- riskB[,2]
### 4. G1 == 1, max risk environmental:
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g1I, e_Ivec=eImax, maf.vec=maf.vec, levelsE=levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g1I_eImax_YR1MB <- riskB[,1]
risk_g1I_eImax_YR0MB <- riskB[,2]
### Create output matrix for method B:
out_methodB <- cbind(c(rep(1, 12), rep(0, 12)),c(rep(0, 6), rep(1,6), rep(0,6), rep(1,6)), rep(mIvec, 8), c(rep(1,3), rep(0,3), rep(1,3), rep(0,3), rep(1,3), rep(0,3), rep(1,3), rep(0,3)), c(risk_g1I_eImin_YR1MB, risk_g1I_eImin_YR0MB, risk_g1I_eImax_YR1MB, risk_g1I_eImax_YR0MB, risk_gImin_eImin_YR1MB, risk_gImin_eImin_YR0MB, risk_gImin_eImax_YR1MB, risk_gImin_eImax_YR0MB))
colnames(out_methodB) <- c("GI", "EI", "mI", "YR", "risk")
### Write output to file. Update file path and name:
# write.csv(out_methodB, file="XXX/methodB_scz_example.csv")
###################################################
### Generate output for method C...
###################################################
### Commented out as calculation of mean and variances
# for individual {R} takes time. Uncomment to run...
### Additional method C inputs:
### Possible combinations of major risk loci
# n_G <- length(maf.vec)
# df_G <- NULL
# for(i in 1:n_G){
	# df_G <- cbind(df_G, 0:2)
# }
# df_G <- data.frame(df_G)
# combosR_G <- expand.grid(df_G)
# ### Means: LR | GI=gI
# # Min risk major loci
# mu_LR_gImin <- mu_LR_gI_f(gIvec=gImin, tau=tau, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, R="P", combosR_E=combosR_E, combosR_G=combosR_G)
# # G1 == 1
# mu_LR_g1I <- mu_LR_gI_f(gIvec=g1I, tau=tau, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, R="P", combosR_E=combosR_E, combosR_G=combosR_G)
# ### Variances: LR | GI=gI
# # Min risk major loci
# sigma2_LR_gImin <- sigma2_LR_gI_f(gIvec=gImin, tau=tau, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, R="P", mu_LR_gI=mu_LR_gImin, VhG=VhG, VhE=VhE, combosR_E=combosR_E, combosR_G=combosR_G)
# # G1 == 1
# sigma2_LR_g1I <- sigma2_LR_gI_f(gIvec=g1I, tau=tau, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, R="P", mu_LR_gI=mu_LR_g1I, VhG=VhG, VhE=VhE, combosR_E=combosR_E, combosR_G=combosR_G)
# ### Means individual I: 
# # using mI=0 regardless of mI value as want mean given GI and EI only here. 
# # Value used for x_R also does not matter- reference used.
# mu_LI_gImin_eImin <- mean_mI_f(tau=tau, x_I=c(0,0,rep(0, 2*9), rep(0,2), 0, rep(0, 2), rep(0,6), 0), x_R=c(rep(0, 20), rep(0,12)), r=0.5, mI=0)
# mu_LI_gImin_eImax <- mean_mI_f(tau=tau, x_I=c(0,0,rep(0, 2*9), 0, 1, 1, 0, 1, rep(0,5), 1, 1), x_R=c(rep(0, 20), rep(0,12)), r=0.5, mI=0)
# mu_LI_g1I_eImin <- mean_mI_f(tau=tau, x_I=c(1,0,rep(0, 2*9), rep(0,2), 0, rep(0, 2), rep(0,6), 0), x_R=c(rep(0, 20), rep(0,12)), r=0.5, mI=0)
# mu_LI_g1I_eImax <- mean_mI_f(tau=tau, x_I=c(1,0,rep(0, 2*9), 0, 1, 1, 0, 1, rep(0,5), 1, 1), x_R=c(rep(0, 20), rep(0,12)), r=0.5, mI=0)
# ### Output from method C
# outC_gImin_eImin_YR1 <- rep(0, length(mIvec))
# outC_gImin_eImax_YR1 <- rep(0, length(mIvec))
# outC_g1I_eImin_YR1 <- rep(0, length(mIvec))
# outC_g1I_eImax_YR1 <- rep(0, length(mIvec))
# outC_gImin_eImin_YR0 <- rep(0, length(mIvec))
# outC_gImin_eImax_YR0 <- rep(0, length(mIvec))
# outC_g1I_eImin_YR0 <- rep(0, length(mIvec))
# outC_g1I_eImax_YR0 <- rep(0, length(mIvec))
# for(i in 1:length(mIvec)){
	# outC_gImin_eImin_YR1[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_gImin_eImin[1], mu_LR_gI= mu_LR_gImin, sigma2_LR_gI= sigma2_LR_gImin, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=1, D_threshold=D_threshold)
	# outC_gImin_eImin_YR0[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_gImin_eImin[1], mu_LR_gI= mu_LR_gImin, sigma2_LR_gI= sigma2_LR_gImin, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=0, D_threshold=D_threshold)
	# outC_gImin_eImax_YR1[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_gImin_eImax[1], mu_LR_gI= mu_LR_gImin, sigma2_LR_gI= sigma2_LR_gImin, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=1, D_threshold=D_threshold)
	# outC_gImin_eImax_YR0[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_gImin_eImax[1], mu_LR_gI= mu_LR_gImin, sigma2_LR_gI= sigma2_LR_gImin, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=0, D_threshold=D_threshold)
	# outC_g1I_eImin_YR1[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_g1I_eImin[1], mu_LR_gI= mu_LR_g1I, sigma2_LR_gI= sigma2_LR_g1I, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=1, D_threshold=D_threshold)
	# outC_g1I_eImin_YR0[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_g1I_eImin[1], mu_LR_gI= mu_LR_g1I, sigma2_LR_gI= sigma2_LR_g1I, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=0, D_threshold=D_threshold)
	# outC_g1I_eImax_YR1[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_g1I_eImax[1], mu_LR_gI= mu_LR_g1I, sigma2_LR_gI= sigma2_LR_g1I, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=1, D_threshold=D_threshold)
	# outC_g1I_eImax_YR0[i] <- riskLTMM_methodC(mu_LI_gIeI= mu_LI_g1I_eImax[1], mu_LR_gI= mu_LR_g1I, sigma2_LR_gI= sigma2_LR_g1I, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=0, D_threshold=D_threshold)
# }
# ### Put into an output matrix
# out_methodC <- cbind(c(rep(1, 12), rep(0, 12)),c(rep(0, 6), rep(1,6), rep(0,6), rep(1,6)), rep(mIvec, 8), c(rep(1,3), rep(0,3), rep(1,3), rep(0,3), rep(1,3), rep(0,3), rep(1,3), rep(0,3)), c(outC_g1I_eImin_YR1, outC_g1I_eImin_YR0, outC_g1I_eImax_YR1, outC_g1I_eImax_YR0, outC_gImin_eImin_YR1, outC_gImin_eImin_YR0, outC_gImin_eImax_YR1, outC_gImin_eImax_YR0))
# colnames(out_methodC) <- c("GI", "EI", "mI", "YR", "risk")
### Write output to file. Update file path and name:
# write.csv(out_methodC, file="XXX/methodC_scz_example.csv")
###################################################
### Output from the 3 computational methods can now
# be compared. Code not given. In thesis methods are
# presented in a table==easy to create within R.
###################################################
### It is shown that all methods produce very similar 
# estimates of risk. Method B is the quickest here and
# we therefore carry on in this chapter using method B
# only for the LTMM risk calculations
###################################################
### STEP 3: Risk estimation #######################
###################################################
### Risk is calculated for an updated, but still 
# restricted set of risk variables. We add:
### x2 additional CNV risk profiles
# 3. carry risk variant at G3
g3I <- c(0,0,1,rep(0, 7))
# 4. carry risk variant at G9
g9I <- c(rep(0, 8),1,0)
### More PRS values:
Mlower <- -1.31
MUpper <- 1.31
mIvec <- c(0, seq(from= Mlower, to= MUpper, length.out = 100))
### RISK via method B
# g1I and eImin
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g1I, e_Ivec=eImin, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g1I_eImin_YR1MB <- riskB[,1]
risk_g1I_eImin_YR0MB <- riskB[,2]
# Create output matrix
scz_results <- cbind(rep(1, 2*length(risk_g1I_eImin_YR1MB)), rep(0, 2*length(risk_g1I_eImin_YR1MB)), rep(0, 2*length(risk_g1I_eImin_YR1MB)), rep(0, 2*length(risk_g1I_eImin_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_g1I_eImin_YR1MB)), rep(0, length(risk_g1I_eImin_YR1MB))), c(risk_g1I_eImin_YR1MB, risk_g1I_eImin_YR0MB))
# g1I and eImax
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g1I, e_Ivec=eImax, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g1I_eImax_YR1MB <- riskB[,1]
risk_g1I_eImax_YR0MB <- riskB[,2]
# Update output matrix
scz_results <- rbind(scz_results, cbind(rep(1, 2*length(risk_g1I_eImax_YR1MB)), rep(0, 2*length(risk_g1I_eImax_YR1MB)), rep(0, 2*length(risk_g1I_eImax_YR1MB)), rep(1, 2*length(risk_g1I_eImax_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_g1I_eImax_YR1MB)), rep(0, length(risk_g1I_eImax_YR1MB))), c(risk_g1I_eImax_YR1MB, risk_g1I_eImax_YR0MB)))
# gImin and eImin
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=gImin, e_Ivec=eImin, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_gImin_eImin_YR1MB <- riskB[,1]
risk_gImin_eImin_YR0MB <- riskB[,2]
# Update output matrix
scz_results <- rbind(scz_results, cbind(rep(0, 2*length(risk_gImin_eImin_YR1MB)), rep(0, 2*length(risk_gImin_eImin_YR1MB)), rep(0, 2*length(risk_gImin_eImin_YR1MB)), rep(0, 2*length(risk_gImin_eImin_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_gImin_eImin_YR1MB)), rep(0, length(risk_gImin_eImin_YR1MB))), c(risk_gImin_eImin_YR1MB, risk_gImin_eImin_YR0MB)))
# gImin and eImax
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=gImin, e_Ivec=eImax, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_gImin_eImax_YR1MB <- riskB[,1]
risk_gImin_eImax_YR0MB <- riskB[,2]
# Update output matrix
scz_results <- rbind(scz_results, cbind(rep(0, 2*length(risk_gImin_eImax_YR1MB)), rep(0, 2*length(risk_gImin_eImax_YR1MB)), rep(0, 2*length(risk_gImin_eImax_YR1MB)), rep(1, 2*length(risk_gImin_eImax_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_gImin_eImax_YR1MB)), rep(0, length(risk_gImin_eImax_YR1MB))), c(risk_gImin_eImax_YR1MB, risk_gImin_eImax_YR0MB)))
# g3I and eImin
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g3I, e_Ivec=eImin, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g3I_eImin_YR1MB <- riskB[,1]
risk_g3I_eImin_YR0MB <- riskB[,2]
# g3I and eImax
riskB <- EXACTrisk_function(mIvec=mIvec, g_Ivec=g3I, e_Ivec=eImax, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g3I_eImax_YR1MB <- riskB[,1]
risk_g3I_eImax_YR0MB <- riskB[,2]
# g9I and eImin
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g9I, e_Ivec=eImin, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g9I_eImin_YR1MB <- riskB[,1]
risk_g9I_eImin_YR0MB <- riskB[,2]
# g9I and eImax
riskB <- riskLTMM_methodB(mIvec=mIvec, g_Ivec=g9I, e_Ivec=eImax, maf.vec=maf.vec, levelsE= levelsE, pEmat= pEmat, D_threshold= D_threshold, tau= tau, VM= VM, VhG= VhG, VhE= VhE, H2= H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25, combosR_E= combosR_E, R="P")
risk_g9I_eImax_YR1MB <- riskB[,1]
risk_g9I_eImax_YR0MB <- riskB[,2]
# Update output matrix
scz_results <- rbind(scz_results, cbind(rep(0, 2*length(risk_g3I_eImin_YR1MB)), rep(1, 2*length(risk_g3I_eImin_YR1MB)), rep(0, 2*length(risk_g3I_eImin_YR1MB)), rep(0, 2*length(risk_g3I_eImin_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_g3I_eImin_YR1MB)), rep(0, length(risk_g3I_eImin_YR1MB))), c(risk_g3I_eImin_YR1MB, risk_g3I_eImin_YR0MB)))
scz_results <- rbind(scz_results, cbind(rep(0, 2*length(risk_g3I_eImax_YR1MB)), rep(1, 2*length(risk_g3I_eImax_YR1MB)), rep(0, 2*length(risk_g3I_eImax_YR1MB)), rep(1, 2*length(risk_g3I_eImax_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_g3I_eImax_YR1MB)), rep(0, length(risk_g3I_eImax_YR1MB))), c(risk_g3I_eImax_YR1MB, risk_g3I_eImax_YR0MB)))
scz_results <- rbind(scz_results, cbind(rep(0, 2*length(risk_g9I_eImin_YR1MB)), rep(0, 2*length(risk_g9I_eImin_YR1MB)), rep(1, 2*length(risk_g9I_eImin_YR1MB)), rep(0, 2*length(risk_g9I_eImin_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_g9I_eImin_YR1MB)), rep(0, length(risk_g9I_eImin_YR1MB))), c(risk_g9I_eImin_YR1MB, risk_g9I_eImin_YR0MB)))
scz_results <- rbind(scz_results, cbind(rep(0, 2*length(risk_g9I_eImax_YR1MB)), rep(0, 2*length(risk_g9I_eImax_YR1MB)), rep(1, 2*length(risk_g9I_eImax_YR1MB)), rep(1, 2*length(risk_g9I_eImax_YR1MB)), c(mIvec, mIvec),c(rep(1, length(risk_g9I_eImax_YR1MB)), rep(0, length(risk_g9I_eImax_YR1MB))), c(risk_g9I_eImax_YR1MB, risk_g9I_eImax_YR0MB)))
colnames(scz_results) <- c("G1I", "G3I","G9I","EI", "mI", "YR", "risk")
### write output to file. UPDATE file path and name
# write.csv(scz_results, file="XXX/scz_results_LTMM.csv")
### Results can then be read in and explored
### We used ggplot2 package to visually explore results
### Please see the thesis for the resulting graphs 
###################################################
###################################################
###################################################
### Log-linear method #############################
###################################################
### Additional and differing inputs:
###################################################
### Relative risk calculated using penetrance function
###################################################
### CNVs:
###################################################
RR_G1nc <- 1
RR_G1c <- pD_G1c/pD_G1nc
RR_G2nc <- 1
RR_G2c <- pD_G2c/pD_G2nc
RR_G3nc <- 1
RR_G3c <- pD_G3c/pD_G3nc
RR_G4nc <- 1
RR_G4c <- pD_G4c/pD_G4nc
RR_G5nc <- 1
RR_G5c <- pD_G5c/pD_G5nc
RR_G6nc <- 1
RR_G6c <- pD_G6c/pD_G6nc
RR_G7nc <- 1
RR_G7c <- pD_G7c/pD_G7nc
RR_G8nc <- 1
RR_G8c <- pD_G8c/pD_G8nc
RR_G9nc <- 1
RR_G9c <- pD_G9c/pD_G9nc
RR_G10nc <- 1
RR_G10c <- pD_G10c/pD_G10nc
RR_Gmat <- rbind(
c(1, RR_G1c, exp(log(RR_G1c)*2)),
c(1, RR_G2c, exp(log(RR_G2c)*2)),
c(1, RR_G3c, exp(log(RR_G3c)*2)),
c(1, RR_G4c, exp(log(RR_G4c)*2)),
c(1, RR_G5c, exp(log(RR_G5c)*2)),
c(1, RR_G6c, exp(log(RR_G6c)*2)),
c(1, RR_G7c, exp(log(RR_G7c)*2)),
c(1, RR_G8c, exp(log(RR_G8c)*2)),
c(1, RR_G9c, exp(log(RR_G9c)*2)),
c(1, RR_G10c, exp(log(RR_G10c)*2))
)
###################################################
### Environmental variables:
###################################################
RR_E1_0 <- 1
RR_E1_1 <- pD_E1_1/pD_E1_0
RR_E1_2 <- pD_E1_2/pD_E1_0
RR_E2_0 <- 1
RR_E2_1 <- pD_E2_1/pD_E2_0
RR_E3_0 <- 1
RR_E3_1 <- pD_E3_1/pD_E3_0
RR_E3_2 <- pD_E3_2/pD_E3_0
RR_E4_0 <- 1
RR_E4_1 <- pD_E4_1/pD_E4_0
RR_E4_2 <- pD_E4_2/pD_E4_0
RR_E4_3 <- pD_E4_3/pD_E4_0
RR_E4_4 <- pD_E4_4/pD_E4_0
RR_E4_5 <- pD_E4_5/pD_E4_0
RR_E4_6 <- pD_E4_6/pD_E4_0
RR_E5_0 <- 1
RR_E5_1 <- pD_E5_1/pD_E5_0
RR_Emat <- matrix(1, nrow=length(levelsE), ncol=max(levelsE))
RR_Emat[1, 2:levels_E1] <- c(RR_E1_1, RR_E1_2)
RR_Emat[2, 2:levels_E2] <- c(RR_E2_1)
RR_Emat[3, 2:levels_E3] <- c(RR_E3_1, RR_E3_2)
RR_Emat[4, 2:levels_E4] <- c(RR_E4_1, RR_E4_2, RR_E4_3, RR_E4_4, RR_E4_5, RR_E4_6)
RR_Emat[5, 2:levels_E5] <- c(RR_E5_1)
###################################################
### Var[PRS] for the log-linear scale
###################################################
### Known inputs
ORprs <- c(1.072, 1.072, 1.069, 0.934, 1.078, 0.914, 1.075, 1.063, 0.904, 1.086, 0.857, 0.944, 0.929, 0.937, 0.929, 1.081, 0.909, 0.939, 0.930, 1.101, 0.940, 1.085, 1.071, 0.941, 0.933, 1.064, 0.857, 1.065, 0.934, 0.926, 0.922, 1.101, 1.076, 1.058, 0.942, 1.205, 0.942, 1.068, 0.849, 0.922, 1.101, 0.904, 1.068, 1.061, 1.083, 0.924, 1.066, 1.073, 0.908, 0.941, 1.079, 0.919, 1.087, 1.069, 1.125, 1.064, 0.906, 1.068, 0.941, 0.927, 1.066, 0.926, 1.076, 1.091, 0.941, 0.943, 1.060, 1.068, 0.915, 0.933, 1.073, 0.939, 1.088, 1.074, 1.067, 1.060, 0.923, 0.922, 1.067, 1.077, 1.073, 0.920, 1.081, 1.071, 0.939, 0.934, 0.931, 1.071, 0.928, 0.937, 0.930, 1.087, 0.947, 1.090, 0.938, 0.846, 1.317, 0.914, 1.069, 1.155, 0.843, 0.912, 0.846, 1.093, 1.188, 1.076, 0.902, 0.873)
fcontrol_prs <- c(0.527, 0.301, 0.296, 0.677, 0.358, 0.164, 0.184, 0.685, 0.101, 0.163, 0.961, 0.593, 0.458, 0.337, 0.643, 0.805, 0.754, 0.326, 0.552, 0.156, 0.48, 0.324, 0.529, 0.615, 0.449, 0.314, 0.922, 0.47, 0.761, 0.802, 0.532, 0.883, 0.213, 0.475, 0.621, 0.85, 0.368, 0.505, 0.0476, 0.423, 0.0959, 0.123, 0.332, 0.628, 0.647, 0.152, 0.642, 0.219, 0.116, 0.6, 0.174, 0.803, 0.424, 0.337, 0.889, 0.499, 0.85, 0.311, 0.334, 0.314, 0.322, 0.335, 0.46, 0.797, 0.337, 0.366, 0.646, 0.474, 0.741, 0.771, 0.249, 0.418, 0.287, 0.274, 0.52, 0.465, 0.257, 0.562, 0.281, 0.223, 0.51, 0.859, 0.162, 0.627, 0.614, 0.628, 0.769, 0.766, 0.322, 0.754, 0.594, 0.232, 0.753, 0.831, 0.759, 0.964, 0.0191, 0.208, 0.252, 0.046, 0.0756, 0.624, 0.0708, 0.725, 0.0292, 0.81, 0.913, 0.0797)
### Var[PRS]
VMRR <- sum((log(ORprs)^2)*2*fcontrol_prs*(1-fcontrol_prs))
mIvec <- mIvec/sqrt(VML)
### Generating results using the log-linear model
### Affected relatives only:
# gImin and eImin
outll <- riskLL(mI_vec=mIvec, eI_vec=eImin, gI_vec=gImin, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Create log-linear output matrix
scz_results_ll <- cbind(rep(0,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), mIvec, rep(1,length(mIvec)), outll)
# g1I and eImin
outll <- riskLL(mI_vec=mIvec, eI_vec=eImin, gI_vec=g1I, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(1,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
# g3I and eImin
outll <- riskLL(mI_vec=mIvec, eI_vec=eImin, gI_vec=g3I, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(0,length(mIvec)), rep(1,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
# g9I and eImin
outll <- riskLL(mI_vec=mIvec, eI_vec=eImin, gI_vec=g9I, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(0,length(mIvec)), rep(0,length(mIvec)), rep(1,length(mIvec)), rep(0,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
# gImin and eImax
outll <- riskLL(mI_vec=mIvec, eI_vec=eImax, gI_vec=gImin, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(0,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), rep(1,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
# g1I and eImax
outll <- riskLL(mI_vec=mIvec, eI_vec=eImax, gI_vec=g1I, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(1,length(mIvec)), rep(0,length(mIvec)), rep(0,length(mIvec)), rep(1,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
# g3I and eImax
outll <- riskLL(mI_vec=mIvec, eI_vec=eImax, gI_vec=g3I, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(0,length(mIvec)), rep(1,length(mIvec)), rep(0,length(mIvec)), rep(1,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
# g9I and eImax
outll <- riskLL(mI_vec=mIvec, eI_vec=eImax, gI_vec=g9I, f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM=VMRR, R="P", r=0.5, lambdaR=8)
# Update output matrix
scz_results_ll <- rbind(scz_results_ll , cbind(rep(0,length(mIvec)), rep(0,length(mIvec)), rep(1,length(mIvec)), rep(1,length(mIvec)), mIvec, rep(1,length(mIvec)), outll))
colnames(scz_results_ll) <- c("G1I", "G3I","G9I","EI", "mI", "YR", "risk")
### Write output to file. Change file path and name
# write.csv(scz_results_ll, file="XXX/scz_results_loglinear.csv")
### Results can then be explored within R, for example
# using the plotting package ggplot2

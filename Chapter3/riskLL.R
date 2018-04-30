#####################################################################
### R code for thesis: Chapter 3 ####################################
#####################################################################
### Risk estimation, log linear model ###############################
### See Section 3.4 (main thesis) for details #######################
#####################################################################
### Additional internal functions: ##################################
#####################################################################
# 1. jointG_I_R_function: outputs the joint genotype probability
# for individual {I} and relative {R},p(G_I = g_I, G_R = g_R), 
# for all G_R values (0, 1, 2). 
# That is outputs vector: 
# c(p(G_I = g_I, G_R = 0), p(G_I = g_I, G_R = 1), p(G_I = g_I, G_R = 2))
#####################################################################
# INPUTS: 
### G_I = numeric. Observed genotype for individual I. 
# Assumes G ~ Bi(2, maf) (G is bi-allelic) and so can equal 0, 1, 2
### R = character. The relative type of individual {R}. 
# Can equal: "S" = sibling, "P" = parent, "O" = offspring, 
# "GP" = grandparent, "Av" = avuncular (aunt/ uncle) (by blood),
# "C" = cousin.
### maf = minor allele frequency (or risk allele frequency) for G.
# OUTPUT:
### vector of length 3 = 
# c(p(G_I = g_I, G_R = 0), p(G_I = g_I, G_R = 1), p(G_I = g_I, G_R = 2))
jointG_I_R_function <- function(G_I, R, maf){
	### Aim to calc p(G_I = g_I, G_R = g_R)
	pG <- c((1-maf)^2, 2*maf*(1-maf), maf^2)
	### Most will rely on conditioning on parents. G_P1, G_P2, G_O, p(G_P1), p(G_P2), p(G_O | G_P1, G_P2)
	po_mat <- cbind(c(rep(0, 3^2), rep(1, 3^2), rep(2, 3^2)), rep(c(rep(0, 3), rep(1, 3), rep(2, 3)), 3), rep(0:2,  3^2))
	po_mat <- cbind(po_mat, pG[po_mat[,1]+1], pG[po_mat[,2]+1], c(1, 0, 0, 0.5, 0.5, 0, 0, 1, 0, 0.5, 0.5, 0, 0.25, 0.5, 0.25, 0, 0.5, 0.5, 0, 1, 0, 0, 0.5, 0.5, 0,0,1))
	colnames(po_mat) <- c("G_P1", "G_P2", "G_O", "pG_P1", "pG_P2", "pG_O_P1P2")
	po_mat <- data.frame(po_mat)
	
	if(R == "S"){
		g_I_mat <- po_mat[po_mat$G_O == G_I,]
		sibling_jointG_f_internal <- function(g_I_mat, G_R, po_mat){
			g_R_mat <- po_mat[po_mat$G_O == G_R,]
			out <- sum(g_I_mat$pG_O_P1P2*g_R_mat$pG_O_P1P2*g_I_mat$pG_P1*g_I_mat$pG_P2)
			out <- sum(g_I_mat$pG_O_P1P2*g_R_mat$pG_O_P1P2*g_I_mat$pG_P1*g_I_mat$pG_P2)
			out
		}
		out_R <- mapply(sibling_jointG_f_internal, G_R=0:2, g_I_mat=list(g_I_mat), po_mat=list(po_mat))		
	}
	if(R == "P"){
		g_I_mat <- po_mat[po_mat$G_O == G_I,]
		parent_jointG_f_internal <- function(g_I_mat, G_R){
			g_I_R_mat <- g_I_mat[g_I_mat$G_P1 == G_R,]
			out <- sum(g_I_R_mat$pG_O_P1P2*g_I_R_mat$pG_P1*g_I_R_mat$pG_P2)
			out
		}
		out_R <- mapply(parent_jointG_f_internal, G_R=0:2, g_I_mat=list(g_I_mat))
	}
	if(R == "O"){
		g_I_mat <- po_mat[po_mat$G_P1 == G_I,]
		offspring_jointG_f_internal <- function(g_I_mat, G_R){
			g_I_R_mat <- g_I_mat[g_I_mat$G_O == G_R, ]
			out <- sum(g_I_R_mat$pG_O_P1P2*g_I_R_mat$pG_P1*g_I_R_mat$pG_P2)
			out
		}
		out_R <- mapply(offspring_jointG_f_internal, G_R=0:2, g_I_mat=list(g_I_mat))
	}
	if(R == "GP"){
		g_I_mat <- po_mat[po_mat$G_O == G_I,]
		grandparent_jointG_f_internal <- function(g_I_mat, G_R, po_mat){
			g_P1_R_mat <- po_mat[po_mat$G_P1 == G_R, ]
			g_P10_R_mat <- g_P1_R_mat[g_P1_R_mat$G_O == 0, ]
			g_P11_R_mat <- g_P1_R_mat[g_P1_R_mat$G_O == 1, ]
			g_P12_R_mat <- g_P1_R_mat[g_P1_R_mat$G_O == 2, ]
			pG_P1_0_R <- sum(g_P10_R_mat$pG_O_P1P2*g_P10_R_mat$pG_P1*g_P10_R_mat$pG_P2)
			pG_P1_1_R <- sum(g_P11_R_mat$pG_O_P1P2*g_P11_R_mat$pG_P1*g_P11_R_mat$pG_P2)
			pG_P1_2_R <- sum(g_P12_R_mat$pG_O_P1P2*g_P12_R_mat$pG_P1*g_P12_R_mat$pG_P2)
			g_I_mat$pG_P1[g_I_mat$G_P1 == 0] <- pG_P1_0_R
			g_I_mat$pG_P1[g_I_mat$G_P1 == 1] <- pG_P1_1_R
			g_I_mat$pG_P1[g_I_mat$G_P1 == 2] <- pG_P1_2_R
			out <- sum(g_I_mat$pG_O_P1P2*g_I_mat$pG_P1*g_I_mat$pG_P2)
			out
		}
		out_R <- mapply(grandparent_jointG_f_internal, G_R=0:2, g_I_mat=list(g_I_mat), po_mat=list(po_mat))
	}
	if(R == "Av"){
		g_I_mat <- po_mat[po_mat$G_O == G_I,]
		avuncular_jointG_f_internal <- function(g_I_mat, G_R, po_mat){
			g_R_mat <- po_mat[po_mat$G_O == G_R,]
			g_P1_0_mat <- po_mat[po_mat$G_O == 0,]
			g_P1_1_mat <- po_mat[po_mat$G_O == 1,]
			g_P1_2_mat <- po_mat[po_mat$G_O == 2,]
			pG_R_G_P1_0 <- sum(g_R_mat$pG_O_P1P2*g_R_mat$pG_P1*g_R_mat$pG_P2*g_P1_0_mat$pG_O_P1P2)
			pG_R_G_P1_1 <- sum(g_R_mat$pG_O_P1P2*g_R_mat$pG_P1*g_R_mat$pG_P2*g_P1_1_mat$pG_O_P1P2)
			pG_R_G_P1_2 <- sum(g_R_mat$pG_O_P1P2*g_R_mat$pG_P1*g_R_mat$pG_P2*g_P1_2_mat$pG_O_P1P2)
			g_I_mat$pG_P1[g_I_mat$G_P1 == 0] <- pG_R_G_P1_0
			g_I_mat$pG_P1[g_I_mat$G_P1 == 1] <- pG_R_G_P1_1
			g_I_mat$pG_P1[g_I_mat$G_P1 == 2] <- pG_R_G_P1_2
			out <- sum(g_I_mat$pG_O_P1P2*g_I_mat$pG_P1*g_I_mat$pG_P2)
			out
		}
		out_R <- mapply(avuncular_jointG_f_internal, G_R=0:2, g_I_mat=list(g_I_mat), po_mat=list(po_mat))
	}
	if(R == "C"){
		g_I_mat <- po_mat[po_mat$G_O == G_I,]
		cousin_jointG_f_internal <- function(g_I_mat, G_R, po_mat){
			g_R_mat <- po_mat[po_mat$G_O == G_R,]

			## P1 of individual I == 0:
			g_I_g_P1_0mat <- g_I_mat[g_I_mat$G_P1 == 0, ]
			g_P1I_0_mat <- po_mat[po_mat$G_O == 0, ]
				## And... parent 1 of R == 0
			g_P1R_0_mat <- po_mat[po_mat$G_O == 0, ]
			pG_P1I0_P1R0 <- sum(g_P1I_0_mat$pG_O_P1P2*g_P1I_0_mat$pG_P1*g_P1I_0_mat$pG_P2*g_P1R_0_mat$pG_O_P1P2)
				## And... parent 1 of R == 1
			g_P1R_1_mat <- po_mat[po_mat$G_O == 1, ]
			pG_P1I0_P1R1 <- sum(g_P1I_0_mat$pG_O_P1P2*g_P1I_0_mat$pG_P1*g_P1I_0_mat$pG_P2*g_P1R_1_mat$pG_O_P1P2)
				## And... parent 2 of R == 2
			g_P1R_2_mat <- po_mat[po_mat$G_O == 2, ]
			pG_P1I0_P1R2 <- sum(g_P1I_0_mat$pG_O_P1P2*g_P1I_0_mat$pG_P1*g_P1I_0_mat$pG_P2*g_P1R_2_mat$pG_O_P1P2)
			## P1 individual I = 0, updating pG_R
			g_R_P1I0_mat <- g_R_mat
			g_R_P1I0_mat$pG_P1[g_R_P1I0_mat$G_P1 == 0] <- pG_P1I0_P1R0
			g_R_P1I0_mat$pG_P1[g_R_P1I0_mat$G_P1 == 1] <- pG_P1I0_P1R1
			g_R_P1I0_mat$pG_P1[g_R_P1I0_mat$G_P1 == 2] <- pG_P1I0_P1R2
			### P1 individual I = 0, updating pG_I
			pG_P1I_0 <- sum(g_R_P1I0_mat$pG_O_P1P2*g_R_P1I0_mat$pG_P1*g_R_P1I0_mat$pG_P2)			
			g_I_mat$pG_P1[g_I_mat$G_P1 == 0] <- pG_P1I_0
			
			## P1 of individual I == 1:
			g_I_g_P1_1mat <- g_I_mat[g_I_mat$G_P1 == 1, ]
			g_P1I_1_mat <- po_mat[po_mat$G_O == 1, ]
				## And... parent 1 of R == 0
			g_P1R_0_mat <- po_mat[po_mat$G_O == 0, ]
			pG_P1I1_P1R0 <- sum(g_P1I_1_mat$pG_O_P1P2*g_P1I_1_mat$pG_P1*g_P1I_1_mat$pG_P2*g_P1R_0_mat$pG_O_P1P2)
				## And... parent 1 of R == 1
			g_P1R_1_mat <- po_mat[po_mat$G_O == 1, ]
			pG_P1I1_P1R1 <- sum(g_P1I_1_mat$pG_O_P1P2*g_P1I_1_mat$pG_P1*g_P1I_1_mat$pG_P2*g_P1R_1_mat$pG_O_P1P2)
				## And... parent 2 of R == 2
			g_P1R_2_mat <- po_mat[po_mat$G_O == 2, ]
			pG_P1I1_P1R2 <- sum(g_P1I_1_mat$pG_O_P1P2*g_P1I_1_mat$pG_P1*g_P1I_1_mat$pG_P2*g_P1R_2_mat$pG_O_P1P2)
			g_R_P1I1_mat <- g_R_mat
			g_R_P1I1_mat$pG_P1[g_R_P1I1_mat$G_P1 == 0] <- pG_P1I1_P1R0
			g_R_P1I1_mat$pG_P1[g_R_P1I1_mat$G_P1 == 1] <- pG_P1I1_P1R1
			g_R_P1I1_mat$pG_P1[g_R_P1I1_mat$G_P1 == 2] <- pG_P1I1_P1R2
			pG_P1I_1 <- sum(g_R_P1I1_mat$pG_O_P1P2*g_R_P1I1_mat$pG_P1*g_R_P1I1_mat$pG_P2)
			g_I_mat$pG_P1[g_I_mat$G_P1 == 1] <- pG_P1I_1
			
			## P1 of individual I == 2:
			g_I_g_P1_2mat <- g_I_mat[g_I_mat$G_P1 == 2, ]
			g_P1I_2_mat <- po_mat[po_mat$G_O == 2, ]
				## And... parent 1 of R == 0
			g_P1R_0_mat <- po_mat[po_mat$G_O == 0, ]
			pG_P1I2_P1R0 <- sum(g_P1I_2_mat$pG_O_P1P2*g_P1I_2_mat$pG_P1*g_P1I_2_mat$pG_P2*g_P1R_0_mat$pG_O_P1P2)
				## And... parent 1 of R == 1
			g_P1R_1_mat <- po_mat[po_mat$G_O == 1, ]
			pG_P1I2_P1R1 <- sum(g_P1I_2_mat$pG_O_P1P2*g_P1I_2_mat$pG_P1*g_P1I_2_mat$pG_P2*g_P1R_1_mat$pG_O_P1P2)
				## And... parent 2 of R == 2
			g_P1R_2_mat <- po_mat[po_mat$G_O == 2, ]
			pG_P1I2_P1R2 <- sum(g_P1I_2_mat$pG_O_P1P2*g_P1I_2_mat$pG_P1*g_P1I_2_mat$pG_P2*g_P1R_2_mat$pG_O_P1P2)
			g_R_P1I2_mat <- g_R_mat
			g_R_P1I2_mat$pG_P1[g_R_P1I2_mat$G_P1 == 0] <- pG_P1I2_P1R0
			g_R_P1I2_mat$pG_P1[g_R_P1I2_mat$G_P1 == 1] <- pG_P1I2_P1R1
			g_R_P1I2_mat$pG_P1[g_R_P1I2_mat$G_P1 == 2] <- pG_P1I2_P1R2
			pG_P1I_2 <- sum(g_R_P1I2_mat$pG_O_P1P2*g_R_P1I2_mat$pG_P1*g_R_P1I2_mat$pG_P2)
			g_I_mat$pG_P1[g_I_mat$G_P1 == 2] <- pG_P1I_2
			
			out <- sum(g_I_mat$pG_O_P1P2*g_I_mat$pG_P1*g_I_mat$pG_P2)
			out
		}
		out_R <- mapply(cousin_jointG_f_internal, G_R=0:2, g_I_mat=list(g_I_mat), po_mat=list(po_mat))
	}
	out_R
}
#####################################################################
### MAIN FUNCTION ###################################################
#####################################################################
### Estimates the risk of disease for an individual {I} given their
# observed polygenic risk score (MI = mI), their observed set of 
# major susceptibility risk genotypes (GI=gI), their observed set of
# environmental risk factors, and their family history of disease.
# Function riskLL can output multiple risks for a vector of inputted
# polygenic risk score values (but the same gI and eI values).
# There is no family history variable- the outputted risk is for an
# observed AFFECTED relative
#####################################################################
# INPUTS
### mI_vec = vector. The polygenic risk score of the individuals for whom
# we wish to estimate risk (all denoted by {I})
### eI_vec = vector. Length = number of included environmental risk
# factors. Vector containing the observed environmental risk factors 
# for individual {I} *Reference category = 1*
### gI_vec = vector. Length = number of included major risk loci. Vector
# containing the observed genotypes at the loci for individual {I}
### f_vec = vector. Length = number of included major risk loci. 
# Vector containing the minor (or risk) allele frequencies for the 
# included major risk loci. Order of major risk loci must be the same as
# in g_Ivec and g_Rvec
### RR_Gmat = matrix. Dim = c(Q, 3); where Q is the number of major risk
# loci (rows) and 3 is the number of possible genotypes. This assumes that
# G_q is bi-allelic, q = 1, ..., Q
### pEmat = matrix. dim = c(n_E, max(levelsE)); n_E = number of included 
# environmental risk factors, max(levelsE) = maximum number of categories 
# across all included environmental risk factors. This matrix contains the 
# probability density function (PDF) for each risk factor. Each row 
# contains the PDF for the i^th risk factor and therefore sums to 1. Column
# 1 contains the probability for the reference categories
### RR_Emat = matrix. Dimensions depends on the number of categories within 
# each of the environmental risk factors. Number of rows = S (number of 
# environmental risk factors). Number of columns = maximum number of categories
# across all environmental risk factors
### K = numeric. The population prevalence of disease
### VM = the variability of the polygenic risk score, calculated using the
# estimated relative risk of disease for each included SNP  
### R = character. The relative type of individual {R}. 
# Can equal: "S" = sibling, "P" = parent, "O" = offspring, 
# "GP" = grandparent, "Av" = avuncular (aunt/ uncle) (by blood),
# "C" = cousin.
### r = numeric. Coefficient of relatedness between {I} and {R}
### lambdaR = recurrence risk ratio. Defined as:
# p(YI = 1 | YR = 1)/K
#####################################################################
# OUTPUT
###
riskLL <- function(mI_vec, eI_vec, gI_vec, f_vec, RR_Gmat, pEmat, RR_Emat, K, VM, R="P", r=0.5, lambdaR){

### Calculating lambda_RM:
lambda_RM_function <- function(VM, r){
	exp(r*VM)
}
lambda_RM <- lambda_RM_function(VM=VM, r=r)
### Calculating lambda_RG:
# We assume all major loci are independent:
# Assume additive within locus...
# test: f <- f1, RR_vec <- c(1, RR_G1c, exp(log(RR_G1c)*2))
lambda_RG_g_function <- function(f, RR_vec, levels_g=3, R){
	G <- matrix(0:(levels_g -1))
	probG <- c((1-f)^2, 2*f*(1-f), f^2)
	denom <- sum(probG*RR_vec)^2
	jointprobs <- apply(G, 1, jointG_I_R_function, maf=f, R=R)
	valg <- matrix(0, nrow=3, ncol=3)
	for(i in 1:3){
		for(j in 1:3){
			valg[i,j] <- RR_vec[i]*RR_vec[j]
		}
	}
	numer <- sum(valg*jointprobs)
	numer/denom
}
lambda_RG <- 1
for(i in 1:length(f_vec)){
	lambda_RG <- lambda_RG*lambda_RG_g_function(f = f_vec[i], RR_vec = RR_Gmat[i,], R=R)
}
### And so lambda_RU is...
lambda_RU <- lambdaR/(lambda_RG*lambda_RM)
### m_I contribution...
mIc <- exp(sqrt(VM)*mI_vec)/exp(0.5*VM)
### g_I contribution:
numer_g_I <- 1
denom_g_I <- 1
for(i in 1:length(f_vec)){
	gi <- gI_vec[i]
	numer_g_I <- numer_g_I*RR_Gmat[i, (gi+1)]
	fi <- f_vec[i]
	denom_g_I <- denom_g_I*sum(c((1-fi)^2, 2*fi*(1-fi), fi^2)*RR_Gmat[i,])	
}
gIc <- numer_g_I/denom_g_I
### environmental contribution...
eIc <- 1
for(i in 1:dim(pEmat)[1]){
	eIc <- eIc*(RR_Emat[i, eI_vec[i]]/sum(pEmat[i,]*RR_Emat[i,]))
}

out <- K*lambda_RU*mIc*gIc*eIc
out
}
#####################################################################
### EXAMPLE OF USAGE ################################################
#####################################################################
### Input definitions taken from Chapter 4 Schizophrenia example 
# (without age and sex). Please see the example R-script located in: 
# https://github.com/alexgillett/THESIS/Chapter4/
# for details
K <- 0.01
H2 <- 0.81
### Recurrence risk ratio for a 1st degree relative
lambdaR1 <- 8
maf.vec <- c(3.964324e-05, 1.879158e-04, 8.171624e-05, 5.558299e-05, 3.568971e-04, 2.710037e-05, 2.825040e-05, 3.923439e-04, 1.247620e-03, 1.515856e-04)
levelsE <- c(3,2,3,7,2)
pEmat <- matrix(0, nrow=5, ncol=7)
pEmat[1,1:3] <- c(0.700, 0.150, 0.150)
pEmat[2,1:2] <- c(0.924, 0.076)
pEmat[3,1:3] <- c(0.250, 0.250, 0.5000)
pEmat[4,1:7] <- c(0.342, 0.204, 0.252, 0.123, 0.052, 0.019, 0.008)
pEmat[5,1:2] <- c(0.730, 0.270)

RR_Gmat <- matrix(1, nrow=10, ncol=3)
RR_Gmat[1, 2:3] <- c(40.6628807, 1653.4698667)
RR_Gmat[2, 2:3] <- c(8.6735194, 75.2299394)
RR_Gmat[3, 2:3] <- c(12.7012768, 161.3224313)
RR_Gmat[4, 2:3] <- c(13.6149950, 185.3680899)
RR_Gmat[5, 2:3] <- c(3.6966913, 13.6655266)
RR_Gmat[6, 2:3] <- c(17.2265579, 296.7542977)
RR_Gmat[7, 2:3] <- c(13.9891828, 195.6972344)
RR_Gmat[8, 2:3] <- c(0.1512868, 0.0228877)
RR_Gmat[9, 2:3] <- c(1.7857421, 3.1888748)
RR_Gmat[10, 2:3] <- c(3.4146951, 11.6601430)

RR_Emat <- matrix(1, nrow=5, ncol=max(levelsE))
RR_Emat[1, 2:levelsE[1]] <- c(1.405642, 2.743074)
RR_Emat[2, 2:levelsE[2]] <- c(2.273057)
RR_Emat[3, 2:levelsE[3]] <- c(1.495378, 1.987713)
RR_Emat[4, 2:levelsE[4]] <- c(1.059402, 1.059402, 1.12862, 1.217481, 1.207615, 1.64976)
RR_Emat[5, 2:levelsE[5]] <- c(2.763157)

### Polygenic risk score (prs) inputs:
ORprs <- c(1.072, 1.072, 1.069, 0.934, 1.078, 0.914, 1.075, 1.063, 0.904, 1.086, 0.857, 0.944, 0.929, 0.937, 0.929, 1.081, 0.909, 0.939, 0.930, 1.101, 0.940, 1.085, 1.071, 0.941, 0.933, 1.064, 0.857, 1.065, 0.934, 0.926, 0.922, 1.101, 1.076, 1.058, 0.942, 1.205, 0.942, 1.068, 0.849, 0.922, 1.101, 0.904, 1.068, 1.061, 1.083, 0.924, 1.066, 1.073, 0.908, 0.941, 1.079, 0.919, 1.087, 1.069, 1.125, 1.064, 0.906, 1.068, 0.941, 0.927, 1.066, 0.926, 1.076, 1.091, 0.941, 0.943, 1.060, 1.068, 0.915, 0.933, 1.073, 0.939, 1.088, 1.074, 1.067, 1.060, 0.923, 0.922, 1.067, 1.077, 1.073, 0.920, 1.081, 1.071, 0.939, 0.934, 0.931, 1.071, 0.928, 0.937, 0.930, 1.087, 0.947, 1.090, 0.938, 0.846, 1.317, 0.914, 1.069, 1.155, 0.843, 0.912, 0.846, 1.093, 1.188, 1.076, 0.902, 0.873)
fcontrol_prs <- c(0.527, 0.301, 0.296, 0.677, 0.358, 0.164, 0.184, 0.685, 0.101, 0.163, 0.961, 0.593, 0.458, 0.337, 0.643, 0.805, 0.754, 0.326, 0.552, 0.156, 0.48, 0.324, 0.529, 0.615, 0.449, 0.314, 0.922, 0.47, 0.761, 0.802, 0.532, 0.883, 0.213, 0.475, 0.621, 0.85, 0.368, 0.505, 0.0476, 0.423, 0.0959, 0.123, 0.332, 0.628, 0.647, 0.152, 0.642, 0.219, 0.116, 0.6, 0.174, 0.803, 0.424, 0.337, 0.889, 0.499, 0.85, 0.311, 0.334, 0.314, 0.322, 0.335, 0.46, 0.797, 0.337, 0.366, 0.646, 0.474, 0.741, 0.771, 0.249, 0.418, 0.287, 0.274, 0.52, 0.465, 0.257, 0.562, 0.281, 0.223, 0.51, 0.859, 0.162, 0.627, 0.614, 0.628, 0.769, 0.766, 0.322, 0.754, 0.594, 0.232, 0.753, 0.831, 0.759, 0.964, 0.0191, 0.208, 0.252, 0.046, 0.0756, 0.624, 0.0708, 0.725, 0.0292, 0.81, 0.913, 0.0797)
VMRR <- sum((log(ORprs)^2)*2*fcontrol_prs*(1-fcontrol_prs))
mIvec <- round(quantile(rnorm(100000, mean=0, sd=sqrt(VMRR)), probs = c(0.25, 0.5, 0.75)), digits=2)

### 'Risk' estimates:
### low risk environmental risk factors, plus carrying risk genotype at major risk locus 1
riskLL(mI_vec = mIvec, eI_vec=rep(1,5), gI_vec=c(1, rep(0,9)), f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM= VMRR, R="P", r=0.5, lambdaR=lambdaR1)
### high risk environmental risk factors, no risk genotypes at major genetic risk loci
riskLL(mI_vec = mIvec, eI_vec=levelsE, gI_vec=rep(0, 10), f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM= VMRR, R="P", r=0.5, lambdaR=lambdaR1)
### high risk environmental risk factors, plus carrying risk genotype at major risk locus 1
riskLL(mI_vec = mIvec, eI_vec=levelsE, gI_vec=c(1, rep(0,9)), f_vec=maf.vec, RR_Gmat=RR_Gmat, pEmat=pEmat, RR_Emat=RR_Emat, K=K, VM= VMRR, R="P", r=0.5, lambdaR=lambdaR1)
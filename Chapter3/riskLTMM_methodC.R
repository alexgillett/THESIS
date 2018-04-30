#####################################################################
### R code for thesis: Chapter 3 ####################################
#####################################################################
### Risk estimation, liability threshold mixture model ##############
### Method C = a method to reduce computational time compared to  ###
### Method A (the original method) ##################################
### Idea: a Pearson-Aitken approach to reduce computation time  #####
### compared to Method A ############################################
### See Section 3.3.5 for details (main thesis) #####################
### EXACT method (multivariate integration) only ####################
#####################################################################
### Required libraries: #############################################
#####################################################################
# Solving (non-)linear equations using iterative methods 
library(BB)
# Multivariate integration for normally distributed variables
library(mvtnorm)
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
# 2. mean_mI_f: calculates the mean vector for the bivariate normal 
# distribution of:
# [L_{I}|M_{I}=m_{I},G_{I}=g_{I},E_{I}=e_{I}, L_{R}|M_{I}=m_{I},G_{R}=g_{R},E_{R}=e_{R}]^{T}
#####################################################################
# INPUTS:
### tau = vector. Contains the effect sizes for the major risk loci and 
# environmental risk factors. Reference categories are assumed to be:
# G_q = 0 for the Q included major risk loci, and,
# E_s = 1 for the S included environmental risk factors.
# Reference categories have effect sizes set to 0 and are not included 
# in tau.
# Major risk are assumed to be bi-allelic.
# The number of categories/ levels for environmental risk factor s is:
# n_E_s.
# The length of tau is therefore:
# (Q*2)+(n_E_1 - 1)+(n_E_2 - 1)+ ... +(n_E_S - 1)
### x_I = vector. length = length(tau). It is the `design' vector for 
# individual {I} corresponding to effect size vector tau
### x_R = vector. length = length(tau). It is the `design' vector for 
# relative {R} corresponding to effect size vector tau
### r = numeric. Coefficient of relatedness between {I} and {R}
### m_I = numeric. Observed polygenic risk score for individual {I}
mean_mI_f <- function(tau, x_I, x_R, r, mI){
	### Generate mean vector for joint conditional distribution L_I, L_R given M_I=m_I, G_I=g_I, G_R=g_R, E_I=e_I, E_R=e_R
	c(t(matrix(tau))%*%matrix(x_I) + mI, t(matrix(tau))%*%matrix(x_R) + r*mI)
}
#####################################################################
# 3. sigma_mI_f: calculates the covariance matrix for the bivariate  
# normal distribution of:
# [L_{I}|M_{I}=m_{I},G_{I}=g_{I},E_{I}=e_{I}, L_{R}|M_{I}=m_{I},G_{R}=g_{R},E_{R}=e_{R}]^{T}
#####################################################################
# INPUTS:
### VM = numeric. The variance on the liability scale of the 
# polygenic risk score
### VhG = numeric. Variance component of the major risk loci on the 
# liability scale
### VhE = numeric. Variance component of the environmental risk factors 
# on the liability scale
### H2 = numeric. Broad-sense heritability on the liability scale.
### h2 = numeric. Narrow-sense heritability on the liability scale.
### VAhG = numeric. Additive variance component of the major risk loci
# on the liability scale. Will equal VhG if no dominance
### VDhG = numeric. Dominance variance component of the major risk loci
# on the liability scale
### r = numeric. The coefficient of relatedness for {I} and {R}
### theta = numeric. Coancestry coefficient for {I} and {R}
sigma_mI_f <- function(VM, VhG, VhE, H2, h2=H2, VAhG=VhG, VDhG=0, r=0.5, theta=0.25){
	### Generate covariance matrix for joint conditional distribution L_I, L_R given M_I=m_I, G_I=g_I, G_R=g_R, E_I=e_I, E_R=e_R
	sigma.out <- matrix(r*(h2-VAhG - VM) + theta*(H2-h2-VDhG), nrow=2, ncol=2)
	sigma.out[1,1] <- 1 - VM - VhG - VhE
	sigma.out[2,2] <- 1 - VhG - VhE - (r^2)*VM
	sigma.out
}
#####################################################################
# 4. mu_LR_gI_f: calculates the mean of:
# L_{R}|G_{I}=g_{I}
#####################################################################
# INPUTS:
### g_Ivec = vector. Length = number of included major risk loci. Vector
# containing the observed genotypes at the loci for individual {I}
### tau = vector. Contains the effect sizes on the liability
# scale for: 1. the genetic major risk loci and, 2. the environmental
# risk factors. Please see description in Function 2. mean_mI_f for 
# details 
### maf.vec = vector. Length = number of included major risk loci. 
# Vector containing the minor (or risk) allele frequencies for the 
# included major risk loci. Order of major risk loci must be the same as
# in g_Ivec
### levelsE = vector. Length = number of included environmental risk
# factors. Vector containing the number of categories (a.k.a levels)
# for each included environmental risk factor. Order of the environmental
# risk factors must be the same as e_Ivec and e_Rvec
### pEmat = matrix. dim = c(n_E, max(levelsE)); n_E = number of included 
# environmental risk factors, max(levelsE) = maximum number of categories 
# across all included environmental risk factors. This matrix contains the 
# probability density function (PDF) for each risk factor. Each row 
# contains the PDF for the i^th risk factor and therefore sums to 1. Column
# 1 contains the probability for the reference categories
### R = character. The relative type of individual {R}. 
# Can equal: "S" = sibling, "P" = parent (DEFAULT), "O" = offspring, 
# "GP" = grandparent, "Av" = avuncular (aunt/ uncle) (by blood),
# "C" = cousin.
### combosR_E = matrix. Number of columns = number of environmental risk factors
# included. Contains every possible combination of risk factors that could be 
# observed. *Reference category coded as 1 for all risk factors*. Can be created
# using the function expand.grid from the base package. Please see
# R code for Chapter 4 examples for details
### combosR_G = matrix. Number of columns = number of major risk loci to be
# included (Q). Contains every possible combination of the major risk loci that 
# could be observed. Here we assume that G_q = bi-alleleic (q = 1, ..., Q). 
# Therefore, number of rows = 3^Q
# Can be created using expand.grid function
#####################################################################
# OUTPUT:
### output is the expected value of L_{R} given the observed major risk 
# loci for individual {I}. This means that for the final risk estimation
# here we are assuming that the environmental risk factors observed for 
# individual {I} are independent from the liability to disease for the 
# relative {R}
mu_LR_gI_f <- function(gIvec, tau, maf.vec, levelsE, pEmat, R="P", combosR_E, combosR_G){
	n_E <- length(levelsE)
	n_G <- length(maf.vec)
	n_tau_G <- n_G*2
	n_tau_E <- sum(levelsE) - length(levelsE)
	mean <- 0
	#pGmat <- matrix(0, nrow=n_G, ncol=3)
	maf.mat <- matrix(maf.vec)
	pGmat <- t(apply(maf.mat,1,function(x){c((1-x)^2, (2*x*(1-x)), x^2)}))
	#gImat <- matrix(gIvec)
	pGI <- rep(0, n_G)
	for(i in 1:n_G){
		pGI[i] <- pGmat[i, gIvec[i] + 1]
	}
	pGI <- prod(pGI)
n_tau_G <- n_G*2
n_tau_E <- sum(levelsE) - length(levelsE)
for(i in 1:dim(combosR_G)[1]){
	g_Rvec <- unlist(combosR_G[i,])
	names(g_Rvec)<-NULL
	pGRGI <- rep(0, n_G)
	for(k in 1:n_G){
		grk <- g_Rvec[k]
		pwork <- jointG_I_R_function(G_I = gIvec[k], R=R, maf=maf.vec[k])
		pGRGI[k] <- pwork[grk + 1]	
	}
	pGRGI <- prod(pGRGI)
	x_Rv <- rep(0, length(tau))
	for(k in 1:n_G){
		x_Rv[(2*k) - 1] <- as.numeric(g_Rvec[k] == 1)
		x_Rv[(2*k)] <- as.numeric(g_Rvec[k] == 2)
	}
	for(j in 1:dim(combosR_E)[1]){
		varf <- (2*n_G)
		e_Rvec <- unlist(combosR_E[j,])
		names(e_Rvec)<- NULL
		for(k in 1:n_E){
			n_tau_Ek <- levelsE[k] - 1
			for(b in 1:n_tau_Ek){
				x_Rv[varf + b] <- as.numeric(e_Rvec[k] == (b+1))
			}
			varf <- varf + n_tau_Ek
		}
		hgrer <- as.vector(t(matrix(tau))%*%matrix(x_Rv))
		pER <- rep(0, n_E)
		for(k in 1:n_E){
			pER[k] <- pEmat[k, e_Rvec[k]]
		}
		pER <- prod(pER)
		meaniter <- (pGRGI/pGI)*pER*hgrer
		mean <- mean+meaniter	
	}
}
mean	
}
#####################################################################
# 5. sigma2_LR_gI_f: calculates the variance of:
# L_{R}|G_{I}=g_{I}
#####################################################################
# INPUTS:
### g_Ivec = vector. Length = number of included major risk loci. Vector
# containing the observed genotypes at the loci for individual {I}
### tau = vector. Contains the effect sizes on the liability
# scale for: 1. the genetic major risk loci and, 2. the environmental
# risk factors. Please see description in Function 2. mean_mI_f for 
# details 
### maf.vec = vector. Length = number of included major risk loci. 
# Vector containing the minor (or risk) allele frequencies for the 
# included major risk loci. Order of major risk loci must be the same as
# in g_Ivec
### levelsE = vector. Length = number of included environmental risk
# factors. Vector containing the number of categories (a.k.a levels)
# for each included environmental risk factor. Order of the environmental
# risk factors must be the same as e_Ivec and e_Rvec
### pEmat = matrix. dim = c(n_E, max(levelsE)); n_E = number of included 
# environmental risk factors, max(levelsE) = maximum number of categories 
# across all included environmental risk factors. This matrix contains the 
# probability density function (PDF) for each risk factor. Each row 
# contains the PDF for the i^th risk factor and therefore sums to 1. Column
# 1 contains the probability for the reference categories
### R = character. The relative type of individual {R}. 
# Can equal: "S" = sibling, "P" = parent (DEFAULT), "O" = offspring, 
# "GP" = grandparent, "Av" = avuncular (aunt/ uncle) (by blood),
# "C" = cousin.
### mu_LR_gI = numeric. Calculated using function number 4 above.
### VhG = numeric. Variance component of the major risk loci on the 
# liability scale
### VhE = numeric. Variance component of the environmental risk factors 
# on the liability scale
### combosR_E = matrix. Number of columns = number of environmental risk factors
# included. Contains every possible combination of risk factors that could be 
# observed. *Reference category coded as 1 for all risk factors*. Can be created
# using the function expand.grid from the base package. Please see
# R code for Chapter 4 examples for details
### combosR_G = matrix. Number of columns = number of major risk loci to be
# included (Q). Contains every possible combination of the major risk loci that 
# could be observed. Here we assume that G_q = bi-alleleic (q = 1, ..., Q). 
# Therefore, number of rows = 3^Q
# Can be created using expand.grid function
#####################################################################
# OUTPUT:
### output is the variance of L_{R} given the observed major risk 
# loci for individual {I}. This means that for the final risk estimation
# here we are assuming that the environmental risk factors observed for 
# individual {I} are independent from the liability to disease for the 
# relative {R}
sigma2_LR_gI_f <- function(gIvec, tau, maf.vec, levelsE, pEmat, R="P", mu_LR_gI, VhG, VhE, combosR_E, combosR_G){
	n_E <- length(levelsE)
	n_G <- length(maf.vec)
	n_tau_G <- n_G*2
	n_tau_E <- sum(levelsE) - length(levelsE)
	maf.mat <- matrix(maf.vec)
	pGmat <- t(apply(maf.mat,1,function(x){c((1-x)^2, (2*x*(1-x)), x^2)}))
	#gImat <- matrix(gIvec)
	pGI <- rep(0, n_G)
	for(i in 1:n_G){
		pGI[i] <- pGmat[i, gIvec[i] + 1]
	}
	pGI <- prod(pGI)
	
	#sigma2GR <- 1 - VhG
	sigma2LRgI <- 1 - VhG - VhE
	
for(i in 1:dim(combosR_G)[1]){
	g_Rvec <- unlist(combosR_G[i,])
	names(g_Rvec)<-NULL
	pGRGI <- rep(0, n_G)
	for(k in 1:n_G){
		grk <- g_Rvec[k]
		pwork <- jointG_I_R_function(G_I = gIvec[k], R=R, maf=maf.vec[k])
		pGRGI[k] <- pwork[grk + 1]	
	}
	pGRGI <- prod(pGRGI)
	x_Rv <- rep(0, length(tau))
	for(k in 1:n_G){
		x_Rv[(2*k) - 1] <- as.numeric(g_Rvec[k] == 1)
		x_Rv[(2*k)] <- as.numeric(g_Rvec[k] == 2)
	}
	for(j in 1:dim(combosR_E)[1]){
		varf <- (2*n_G)
		e_Rvec <- unlist(combosR_E[j,])
		names(e_Rvec)<- NULL
		for(k in 1:n_E){
			n_tau_Ek <- levelsE[k] - 1
			for(b in 1:n_tau_Ek){
				x_Rv[varf + b] <- as.numeric(e_Rvec[k] == (b+1))
			}
			varf <- varf + n_tau_Ek
		}
		hgrer <- as.vector(t(matrix(tau))%*%matrix(x_Rv))
		pER <- rep(0, n_E)
		for(k in 1:n_E){
			pER[k] <- pEmat[k, e_Rvec[k]]
		}
		pER <- prod(pER)
		siter <- (pGRGI/pGI)*pER*((hgrer - mu_LR_gI)^2)
		sigma2LRgI <- sigma2LRgI +siter	
	}
}
sigma2LRgI	
}
### We note here that functions 4 and 5 (mu_LR_gI_f and sigma2_LR_gI_f) 
### can take a long time to run
#####################################################################
### MAIN RISK ESTIMATION FUNCTION ##################################
#####################################################################
### This code calculates the risk of disease for an individual {I} 
# given their observed polygenic risk score, major genetic risk loci, 
# and environmental risk factors, if they: 1. have a relative {R} observed to 
# be affected or, 2. have a relative {R} observed to be unaffected with
# the disease of interest.
### The polygenic risk score variable is numeric here. 
# To estimate risk for multiple polygenic risk score values a loop or
# apply function ca be used.
### Output presented is for Scenario B, Chapter 3, from the thesis by Alexandra
# Gillett.
#####################################################################
# INPUTS:
### mu_LI_gIeI = numeric. Expected liability to disease for individual
# {I} given their observed major genetic suspectibility loci and their
# observed environmental risk factors. Can be found using function
# mean_mI_f with mI = 0. It will be the first vector entry: output[1]
### mu_LR_gI = numeric. Expected liability disease for relative {R}
# given the observed major genetic susceptibility loci of individual
# {I}. Found using function mu_LR_gI_f
### sigma2_LR_gI = numeric. Variance of the liability to disease for
# relative {R} given the observed major genetic susceptibility loci of 
# individual {I}. Found using function sigma2_LR_gI_f
### mI = numeric. The polygenic risk score for individual {I}
### r = numeric. Coefficient of relatedness between {I} and {R}
### theta = numeric. Coefficient of co-ancestry between {I} and {R} 
### VM = numeric. Variability in liability to disease attributable to M,
# the polygenic risk score (Var[M] = VM)
### H2 = numeric. Broad-sense heritability
### h2 = numeric. Narrow-sense heritability
### VhG = numeric. The variability in liability to disease attributable to
# the major risk loci
### VhE = numeric. The variability in liability to disease attributable to
# the environmental risk factors
### VAhG = numeric. The variability in liability to disease attributable to
# the \emph{additive} contribution of the major risk loci
### VDhG = numeric. The variability in liability to disease attributable to
# the \emph{quasi-dominant} contribution of the major risk loci
### YR = numeric. Observed disease status of relative {R}. Can be 0 (unaffected)
# or 1 (affected). Default = 1
### D_threshold = numeric. Disease threshold
#####################################################################
# OUTPUT:
### Function outputs the risk of disease using method C. This main
# function itself should be very quick. However, the functions to calculate
# mu_LR_gI and sigma2_LR_gI will take time
riskLTMM_methodC <- function(mu_LI_gIeI, mu_LR_gI, sigma2_LR_gI, mI, r, theta, VM, H2, h2, VhG, VhE, VAhG, VDhG, YR=1, D_threshold){
	mu_C <- c(mu_LI_gIeI + mI, mu_LR_gI + (r*mI))
	sigma2_C <- matrix((r*(h2 - VAhG - VM) + theta*(H2-h2 - VDhG)), ncol=2, nrow=2)
	sigma2_C[1,1] <- 1 - VhG - VhE - VM
	sigma2_C[2,2] <- sigma2_LR_gI - ((r^2)*VM)
	alpha <- (D_threshold - mu_C[2])/sqrt(sigma2_C[2,2])
	if(YR == 1){
		denom <- 1 - pnorm(alpha)
		lower.int <- rep(D_threshold, 2)
		upper.int <- rep(Inf, 2)
		numer <- pmvnorm(lower=lower.int, upper=upper.int, mean = mu_C, sigma=sigma2_C)[1]
	}else{
		denom <- pnorm(alpha)
		lower.int <- c(D_threshold, -Inf)
		upper.int <- c(Inf, D_threshold)
		numer <- pmvnorm(lower=lower.int, upper=upper.int, mean = mu_C, sigma=sigma2_C)[1]
	}
	numer/denom
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
maf.vec <- c(3.964324e-05, 1.879158e-04, 8.171624e-05, 5.558299e-05, 3.568971e-04, 2.710037e-05, 2.825040e-05, 3.923439e-04, 1.247620e-03, 1.515856e-04)
levelsE <- c(3,2,3,7,2)
pEmat <- matrix(0, nrow=5, ncol=7)
pEmat[1,1:3] <- c(0.700, 0.150, 0.150)
pEmat[2,1:2] <- c(0.924, 0.076)
pEmat[3,1:3] <- c(0.250, 0.250, 0.5000)
pEmat[4,1:7] <- c(0.342, 0.204, 0.252, 0.123, 0.052, 0.019, 0.008)
pEmat[5,1:2] <- c(0.730, 0.270)
tau <- c(2.08765243, 4.17530487, 0.96454210, 1.92908420, 1.18513708, 2.37027416, 1.22811566, 2.45623132, 0.53911125, 1.07822251, 1.38076307, 2.76152614, 1.24530173, 2.49060346, -0.63870613, -1.27741226, 0.22611380, 0.45222760, 0.50321049, 1.00642097, 0.12484618, 0.38691427, 0.32064091, 0.14512302, 0.25144700, 0.02155491, 0.02155491, 0.04537918, 0.07416491, 0.07106093, 0.19254779, 0.38256397)
D_threshold <- 2.712044
VM <- 0.07
VhG <- 0.002014081
VhE <- 0.06625238
combosR_E <- expand.grid(1:(levelsE[1]), 1:(levelsE[2]), 1:(levelsE[3]), 1:(levelsE[4]), 1:(levelsE[5]))
n_G <- length(maf.vec)
df_G <- NULL
for(i in 1:n_G){
	df_G <- cbind(df_G, 0:2)
}
df_G <- data.frame(df_G)
combosR_G <- expand.grid(df_G)
### Estimates the risk of disease for an individual with low risk (reference category)
### environmental risk factors and carrying one major risk genotype at risk locus 1
mu_LR_gI <- mu_LR_gI_f(gIvec=c(1, rep(0, 9)), tau=tau, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, R="P", combosR_E=combosR_E, combosR_G=combosR_G)
sigma2_LR_gI <- sigma2_LR_gI_f(gIvec=c(1, rep(0, 9)), tau=tau, maf.vec=maf.vec, levelsE=levelsE, pEmat=pEmat, R="P", mu_LR_gI=mu_LR_gI, VhG=VhG, VhE=VhE, combosR_E=combosR_E, combosR_G=combosR_G)
mu_LI_gIeI <- mean_mI_f(tau=tau, x_I=c(1,0,rep(0, 2*9), rep(0,2), 0, rep(0, 2), rep(0,6), 0), x_R=c(rep(0, 20), rep(0,12)), r=0.5, mI=0)
mu_LI_gIeI <- mu_LI_gIeI[1]

example.output <- riskLTMM_methodC(mu_LI_gIeI=mu_LI_gIeI, mu_LR_gI=mu_LR_gI, sigma2_LR_gI=sigma2_LR_gI, mI=0, r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=1, D_threshold=D_threshold)

### multiple mI:
mIvec <- c(-0.18, 0, 0.18)
example.output.mI <- rep(0, length(mIvec))
for(i in 1:length(mIvec)){
	example.output.mI[i] <- riskLTMM_methodC(mu_LI_gIeI=mu_LI_gIeI, mu_LR_gI=mu_LR_gI, sigma2_LR_gI=sigma2_LR_gI, mI=mIvec[i], r=0.5, theta=0, VM=VM, H2=H2, h2=H2, VhG=VhG, VhE=VhE, VAhG=VhG, VDhG=0, YR=1, D_threshold=D_threshold)
}
# ================================== CODE METADATA =============================
# AUTHOR: Fernando Colchero
# DATE CREATED: 2024-08-26
# DATE MODIFIED: 
# DESCRIPTION: Functions for reproducibility code for Staerk et al. (2024)
# NOTES: Uses packages:
#           - paramDemo (parametric and non-parametric demographic functions)
#           - BaSTA (Bayesian survival trajectory analysis)
#           - snowfall (parallel computing)
#           - BayesPGLS (Bayesian Phylogenetic Generalized Least Squares)
#           - phytools ()
#           - caper ()
# ================================ CODE START ==================================
# ==================== #
# ==== LIBRARIES: ====
# ==================== #
# --------------------------- #
# ---- Install packages: ----
# --------------------------- #
# Vector of installed packages:
instPacks <- installed.packages()[, 1]

# Install and load 'devtools':
if (!"devtools" %in% instPacks) {
  install.packages("devtools")
}
library(devtools)

# snowfall for parallel computing:
if (!"snowfall" %in% instPacks) {
  install.packages("snowfall")
}
library(snowfall)

# Install and load latest version of 'BaSTA':
if (!"BaSTA" %in% instPacks) {
  install_git("https://github.com/fercol/basta2.0", subdir = "pkg/")
}
library(BaSTA)

# Install and load 'paramDemo':
if (!"paramDemo" %in% instPacks) {
  install.packages("paramDemo")
}
library(paramDemo)

# Install and load 'BayesPGLS':
if (!"BayesPGLS" %in% instPacks) {
  install_git("https://github.com/fercol/BayesPGLS", subdir = "pkg/")
}
library(BayesPGLS)

# Install and load 'phytools':
if (!"phytools" %in% instPacks) {
  install.packages("phytools")
}
library(phytools)

# Install and load 'caper':
if (!"caper" %in% instPacks) {
  install.packages("caper")
}
library(caper)


# =============================== #
# ==== ADDITIONAL FUNCTIONS: ==== 
# =============================== #
# Demographic rates:
surv <- function(th, x) {
  exp(exp(th["a0"]) / th["a1"] * (exp(-th["a1"] * x) - 1) - th["c"] * x +
        exp(th["b0"]) / th["b1"] * (1 - exp(th["b1"] * x)))
}

mort <- function(th, x) {
  exp(th["a0"] - th["a1"] * x) + th["c"] + exp(th["b0"] + th["b1"] * x)
}

# Parallel function for life expectancy:
exParal <- function(isim, theFem, theMal, surv, mort, xv, dx, ncpus) {
  nthe <- nrow(theFem)
  idsim <- seq(0, nthe, nthe / ncpus)
  isim <- 1
  
  id <- (idsim[isim] + 1):idsim[isim + 1]
  # Remaining life expectancy:
  exFull <- t(sapply(id, function(ith) {
    thef <- theFem[ith, ]
    them <- theMal[ith, ]
    Sxf <- surv(th = thef, x = xv) / 
      surv(th = thef, x = xv[1])
    Sxm <- surv(th = them, x = xv) / 
      surv(th = them, x = xv[1])
    exf <- sum(Sxf * dx)
    exm <- sum(Sxm * dx)
    maxex <- max(c(exf, exm))
    exdiff <- (exf - exm) / maxex
    return(c(exFemale = exf, exMale = exm, exDiff = exdiff))
  }))
  return(exFull)
}

# Calculate posterior life expectancy and sex differences:
CalcPostALE <- function(out, ncpus = 4) {
  # Age vector for integration:
  dx <- 0.005
  x <- seq(0, 150, dx)
  
  # Extract Siler mortality parameters (theta) per sex:
  theFem <- out$params[, grep("Female", colnames(out$params))]
  theMal <- out$params[, grep("Male", colnames(out$params))]
  nthe <- nrow(theFem)
  colnames(theFem) <- colnames(theMal) <- c("a0", "a1", "c", "b0", "b1")
  
  # Remaining life expectancy from each converged MCMC iteration:
  sfInit(parallel = TRUE, cpus = ncpus)
  outparal <- sfClusterApplyLB(x = 1:ncpus, exParal, theFem = theFem, 
                               theMal = theMal, mort = mort, surv = surv, 
                               xv = x, dx = dx, ncpus = ncpus)
  sfStop()
  
  # Extract parallel runs into a single data frame:
  for (isim in 1:ncpus) {
    if (isim == 1) {
      exFull <- outparal[[isim]]
    } else {
      exFull <- rbind(exFull, outparal[[isim]])
    }
  }
  return(exFull)
}

# Calculate quantiles of ALE differences:
CalcQuants <- function(object) {
  # Summary statistics names:
  statsnames <- c("Mean", "SD", "Lower", "Upper", "zeroOverlap")
  nstats <- length(statsnames)
  
  # Variables:
  varNames <- c("ALE_Female", "ALE_Male", "ALEdiff") 
  nvar <- length(varNames)
  
  # Calculate summary statistics (posterior mean, SD, and quantiles):
  quants <- matrix(NA, nrow = nvar, ncol = nstats,
                   dimnames = list(varNames, statsnames))
  for (ii in 1:ncol(object)) {
    xi <- object[, ii]
    qi <- c(Mean = mean(xi, na.rm = TRUE), 
            SD = sd(xi, na.rm = TRUE), 
            Lower = quantile(xi, 0.025, na.rm = TRUE),
            Upper = quantile(xi, 0.975, na.rm = TRUE), NA)
    names(qi) <- sprintf("%s_%s", colnames(object)[ii], statsnames)
    quants[ii, ] <- qi
  }
  
  # Calculate zero overlap for ALE differences:
  pv <- pnorm(q = 0, mean = quants["ALEdiff", "Mean"],
              sd = quants["ALEdiff", "SD"])
  if (0 < quants["ALEdiff", "Mean"]) {
    peq <- 2 * pv
  } else {
    peq <- 2 * (1 - pv)
  }
  
  # Store results:
  quants[3, nstats] <- peq
  
  return(quants)
}

# ================================== CODE METADATA =========================== #
# AUTHORS: Fernando Colchero
# DATE CREATED: 2024-08-22
# DESCRIPTION: Reproducibility code for Staerk et al. (2024).
# NOTES: 
# ================================ CODE START ================================ #
# ======================== #
# ==== GENERAL SETUP: ====
# ======================== #
# ----------------------------------------------------- #
# ---- Libraries, functions and working directory: ----
# ----------------------------------------------------- #
# Working directory:
setwd("Your Path/StaerkEtal2024ALEdifferences-main/")

# Libraries and functions:
source("01code/ALEdifferencesFunctions.R")

# -------------------- #
# ---- Load data: ----
# -------------------- #
# Load life tables:
lifeTabs <- read.csv(file = "02data/tables/LifeTables.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Load Siler mortality parameters from maturity:
thetaMat <- read.csv(file = "02data/tables/SilerParamsFromMatur.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Load estimated ALE differences:
leDiff <- read.csv(file = "02data/tables/ALEdifferences.csv",
                   header = TRUE, stringsAsFactors = FALSE)

# Load data for BPGLS:
pglsDat <- read.csv(file = "02data/tables/PGLSdata.csv",
                    header = TRUE, stringsAsFactors = FALSE)

# Test data for BaSTA analysis:
survDat <- read.csv(file = "02data/tables/indivTestDat.csv",
                    header = TRUE, stringsAsFactors = FALSE)

# Load phylogeny:
load("02data/rdata/fullphylo.RData")

# -------------------------- #
# ---- Other variables: ----
# -------------------------- #
# Sexes:
sexes <- c("Female", "Male")
nsex <- length(sexes)

# Plotting colors for sexes:
sexCol <- c(Female = "#BD0026", Male = "#377EB8")

# ================================ #
# ==== VISUALIZE LIFE TABLES: ====
# ================================ #
# Chose species:
species <- "Pan troglodytes"

# Find species on life table data:
idSps <- which(lifeTabs$Species == species)

# Subset life table:
ltSps <- lifeTabs[idSps, ]

# find X limits:
xlim <- c(0, max(ltSps$Ages))

# Produce plot of survival, lx:
par(mfrow = c(1, 1), mar = c(4, 4, 4, 1))
plot(xlim, c(0, 1), col = NA, xlab = "Age", ylab = "Survival")
mtext(text = species, side = 3, line = 2, font = 3)
for (sx in sexes) {
  idsx <- which(ltSps$Sex == sx)
  lines(ltSps$Ages[idsx], ltSps$lx[idsx], col = sexCol[sx], type = 's',
        lwd = 4)
}
legend("topright", legend = sexes, lwd = 4, col = sexCol, bty = 'n')


# ========================================= #
# ==== AGE-SPECIFIC SURVIVAL ANALYSIS: ====
# ========================================= #
# Find life history data for chosen species:
idSps <- which(pglsDat$species == species)

# Find age at maturity:
alpha <- floor(max(c(pglsDat$FemaleAFB[idSps], pglsDat$MaleAFB[idSps])))

# Run BaSTA on survival data:
out <- basta(object = survDat, dataType = "census", model = "GO",
             shape = "bathtub", formulaMort = ~ Sex - 1, minAge = alpha,
             nsim = 4, parallel = TRUE, ncpus = 4)

# Print output to the console:
summary(out)

# Plot goodness of fit plots:
plot(out, plot.type = 'gof')

# ================================================================ #
# ==== CALCULATE POSTERIOR LIFE EXPECTANCY FROM BaSTA OUTPUT: ==== 
# ================================================================ #
# Calculate posterior for ALE per sex and of ALE differences:
exFull <- CalcPostALE(out, ncpus = 4)

# Calculate summary statistics (posterior mean, SD, and quantiles):
quants <- CalcQuants(object = exFull)

# Print results to the console:
print(quants)

# ================================================ #
# ==== EVOLUTIONARY PREDICTORS OF ALE DIFFS.: ====
# ================================================ #
# --------------------------- #
# ---- Test on all data: ----
# --------------------------- #
# Example with cost of reproduction:

# Run BPGLS (note that this may take around one hour):
outPGLSclass <- RunBayesPGLS(formula = exDiff ~ Class + log(MaleBM) + 
                               log(FemaleBM) + Monogamy, data = pglsDat, 
                             weights = "weights", phylo = fullphylo, 
                             ncpus = 6, nsim = 6)

# Print summary:
outPGLSclass

# Plot parameter posterior densities:
plot(outPGLSclass, plot.type = "density")

# ------------------------------- #
# ---- Test on single order: ----
# ------------------------------- #
# Subset data to only include Artiodactyls (for faster computing):
idOrder <- which(pglsDat$Order == "Artiodactyla")
pglsDatOrder <- pglsDat[idOrder, ]

# Run BPGLS:
outPGLSorder <- RunBayesPGLS(formula = exDiff ~ log(MaleBM) + log(FemaleBM) + 
                          Monogamy, data = pglsDatOrder, weights = "weights",
                        phylo = fullphylo, ncpus = 6, nsim = 6)

# Print summary:
outPGLSorder

# Plot parameter posterior densities:
plot(outPGLSorder, plot.type = "density")

# ============================================ #
# ==== EFFECT OF ALLOMETRY ON REGRESSION: ==== 
# ============================================ #
# ---------------------------- #
# Allometric relation (e.g., testes vs male body mass): Y = alpha * X^gamma
# Predictor (e.g., RL): XY = log(Y) - gamma * log(X)
# ---------------------------- #
# Number of observations:
n <- 1000

# Alpha:
al <- 1

# Gamma:
gam <- 0.8

# Simulate X (e.g., male body mass)
x <- runif(n = n, 1, 100)

# Variability in allometry:
del <- 2

# Simulate Y (e.g., testes mass)
y <- exp(log(al) + gam * log(x) + rnorm(n = n, mean = 0, sd = del))

# Simulate predictor (e.g., relative testes mass, RL):
xy <- log(y) - gam * log(x)

# Regression coefficients between response and predictor (e.g., ALE diffs ~ RL)
sig <- 0.5
eta <- 0.5
bet <- 1

# Simulate response Z (e.g., ALE diffs.):
z <- eta + bet * xy + rnorm(n = 100, mean = 0, sd = sig)

# Run regression between Z and log Y, log X:
reg1 <- lm(z ~ log(y) + log(x))

# Run regression between Z and XY (e.g., RL):
reg2 <- lm(z ~ xy)

# Verify that beta in regression 3 is equal to beta1 in regression 2 and equal 
# to beta used for simulation:
bet
reg2$coefficients[2]
reg1$coefficients[2]

# Verify that gamma can be retrieved from the ratio -beta_2 / beta_1 
# in regression 2:
gam
-reg1$coefficients[3] / reg1$coefficients[2]

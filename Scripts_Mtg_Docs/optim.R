########################### R Club August 14, 2018 ############################
########################## optim(ize) Your R Usage! ###########################
# Function descriptions largely from R Documentation.
# Example 1 largely from Ben Bolker's Ecological Models and Data in R (2008).
# This is an excellent book that I highly recommend. If you don't want to shell
# out $57 for the book, PDFs of a draft version are available on his website 
# for free (https://ms.mcmaster.ca/~bolker/emdbook/index.html).
#
############################### optim Function ################################
# General-purpose optimization
# Finds the minimum (by default) of a function, along with the parameter values
# that lead to that minimum. Use "method" to choose optimization algorithm 
# that works best for your situation.
# Arguments:
# par	
#   Initial values for the parameters to be optimized over.
# fn	
#   A function to be minimized (or maximized), with first argument the vector
#   of parameters over which minimization is to take place. It should return a 
#   scalar result.
# gr	
#   A function to return the gradient for derivative-based methods ("BFGS", 
#   "CG", and "L-BFGS-B"). If it is NULL (default), a finite-difference 
#   approximation is used. For the "SANN" method it specifies a function to
#   generate a new candidate point. If it is NULL a default Gaussian Markov
#   kernel is used.
# ...	
#   Further arguments to be passed to fn and gr (set, not optimized over). 
# method	
#   The method to be used. Can be abbreviated.
#   "Nelder-Mead": default. Robust simplex-based approach, slower than
#     derivative-based approaches. Unreliable for single parameter optimization 
#     (use "BFGS" or "Brent" instead).
#   "BFGS": derivative-based (quasi-Newton) approach, generally less robust but
#     faster than Nelder-Mead.
#   "CG": conjugate gradients method, more fragile than BFGS but may be more 
#     successful in really large optimization problems.
#   "L-BFGS-B": Variant on BFGS that allows each parameter to have lower and/or
#     upper bounds.
#   "SANN": A variant of simulated annealing. Should be even more robust (but 
#     slower) than Nelder-Mead.
#   "Brent": derivative-free one-dimensional bounded optimization method. Calls
#     optimize (see below).
# lower, upper	
#   Bounds on the variables for the "L-BFGS-B" method, or bounds in which to 
#   search for method "Brent".
# control	
#   a list of control parameters. If control$fnscale is negative, does 
#   maximization. Lots of other potentially useful stuff too, see Details of 
#   function help.
# hessian	
#   Logical (default FALSE). Should a numerically differentiated Hessian matrix
#   be returned? Can be used to obtain approximate (co)variances and confidence
#   intervals.
#
# Value:
# A list with components:
# par	
#   The best set of parameters found.
# value	
#   The value of fn corresponding to par.
# counts	
#   A two-element integer vector giving the number of calls to fn and gr,
#   respectively. 
# convergence	
#   An integer code. 0 indicates successful completion. See optim help for 
#   error codes.
# message	
#   A character string giving any additional information returned by the 
#   optimizer, or NULL.
# hessian	
#   Only if argument hessian is true. A symmetric matrix giving an estimate of 
#   the Hessian at the solution found. 
#
########################## optimize/optimise Function #########################
# The function optimize searches the interval from lower to upper for a minimum
# or maximum of the function f with respect to its first argument.
# Arguments:
# f	
#   the function to be optimized. The function is either minimized or maximized
#   over its first argument depending on the value of maximum.
# interval	
#   a vector containing the end-points of the interval to be searched.
# ...	
#   additional arguments to be passed to f (not optimized over).
# lower	
#   the lower end of the interval to be searched, set here or in interval.
# upper	
#   the upper end of the interval to be searched, set here or in interval.
# maximum	
#   logical. Should we maximize or minimize (the default)?
# tol	
#   the desired accuracy, affects how close together algorithm searches. 
#   Default .Machine$double.eps^0.25 (0.0001220703 for me). Make smaller for 
#   more precise, slower searches.
#
# Value:
# A list with components minimum (or maximum) and objective, which give the
# location of the minimum (or maximum) and the value of the function at that 
# point, respectively.
#
################################## Examples ###################################
# Example 1: Maximum Likelihood Estimation of Binomial Model
# (tadpole predation experiment)
###############################################################################
# Some definitions:
# Likelihood: probability of seeing the data collected given a particular model
# (either general form or a specific set of parameter values).
# Maximum likelihood estimate (MLE): parameter values that make the observed
# data most likely to have happened.
# Binomial distribution: Discrete probability distribution for number of 
# successes in fixed number of trials, each of which has binary outcome 
# (heads/tails, alive/dead, etc.)
library(emdbook); library(tidyverse)
# Data from field predation trials (Vonesh and Bolker 2005) of African treefrog
# tadpoles. 
data("ReedfrogFuncresp")
# This dataset comes from experiment to determine effects of tadpole density on
# number which are predated.
summary(ReedfrogFuncresp)
head(ReedfrogFuncresp)
tad <- mutate(ReedfrogFuncresp,
              Survived = Initial - Killed,
              PropSurv = Survived / Initial)
summary(tad)
colSums(tad)

# Example 1a: Simple binomial model (single parameter)
# Attempting to estimate per tadpole probability of not being eaten (assuming
# same for all initial densities for now). Uses binomial distribution.
# This problem is simple enough that it can be solved analytically, and it 
# turns out the maximum likelihood estimate of the survival probability is also
# the intuitive answer: total number that survived divided by total number at 
# risk:
with(tad, sum(Survived)/sum(Initial))

# In most cases, have to solve numerically, using optim or something similar.
# We maximize the log-likelihood instead the likelihood for mathematical
# accuracy. Because optim (by default) minimizes a function, we actually
# minimize the negative log-likelihood to get the MLE. We're assuming 
# experiments are independent and identically distributed (iid), so likelihood
# of whole shebang is product of likelihoods of individual experiments. 
# Therefore, log-likelihood of whole thing is sum of log-likelihoods of 
# individual experiments, which can be given using dbinom function:
binomNLL1 <- function(p, k, N) {
  -sum(dbinom(k, prob = p, size = N, log = TRUE))
}

opt1 <- optim(fn = binomNLL1, par = c(p = 0.5), N = tad$Initial, 
              k = tad$Survived, method = "BFGS", hessian = TRUE)
opt1

# Approximate confidence intervals function
confInt <- function(opt) {
  fisherInfo <- solve(opt$hessian)
  sigma <- sqrt(diag(fisherInfo))
  return(opt$par + c(-1.96 * sigma, 1.96 * sigma))
}

confInt(opt1)
m1 <- glm(cbind(Survived, Killed) ~ 1, binomial, tad)
m1
m1.resp <- predict(m1, newdata = data.frame(Initial = 10), type = 'response', 
                   se.fit = TRUE)
m1.resp$fit
m1.resp$fit + c(-1.96, 1.96) * m1.resp$se.fit

# Example 1b: Functional response model (two parameters)
# Supose tadpole predators have a Holling type II function response
# (predation rate = aN/(1 + ahN)) Per capita predation rate of tadpoles
# decreases hyperbolically with density = a/(1 + ahN). Distribution of number
# of tadpoles eaten in a tank likely to be binomial with this probability:
# p = a / (1 + ahN)
# k ~ Binom(p, N)
# where N is number of tadpoles in a tank, a is attack rate, and h is handling
# time.
binomNLL2 <- function(params, k, N) {
  a <- params[1]
  h <- params[2]
  predprob <- a / (1 + a * h * N)
  -sum(dbinom(k, prob = predprob, size = N, log = TRUE))
}
plot(Killed ~ Initial, data = tad, xlim = c(0, 100), ylim = c(0, 100))

opt2 <- optim(fn = binomNLL2, par = c(a = 0.5, h = 0.0125), N = tad$Initial, 
              k = tad$Killed, lower = c(a = 0, h = 0), method = "L-BFGS-B",
              hessian = TRUE)
opt2
opt2a <- optim(fn = binomNLL2, par = c(a = 0.5, h = 0.0125), N = tad$Initial, 
               k = tad$Killed, hessian = TRUE)
opt2a$par
opt2$par
confInt(opt2)
confInt(opt2a)
# Also check out mle2 function in bbmle package, a nice "wrapper" for optim
# specifically for maximum likelihood estimation. To get even more specific,
# check out frair package (functional response analysis). Another option in 
# some simple cases is fitdistr function of MASS package.

###############################################################################
# Example 2: Inverting Matrix Projection (manatee retrospective model)
###############################################################################
# Problem: we have manatee abundance estimate for one year (2011) and want to
# obtain estimates for years before and after that (with uncertainty measures), 
# using estimated survival and reproductive rates each year. For years after 
# 2011, this is relatively easy with a resampling approach:
# 1) Generate random values from distributions for 2011 abundance and survival/
# reproductive rates.
# 2) Construct population projection matrix for each year from survival and 
# reproductive rates.
# 3) Pick a starting stage distribution for 2011, multiply by abundance to get
# number in each stage that year.
# 4) Multiply that vector by population matrix to get population vector next
# year, etc.
# 5) Repeat steps 1 - 4 a large number of times, report distribution of results
#
# In theory, going backwards (hindcast) should be easy too, by inverting the 
# population projection matrices you can get population vector for year before.
# But when you try this in R, small errors creep in, get magnified, and blow 
# things up. So we needed another solution.
library(foreign)
## Functions
# Population projection matrix
fem.manatee.matrix <- function(s, g) {
  return(matrix(c(0,    0,    0,             0,             0.5*s[1], 0,
                  s[2], 0,    0,             0,             0,        0,
                  0,    s[3], 0,             0,             0,        0,
                  0,    0,    s[4]*(1-g[1]), s[5]*(1-g[2]), 0,        0,
                  0,    0,    s[4]*g[1],     s[5]*g[2],     0,        s[6]*g[3],
                  0,    0,    0,             0,             s[6],     s[6]*(1-g[3])),
                nrow = 6, ncol = 6, byrow = TRUE))
}

# Function to optimize:
# Feed in possible starting population size, "known" ending population size,
# starting stage structure, survival rates (6xnumber of years matrix),
# and reproductive rates (3xnumber of years matrix)
# Returns squared difference between resulting ending population size
# and input ending population size.  optimize will find best starting
# population size to minimize this squared difference.
project.to.goal <- function(N.start, N.end, prop.start, s, g) {
  nYears <- ncol(s)
  nStages <- length(prop.start)
  N <- matrix(NA, nStages, nYears + 1)
  N[,1] <- N.start * prop.start
  for (jj in 1:nYears) {
    N[,jj+1] <- fem.manatee.matrix(s[,jj], g[,jj]) %*% N[,jj]
  }
  tot <- sum(N[,nYears + 1])
  return((tot - N.end) ^ 2)
}

# Transforms estimates from sin link scale to response (0 - 1) scale.
sin.link <- function(betas) {
  reals <- (sin(betas)+1)/2
  return(reals)
}

# Generates random values on response scale given estimates and covariance
# matrix on the sin link scale.
rand.sin.mnorm <- function(mu, Sigma, n = 1) {
  require(MASS)
  rand.betas <- mvrnorm(n, mu, Sigma)
  rand.reals <- sin.link(rand.betas)
  return(rand.reals)
}

## Settings
it <- 10000 ##number iterations in simulation  
yrs.future <- 3 # 2010/2011 - 2013/2014
yrs.past <- 14 # 1996/1997 - 2010/2011
yrs <- yrs.future + yrs.past 
#t0 <- 3 # Where to end hindcast 
#yrs.hin <- t0 + yrs.past
years <- 1996 + 1:yrs
Time <- 1997 + 0:yrs

stages <- c(2:4, 'P', 'C', 'B')
nStages <- length(stages)

## Inputs
# Read in betas for subadult/adult survival estimates by year
beta.surv.rand <- read.csv('MARK/beta shrinkage estimates 97-13.csv')
# Rand VCV matrix actually has correlations above diagonal.  Since it rounded off covariances too much, 
# use correlations to get more precise covariances
beta.vcv.surv.rand <- as.matrix(read.csv('MARK/S-tilde CV matrix beta parameters.csv'))
for (i in 2:yrs)
  for (j in 1:(i - 1))
    beta.vcv.surv.rand[i, j] <- beta.vcv.surv.rand[j, i] <- beta.vcv.surv.rand[j, i] * 
  sqrt(beta.vcv.surv.rand[i, i]) * sqrt(beta.vcv.surv.rand[j, j])
# Read in betas for adult reproduction estimates by year
beta.repro <- read.csv('MARK/sw repro beta estimates.csv')
beta.vcv.repro <- as.matrix(read.dbf('MARK/SW REPRO BETAVCMATRIX.DBF'))
which.beta.surv <- 1:yrs
which.beta.surv.vcv <- 1:yrs
which.beta.repro <- 18 + 1:yrs
which.beta.surv.av <- 1 + 1:(yrs-1) 
beta.repro.rand <- read.csv('MARK/repro shrinkage est - beta parameters.csv')
beta.vcv.repro.rand <- as.matrix(read.csv('MARK/repro S-tilde CV matrix beta parameters.csv'))
for (i in 2:yrs)
  for (j in 1:(i - 1))
    beta.vcv.repro.rand[i, j] <- beta.vcv.repro.rand[j, i] <- beta.vcv.repro.rand[j, i] * 
  sqrt(beta.vcv.repro.rand[i, i]) * sqrt(beta.vcv.repro.rand[j, j])

beta.surv.usj <- read.csv('MARK/usj survival beta estimates age.csv')
beta.vcv.surv.usj <- as.matrix(read.dbf('MARK/BETAVCMATRIX.DBF'))
which.beta.usj <- 1:3

umk <- 7.675##mean initial population lognormal
sdmk <- 0.144##uncertainty initial population

l.g4.mean <- -6.907 ##mean 4th year reproduction rate, on logit scale, from appendix B
l.g4.sd <- 3.055 ##uncertainty in reproduction rate, on logit scale, from appendix B

l.gP.mean <- -0.8283 ##mean prebreeder reproduction rate, on logit scale, from appendix B and discussion with Mike
l.gP.sd <- 0.4532 ##uncertainty in reproduction rate, on logit scale, from appendix B and discussion with Mike

mu.usj <- beta.surv.usj$Beta[which.beta.usj]
Sigma.usj <- beta.vcv.surv.usj[which.beta.usj,which.beta.usj]

# stable stage structure 
prop.start0 <- c(0.06923249, 0.06049673, 0.05650075, 0.14519273, 0.18396794, 0.48460936)

#Get mean and uncertainty for adult survival rates
mu.adult <- beta.surv.rand$S.tilde[which.beta.surv]
Sigma.adult <- beta.vcv.surv.rand[which.beta.surv.vcv,which.beta.surv.vcv]

#Get mean and uncertainty for adult reproductive rates
mu.repro.adult <- beta.repro.rand$S.tilde
Sigma.repro.adult <- beta.vcv.repro.rand

## Results storage
N <- array(0, c(it, nStages, yrs+1)) # Storing abundance by stage in one big array
Ntot <- array(0, c(it, yrs+1)) # And total abundance in somewhat smaller array
dimnames(N) <- list(1:it, stages, Time)
dimnames(Ntot) <- list(1:it, Time)

ss <- array(0, c(it, nStages, yrs))##"true" survival rates for all stages in one big array
gg <- array(0, c(it, 3, yrs))##"true" reproduction rate

## Run resampling
set.seed(456)
for (ii in 1:it) ##initiates iteration loop ## statistical uncertainty
{
  if (ii %% (it / 10) == 0)
    cat('\tSimulation', ii, 'of', it, '\n')
  # Resample from statistical distribution for adult survival for each year
  ss.adult <- rand.sin.mnorm(mu.adult, Sigma.adult) 
  # Resample from statistical distribution for 4th year and prebreeder reproductive rates
  gg4 <- rep(plogis(rnorm(1, l.g4.mean, l.g4.sd)), yrs)
  ggP <- rep(plogis(rnorm(1, l.gP.mean, l.gP.sd)), yrs) 
  # Force 4th year reproductive rate to be lower than prebreeder
  while (gg4[1] >= mean(ggP)) {
    gg4 <- rep(plogis(rnorm(1, l.g4.mean, l.g4.sd)), yrs)
  }
  # Resample from statistical distribution for adult reproductive rates each year
  gg.adult <- rand.sin.mnorm(mu.repro.adult, Sigma.repro.adult)
  # Resample from distribution of calf and subadult mortality ratios (compared to adult)
  survUSJ <- rand.sin.mnorm(mu.usj, Sigma.usj)   
  mortRatioCalf <- (1 - survUSJ[1:2]) / (1 - survUSJ[3])
  ss.calf <- 1 - c(mortRatioCalf[1] * (1 - ss.adult), mortRatioCalf[2] * (1 - ss.adult))
  subAdSurvUSJ <- sin.link(rnorm(2, mu.usj[3], sqrt(Sigma.usj[3, 3])))
  mortRatioSub <- (1 - subAdSurvUSJ) / (1 - survUSJ[3])
  ss.sub <- 1 - c(mortRatioSub[1] * (1-ss.adult), mortRatioSub[2] * (1-ss.adult))
  ss.calf[ss.calf < 0] <- 0
  # Store survival and reproductive rates for this iteration
  ss[ii,,] <- matrix(c(ss.calf, ss.sub, rep(ss.adult, 2)), 6, yrs, byrow = T)
  gg[ii,,] <- matrix(c(gg4, ggP, gg.adult), 3, yrs, byrow = T)
  ss.past <- ss[ii,,1:yrs.past]
  gg.past <- gg[ii,,1:yrs.past]
  #sample lognormal distribution for 2011 population estimate
  xxx <- rlnorm(1,umk,sdmk)
  opt.N.start <- optimize(project.to.goal, c(1, 60000), N.end = xxx, 
                          prop.start=prop.start0, s=ss.past, g=gg.past)
  N.start <- opt.N.start$minimum
  N[ii, , 1] <- N.start * prop.start0 ##initial population
  for(jj in 1:yrs)  #initiate annual loop ## environmental variance
  {
    N[ii,,jj+1] <- fem.manatee.matrix(ss[ii,,jj], gg[ii,,jj]) %*% N[ii,,jj]
  } # Close the annual loop
  Ntot[ii,] <- apply(N[ii,,], 2, sum)
} # Close the iteration loop


## Summarize and plot some of the outputs
# 2011 initial abundance distribution
x2011 <- rlnorm(it*10,umk,sdmk)

# Total population size, summarized across iterations
Ntot.sum <- t(apply(Ntot, 2, quantile, probs = c(0.025, 0.5, 0.975)))
Ntot.sum2 <- data.frame(Year = as.integer(row.names(Ntot.sum)), as.data.frame(Ntot.sum))
summary(Ntot.sum2)
colnames(Ntot.sum2)[2:4] <- c('Low', 'Median', 'High')
Ntot.sum2 <- rbind(Ntot.sum2, 
                   data.frame(Year = 2011, Low = quantile(x2011, probs = 0.025),
                              Median = quantile(x2011, probs = 0.5),
                              High = quantile(x2011, probs = 0.975)))
Ntot.sum2$Type <- factor(c(rep('Hindcast', yrs.past), 'Starting', 
                           rep('Projection', yrs.future), 'Original'), 
                         levels = c('Original', 'Hindcast', 'Starting', 'Projection'))

# Plot total abundance
plot.Ntot<-ggplot(Ntot.sum2, 
                  aes(x = Year, y = Median))+
  geom_pointrange(aes(ymin = Low, ymax = High, color = Type)) + 
  labs(x="Year",y="Estimated Abundance")+theme_bw() +
  scale_color_brewer(type = 'qual', palette = 2) +
  ylim(0, NA)
plot.Ntot

# Add dodge so can see both estimates from 2011
plot.Ntot<-ggplot(Ntot.sum2, 
                  aes(x = Year, y = Median))+
  geom_pointrange(aes(ymin = Low, ymax = High, color = Type), position = position_dodge(width=0.5)) + 
  labs(x="Year",y="Estimated Abundance")+theme_bw() +
  scale_color_brewer(type = 'qual', palette = 2) +
  ylim(0, NA)
plot.Ntot

require("deming")
library(deming)
library(here)

#TESTESTTEST

## Resultaten van regressie met package 'deming' geven confidence intervals die niet convergeren
## naar die van standaard lm, als de sd van 'x' klein wordt.  Ik gebruik een empirische methode
## om CIs te schatten.

sample_data <- function(N=25, lambda=200, alpha=-0.25, beta=0.0, cer.shape=c(0.75,2), i2=0, cer0 = 0.5) {
  ##
  ## N = number of studies
  ## lambda = Poisson rate for participants per study
  ## alpha = scaled relative risk (RR intercept - 1 at CER=infty)
  ## beta = scaled relative risk slope ( 0 < beta < 1)
  ## cer.shape = shape parameters for beta distribution of CER
  ## i2 = I squared parameter - fraction of variance due to random effects
  ##
  ## use beta = 0.25 for strong CER dependence, beta = 0 for none
  ##
  ## cer.shape was (1,2)
  nc = 1 + rpois(N, lambda)
  nt = 1 + rpois(N, lambda)
  cer = rbeta(N, cer.shape[1], cer.shape[2])
  
  ## add random effect
  cer.sd = sqrt( cer * (1-cer) / nc )
  ter0 = cer * (1 + alpha - alpha * beta / cer)    ## approximation of ter - before adding random effects
  ter0.sd = sqrt( ter0 * (1-ter0) / nt )
  rr.sd = (ter0/cer) * sqrt( (ter0.sd/ter0)**2 + (cer.sd/cer)**2 )
  random.sd = sqrt(1/(1 - i2) - 1) * rr.sd
  epsilon = rnorm( N, 0, sd = random.sd )    ## this amount of variation corresponds to the I^2 value provided
  rr = (1 + alpha - alpha * beta / cer) + epsilon
  
  ## calculate ter
  ter = rr * cer
  
  ## avoid occasional ter > 1, ter < 0
  cer[ ter > 1 ] = cer[ ter > 1 ] / ter[ ter > 1 ]
  ter[ ter > 1 ] = 1.0
  ter[ ter < 0 ] = 0.0
  
  ## simulate data
  ec = nc
  et = nt
  for(n in 1:N) {
    ec[n] = rbinom(1, nc[n], cer[n])
    et[n] = rbinom(1, nt[n], ter[n])
  }
  data <- list(
    N = N,
    nc = nc,
    ec = ec,
    nt = nt,
    et = et,
    cer0 = cer0,
    true_cer.sd = sqrt(cer * (1-cer) / nc),
    true_ter.sd = sqrt(ter * (1-ter) / nt),
    true_rr0 = 1 + alpha - alpha * beta / cer0,
    true_rr1 = 1 + alpha - alpha * beta)
  return(data)
}

run.deming <- function(data, cer0 = 0.5, replicates = 100, use.new.var=TRUE, use.true.sd=FALSE) {
  ##
  ## cer0 : control event rate for which the estimated RR and confidence interval is returned
  ## replicates : number of simulated data sets to use for estimation of vcov
  ##
  
  cer <- (0.5 + data$ec) / (1 + data$nc)   ## Expected posterior CER under Jeffrey's prior
  ter <- (0.5 + data$et) / (1 + data$nt)   ## Expected posterior TER under Jeffrey's prior
  ## Expectation of p(1-p) under Jeffrey's prior is P(1-P)(N+1/(N+2)), where N = nc or nt
  #cer.sd <- sqrt( (0.5 + data$ec) * (0.5 + data$nc - data$ec) / (2 + data$nc) ) / (1 + data$nc)
  #ter.sd <- sqrt( (0.5 + data$et) * (0.5 + data$nt - data$et) / (2 + data$nt) ) / (1 + data$nt)
  ## Other attempt - sample from the Dirichlet, and use empirical estimate
  ## This one gives better calibrated CI's
  cer.sd <- cer
  ter.sd <- ter
  for (i in 1:length(cer)) {
    cer.post <- rbeta(1000, 0.5 + data$ec[i], 0.5 + data$nc[i] - data$ec[i])
    ter.post <- rbeta(1000, 0.5 + data$et[i], 0.5 + data$nt[i] - data$et[i])
    cer.sd[i] <- mean(sqrt(cer.post * (1-cer.post) / (1+data$nc[i])))
    ter.sd[i] <- mean(sqrt(ter.post * (1-ter.post) / (1+data$nt[i])))
    #cer.sd[i] <- 1/mean(1/sqrt(cer.post * (1-cer.post) / (1+data$nc[i])))
    #ter.sd[i] <- 1/mean(1/sqrt(ter.post * (1-ter.post) / (1+data$nt[i])))
  }
  if (use.true.sd) {
    ## use actual SD (based on true parameter) rather than SD estimated from counts
    ## (for use in simulation)
    cer.sd <- data$true_cer.sd
    ter.sd <- data$true_ter.sd
  }
  
  ## regress TER against CER, and take account of errors in both
  ## (note - these are the expected errors under a binomial event count model, and do not 
  ##  include any additional random variation in RR)
  
  result <- deming( ter ~ cer, xstd = cer.sd, ystd = ter.sd, x=TRUE, y=TRUE, model=TRUE)
  
  ## calculate semi-empirical variance-covariance matrix - the values provided by the package
  ## are not correct
  
  result$variance.old <- result$variance
  result$ci.old <- result$ci
  if (use.new.var) {
    xx = cbind(1, cer)
    a = result$coef[2]  ## slope estimate
    lambda = (ter - xx %*% result$coef) / (a * a * cer.sd + ter.sd)
    xx0 <- cer + lambda * a * cer.sd
    yy0 <- ter - lambda * ter.sd          ## points (xx0, yy0) now fall on the estimated line
    nn.gen = length(cer) * replicates     ## generate data, as a matrix with data in rows: cer.gen[row,]
    cer.gen <- matrix(rnorm( n = nn.gen, mean = xx0, sd = cer.sd ), nrow = replicates, byrow=TRUE)
    ter.gen <- matrix(rnorm( n = nn.gen, mean = yy0, sd = ter.sd ), nrow = replicates, byrow=TRUE)
    ## fit curves to generated data, and collect intercept and slope estimates
    results.gen <- matrix( 0, nrow = replicates, ncol = 2)  ## intercept and slope
    for (i in 1:replicates) {
      result.gen <- deming( ter.gen[i,] ~ cer.gen[i,], xstd = cer.sd, ystd=ter.sd, x=TRUE, y=TRUE, model=TRUE)
      results.gen[i,] <- result.gen$coef
    }
    ## calculate variance-covariance matrix and 95% CI's, and assign to results
    result$variance <- cov(results.gen)
    result$ci[1,1] <- c(result$coef[1] - 1.96*sqrt(result$variance[1,1]))
    result$ci[1,2] <- c(result$coef[1] + 1.96*sqrt(result$variance[1,1]))
    result$ci[2,1] <- c(result$coef[2] - 1.96*sqrt(result$variance[2,2]))
    result$ci[2,2] <- c(result$coef[2] + 1.96*sqrt(result$variance[2,2]))
  }
  
  ## calculate RR at CER0.  First coefficient is intercept, second is coefficient of CER
  x = c(1, cer0)
  ter.pred <- result$coef %*% x
  ter.var <- x %*% result$variance %*% x
  rr <- ter.pred / cer0
  rr.sd <- sqrt(ter.var) / cer0
  result$RR.at.CER0 <- list(rr = rr, sd = rr.sd, lower = rr - 1.96*rr.sd, upper = rr + 1.96*rr.sd)    
  return(result)
}

assessment <- function(params, n, use.new.var=TRUE, use.true.sd=FALSE) {
  ## how often do we get a significant effect of 1/CER on the RR?
  signif_nonzero_count = 0
  param = 1  ## 1 (intercept) for new (ter ~ cer) model, 2 (slope) for old (rr ~ 1/cer) model
  for (i in 1:n) {
    data <- do.call(sample_data, params)
    dem <- run.deming(data, use.new.var=use.new.var, use.true.sd = use.true.sd)
    signif <- dem$ci[param,2]< 0 | dem$ci[param,1] > 0  ## does CI for intercept parameter exclude 0?
    signif_nonzero_count <- signif_nonzero_count + signif
  }
  return(signif_nonzero_count / n)
}

##############################################################
## SIMULATION - wordt niet gerund, neemt veel tijd
##############################################################

##
## make data.frame with all combinations of the listed values
##
tests <- expand.grid(  i2 = c(0,0.25,0.5), N=c(10,50),  beta=c(0,0.1,0.25) )

## apply all tests in a data.frame
apply.tests <- function(tests, n) {
  powers <- c()
  for (i in 1:dim(tests)[1]) {
    power <- assessment( tests[i,], n=n, use.new.var=TRUE, use.true.sd=FALSE )
    powers <- c(powers, power)
  }
  tests$power <- powers
  return(tests)
}

## 1: 0.03 (200); 0.069 (1000); 0.058 (1000) -- using true sd
## 1: 0.025 (1000) -- using sd estimated from expected variance
## 1: 0.032 (1000) -- using sd estimated as E(sd | posterior)     <-- deze gekozen
## 1: 0.056 0.044 (1000) -- using sd estimated as 1/E(1/sd | posterior)
## 1: 0.211 (2000) - using sd estimated as 1/sqrt(E(1/var | posterior))


## n=200 takes 60 seconds per test
if (FALSE) {
  tests.results <- apply.tests(tests[1,], n=2000)
  tests.results
}

## Visualisatie hiervan?  Voeg i2 = 0.25 toe?

analyse.thijsdata <- function(fname, cer0=0.5) {
  data <- read.csv(paste("csv files/",fname,sep=""))
  colnames(data) <- c("RCT","et","nt","ec","nc")
  out <- run.deming(data, cer0=cer0)
  return(out)
}

summary.thijsdata <- function(fname) {
  data <- read.csv(paste("csv files/",fname,sep=""))
  colnames(data) <- c("RCT","et","nt","ec","nc")
  data$cer <- (data$ec + 1) / (data$nc + 1)
  return(data$cer)
}

studies <- list.files(path = "csv files")

param = 1
cer0 = 0.5       ## verander dit indien nodig
results <- NULL
for (study in studies) {
  out <- analyse.thijsdata( study, cer0=cer0 )
  result <- data.frame(study,
                       out$coeff[1],  ## intercept
                       out$ci[1,1],   ## lower CI
                       out$ci[1,2],   ## upper CI
                       out$RR.at.CER0$rr,    ## RR estimate, at CER = 0.5 (of zoals boven gekozen)
                       out$RR.at.CER0$lower, ## lower CI
                       out$RR.at.CER0$upper, ## upper CI
                       out$coeff[2],  ## slope
                       out$ci[2,1],   ## lower
                       out$ci[2,2])   ## upper
  row.names(result) <- NULL
  names(result) <- c("study","b", "b.lower95", "b.upper95",
                     "RR", "RR.lower95", "RR.upper95",
                     "a", "a.lower95", "a.upper95")
  results <- rbind(results, result)
  print(study)
}

results

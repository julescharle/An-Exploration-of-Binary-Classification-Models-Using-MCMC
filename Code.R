library(ggplot2)
library(gridExtra)
library(mvtnorm)
library(coda)
set.seed(9323) #student id
############################################################
############################################################
######################## Task 1 ############################
############################################################
############################################################

 

n <- 150 # number of observations
d <- 10 # number of beta parameters

# create matrix to populate with covariates
X <- matrix(nrow = n, ncol = d)
X[,1] <- rep(1, n) # first column is an intercept

# create base uncorrelated random numbers to turn into x_i's
z <- matrix(rnorm(n*(d-1)), nrow = n, ncol = d-1)

# create x_i's (ith column of matrix corresponds to variable x_i)
X[,2] <- z[,1]
X[,3] <- z[,1] + 0.2*z[,2]
X[,4] <- 0.5*z[,3]
X[,5] <- z[,4]
X[,6] <- 2*z[,5] + 20*z[,6]
X[,7] <- z[,6]
X[,8] <- 0.5 * (z[,7] + z[,4] + z[,8] + z[,1])
X[,9] <- z[,8] + 10*z[,4]
X[,10] <- z[,5] + 0.5*z[,9]

# create true beta values
beta <- seq(-2,2, length = 10)

############################################################
############################################################
######################## Task 2 ############################
############################################################
############################################################

# using the latent variable interpretation as instructed
epsilon_C <- rt(n,1)
Y_star <- apply(X, 1, function(row) sum(row * beta)) + epsilon_C
Y_C <- as.integer(Y_star > 0)



############################################################
############################################################
######################## Task 3 ############################
############################################################
############################################################



rejection_sampler <- function(n,M) {
  # logistic distribution
  Pi <- function(x) {exp(-x)/(1 + exp(-x))^2}
  # cauchy distribution 
  q <- function(x) {1 / (pi * (1 + x^2))}
  
  epsilon <- numeric(n)  
  iter <- 0           
  accept <- 0
  while (accept < n) {
    
    candidate <- rt(1,1)
    
    u <- runif(1)
    
    ratio <- Pi(candidate) / (M * q(candidate))
    
    # due to the candidate distribution being cauchy, we can sometimes get extremely small values.
    # these small values will cause the ratio to be incalculable by R
    # to fix this we have to add an !is.nan argument. 
    # this still works as if the ratio was incalcuable it would be so small that it would be rejected anyway
    # this will be carried forwards into the next task as well.
    
    if (!is.nan(ratio) && u < ratio) {
      accept <- accept + 1
      epsilon[accept] <- candidate
    }
    iter <- iter + 1
  }
  
  list(epsilon,iter)
}


# what value of M? we want the smallest M that is larger than the ratio of pi(x)/q(x):
# look at the maxima of the ratio graphically


x_values <- seq(-9, 9, by = 0.1)
cauchy <- 1 / (pi * (1 + x_values^2))
logistic <- exp(-x_values) / (1 + exp(-x_values))^2
ratio <- exp(-x_values) * pi * (1 + x_values^(2)) / (1 + exp(-x_values))^2
ratio_as_function <- function(x){
  exp(-x) * pi * (1 + x^(2)) / (1 + exp(-x))^2
}
# finding the maxima of M
optimal_m <- optimize(ratio_as_function, interval = c(-10, 10), maximum = TRUE)$objective

# show this graphically
rejection_data <- data.frame(x = x_values, pi = logistic, q = cauchy, ratio = ratio)
ggplot(rejection_data, aes(x = x)) +
  geom_line(aes(y = pi, colour = "pi"), size = 1.25) +
  geom_line(aes(y = q, colour = "q"), size = 1.25) +
  geom_line(aes(y = ratio, colour = "ratio"), size = 1.25) +
  scale_colour_manual("", 
                      breaks = c("pi", "q", "ratio"),
                      values = c("gold", "firebrick", "royalblue")) +
  labs(y = "density or ratio") +
  geom_hline(yintercept = optimal_m, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = 3, y = 1.7, label = "maxima = 1.652", hjust = 1, size = 5) +
  theme_classic() +
  ylim(0, 1.7)  + 
  xlim(-9,9)
  
# now, to generate n samples of Y^L
samples <- n
epsilon_L_output <- rejection_sampler(samples,M = optimal_m)
epsilon_L <- epsilon_L_output[[1]]

# theoretical vs actual samples
cat("Theoretical number of samples: ", optimal_m*n)
cat("Actual number of samples: ", epsilon_L_output[[2]])


Y_star_L <- apply(X, 1, function(row) sum(row * beta)) 
Y_star_L <- Y_star_L + epsilon_L

# Y^L values
Y_L <- as.integer(Y_star_L > 0)

############################################################
############################################################
######################## Task 4 ############################
############################################################
############################################################
# pre-conditioned RWM

RWM <- function(logpi, nits, h, beta_curr, X, Y, covariancemat)  {
  # calculate the current density
  logpi_curr <- logpi(beta_curr, X, Y)
  
  # initialise
  accepted <- 0
  beta_store <- matrix(nrow = nits, ncol = 10)
  
  for (i in 1:nits) {
    #proposal move
    beta_prop <- beta_curr + h * as.vector(rmvnorm(1, sigma = covariancemat))
    logpi_prop <- logpi(beta_prop, X, Y)
    
    #log acceptance rate
    loga <- logpi_prop - logpi_curr
    u <- runif(1)
    if (!is.na(loga)){
      if (log(u) < loga) { #if accepted, change the current state
        beta_curr <- beta_prop
        logpi_curr <- logpi_prop
        accepted <- accepted + 1
      }
      #if rejected just keep chain where it is
    }
    #store this state
    beta_store[i,] <- beta_curr
    
  }
  # return both the chain history and the acceptance rate, we need the acceptance rate output to be tuned to ~ 0.23 via the step size
  return(list(beta_store = beta_store, a_rate = accepted/nits))
}





### standard prior - cauchit and logistic


logpi_C <- function(beta,X,Y){
  w <- X %*% beta # this is F(x_i^T \beta)
  log_y_givenbeta <- sum(dbinom(Y, size = 1, prob = pcauchy(w), log = TRUE)) # log likelihood function
  return(log_y_givenbeta - 0.5 * sum((beta)^2)) # log likelihood + log prior = log posterior
}
# same logic for logistic
logpi_L <- function(beta,X,Y){
  w <- X %*% beta
  log_y_givenbeta <- sum(dbinom(Y, size = 1, prob = plogis(w), log = TRUE))
  return(log_y_givenbeta - 0.5 * sum((beta)^2))
}

# generic optimisation function. start at 0 if not logit sometimes has some problems converging.
optimize_log_posterior <- function(log_posterior, X, Y) {
  opt_result <- optim(c(rep(0,10)), log_posterior, X = X, Y = Y, control = list(fnscale = -1), hessian = TRUE)
  return(opt_result)
}


modes_L <- optimize_log_posterior(logpi_L, X, Y_L)$par # check if these look reasonable
modes_L_hessian <- - (optimize_log_posterior(logpi_L, X, Y_L)$hessian) # take the negative hessian
eigen(modes_L_hessian)$values #positive definite! nice 

modes_C <- optimize_log_posterior(logpi_C, X, Y_C)$par # check if these look reasonable
modes_C_hessian <- - (optimize_log_posterior(logpi_C, X, Y_C)$hessian) # take the negative hessian
eigen(modes_C_hessian)$values #positive definite ! nice




covariance_matrix_C <- solve(modes_C_hessian) # take the inverse of negative hessian
covariance_matrix_L <- solve(modes_L_hessian) # take the inverse of negative hessian

output_id_L<- RWM(logpi = logpi_L, nits = 40000, h = 0.84, beta_curr = c(rmvnorm(1, rep(0,10), diag(10))),X = X, Y = Y_L, covariancemat = covariance_matrix_L)
output_id_L$a_rate # approx 0.23

output_id_C<- RWM(logpi = logpi_C, nits = 40000, h = 0.91, beta_curr = c(rmvnorm(1, rep(0,10), diag(10))),X = X, Y = Y_C, covariancemat = covariance_matrix_C) 
output_id_C$a_rate # approx 0.23

### UIP

UIP_sigma <- solve(t(X) %*% X) * n # this is our UIP covariance as stated

# same as previous for the log posteriors. no summation needed on the prior this time.
logpi_L_UIP <- function(beta,X,Y){
  w <- X %*% beta
  log_y_givenbeta <- sum(dbinom(Y, size = 1, prob = plogis(w), log = TRUE))
  return(log_y_givenbeta + dmvnorm(beta,mean = rep(0,10), sigma = UIP_sigma, log = TRUE))
}

logpi_C_UIP <- function(beta,X,Y){
  w <- X %*% beta
  log_y_givenbeta <- sum(dbinom(Y, size = 1, prob = pcauchy(w), log = TRUE))
  return(log_y_givenbeta + dmvnorm(beta,mean = rep(0,10), sigma = UIP_sigma, log = TRUE))
}

modes_L_UIP <- optimize_log_posterior(logpi_L_UIP, X, Y_L)$par # check if these look reasonable
modes_L_hessian_UIP <- - (optimize_log_posterior(logpi_L_UIP, X, Y_L)$hessian) # take the negative hessian
eigen(modes_L_hessian_UIP)$values #positive definite ! nice


modes_C_UIP <- optimize_log_posterior(logpi_C_UIP, X, Y_C)$par # check if these look reasonable
modes_C_hessian_UIP <- - (optimize_log_posterior(logpi_C_UIP, X, Y_C)$hessian) # take the negative hessian
eigen(modes_C_hessian_UIP)$values #positive definite ! nice


covariance_matrix_C_UIP <- solve(modes_C_hessian_UIP) # take the inverse of the negative hessian
covariance_matrix_L_UIP <- solve(modes_L_hessian_UIP) # take the inverse of the negative hessian

output_id_L_UIP<- RWM(logpi = logpi_L_UIP, nits = 40000, h = 0.85, beta_curr = c(rmvnorm(1, rep(0,10), covariance_matrix_L_UIP)),X = X, Y = Y_L, covariancemat = covariance_matrix_L_UIP)
output_id_L_UIP$a_rate # approx 0.23
output_id_C_UIP<- RWM(logpi = logpi_C_UIP, nits = 40000, h = 0.9, beta_curr = c(rmvnorm(1, rep(0,10), covariance_matrix_C_UIP)),X = X, Y = Y_C, covariancemat = covariance_matrix_C_UIP)
output_id_C_UIP$a_rate # approx 0.23


# at the end of the code, i also go through adaptive mcmc as it was referenced in the report.

################################################################################################
################################################################################################
################################################################################################

### diagnostic testing

# investigating mixing - just looking at beta_1 for a visual representation, multiple betas is too messy to infer anything from.
# in this order: logistic standard, cauchit standard, logistic UIP, cauchit UIP

beta_1_L <- ggplot(data = data.frame(nits = 1:length(output_id_L$beta_store[,1]), beta_1 = output_id_L$beta_store[,1]), aes(x = nits, y = beta_1)) +
  geom_line() +
  xlab("nits") +
  ylab("beta - L") +
  theme_classic()
beta_1_C <- ggplot(data = data.frame(nits = 1:length(output_id_C$beta_store[,1]), beta_1 = output_id_C$beta_store[,1]), aes(x = nits, y = beta_1)) +
  geom_line() +
  xlab("nits") +
  ylab("beta - C") +
  theme_classic()
beta_1_L_UIP <- ggplot(data = data.frame(nits = 1:length(output_id_L_UIP$beta_store[,1]), beta_1 = output_id_L_UIP$beta_store[,1]), aes(x = nits, y = beta_1)) +
  geom_line() +
  xlab("nits") +
  ylab("beta - L UIP") +
  theme_classic()
beta_1_C_UIP <- ggplot(data = data.frame(nits = 1:length(output_id_C_UIP$beta_store[,1]), beta_1 = output_id_C_UIP$beta_store[,1]), aes(x = nits, y = beta_1)) +
  geom_line() +
  xlab("nits") +
  ylab("beta - C UIP") +
  theme_classic()
grid.arrange(beta_1_L, beta_1_C, beta_1_L_UIP, beta_1_C_UIP, ncol = 4, nrow = 1)


# function which takes 4 chains and compares its density estimates of beta to the true beta
plot_beta_densities <- function(beta_store_logit_std, beta_store_cauchit_std, beta_store_logit_UIP, beta_store_cauchit_UIP, i) {
  
  beta_logit_std <- mcmc(beta_store_logit_std[,i])
  beta_cauchit_std <- mcmc(beta_store_cauchit_std[,i])
  beta_logit_UIP <- mcmc(beta_store_logit_UIP[,i])
  beta_cauchit_UIP <- mcmc(beta_store_cauchit_UIP[,i])
  
  dens1 <- density(beta_logit_std)
  dens2 <- density(beta_cauchit_std)
  dens3 <- density(beta_logit_UIP)
  dens4 <- density(beta_cauchit_UIP)
  
  df1 <- data.frame(x = dens1$x, y = dens1$y, group = "Logit Standard")
  df2 <- data.frame(x = dens2$x, y = dens2$y, group = "Cauchit Standard")
  df3 <- data.frame(x = dens3$x, y = dens3$y, group = "Logit UIP")
  df4 <- data.frame(x = dens4$x, y = dens4$y, group = "Cauchit UIP")
  
  combined_df <- rbind(df1, df2, df3, df4)
  
  ggplot(combined_df, aes(x = x, y = y, colour = group)) +
    geom_line(size = 1.4) +
    labs(x = paste("beta", i), y = "Density") +
    scale_colour_manual(values = c("Logit Standard" = "lightblue", "Cauchit Standard" = "peachpuff", "Logit UIP" = "dodgerblue", "Cauchit UIP" = "darkorange")) +
    geom_vline(xintercept = beta[i], linetype = "dashed", color = "black", size = 0.5) +
    xlim(-3.5, 3.5) +
    theme_classic() +
    theme(legend.position = "none")
  
}

# make a plot for each \beta_i
plots <- list()
for (i in 1:10) {
  plots[[i]] <- plot_beta_densities(output_id_L$beta_store, output_id_C$beta_store, output_id_L_UIP$beta_store, output_id_C_UIP$beta_store, i)
}
# arrange in a grid
grid.arrange(grobs = plots, ncol = 5, nrow = 2)


# legend in case of any confusion
plot(1, type="n", axes=FALSE, xlab="", ylab="", main="")
legend("topright", legend=c("Logit, standard prior", "Cauchit, standard prior", "Logit, UIP", "Cauchit, UIP"), fill=c("lightblue", "dodgerblue", "peachpuff", "darkorange"), bty="n", cex=1.2)

# brier scores for each model/prior combination
brier_score <- function(beta_est, G){
  w <- X %*% colMeans(beta_est)
  p_y_est = G(w) #F(x_i^T \beta)
  return(1/n * (sum((Y_L-p_y_est)^2))) #brier score
}

brier_score(output_id_L$beta_store, plogis) #logistic - standard
brier_score(output_id_C$beta_store, pcauchy) #cauchit - standard
brier_score(output_id_L_UIP$beta_store, plogis) #logistic - UIP
brier_score(output_id_C_UIP$beta_store, pcauchy) #cauchit - UIP

# Gelman-Rubin

# we need two chains of each model/prior combination for GR
output_id_L_2<- RWM(logpi = logpi_L, nits = 40000, h = 0.84, beta_curr = c(rmvnorm(1, rep(0,10), diag(10))),X = X, Y = Y_L, covariancemat = covariance_matrix_L)
output_id_C_2<- RWM(logpi = logpi_C, nits = 40000, h = 0.91, beta_curr = c(rmvnorm(1, rep(0,10), diag(10))),X = X, Y = Y_C, covariancemat = covariance_matrix_C)
output_id_L_UIP_2<- RWM(logpi = logpi_L_UIP, nits = 40000, h = 0.85, beta_curr = c(rmvnorm(1, rep(0,10), covariance_matrix_L_UIP)),X = X, Y = Y_L, covariancemat = covariance_matrix_L_UIP)
output_id_C_UIP_2<- RWM(logpi = logpi_C_UIP, nits = 40000, h = 0.9, beta_curr = c(rmvnorm(1, rep(0,10), covariance_matrix_C_UIP)),X = X, Y = Y_C, covariancemat = covariance_matrix_C_UIP)


# apologies for the slightly ugly code, just need to convert each beta store to mcmc and compute GR statistic
mcmc1_L <- mcmc(output_id_L$beta_store)
mcmc2_L <- mcmc(output_id_L_2$beta_store)
mcmc1_C <- mcmc(output_id_C$beta_store)
mcmc2_C <- mcmc(output_id_C_2$beta_store)
mcmc1_L_UIP <- mcmc(output_id_L_UIP$beta_store)
mcmc2_L_UIP <- mcmc(output_id_L_UIP_2$beta_store)
mcmc1_C_UIP <- mcmc(output_id_C_UIP$beta_store)
mcmc2_C_UIP <- mcmc(output_id_C_UIP_2$beta_store)
combinedchains_L = mcmc.list(mcmc1_L, mcmc2_L)
combinedchains_C = mcmc.list(mcmc1_C, mcmc2_C)
combinedchains_L_UIP = mcmc.list(mcmc1_L_UIP, mcmc2_L_UIP)
combinedchains_C_UIP = mcmc.list(mcmc1_C_UIP, mcmc2_C_UIP)
gelman.diag(combinedchains_L)
gelman.diag(combinedchains_C)
gelman.diag(combinedchains_L_UIP)
gelman.diag(combinedchains_C_UIP)
# they all seem approx= 1, very good sign.

################################################################################################  
################################################################################################
################################################################################################

# extra code, testing with adaptive MCMC

# begin with the negative inverse hessian
sig_L <- covariance_matrix_L
sig_C <- covariance_matrix_C

for (i in 1:25) {
  run_L <- RWM(logpi = logpi_L, nits = 10000, h = 0.8, beta_curr = colMeans(output_id_L$beta_store), X = X, Y = Y_L, covariancemat = sig_L)
  run_C <- RWM(logpi = logpi_C, nits = 10000, h = 0.8, beta_curr = colMeans(output_id_C$beta_store), X = X, Y = Y_C, covariancemat = sig_C)
  sig_L <- (2.38^2/d) * cov(run_L$beta_store)
  sig_C <- (2.38^2/d) * cov(run_C$beta_store)
  # continually updates sigma to the sample covariance found in the chain
}
# again, need two runs for each model. org is just what i put for "original" i.e. not adaptive

run_1_L <- RWM(logpi = logpi_L, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_L, covariancemat = sig_L)
run_2_L <- RWM(logpi = logpi_L, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_L, covariancemat = sig_L)
run_1_C <- RWM(logpi = logpi_C, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_C, covariancemat = sig_C)
run_2_C <- RWM(logpi = logpi_C, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_C, covariancemat = sig_C)

run_1_L_org <- RWM(logpi = logpi_L, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_L, covariancemat = covariance_matrix_L)
run_2_L_org <- RWM(logpi = logpi_L, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_L, covariancemat = covariance_matrix_L)
run_1_C_org <- RWM(logpi = logpi_C, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_C, covariancemat = covariance_matrix_C)
run_2_C_org <- RWM(logpi = logpi_C, nits = 10000, h = 0.8, beta_curr = rep(0,10), X = X, Y = Y_C, covariancemat = covariance_matrix_C)

l_combined = mcmc.list(mcmc(run_1_L$beta_store), mcmc(run_2_L$beta_store))
gelman.diag(l_combined)

c_combined = mcmc.list(mcmc(run_1_C$beta_store), mcmc(run_2_C$beta_store))
gelman.diag(c_combined)

l_org_combined = mcmc.list(mcmc(run_1_L_org$beta_store), mcmc(run_2_L_org$beta_store))
gelman.diag(l_combined)

c_org_combined = mcmc.list(mcmc(run_1_C_org$beta_store), mcmc(run_2_C_org$beta_store))
gelman.diag(c_combined)

# no real difference in the gelman-rubin statistic, no real reason to change the covariance as the main point of this is to improve mixing and seeing no improvement from the NIH, so continue with it.





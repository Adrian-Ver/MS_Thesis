#####################################
####### parameters, packages ######## ---------
# n      = number of subjects
# time   = number of time points
# rho    = common autoregressive parameter
# p      = proportion of 1s in the data  
# gamma  = vector of individual-specific effects (under FE)
# delta  = vector of time-specific effects (under FE)
# gsd    = gamma standard deviation (under RE)
# dsd    = delta standard deviation (under RE)
library(plm) # for FDGMM
library(mvtnorm) # DGP
library(lme4) # for REML (not used)
library(nlme) # for REML (Elfred)
#####################################
################ DGP ################ -------
# p <- 0.4;n <- 10;time <- 20;gsd <- 2; dsd <- 2
# corrx <- 0.9;rho <- 0.5; distx <- 3
make_ldpd <- function(n,time,rho,p,corrx,distx,gsd,dsd){
    n <<- n
    tt <<- time
    # beta, and generating correlated random variables for X ---------
    set.seed(2018)
    beta <- c(0,2,-2)
    Sigma <- matrix(c(1,corrx,corrx,1),nrow = 2)
    mu <- c(0,0)
    cholS <- chol(Sigma)
    x <- matrix(runif(2*n*time), ncol = 2)
    
    # determining the distribution of X ---------
    switch(distx,
        X <- qnorm(x,0,1),                      # for dist = 1
        X <- qnorm(x,0,0.5),                    # for dist = 2
        X <- qweibull(x,shape = 1,scale = 0.5)) # for dist = 3
    X  <- X %*% cholS + mu
    X1 <- matrix(X[,1],ncol = time)
    X2 <- matrix(X[,2],ncol = time)
    cor(X1,X2) 
    
    # setting up the random effects ----------
    set.seed(2018)
    re_gamma <- rnorm(n, mean = 0, sd = gsd)     # is this the correct way 
    re_delta <- rnorm(time+1, mean = 0, sd = dsd)   # of specifying REs :/
    
    # setting up matrix of outcomes: y, pi, and L -----------
    y   <- matrix(0, nrow = n, ncol = time)
    L   <- matrix(0, nrow = n, ncol = time)
    pi  <- matrix(0, nrow = n, ncol = time)
    L0 <<- log(p/(1-p))
    
    # setting up model, computing for L --------------
    L[,1] <- rho*L0 + beta[1] + beta[2]*X1[,1] + beta[3]*X2[,1]
                + re_gamma + re_delta[1]
    for(t in 2:time){ # random effects model
        L[,t] <- rho*L[,t-1] + beta[1] + beta[2]*X1[,t] + beta[3]*X2[,t] 
                + re_gamma + re_delta[t]
    }
    
    # computing pi; return X, y_{it} and y_{i,t-1} --------------
    pi <- exp(L)/(1+exp(L))
    for(t in 1:time){
        y[,t] <- rbinom(y[,t],1,pi[,t])   
    }
    
    # restructuring to a "proper" data frame
    ldpd <- data.frame(i = rep(1:n, each = time), t = rep(1:time, n), y = c(y), 
                       L = c(L), X1 = c(X1), X2 = c(X2))
    return(ldpd)
}
# (pd <- make_ldpd(n=10,time=5,rho=.5,p=0.9,corrx = 0.1,distx = 1,dsd=1,gsd=1))
# remark: hindi stable yung proportion of 1s pagdating sa time 'time' :<

#####################################
################ MSP ################ -------------
msp <- function(pd){
    pd <- data.frame(i = rep(1:n, each = time), t = rep(1:tt, n), 
                       y = c(pd$y), L = c(pd$L), 
                       X1 = c(pd$X1), X2 = c(pd$X2))
    spacing <- function(beta){
        b0  <- beta[1]; b1 <- beta[2]; b2 <- beta[3]
        # specify the model to the location parameter
        Lhat <- b0+pd$X1*b1+pd$X2*b2
        # sort the obs
        Lhat_t <- c(sort(Lhat),Inf)
        # vector containing the y_(i-1) obs
        Lhat_t1 <- c(-Inf,sort(Lhat))
        # input them to the logistic cdf
        FLhat_t <- plogis(Lhat_t, location = beta[1]+pd$X1*beta[2]+pd$X2*beta[3],scale = 1)
        FLhat_t1 <- plogis(Lhat_t1, location = beta[1]+pd$X1*beta[2]+pd$X2*beta[3],scale = 1)
        # calculate S_n
        S_n <- mean(log(FLhat_t-FLhat_t1))
        # return(S_n)
    }
    lm_est <- lm(L ~ X1 + X2, data = pd)$coeff
    betahat <- optim(lm_est,spacing) 
    # probably not right to have ls estimates as initial values for MSP ?
    return(betahat$par)
}

testing <- make_ldpd(n=100,time=10,rho=0.5,p=.6,corrx = 0.5,distx = 3,dsd=5,gsd=1)
msp(testing) # ampangit ng estimates :(((
lm(L ~ X1 + X2,data = testing)$coeff

#####################################
######## (bootstrap) FD-GMM ######### ------------

# n <- nrow(pd$L); time <- ncol(pd$L)
# pdata <- data.frame(i = rep(1:n, each = time), t = rep(1:time, n), 
#                        y = c(pd$y), L = c(pd$L), 
#                        X1 = c(pd$X1), X2 = c(pd$X2))
# pdata <- pdata.frame(pdata, index = c("i","t"))

###### iid bootstrap procedure
boot_rho <- rep(0,10)
for(b in 1:10){
    sub <- sample(1:n*time, size = n*time, replace = T)
    panel_boot <- pdata[sub,]
    boot_rho[b] <- pgmm(L ~ lag(L,1)|lag(L,1:99), data = panel_boot, 
                    effect = "individual", model = "onestep", 
                    transformation = "d")$coefficient
    return(boot_rho)
} # problem: same y_{it} may be obtained in a bootstrap sample (use sieve? [not yet performed])

###### unreplicated FDGMM #################
fdgmm_unrep <- function(pd){
    rho_hat <- abs(pgmm(L ~ lag(L,1)|lag(L, 1:99), data = pd, 
                    effect = "individual", model = "onestep", 
                    transformation = "d")$coefficient) 
    return(rho_hat)
}
# rho_hat <- abs(pgmm(L ~ lag(L,1)|lag(L, 1:99), data = testing, 
#                     effect = "individual", model = "onestep", 
#                     transformation = "d")$coefficient) 
# is taking the absolute value of a negative coefficient allowed?
# aralin pa ang specifications ng pgmm function

###### manual implementation of the FD-GMM estimator (gg)
# n = 10, time = 20
# time <- (ncol(test)-4)/2
# yit <- t(test[,-c(1:(time+4))]) # replace test with df later on!
# xit <- t(test[,-c(1:(time+3),ncol(test))])
# oneT <- matrix(rep(1,time),ncol=1)
# Q <- diag(rep(1,time))- (1/time)* oneT %*% t(oneT) 
# A <- chol(Q) # chol in R has form A'A = Q
# xstar <- A %*% as.matrix(xit); xstar <- as.matrix(c(xstar),ncol=1)
# ystar <- A %*% as.matrix(yit); ystar <- as.matrix(c(xstar),ncol=1)

#####################################
################ REML ############### ------------
do_reml <- function(pd){
    panel_reml  <- lmer(L ~ (1|i) + (1|t), data = pd)
    # gamma_delta <- panel_reml@theta # gives the relative standard deviation from the residuals 
    # return(gamma_delta)
    return(panel_reml)
}
# remltest <- do_reml(testing)
# mean(fitted(remltest))
# hindi pa kasi nacocompute yung effects ng rho at beta so hindi pa accurate yung shits
# do we need to add the error terms to the model, so that the relative std dev from the residuals can be useful?

#####################################
#### estimating initial L from y #### ------------
estimate_L_from_y <- function(pd){
    pi_hat_i <- rep(0,tt) # estimates of t for a given i
    for(l in 1:tt){
        pi_hat_i[l] <- mean(subset(pd, t == l)$y)
    }
    pi_hat <- rep(pi_hat_i, n) # since initial p (across t) ...
    # ... estimates are the same for every i, we just replicate it for all is
    L_hat <- log(pi_hat/(1-pi_hat))
    return(L_hat)
}
# estimate_L_from_y(testing)
lag(testing$L,1)
#####################################
############ backfitting ############ -----------
backfit <- function(pd, max_iter = 100){
    pi_hat <- matrix(0, ncol = max_iter, nrow = n*tt)
    L_hat  <- matrix(0, ncol = max_iter, nrow = n*tt)
    par_ests <- data.frame(rho = rep(0,max_iter), b0 = rep(0,max_iter),
                           b1 = rep(0,max_iter), b2 = rep(0,max_iter),
                           gamma = rep(0,max_iter), delta = rep(0,max_iter))
    # maybe use nested lists for par_ests? then use apply?
    while(!conv && max_iter){
        iter <- iter + 1
        iter <- 1
        if(iter == 1){
# ______________________________ during 1st loop __________________________
            L_hat[,iter] <- estimate_L_from_y(testing)
            # (1) MSP, and obtaining the residual from the method
            beta_hats <- msp(testing)
            resid_a <- L_hat[,iter]-(beta_hats[1]+beta_hats[2]*testing$X1 + 
                                         beta_hats[3]*testing$X2)
            r_df_a <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                            L = resid_a, X1 = testing$X1, X2 = testing$X2)
            # dapat bang kasama yung L0 sa estimation?
            # (2) we now treat resid_a as our L, and perform FDGMM
            rho_hat <- fdgmm_unrep(r_df_a)
            resid_b <- resid_a - rho_hat*lag(resid_a,1) # mali to pero for the sake of running lang muna
            r_df_b  <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                                L = resid_b, X1 = testing$X1, X2 = testing$X2)
            # (3) estimation of the random effects using REML
            rand_effs <- do_reml(r_df_b)
            gamma_hat <- as.numeric(unlist(coef(rand_effs)$i))
            delta_hat <- as.numeric(unlist(coef(rand_effs)$t))
            resid_c <- resid_b - fitted(rand_effs)
            r_df_c  <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                                L = resid_c, X1 = testing$X1, X2 = testing$X2)
            L_hat[,iter] <- r_df_c$L
# ________________________ after 1st loop (di pa naaayos!) ____________________
        }else{
            # (1) MSP, and obtaining the residual from the method
            iter_data <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                            L = L_hat[,iter-1], X1 = testing$X1, X2 = testing$X2)
            beta_hats <- msp(iter_data)
            resid_a <- L_hat[,iter]-(beta_hats[1]+beta_hats[2]*testing$X1 + 
                                         beta_hats[3]*testing$X2)
            r_df_a <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                            L = resid_a, X1 = testing$X1, X2 = testing$X2)
            # dapat bang kasama yung L0 sa estimation?
            # (2) we now treat resid_a as our L, and perform FDGMM
            rho_hat <- fdgmm_unrep(r_df_a)
            resid_b <- resid_a - rho_hat*lag(resid_a,1) # mali to pero for the sake of running lang muna
            r_df_b  <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                                L = resid_b, X1 = testing$X1, X2 = testing$X2)
            # (3) estimation of the random effects using REML
            rand_effs <- do_reml(r_df_b)
            gamma_hat <- as.numeric(unlist(coef(rand_effs)$i))
            delta_hat <- as.numeric(unlist(coef(rand_effs)$t))
            resid_c <- resid_b - fitted(rand_effs)
            r_df_c  <- data.frame(i = testing$i, t = testing$t, y = testing$y, 
                                L = resid_c, X1 = testing$X1, X2 = testing$X2)
            L_hat[,iter] <- r_df_c$L
        }
    }
}

#####################################
######## MLE (benchmarking) ######### -------------
library(pglm)
#####################################
######### model evaluation ########## ------------
# bias (for each of the params?)
# se(bias)
# sensitivity and specificity
# MAPE
# MSE


# note: gumawa ng error-message for certain conditions!
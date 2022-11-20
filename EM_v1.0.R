library(tidyverse)
library(numDeriv)

#source("density.R")

# E step

Qfunc <- function(x,w1,loc,scale1,scale2,w1t,loct,scale1t,scale2t){
    T1 <- w1t * rgumbel_pdf(x=x,loct,scale1t)
    T2 <- (1-w1t) * lgumbel_pdf(x=x,loct,scale2t)
    Tp <- T1/(T1+T2)
    p1 <- Tp*log(w1)
    p2 <- Tp*rgumbel_pdf(x=x,loc,scale1,log=TRUE)
    p3 <- (1-Tp)*log(1-w1)
    p4 <- (1-Tp)*lgumbel_pdf(x=x,loc,scale2,log=TRUE)
    output <- sum(p1+p2+p3+p4)
    return(output)
}

em <- function(x,w1,loc,scale1,scale2,silent=TRUE,
               lower_loc=-1000,upper_loc=1000,
               lower_scale1=0,upper_scale1=1000,
               lower_scale2=0,upper_scale2=1000){
    w1t <- w1 #_t for current stage value
    loct <- loc
    scale1t <- scale1
    scale2t <- scale2
    delta <- 10
    iteration <- 1
    while (delta >= 1e-10) {
        T1 <- w1t * rgumbel_pdf(x=x,loct,scale1t)
        T2 <- (1-w1t) * lgumbel_pdf(x=x,loct,scale2t)
        Tp <- T1/(T1+T2)
        
        # update w1
        w1n <- mean(Tp)
        
        # update loc
        efunc <- function(locn){
            output <- Qfunc(x=x,
                            w1=w1n,
                            loc=locn,
                            scale1=scale1t,
                            scale2=scale2t,
                            w1t=w1t,
                            loct=loct,
                            scale1t=scale1t,
                            scale2t=scale2t)
            return(-2*output)
        }
        estout <- optim(par=loct,
                        fn=efunc,
                        method="Brent",
                        lower=lower_loc,
                        upper=upper_loc)
        locn <- estout$par
 
        # update scale1
        efunc <- function(scale1n){
            output <- Qfunc(x=x,
                            w1=w1n,
                            loc=locn,
                            scale1=scale1n,
                            scale2=scale2t,
                            w1t=w1t,
                            loct=loct,
                            scale1t=scale1t,
                            scale2t=scale2t)
            return(-2*output)
        }
        estout <- optim(par=scale1t,
                        fn=efunc,
                        method="Brent",
                        lower=lower_scale1,
                        upper=upper_scale1)
        scale1n <- estout$par
        
        # update scale2
        efunc <- function(scale2n){
            output <- Qfunc(x=x,
                            w1=w1n,
                            loc=locn,
                            scale1=scale1n,
                            scale2=scale2n,
                            w1t=w1t,
                            loct=loct,
                            scale1t=scale1t,
                            scale2t=scale2t)
            return(-2*output)
        }
        estout <- optim(par=scale2t,
                        fn=efunc,
                        method="Brent",
                        lower=lower_scale2,
                        upper=upper_scale2)
        scale2n <- estout$par

        # convergence criteria
        delta1 <- Qfunc(x=x,
                        w1=w1t,
                        loc=loct,
                        scale1=scale1t,
                        scale2=scale2t,
                        w1t=w1t,
                        loct=loct,
                        scale1t=scale1t,
                        scale2t=scale2t)
        delta2 <- Qfunc(x=x,
                        w1=w1n,
                        loc=locn,
                        scale1=scale1n,
                        scale2=scale2n,
                        w1t=w1t,
                        loct=loct,
                        scale1t=scale1t,
                        scale2t=scale2t)
        delta <- delta2 - delta1
        # update t stage value
        w1t <- w1n
        loct <- locn
        scale1t <- scale1n
        scale2t <- scale2n
        if (silent==FALSE){
            loglikelihood <- sum(mixture_gumbel_pdf(x,
                                                    w1t,
                                                    loct,
                                                    scale1t,
                                                    scale2t,
                                                    log=TRUE))
            print(data.frame(iteration=iteration,
                             w1=w1n,
                             loc=locn,
                             scale1=scale1n,
                             scale2=scale2n,
                             loglikelihood=loglikelihood))
        }
        iteration <- iteration + 1
    }
    loglikelihood <- sum(mixture_gumbel_pdf(x,w1t,loct,scale1t,scale2t,log=TRUE))
    output <- data.frame(w1=w1t,
                         loc=loct,
                         scale1=scale1t,
                         scale2=scale2t,
                         loglikelihood=loglikelihood)
    return(output)
}

multiple_initial_em <- function(x,w1,loc,scale1,scale2,n_multiple=50,seed=100,...){
    set.seed(100)
    w1m <- seq(0.001,1-0.001,length.out=n_multiple)
    locm <- rnorm(n_multiple,loc,1)
    scale1m <- rgamma(n_multiple,scale1,1)
    scale2m <- rgamma(n_multiple,scale2,1)
    output <- list()
    for (i in 1:n_multiple) {
        paste0("iteration:",i)
        temp <- tryCatch(em(x=x,
                            w1=w1m[i],
                            loc=locm[i],
                            scale1=scale1m[i],
                            scale2=scale2m[i],...),
                         error=function(cond){
                             output <- data.frame(w1=NA,
                                                  loc=NA,
                                                  scale1=NA,
                                                  scale2=NA,
                                                  loglikelihood=NA)
                             return(output)
                         },
                         warning=function(cond){
                             output <- data.frame(w1=NA,
                                                  loc=NA,
                                                  scale1=NA,
                                                  scale2=NA,
                                                  loglikelihood=NA)
                             return(output)
                         })
        output[[i]] <- temp
    }
    output <- bind_rows(output)
    output <- output%>%
        drop_na()%>%
        arrange(loglikelihood)
    output <- output[1,]
    return(output)
}
##sandwich variance estimation
matrixA <- function(x,w1,loc,scale1,scale2){
    margin_loglikelihood <- function(parameters){
        w1 <- parameters[1]
        loc <- parameters[2]
        scale1 <- parameters[3]
        scale2 <- parameters[4]
        output <- sum(mixture_gumbel_pdf(x,w1,loc,scale1,scale2,log=TRUE))
        return(output)
    }
    n <- length(x)
    A <- -hessian(margin_loglikelihood,c(w1,loc,scale1,scale2))*(1/n)
    return(A)
}

matrixB <- function(x,w1,loc,scale1,scale2){
    n <- length(x)
    output <- list()
    for (i in 1:n) {
        B <- grad(function(parameters){
            w1 <- parameters[1]
            loc <- parameters[2]
            scale1 <- parameters[3]
            scale2 <- parameters[4]
            output <- sum(mixture_gumbel_pdf(x[i],w1,loc,scale1,scale2,log=TRUE))
            return(output)
        },
        c(w1,loc,scale1,scale2))
        B <- matrix(B,ncol=1)
        output[[i]] <- B %*% t(B)
    }
    B <- Reduce('+',output)*(1/n)
    return(B)
}

sandwichvar_se <- function(x,w1,loc,scale1,scale2){
    n <- length(x)
    A <- matrixA(x,w1,loc,scale1,scale2)
    A_inv <- solve(A)
    B <- matrixB(x,w1,loc,scale1,scale2)
    output <- sqrt(diag(A_inv%*%B%*%t(A_inv))*(1/n))
    output <- data.frame(w1_se=output[1],
                         loc_se=output[2],
                         scale1_se=output[3],
                         scale2_se=output[4])
    return(output)
}

##test
# df1 <-  read.csv("Gumbel_mixture.csv")
# est <- em(df1$Y,w1=0.5,loc=0,scale1=1,scale2=1,silent=FALSE)
# sandwichvar_se(df1$Y,est$w1,est$loc,est$scale1,est$scale2)
#output <- multiple_initial_em(df1$Y,w1=0.5,loc=0,scale1=1,scale2=1)



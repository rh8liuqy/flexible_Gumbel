##Right Skewed Gumbel pdf
rgumbel_pdf <- function(x,loc,scale,log=FALSE){
  if (log) {
    if (scale <= 0 ){
      print("scale parameter must be > 0!")
      return(-Inf)
    }
    else {
      p1 <- -log(scale)
      z <- (x-loc)/scale
      p2 <- -(z+exp(-z))
      output <- p1+p2
      return(output)
    }
  }
  else {
    if (scale <= 0) {
      print("scale parameter must be > 0!")
      return(0)
    }
    else {
      p1 <- 1/scale
      z <- (x-loc)/scale
      p2 <- exp(-(z+exp(-z)))
      output <- p1*p2
      return(output)
    }
  }
}

##Left Skewed Gumbel pdf
lgumbel_pdf <- function(x,loc,scale,log=FALSE){
  if (log) {
    if (scale <= 0 ){
      print("scale parameter must be > 0!")
      return(-Inf)
    }
    else {
      p1 <- -log(scale)
      z <- (x-loc)/scale
      p2 <- z-exp(z)
      output <- p1+p2
      return(output)
    }
  }
  else {
    if (scale <= 0) {
      print("scale parameter must be > 0!")
      return(0)
    }
    else {
      p1 <- 1/scale
      z <- (x-loc)/scale
      p2 <- exp(z-exp(z))
      output <- p1*p2
      return(output)
    }
  }
}

##Log-sum-exp

log_sum_exp <- function(u,v){
  m <- max(u,v)
  p1 <- m
  p2 <- log(exp(u-m)+exp(v-m))
  output <- p1+p2
  return(p1+p2)
}

##Mixture Gumbel pdf
mixture_gumbel_pdf <- function(x,w1,loc,scale1,scale2,log=FALSE){
  if (log) {
    if ((w1<0) | (w1>1)) {
      print("weight parameter must fall within 0 and 1")
      return(-Inf)
    }
    else if ((w1>=0)*(w1<=1)){
      p1 <- log(w1)+rgumbel_pdf(x,loc,scale1,log=TRUE)
      p2 <- log(1-w1)+lgumbel_pdf(x,loc,scale2,log=TRUE)
      output <- log_sum_exp(p1,p2)
      return(output)
    }
  }
  else {
    if ((w1<0) | (w1>1)) {
      print("weight parameter must fall within 0 and 1")
      return(0)
    }
    else if ((w1>=0)*(w1<=1)){
      p1 <- w1*rgumbel_pdf(x,loc,scale1)
      p2 <- (1-w1)*lgumbel_pdf(x,loc,scale2)
      output <- p1+p2
      return(output)
    }
  }
}

##Right Skewed Gumbel Generation
rgumbel_RNG <- function(n,loc,scale) {
  if (scale <= 0) {
    print("scale parameter must be > 0")
    return(NULL)
  }
  else {
    U <- runif(n=n)
    output <- loc - scale*log(-log(U))
    return(output)
  }
}

##pdf of inverse gamma
dinvgamma <- function(x,shape,scale,log=FALSE){
  alpha <- shape
  beta <- scale
  if (log) {
    if ((shape <= 0)|(scale <= 0)) {
      return(-Inf)
    }
    else if((shape>0)*(scale>0)){
      p1 <- alpha*log(beta)
      p2 <- -lgamma(alpha)
      p3 <- (-alpha-1)*log(x)
      p4 <- -beta/x
      return(p1+p2+p3+p4)
    }

  }
  else {
    if ((shape <= 0)|(scale <= 0)) {
      return(0)
    }
    else if((shape>0)*(scale>0)){
      p1 <- (beta^alpha)/gamma(alpha)
      p2 <- x^(-alpha-1)*exp(-beta/x)
      return(p1*p2)
    }
  }
}

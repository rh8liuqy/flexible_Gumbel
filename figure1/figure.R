rm(list=ls())
library(latticeExtra)
##Right Skewed Gumbel cdf
rgumbel_cdf <- function(x,loc,scale) {
  if (scale > 0) {
    z <- (x-loc)/scale
    output <- exp(-exp(-z))
    return(output)
  }
  else if (scale <= 0){
    print("scale parameter must be > 0!")
    return(0)
  }
}
##Left skewed Gumbel cdf
lgumbel_cdf <- function(x,loc,scale) {
  if (scale > 0) {
    z <- (x-loc)/scale
    output <- 1 - exp(-exp(z))
    return(output)
  }
  else if (scale <= 0){
    print("scale parameter must be > 0!")
    return(0)
  }
}

detfunction <- function(x1,x2,loc1,loc2,scale1,scale2) {
  p11 <- rgumbel_cdf(x = x1,loc = loc1, scale = scale1)
  p12 <- lgumbel_cdf(x = x1,loc = loc2, scale = scale2)
  p21 <- rgumbel_cdf(x = x2,loc = loc1, scale = scale1)
  p22 <- lgumbel_cdf(x = x2,loc = loc2, scale = scale2)
  dmatrix <- matrix(c(p11,p12,p21,p22), nrow = 2, ncol = 2)
  output <- det(dmatrix)
  return(output)
}


df <- data.frame(x = rep(seq(-20,20,1),each=length(seq(-20,20,1))),
                 y = rep(seq(-20,20,1),times=length(seq(-20,20,1))))
df$z <- NA
n <- 1
for (i in seq(-20,20,1)) {
  for (j in seq(-20,20,1)) {
    df$z[n]  <- detfunction(x1 = i, 
                            x2 = j, 
                            loc1 = 0,
                            loc2 = 0,
                            scale1 = 1,
                            scale2 = 5)
    n <- n + 1
  }
  
}
pdf("determinant.pdf",width = 6,height = 6)
wireframe(z ~ x * y, data = df,
          scales = list(arrows = FALSE),
          drape = TRUE, colorkey = TRUE,
          main = 'Determinant of matrix for model identifiability',
          xlab = 'x1',
          ylab = 'x2',
          zlab = "Value")
dev.off()
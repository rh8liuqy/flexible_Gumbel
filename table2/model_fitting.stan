functions {
  real lgumbel_lpdf(real y, real loc, real scale) {
    real z;
    real p1;
    real p2;
    real output;
    z = (y-loc)/scale;
    p1 = -log(scale);
    p2 = z-exp(z);
    output = p1+p2;
    return output;
  }
}

data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real<lower=0,upper=1> w1;
  real theta;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

model {
  w1 ~ beta(1,1);
  sigma1 ~ inv_gamma(1,1);
  sigma2 ~ inv_gamma(1,1);
  for (n in 1:N) {
    target += log_sum_exp(gumbel_lpdf(y[n]|theta,sigma1)+log(w1),
    lgumbel_lpdf(y[n]|theta,sigma2)+log(1-w1));
  }
}

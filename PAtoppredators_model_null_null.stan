data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 int b[res];  // response variable 
 int b2[fis];  // response variable 
}
parameters {
  real I_reserves;
 real I_fished;
}
transformed parameters {
 vector[res] mu;//mean 
 vector[fis] mu2;//mean  
//reserve component
for (i in 1:res){ 
  mu[i] = I_reserves ;
}
//fished component
for (i in 1:fis){ 
  mu2[i] = I_fished;
  }
}
model {
 //priors
 I_fished ~ normal(-5,10);
 I_reserves ~ normal(-5,10);
 //likelihoods  
 for(n in 1:res){
      b[n] ~ bernoulli_logit(mu[n]);
}
 for(n in 1:fis){
      b2[n] ~ bernoulli_logit(mu2[n]);
}
}
generated quantities {
vector[res+fis] log_lik; //log-likelihood (for loo)
vector[res] y_rep;
vector[fis] y_rep2;
 for (n in 1:res) {
 log_lik[n] = bernoulli_logit_lpmf(b[n]| mu[n]);
 y_rep[n] =bernoulli_logit_rng(mu[n]);
}
for (n in 1:fis) {
 log_lik[n+res] = bernoulli_logit_lpmf(b2[n]| mu2[n]);
 y_rep2[n] =bernoulli_logit_rng(mu2[n]);
}
}

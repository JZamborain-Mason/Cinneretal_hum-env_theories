data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 vector[res] b;  // response variable (log-transformed)
 vector[fis] b2;  // response variable (log-transformed)
}
parameters {
 real I_reserves;
 real I_fished;
 real<lower=0> sigma_e; //error sd for biomass reserves
 real<lower=0> sigma_f; //error sd for biomass fished
}
transformed parameters {
 vector[res] mu;//mean biomass  reserves
 vector[fis] mu2;//mean biomass fished
//reserve component
for (i in 1:res){ 
  mu[i] = I_reserves;
}
//fished component
for (i in 1:fis){ 
  mu2[i] = I_fished;
  }
}
model {
 //priors
  sigma_e ~ cauchy(0,2.5); //uninformative prior sd
 sigma_f ~ cauchy(0,2.5); //uninformative prior sd
 I_fished ~ normal(5,10);
 I_reserves ~ normal(5,10);
 //likelihoods  
 for(n in 1:res){
      b[n] ~ normal(mu[n],sigma_e);
}
 for(n in 1:fis){
      b2[n] ~ normal(mu2[n],sigma_f);
}
}
generated quantities {
vector[res+fis] log_lik; //log-likelihood (for loo)
 for (n in 1:res) {
 log_lik[n] = normal_lpdf(b[n]| mu[n], sigma_e);
}
for (n in 1:fis) {
 log_lik[n+res] = normal_lpdf(b2[n]| mu2[n], sigma_f);
}
}

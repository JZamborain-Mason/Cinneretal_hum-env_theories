data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 vector[res] b;  // response variable 
 vector[fis] b2;  // response variable 
}
parameters {
real<lower=0,upper=1> theta;// probability of observing a 0
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
  mu[i] = I_reserves ;
}
//fished component
for (i in 1:fis){ 
  mu2[i] = I_fished;
  }
}
model {
 //priors
  theta~ uniform (0,1); // prior 
  sigma_e ~ cauchy(0,2.5); //uninformative prior sd
 sigma_f ~ cauchy(0,2.5); //uninformative prior sd
 I_fished ~ normal(5,10);
 I_reserves ~ normal(5,10);
 //likelihoods  
 for(n in 1:res){
     if (b[n] == 0)
      1 ~ bernoulli(theta);
        else {
      0 ~ bernoulli(theta);
      b[n] ~ lognormal(mu[n],sigma_e);
  }
}
 for(n in 1:fis){
  if (b2[n] == 0)
      1 ~ bernoulli(theta);
        else {
      0 ~ bernoulli(theta);
      b2[n] ~ lognormal(mu2[n],sigma_f);
  }
}
}
generated quantities {
vector[res+fis] log_lik; //log-likelihood (for loo)
real y_rep[res];
real y_rep2[fis];
for (n in 1:res) {
 if (bernoulli_rng(theta)) 
 y_rep[n] = 0;
 else
 y_rep[n]= lognormal_rng(mu[n],sigma_e);
 if (b[n] == 0)
   log_lik[n] = bernoulli_lpmf(1|theta);
        else {
    log_lik[n] = bernoulli_lpmf(0 | theta)+
                lognormal_lpdf(b[n]| mu[n], sigma_e);
}
}
for (n in 1:fis) {
 if (bernoulli_rng(theta)) 
 y_rep2[n] = 0;
 else
 y_rep2[n]= lognormal_rng(mu2[n],sigma_f);
 if (b2[n] == 0)
   log_lik[n+res] = bernoulli_lpmf(1|theta);
        else {
    log_lik[n+res] = bernoulli_lpmf(0 | theta)+
                lognormal_lpdf(b2[n]| mu2[n], sigma_f);
}
}
}
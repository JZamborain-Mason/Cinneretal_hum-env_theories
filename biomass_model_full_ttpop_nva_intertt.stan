data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 vector[res] b;  // response variable (log-transformed)
 vector[fis] b2;  // response variable (log-transformed)
//explanatory variables for each component
 real ag[res]; //age reserve (only for reserves)
 real si[res]; //predictor reserve size (only for reserves)
 real sa[res]; //predictor sampling area
 real op[res]; //predictor ocean productivity 
 real cs[res]; //sst anom
 real pg[res]; //predictor popgrowth 
 real ttNC[res]; //predictor ttNC
 real ttNP[res]; //predictor ttNP
 real popNC[res]; //predictor popNC
 real popNP[res]; //predictor popNP
 real rl[res]; //predictor reef landings
 real ps[res]; //predictor pop size
 real hdi[res]; //predictor hdi
 int at[res];//predictor atoll
 int rh_b[res];//predictor backreef/lagoon habitat
 int rh_f[res];//predictor flat habitat
 int rh_c[res];//predictor crest habitat
 int cm_pc[res];//predictor census method point count
 int ds[res]; //predictor depth shallow
 int dd[res]; //predictor depth deep
 real sa2[fis]; //predictor sampling area
 real op2[fis]; //predictor ocean productivity 
 real cs2[fis]; //sst anom
 real pg2[fis]; //predictor popgrowth 
 real ttNC2[fis]; //predictor ttNC
 real ttNP2[fis]; //predictor ttNP
 real popNC2[fis]; //predictor popNC
 real popNP2[fis]; //predictor popNP
 real rl2[fis]; //predictor reef landings
 real ps2[fis]; //predictor pop size
 real hdi2[fis]; //predictor hdi
 int at2[fis];//predictor atoll
 int rh_b2[fis];//predictor backreef/lagoon habitat
 int rh_f2[fis];//predictor flat habitat
 int rh_c2[fis];//predictor crest habitat
 int cm_pc2[fis];//predictor census method point count
 int cm_ds[fis]; //predictor census method distance sampling
 int ds2[fis]; //predictor depth shallow
 int dd2[fis]; //predictor depth deep
 int mr[fis]; //predictor restrictions
//random effects level 1: reef clsuter (nested within larger)
 int<lower=1> RC; //number of reef clusters for reserves (groups)
 int<lower=1, upper=RC> prc[res]; //reef clusters id 
 int<lower=1> RC2; //number of reef clusters fished(groups)
 int<lower=1, upper=RC2> prc2[fis]; //reef clusters id 
//random effects level 2:larger
 int<lower=1> R; //number of data regions reserves (groups)
 int<lower=1, upper=R> pr[RC]; //region id for each reef cluster
 int<lower=1> R2; //number of data regions fished(groups)
 int<lower=1, upper=R2> pr2[RC2]; //region id for each reef cluster
}
parameters {
 vector[24] beta; //effect sizes
 real I_reserves;
 real I_fished;
 real<lower=0> sigma_e; //error sd for biomass reserves
 real<lower=0> sigma_f; //error sd for biomass fished
 vector[R] u; // random effects for reserves (deviation between regions)
 vector[R2] u2; // random effects for fished (deviation between regions)
 real<lower=0> sigma_u; //deviation sd 
 real<lower=0> sigma_u2; //deviation sd 
 vector[RC] urc; // random effects for reserve (deviation between reef clusters within regions)
 vector[RC2] urc2; // random effects for remote (deviation between reef clusters within regions)
 real<lower=0> sigma_urc; //deviation sd 
 real<lower=0> sigma_urc2; //deviation sd 
}
transformed parameters {
 vector[res] mu;//mean biomass  reserves
 vector[fis] mu2;//mean biomass fished
 //varying intercepts by reef clusters and regions
  real intercept_reserves_rc[RC];
  real intercept_reserves_r[R];
  real intercept_fished_rc[RC2];
  real intercept_fished_r[R2];
  //compute the varying intercept at the region level
  for(i in 1:R){
    intercept_reserves_r[i] = I_reserves + (sigma_u*u[i]);
   }
  for(i in 1:R2){
    intercept_fished_r[i] = I_fished + (sigma_u2*u2[i]);
   }
  //compute varying intercept at the reef cluster within region level
  for(i in 1:RC){
     intercept_reserves_rc[i] = intercept_reserves_r[pr[i]] + (sigma_urc*urc[i]);
  }
  for(i in 1:RC2){
     intercept_fished_rc[i] = intercept_fished_r[pr2[i]] + (sigma_urc2*urc2[i]);
  }
 
//reserve component
for (i in 1:res){ 
  mu[i] = intercept_reserves_rc[prc[i]] +beta[1]*dd[i]+beta[2]*ds[i]+beta[3]*rh_c[i]+beta[4]*rh_f[i]+beta[5]*rh_b[i]+beta[6]*cm_pc[i]+beta[8]*sa[i]+ beta[9]*si[i]+ beta[10]*ag[i]+ beta[11]*at[i]+ beta[12]*op[i]+ beta[13]*cs[i]+ beta[14]*pg[i]+ beta[15]*ttNC[i]+ beta[16]*ttNP[i]+ beta[17]*rl[i]+ beta[18]*ps[i]+ beta[19]*hdi[i]+ beta[21]*popNC[i]+ beta[22]*popNP[i];
}
//fished component
for (i in 1:fis){ 
  mu2[i] = intercept_fished_rc[prc2[i]]+beta[1]*dd2[i]+beta[2]*ds2[i]+beta[3]*rh_c2[i]+beta[4]*rh_f2[i]+beta[5]*rh_b2[i]+beta[6]*cm_pc2[i]+beta[7]*cm_ds[i]+beta[8]*sa2[i]+ beta[11]*at2[i]+ beta[12]*op2[i]+ beta[13]*cs2[i]+ beta[14]*pg2[i]+ beta[23]*ttNC2[i]+ beta[24]*ttNP2[i]+ beta[17]*rl2[i]+ beta[18]*ps2[i]+ beta[19]*hdi2[i]+ beta[20]*mr[i]+ beta[21]*popNC2[i]+ beta[22]*popNP2[i];
  }
}
model {
 //priors
 beta[1] ~ normal (0,10); //prior slope
 beta[2] ~ normal (0,10); //prior slope
 beta[3] ~ normal (0,10); //prior slope
 beta[4] ~ normal (0,10); //prior slope
 beta[5] ~ normal (0,10); //prior slope
 beta[6] ~ normal (0,10); //prior slope
 beta[7] ~ normal (0,10); //prior slope
 beta[8] ~ normal (0,10); //prior slope
 beta[9] ~ normal (0,10); //prior slope
 beta[10] ~ normal (0,10); //prior slope
 beta[11] ~ normal (0,10); //prior slope
 beta[12] ~ normal (0,10); //prior slope
 beta[13] ~ normal (0,10); //prior slope
 beta[14] ~ normal (0,10); //prior slope
 beta[15] ~ normal (0,10); //prior slope
 beta[16] ~ normal (0,10); //prior slope
 beta[17] ~ normal (0,10); //prior slope
 beta[18] ~ normal (0,10); //prior slope
 beta[19] ~ normal (0,10); //prior slope
 beta[20] ~ normal (0,10); //prior slope
 beta[21] ~ normal (0,10); //prior slope
 beta[22] ~ normal (0,10); //prior slope
 beta[23] ~ normal (0,10); //prior slope
 beta[24] ~ normal (0,10); //prior slope
 sigma_u ~ cauchy (0,2.5); //prior sd for group varying intercept
 sigma_u2 ~ cauchy (0,2.5); //prior sd for group varying intercept
 u ~ normal(0,  1); //prior re (distribution of the varying intercept)
 u2 ~ normal(0,  1); //prior re
 sigma_urc ~ cauchy (0,2.5); //prior sd for group varying intercept
 sigma_urc2 ~ cauchy (0,2.5); //prior sd for group varying intercept
 urc ~ normal(0,  1); //prior re
 urc2 ~ normal(0,  1); //prior re
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

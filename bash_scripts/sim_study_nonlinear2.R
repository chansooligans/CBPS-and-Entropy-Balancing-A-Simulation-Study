library(WeightIt)
library(CBPS)
library(ebal)
library(dplyr)
library(ggplot2)


# Samle size
n = 2000
n_sims = 1000
nlo=TRUE

# Models
model_names = c('Baseline Logistic','Baseline Mahalanobis','CBPS - Over','CBPS - Just','EB - 1','EB - 2')
dgp_names = c('std normally dist covs',
              'std normally dist covs (non-linear ps model)',
              'normally dist covs with unequal variances',
              'normally dist covs with unequal variances (non-linear ps model)',
              'std normally dist covs with 3 count covariates',
              'std normally dist covs with 3 count covariates (non-linear ps model)',
              'normally dist covs with unequal var and 3 count covs',
              'normally dist covs with unequal var and 3 count covs (non-linear ps model)')
wts = c('wt','wt_mh','wt.cbps_over','wt.cbps_just','wt.eb','wt.eb2')

# Inverse Logit Function
inv.logit = function (x) {
  y = 1/(1+exp(-x))
  return(y) 
}


df = list()
SATE = list()
confounders = list()
for(i in 1:4) confounders[[i]] = paste('z',seq(1,4),sep='')
for(i in 5:8) confounders[[i]] = paste('z',seq(1,7),sep='')

#########################
##### DGP
#########################

# Function to Generate Pre-Treatment Covariates
dgp = function(n){
  Z_ = list()
  X_ = list()
  
  # Pre-Treatment Covariates
  z1 = rnorm(n,  0, 1)
  z2 = rnorm(n,  0, 1)
  z3 = rnorm(n,  0, 1)
  z4 = rnorm(n,  0, 1)
  
  # DGP 1: std normally distributed covariates
  Z_[[1]] = cbind(z1, z2, z3, z4)
  X_[[1]] = z1 + z2 + z3 + z4
  
  # DGP 2: std normally distributed covariates (non-linear ps model)
  Z_[[2]] = cbind(z1, z2, z3, z4)
  X_[[2]] = -exp(z1/2) - z2/(1+exp(z1)) + z3 - sqrt(z4^2)
  
  # Pre-Treatment Covariates
  z1 = rnorm(n,  0, 0.5)
  z2 = rnorm(n,  0, 0.9)
  z3 = rnorm(n,  0, 1.1)
  z4 = rnorm(n,  0, 1.2)
  
  # DGP 3: normally distributed covariates with unequal variances
  Z_[[3]] = cbind(z1, z2, z3, z4)
  X_[[3]] = z1 + z2 + z3 + z4
  
  # DGP 4: normally distributed covariates with unequal variances (non-linear ps model)
  Z_[[4]] = cbind(z1, z2, z3, z4)
  X_[[4]] = -exp(z1/2) - z2/(1+exp(z1)) + z3 - sqrt(z4^2)
  
  # Pre-Treatment Covariates
  z1 = rnorm(n,  0, 1)
  z2 = rnorm(n,  0, 1)
  z3 = rnorm(n,  0, 1)
  z4 = rnorm(n,  0, 1)
  z5 = rpois(n, 1)
  z6 = -rbinom(n, 3, 0.8)
  z7 = rchisq(n, df=1)
  
  # DGP 5: std normally distributed covariates with 3 count covariates
  Z_[[5]] = cbind(z1, z2, z3, z4, z5, z6, z7)
  X_[[5]] = z1 + z2 + z3 + z4 + z5 + z6 + z7
  
  # DGP 6: std normally distributed covariates with 3 count covariates (non-linear ps model)
  Z_[[6]] = cbind(z1, z2, z3, z4, z5, z6, z7)
  X_[[6]] = 0.5*exp(z1/2) + z2/(1+exp(z1)) + -.2*z3^2 + z1*z4 - 0.4*sqrt(z5-z6) - 0.2*(z1+1.2*z6)^2 + 0.5*z7
  
  # Pre-Treatment Covariates
  z1 = rnorm(n,  0, 0.5)
  z2 = rnorm(n,  0, 0.9)
  z3 = rnorm(n,  0, 1)
  z4 = rnorm(n,  0, 1.1)
  z5 = rpois(n, 1)
  z6 = -rbinom(n, 2, 0.8)
  z7 = rchisq(n, df=1)
  
  # DGP 7: normally distributed covariates with unequal variances and 3 count covariates
  Z_[[7]] = cbind(z1, z2, z3, z4, z5, z6, z7)
  X_[[7]] = z1 + z2 + z3 + z4 + z5 + z6 + z7
  
  # DGP 8: normally distributed covariates with unequal variances and 3 count covariates (non-linear ps model)
  Z_[[8]] = cbind(z1, z2, z3, z4, z5, z6, z7)
  X_[[8]] = 0.5*exp(z1/2) + z2/(1+exp(z1)) -.2*z3^2 + z1*z4 - 0.4*sqrt(z5-z6) - 0.2*(z1+1.2*z6)^2 + 0.5*z7
  
  return(list(Z_,X_))
}

# Function to simulate DGP given pre-treatment covariates
sim_dgp = function(Z_=Z_[[1]],X_=X_[[1]],sim_no=1,n,B_=c(0.4,0.5,0.3,0.9),plot=TRUE,nlo=FALSE,dgp_name='dgp name here'){
  
  # Propensity Score
  e_Z = inv.logit(X_)
  
  # Generate treatment vector
  # treat = as.integer(e_Z>runif(n,0,1))
  treat = rbinom(n,1,prob=e_Z)
  
  # Treatment Effect
  eff = 4
  
  if(nlo==FALSE){
    Z_2 = Z_  
  }
  
  if(nlo==TRUE){
    Z_2 = Z_
    Z_2[,1] = Z_[,1]^2
    Z_2[,2] = Z_[,2]*Z_[,1]
    Z_2[,3] = Z_[,3]^2
    Z_2[,4] = sqrt(Z_[,4])
  }
  
  # Generate Potential Outcomes
  y_0 = Z_2 %*% B_ + rnorm(n,0,1)
  y_1 = Z_2 %*% B_ + eff + rnorm(n,0,1)
  y = y_0*(1-treat) + y_1*treat
  
  
  # Generate Researcher Dataset
  df.temp = as.data.frame(cbind(y,y_1,y_0, e_Z, treat, Z_))
  colnames(df.temp)[1:3] = c('y','y_1','y_0')
  
  return(df.temp)
}

#########################
##### Generate Data
#########################

dgp_data = dgp(n)
Z_ = dgp_data[[1]]
X_ = dgp_data[[2]]

b_4 = c(0.4,0.5,0.3,0.9)
b_6 = c(0.4,0.5,0.3,0.9,0.1,1,0.3)
B_ = list(b_4,b_4,b_4,b_4,b_6,b_6,b_6,b_6)

for(k in 1:8){
  df[[k]] = sim_dgp(Z=Z_[[k]], X=X_[[k]], sim_no=k ,n=n, B_=B_[[k]], dgp_name=dgp_names[k])
  SATE[[k]] = mean(df[[k]]$y_1[df[[k]]$treat==1] - df[[k]]$y_0[df[[k]]$treat==1])
}

#########################
##### Simulations
#########################


k=2

# For Debugging
print(k)
print('----------------------------------')

# Initialize Matrices to Store Results for Mean Difference and Linear Regression
ci_results  = matrix(nrow=n_sims,
                     ncol=length(model_names),
                     dimnames = list(1:n_sims,model_names))

lin_results = matrix(nrow=n_sims,
                     ncol=length(model_names),
                     dimnames = list(1:n_sims,model_names))

df.temp = df[[k]]

i=1
while(i <= n_sims){
  if(i%%5==0){print(i)}
  formula = as.formula(paste('treat ~ ',paste(confounders[[k]],collapse=' + '),sep=''))
  
  # Generate treatment vector
  df.temp$treat = as.integer(runif(n) <= df.temp$e_Z)
  df.temp$y = df.temp$y_0*(1-df.temp$treat) + df.temp$y_1*df.temp$treat
  
  # Logistic Regression, 1-1 Matching
  mod_match = matchit(formula, 
                      data = df.temp, 
                      method = 'nearest', 
                      distance = 'logit', 
                      replace = TRUE)
  df.temp$wt = mod_match$weights
  
  # Mahalanobis Distance Matching
  mod_mahalo = matchit(formula, 
                       data = df.temp, 
                       method = "nearest", 
                       distance = "mahalanobis", 
                       replace = TRUE)
  df.temp$wt_mh = mod_mahalo$weights
  
  # CBPS Weights (OVER)
  cbps.over = weightit(formula, 
                       data = df.temp, 
                       method = "cbps", 
                       estimand = "ATT", 
                       over=TRUE)
  df.temp$wt.cbps_over = cbps.over$weights
  
  # CBPS Weights (JUST)
  cbps.just = weightit(formula, 
                       data = df.temp, 
                       method = "cbps", 
                       estimand = "ATT", 
                       over=FALSE)
  df.temp$wt.cbps_just = cbps.just$weights
  
  # Entropy Balancing Weights 1
  r1 = NULL
  test1 = try(r1 <- ebalance(Treatment = df.temp$treat, 
                             X = df.temp[,confounders[[k]]]))
  if(class(test1) %in% 'try-error') {
    i = i-1
    next
  } else {
    eb.out = r1
  }
  df.temp$wt.eb = 1
  df.temp[df.temp$treat == 0,'wt.eb'] = eb.out$w  
  
  # Entropy Balancing Weights 2
  r2 = NULL
  test2 = try(r2 <- ebalance(Treatment = df.temp$treat, 
                             X = cbind(df.temp[,confounders[[k]]],
                                       df.temp[,confounders[[k]]]^2)))
  if(class(test2) %in% 'try-error') {next} else {eb.out2 = r2}
  df.temp$wt.eb2 = 1
  df.temp[df.temp$treat == 0,'wt.eb2'] = eb.out2$w   
  
  # Linear Regression
  ##############################
  
  for(j in 1:length(wts)){
    lm.mod = lm(y ~ ., 
                data = df.temp[,c('y','treat',confounders[[k]])], 
                weights = df.temp[,wts[j]])  
    trt.est = summary(lm.mod)$coefficients['treat','Estimate']
    trt.se = summary(lm.mod)$coefficients['treat','Std. Error']
    lin_results[i,j] = trt.est
    ci_results[i,j] = dplyr::between(SATE[[k]],confint(lm.mod)['treat',1],confint(lm.mod)['treat',2])
  }
  if(i%%5==0) print(apply(ci_results[1:i,],2,mean))
  i = i+1
}

save(lin_results, ci_results, SATE, file = paste('sim',k,'_n',n,'_nsims',n_sims,'_nlo',nlo,'.RDATA',sep=''))




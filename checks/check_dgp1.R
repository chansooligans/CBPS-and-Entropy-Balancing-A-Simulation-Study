library(plotly)
packageVersion('plotly')

# p <- plot_ly(df, x = ~z1, y = ~z2, z=~z3, color = ~treat, colors = c('#BF382A', '#0C4B8E')) %>%
#   add_markers() 
# p

nlo=TRUE
ci_results1 = lin_results1 = c()
ci_results2 = lin_results2 = c()

#########################
n = 1000

# Inverse Logit Function
inv.logit = function (x) {
  y = 1/(1+exp(-x))
  return(y) 
}

# Pre-Treatment Covariates
z1 = rnorm(n,  0, 1)
z2 = rnorm(n,  0, 1)
z3 = rnorm(n,  0, 1)
z4 = rnorm(n,  0, 1)

# DGP 1: std normally distributed covariates
Z_ = cbind(z1, z2, z3, z4)
X_ = -exp(z1/2) - z2/(1+exp(z1)) + z3 - sqrt(z4^2)

if(nlo==FALSE){
  Z_2 = Z_  
}

if(nlo==TRUE){
  Z_2 = Z_
  Z_2[,1] = Z_[,1]^2
  Z_2[,2] = Z_[,2]*Z_[,1]
  Z_2[,3] = Z_[,3]^2
  Z_2[,4] = sqrt(Z_[,4]^2)
}

#########################
# Propensity Score
e_Z = inv.logit(X_)

# Generate treatment vector
# treat = as.integer(e_Z>runif(n,0,1))
treat = rbinom(n,1,prob=e_Z)

eff=4
b_4 = c(0.4,0.5,0.3,0.9)
y_0 = Z_2 %*% b_4 + rnorm(n,0,1)
y_1 = Z_2 %*% b_4 + eff + rnorm(n,0,1)
y = y_0*(1-treat) + y_1*treat

df.temp = as.data.frame(cbind(y,y_1,y_0, e_Z, treat, Z_))
colnames(df.temp)[1:3] = c('y','y_1','y_0')

SATE = mean(df.temp$y_1-df.temp$y_0)


ggplot(df.temp, aes(x=e_Z, group=treat, fill = as.factor(treat))) +
  geom_density(alpha=.3) +
  theme_gray() + 
  labs(fill = "Treat")


#########################
# Simulations
#########################
confounders = paste('z',seq(1,4),sep='')
formula = as.formula(paste('treat ~ ',paste(confounders,collapse=' + '),sep=''))

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

lm.mod = lm(y ~ ., 
            data = df.temp[,c('y','treat',confounders)], 
            weights = df.temp[,'wt'])  
trt.est = summary(lm.mod)$coefficients['treat','Estimate']
trt.se = summary(lm.mod)$coefficients['treat','Std. Error']
lin_results1[i] = trt.est
ci_results1[i] = dplyr::between(SATE,confint(lm.mod)['treat',1],confint(lm.mod)['treat',2])


# 
z_tr = cbind(-exp(z1/2),-z2/(1+exp(z1)),z3,sqrt(z4^2))
colnames(z_tr) = c('z1t','z2t','z3t','z4t')
df.temp = cbind(df.temp,z_tr)

lm.mod2 = lm(y ~ ., 
            data = df.temp[,c('y','treat',colnames(z_tr))], 
            weights = df.temp[,'wt'])  
trt.est = summary(lm.mod2)$coefficients['treat','Estimate']
trt.se = summary(lm.mod2)$coefficients['treat','Std. Error']
lin_results2[i] = trt.est
ci_results2[i] = dplyr::between(SATE,confint(lm.mod2)['treat',1],confint(lm.mod2)['treat',2])



rm(list=ls())

confounders = paste('z',seq(1,7),sep='')

inv.logit = function (x) {
  y = 1/(1+exp(-x))
  return(y) 
}

for(i in 1:100){
  print(i)
  n=2000
  z1 = rnorm(n,  0, 1)
  z2 = rnorm(n,  0, 1)
  z3 = rnorm(n,  0, 1)
  z4 = rnorm(n,  0, 1)
  z5 = rpois(n, 1)
  z6 = -rbinom(n, 3, 0.8)
  z7 = rchisq(n, df=1)
  
  # DGP 5: std normally distributed covariates with 3 count covariates
  Z_ = cbind(z1, z2, z3, z4, z5, z6, z7)
  X_ = z1 + z2 + z3 + z4 + z5 + z6 + z7
  
  X_ = z1 + z2 + z3 + z4 + z5 + z6 + z7
  e_Z = inv.logit(X_)
  treat = rbinom(n,1,prob=e_Z)
  
  df.temp = as.data.frame(cbind(e_Z, treat, Z_))
  
  eb.out = ebalance(Treatment = df.temp$treat, X = df.temp[,confounders])
}

print(ggplot(df.temp, aes(x=e_Z, group=treat, fill = as.factor(treat))) +
        geom_density(alpha=.3) +
        theme_gray() + 
        labs(fill = "Treat"))



love.plot(bal.tab(formula, 
                  data = df[[k]][,c('treat',confounders[[k]])], 
                  weights = data.frame(Logistic = get.w(mod_match[[k]]),
                                       Mahalanobis = get.w(mod_mahalo[[k]]),
                                       CBPS.over = get.w(cbps.over[[k]]),
                                       CBPS.just = get.w(cbps.just[[k]]),
                                       EB = df[[k]]$wt.eb,
                                       EB2 = df[[k]]$wt.eb2),
                  m.threshold = 0.1), 
          var.order = "unadjusted",
          abs = TRUE, 
          colors = c("red", "brown", "purple", "cyan", "blue", "orange", "pink"), 
          shapes = c("circle open", "square open", "square open", "triangle open","triangle open", "circle open", "diamond open"),
          stat = c('mean.diffs'),
          line = TRUE) +
  theme_gray() + 
  labs(title=paste('covariate balance: simulation',k),
       subtitle=dgp_names[k]) +
  xlim(0,0.5) 
  




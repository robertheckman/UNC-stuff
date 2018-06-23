## ---- Bayesian hierarchical model - 2014 seed addition ----


source('./UNC Data/Statistics codes/ggplot theme.R')
source('./UNC Data/Seed Addition Experiment/Seed addition 2014 data manipulation.R')

require(tidyverse)
require(R2jags) 
require(ggmcmc)




## ---- Simple Bayesian model, treatments only ----
sink('mod1.R')
cat("
model {
    
  for(i in 1:length(y)) {
    lambda[i] <- exp(alpha + b.fert * Nut[i] + b.spray * Enemy.Exclusion[i] + b.rake * Litter.Removal[i] +
      b.fert.spray * Nut[i] * Enemy.Exclusion[i] + b.fert.rake * Nut[i] * Litter.Removal[i] + 
      b.spray.rake * Enemy.Exclusion[i] * Litter.Removal[i] + 
      b.fert.spray.rake * Nut[i] * Enemy.Exclusion[i] * Litter.Removal[i])
    y[i] ~ dpois(lambda[i])
  }
    

# priors for fixed effects
  alpha ~ dnorm(0, 0.00001)
  b.fert ~ dnorm(0, 0.00001)
  b.spray ~ dnorm(0, 0.00001)
  b.rake ~ dnorm(0, 0.00001)
  b.fert.spray ~ dnorm(0, 0.00001)
  b.fert.rake ~ dnorm(0, 0.00001)
  b.spray.rake ~ dnorm(0, 0.00001)
  b.fert.spray.rake ~ dnorm(0, 0.00001)

} # end model
    
", fill = TRUE)
sink()




# Input lists for JAGs
germ.data.1 <- list(Nut = as.numeric(germ$Nut) - 1,
                    Enemy.Exclusion = as.numeric(germ$Enemy.Exclusion) - 1,
                    Litter.Removal = as.numeric(germ$Litter.Removal) - 1,
                    y = germ$tot.germ.sept)
                    

germ.1 <- jags.model('mod1.R', data = germ.data.1,
                      n.chains = 3)


germ.1.coda <- coda.samples(germ.1, 
                            variable.names = c('alpha', "b.fert", 'b.spray', 'b.rake', 
                                               'b.fert.spray', 'b.fert.rake', 'b.spray.rake', 
                                               'b.fert.spray.rake'), 
                             n.iter = 10000, n.thin = 1)
summary(germ.1.coda)


## ---- Simple Bayesian model, environmental covariates only ---- 
sink('mod2.R')
cat("
model {
    
  for(i in 1:length(y)) {
    lambda[i] <- exp(alpha + b.damage * tot.damage.std[i] + b.light * lightavg.std[i] + 
      b.soil.mois * soil.mois.std[i])
    y[i] ~ dpois(lambda[i])
  }
    
    
# priors for fixed effects
  alpha ~ dnorm(0, 0.01)
  b.damage ~ dnorm(0, 0.01)
  b.light ~ dnorm(0, 0.01)
  b.soil.mois ~ dnorm(0, 0.01)

} # end model
    
    ", fill = TRUE)
sink()


# Input lists for JAGs
germ.data.2 <- list(tot.damage.std = germ$totdamage.standard,
                    lightavg.std = germ$lightavg.standard,
                    soil.mois.std = germ$moisture.standard,
                      y = germ$tot.germ.sept)


germ.2 <- jags.model('mod2.R', data = germ.data.2,
                     n.chains = 3)


germ.2.coda <- coda.samples(germ.2, 
                            variable.names = c('alpha', "b.damage", 'b.light', 'b.soil.mois'), 
                            n.iter = 10000, n.thin = 1)
summary(germ.2.coda)


## ---- Bayesian single-level model with environmental mediators ----

sink('mod3.R')
cat("
model {
    
  for(i in 1:length(y)) {
  lambda[i] <- exp(alpha.germ + b.damage * tot.damage.std[i] + b.light * lightavg.std[i] + 
    b.soil.mois * soil.mois.std[i] + b.fert * Nut[i])
  y[i] ~ dpois(lambda[i])

  tot.damage.std[i] ~ dnorm(mu.damage[i], tau.damage) 
    mu.damage[i] <- alpha.dam + c.fert * Nut[i] + c.spray * Enemy.Exclusion[i] + c.rake * Litter.Removal[i]

  soil.mois.std[i] ~ dnorm(mu.soil.mois[i], tau.soil.mois)
    mu.soil.mois[i] <- alpha.mois + d.fert * Nut[i] + d.spray * Enemy.Exclusion[i] + d.rake * Litter.Removal[i]

  lightavg.std[i] ~ dnorm(mu.light[i], tau.light)
    mu.light[i] <- alpha.light + e.fert * Nut[i] + e.spray * Enemy.Exclusion[i] + e.rake * Litter.Removal[i]
  }
    
    
# priors for seedling response
  alpha.germ ~ dnorm(0, 0.01)
  b.damage ~ dnorm(0, 0.01)
  b.light ~ dnorm(0, 0.01)
  b.soil.mois ~ dnorm(0, 0.01)
  b.fert ~ dnorm(0, 0.01)

# priors for damage response
  alpha.dam ~ dnorm(0, 0.01)
  c.fert ~ dnorm(0, 0.01)
  c.spray ~ dnorm(0, 0.01)
  c.rake ~ dnorm(0, 0.01)
  tau.damage <- 1 / sigma.damage^2
  sigma.damage ~ dunif(0, 50)

# priors for soil moisture response
  alpha.mois ~ dnorm(0, 0.01)
  d.fert ~ dnorm(0, 0.01)
  d.spray ~ dnorm(0, 0.01)
  d.rake ~ dnorm(0, 0.01)
  tau.soil.mois <- 1 / sigma.soil.mois^2
  sigma.soil.mois ~ dunif(0, 50)

# priors for light response
  alpha.light ~ dnorm(0, 0.01)
  e.fert ~ dnorm(0, 0.01)
  e.spray ~ dnorm(0, 0.01)
  e.rake ~ dnorm(0, 0.01)
  tau.light <- 1 / sigma.light^2
  sigma.light ~ dunif(0, 50)
  
} # end model
    
", fill = TRUE)
sink()



# Input lists for JAGs
germ.data.3 <- list(tot.damage.std = germ$totdamage.standard,
                    lightavg.std = germ$lightavg.standard,
                    soil.mois.std = germ$moisture.standard,
                    Nut = as.numeric(germ$Nut) - 1,
                    Enemy.Exclusion = as.numeric(germ$Enemy.Exclusion) - 1,
                    Litter.Removal = as.numeric(germ$Litter.Removal) - 1,
                    y = germ$tot.germ.sept)


germ.3 <- jags.model('mod3.R', data = germ.data.3,
                     n.chains = 3)


germ.3.coda <- coda.samples(germ.3, 
                            variable.names = c('alpha.germ', "b.damage", 'b.light', 'b.soil.mois', 'b.fert', 
                                               'alpha.dam', 'c.fert', 'c.spray', 'c.rake', 
                                               'alpha.mois', 'd.fert', 'd.spray', 'd.rake', 
                                               'alpha.light', 'e.fert', 'e.spray', 'e.rake'), 
                            n.iter = 10000, n.thin = 1)
summary(germ.3.coda)



## ---- Mass per seed Bernoulli-gamma mixture ----
sink('mod.mix1.R')
cat("
  model {

# Likelihood  
  for(i in 1 : length(w)) {
  w[i] ~ dbern(p[i])
    logit(p[i]) <- c0 + c.fert * Nut[i] + c.spray * Enemy.Exclusion[i] + c.rake * Litter.Removal[i]

  y[i] ~ dgamma(mu[i]^2 / sigma^2, mu[i] / sigma^2)
    mu[i] <- exp(b0 + b.fert * Nut[i] + b.spray * Enemy.Exclusion[i] + b.rake * Litter.Removal[i])
  }

# Priors
  b0 ~ dnorm(0, 1)
  b.fert ~ dnorm(0, 0.01)
  b.spray ~ dnorm(0, 0.01)
  b.rake ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 50)

  c0 ~ dnorm(0, 0.01)
  c.fert ~ dnorm(0, 0.01)
  c.spray ~ dnorm(0, 0.01)
  c.rake ~ dnorm(0, 0.01)


} # end model
    
    ", fill = TRUE)
sink()


germ.mass.data.1 <- list(Nut = as.numeric(germ$Nut) - 1,
                         Enemy.Exclusion = as.numeric(germ$Enemy.Exclusion) - 1,
                         Litter.Removal = as.numeric(germ$Litter.Removal) - 1,
                         y = 10 * germ$massperseed, 
                         w = ifelse(germ$massperseed == 0, 0, 1))


germ.mass.1 <- jags.model('mod.mix1.R', data = germ.mass.data.1,
                          n.chains = 3)

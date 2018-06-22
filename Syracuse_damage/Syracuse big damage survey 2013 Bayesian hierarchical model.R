###############################################################################################
### July 2014 Syracuse big damage survey, 25 leaves per plant, ~ 200 plants total.          ###
### Chlorophyll data collected summer 2014 on leaves that I had tagged, several issues.     ###
### SLA and tissue N were from leaves I collected, 3 per plant.                             ###
### We didn't measure wet area, so SLA may be a bit off. Best to use N by mass.             ###
### Lignin is from Insu and comes from senesced leaves, at the species level.               ###
### According to Insu, lignin.area is likely more indicative of lignin in fresh leaves.     ###
### US range is based on occurence in USDA Plants database. I included the area of a        ###
### state if there was a record of the species in that state.                               ###
### Congeners is the number of species of that genus in an area - NYS, Onandoga Co, or      ###
### Onandoga Co and directly adjacent counties.                                             ###
### Factor analysis included 3 traits that could represent the latent variable, phenotype.  ###
### The factors don't really load in a clean manner, not sure it's useful.                  ###
###############################################################################################



#must have R2jags package installed in R, and ggmcmc if you want to run some chain diagnostics (not essential)
require(R2jags) 
require(ggmcmc)
require(ape) #for phylogeny
require(geiger) #for matching tree and traits
require(tidyverse)


## -----------------------------------------------------*
## ---- Data Manipulation ----
## -----------------------------------------------------*

damage <- read.csv('./UNC Data/Syracuse Damage Surveys/Syracuse big damage survey 8_12_2013_w zeros.csv', header=T) %>%
  dplyr::mutate(fungal.damage = Powdery.mildew + Black.leaf.spot + Black.residue + Red.leaf.spot + Grey.residue + 
                  Black.blotch + Brown.leaf.spot + Brown.senesced.edges + Rust + Cercospora.leaf.spot + Yellow.leaf.spot + 
                  Brown.residue) %>%
  dplyr::mutate(insect.damage = Mining + Chewing + Scraping + Thrips + Leaf.curl + Skeleton + Leaf.roll + Leaf.gall + Yellow.leaf.spot.sucking) %>%
  dplyr::mutate(chewing.damage = Chewing + Skeleton) %>%
  dplyr::mutate(insect.sub.damage = Mining + Scraping + Thrips + Leaf.curl + Leaf.roll + Leaf.gall + Yellow.leaf.spot.sucking) %>%
  dplyr::group_by(Plot, Species, ID) %>%
  dplyr::summarise(obs = length(fungal.damage), fungal = mean(fungal.damage), insect = mean(insect.damage),
                   chewing = mean(chewing.damage), sub.insect = mean(insect.sub.damage)) %>%
  dplyr::arrange(Species, ID) %>%
  merge(read.csv('./UNC Data/Syracuse Damage Surveys/Syracuse shrub provenance info.csv', header = T), ., by = c("Plot", "Species", "ID"), all = TRUE) %>%
  merge(., read.csv('./UNC Data/Syracuse Damage Surveys/Syracuse SLA CN 2013.csv', header = T), by = 'ID', all = TRUE) %>%
  dplyr::select( - X) %>%
  dplyr::mutate(total.damage = fungal + insect) %>%
  merge(., read.table('./UNC Data/Syracuse Damage Surveys/Syracuse Lignin 2013.csv', header=T, sep=','), by = 'Species', all = TRUE) %>%
  dplyr::select(-Shade.tolerant) %>%
  merge(., read.csv('./UNC Data/Syracuse Damage Surveys/Syracuse ranges.csv', header=T), by = c('Species', 'Group', 'Provenance', 'Genus'), all = TRUE) %>%
  dplyr::mutate(Plot = factor(Plot)) %>%
  dplyr::select(-growthform, -taxon, -obs) %>%
  merge(., read.csv('./UNC Data/Syracuse Damage Surveys/Syracuse 2013_shrub_chlorophyll.csv', header=T) %>%
          dplyr::mutate(Plot = factor(Plot)) %>%
          tidyr::gather(key = measure, value = chloro, measure.1:measure.5) %>%
          dplyr::mutate(measure = as.numeric(measure)) %>%
          dplyr::group_by(ID, Plot, Date.Collected, Tag.Color) %>%
          dplyr::summarise(mean.chloro = mean(chloro, na.rm = TRUE)) %>%
          dplyr::group_by(ID, Plot, Tag.Color) %>%
          dplyr::summarise(peak.chloro = max(mean.chloro)) %>%
          dplyr::group_by(ID, Plot) %>%
          dplyr::summarise(mean.peak.chloro = mean(peak.chloro, na.rm = TRUE)),
        by = c('ID', 'Plot'), all = TRUE) %>%
  na.omit(.)


# Performs correlations amongs all numerical variables for which there are pairwise observations
# No surprises, everything looks reasonably uncorrelated, including leaf functional traits
cor(damage[,unlist(lapply(damage, is.numeric))], use = 'pairwise.complete.obs')


## -----------------------------------------------------*
## ---- Fixed Effects Bayesian model ----
## -----------------------------------------------------*


# standardize indepedent continuous variables
standard <- function(x) (x - mean(x, na.rm = T)) / (sd(x, na.rm = T))


sink('mod1.R')
cat("
model {

  for(i in 1:length(y)) {
    mu[i] <- b0 + b.prov * prov[i] + b.range * range[i] + b.leafn * leafn[i] + b.lignin * lignin[i] +
    b.provleafn * prov[i] * leafn[i] + b.provlignin * prov[i] * lignin[i] + 
    b.provleafnlignin * prov[i] * leafn[i] * lignin[i]
  # Likelihood
    y[i] ~ dnorm(mu[i], tau)
  }


#priors
  b0 ~ dnorm(0, 0.00001)
  b.prov ~ dnorm(0, 0.00001)
  b.range ~ dnorm(0, 0.00001)
  b.leafn ~ dnorm(0, 0.00001)
  b.lignin ~ dnorm(0, 0.00001)
  b.provleafn ~ dnorm(0, 0.00001)
  b.provlignin ~ dnorm(0, 0.00001)
  b.provleafnlignin ~ dnorm(0, 0.00001)
  tau <- 1 / sigma^2
  sigma ~ dunif(0, 50)

} #end model

", fill = TRUE)
sink()


fungal.data <- list(prov = as.numeric(damage$Provenance) - 1,
                    range = standard(damage$range),
                    leafn = standard(damage$tissue.N),
                    lignin = standard(damage$lignin.mass),
                    y = car::logit(damage$fungal)) 
                    

fungal <- jags.model('mod1.R', data = fungal.data, 
                     n.chains = 3)

fungal.coda <- coda.samples(fungal, 
                            variable.names = c("b0", "b.prov", 'b.range', 'b.leafn', 'b.lignin', 'b.provleafn', 'b.provlignin', 'b.provleafnlignin', "sigma"), 
                            n.iter = 10000, n.thin = 1)
  summary(fungal.coda)
  



## ----------------------------------------------------------------*
## ---- Mixed Effects Bayesian model, Species-level intercepts ----
## ----------------------------------------------------------------*


sink('mod2.R')
cat("
model {

  for(j in 1:N.species) {
    alpha[j] ~ dnorm(mu.species, tau.species) # species effect (mean and variance)
  }


  for(i in 1:length(y)) {
  	mu[i] <-  alpha[species[i]] + b.prov * prov[i] + b.range * range[i] + b.leafn * leafn[i] + b.lignin * lignin[i] +
      b.provleafn * prov[i] * leafn[i] + b.provlignin * prov[i] * lignin[i] + 
      b.provleafnlignin * prov[i] * leafn[i] * lignin[i] # fixed effects
    y[i] ~ dnorm(mu[i], tau)
  }	

  
# priors for variance components
  tau.species <- sigma.species^-2 
  sigma.species ~ dunif(0, 100) 
  mu.species ~ dnorm(5, 1)


# priors for fixed effects
    b.prov ~ dnorm(0, 0.00001)
    b.range ~ dnorm(0, 0.00001)
    b.leafn ~ dnorm(0, 0.00001)
    b.lignin ~ dnorm(0, 0.00001)
    b.provleafn ~ dnorm(0, 0.00001)
    b.provlignin ~ dnorm(0, 0.00001)
    b.provleafnlignin ~ dnorm(0, 0.00001)
    tau <- 1 / sigma^2
    sigma ~ dunif(0, 50)

} # end model

", fill = TRUE)
sink()


  

# Input lists for JAGs
  
fungal2.data <- list(prov = as.numeric(damage$Provenance) - 1,
                     range = standard(damage$range),
                     leafn = standard(damage$tissue.N),
                     lignin = standard(damage$lignin.mass),
                     y = car::logit(damage$fungal), 
                     # Turn random effects (family, genus, species) into numeric vectors
                     species = as.numeric(factor(damage$Species)),
                     # Get length for each group
                     N.species = length(unique(damage$Species)))


fungal2 <- jags.model('mod2.R', data = fungal2.data,
                     n.chains = 3)

fungal2.coda <- coda.samples(fungal2, 
                             variable.names = c('alpha', "b.prov", 'b.range', 'b.leafn', 'b.lignin', 'b.provleafn', 'b.provlignin', 'b.provleafnlignin', "sigma"), 
                             n.iter = 10000, n.thin = 1)
summary(fungal2.coda)

  

## ----------------------------------------------------------------------------*
## ---- Mixed Effects Bayesian model, Species-level intercepts and effects ----
## ----------------------------------------------------------------------------*


sink('mod3.R')
cat("
model {

  for(j in 1:N.species) {
    alpha[j] ~ dnorm(mu.species[j], tau.species) # species effect (mean and variance)
    mu.species[j] <- kappa + c.prov * prov[j] + c.range * range[j]
  }
    
  for(i in 1:length(y)) {
    mu[i] <-  alpha[species[i]] + b.leafn * leafn[i] + b.lignin * lignin[i] +
    b.provleafn * prov[i] * leafn[i] + b.provlignin * prov[i] * lignin[i] + 
    b.provleafnlignin * prov[i] * leafn[i] * lignin[i] # fixed effects
    y[i] ~ dnorm(mu[i], tau)
  }	
    

# priors for variance components
  tau.species <- sigma.species^-2 
  sigma.species ~ dunif(0, 100) 
  kappa ~ dnorm(0, 0.00001)
  
    
# priors for fixed effects
  c.prov ~ dnorm(0, 0.00001)
  c.range ~ dnorm(0, 0.00001)
  b.leafn ~ dnorm(0, 0.00001)
  b.lignin ~ dnorm(0, 0.00001)
  b.provleafn ~ dnorm(0, 0.00001)
  b.provlignin ~ dnorm(0, 0.00001)
  b.provleafnlignin ~ dnorm(0, 0.00001)
  tau <- 1 / sigma^2
  sigma ~ dunif(0, 50)
    
} # end model
    
", fill = TRUE)
sink()


# Input lists for JAGs

fungal3.data <- list(prov = as.numeric(damage$Provenance) - 1,
                     range = standard(damage$range),
                     leafn = standard(damage$tissue.N),
                     lignin = standard(damage$lignin.mass),
                     y = car::logit(damage$fungal), 
                     # Turn random effects (family, genus, species) into numeric vectors
                     species = as.numeric(factor(damage$Species)),
                     # Get length for each group
                     N.species = length(unique(damage$Species)))


fungal3 <- jags.model('mod3.R', data = fungal3.data, 
                      n.chains = 3)

fungal3.coda <- coda.samples(fungal3, 
                             variable.names = c('alpha', "c.prov", 'c.range', 'b.leafn', 'b.lignin', 'b.provleafn', 'b.provlignin', 'b.provleafnlignin', "sigma"), 
                             n.iter = 10000, n.thin = 1)
summary(fungal3.coda)


## ----------------------------------------------------------------------------------------*
## ---- Mixed Effects Bayesian model, Genus-level intercept plus Species-level effects ----
## ----------------------------------------------------------------------------------------*


sink('mod4.R')
cat("
model {

  for(k in 1:N.genus) {
    kappa[k] ~ dnorm(mu.genus, tau.genus) # genus random effect
  }	
    
    
  for(j in 1:N.species) {
    alpha[j] ~ dnorm(mu.species[j], tau.species) # species effect (mean and variance)
    mu.species[j] <- kappa[genus[j]] + c.prov * prov[j] + c.range * range[j]
  }

    
  for(i in 1:length(y)) {
    mu[i] <-  alpha[species[i]] + b.leafn * leafn[i] + b.lignin * lignin[i] +
    b.provleafn * prov[i] * leafn[i] + b.provlignin * prov[i] * lignin[i] + 
    b.provleafnlignin * prov[i] * leafn[i] * lignin[i] # fixed effects
    y[i] ~ dnorm(mu[i], tau)
  }	


# Derived quantities
  for (p in 1:N.species) {
    for (m in 1:length(N)) {
    native.leafn[m] <- alpha[p] + c.prov * 1 + b.leafn * N[m] + b.provleafn * 1 * N[m]
    exotic.leafn[m] <- alpha[p] + c.prov * 0 + b.leafn * N[m] + b.provleafn * 0 * N[m]
    }
  }
    
    
# priors for variance components
  mu.genus ~ dnorm(0, 0.00001)
  tau.genus <- sigma.genus^-2
  sigma.genus ~ dunif(0, 100)
  tau.species <- sigma.species^-2 
  sigma.species ~ dunif(0, 100) 

    
# priors for fixed effects
  c.prov ~ dnorm(0, 0.00001)
  c.range ~ dnorm(0, 0.00001)
  b.leafn ~ dnorm(0, 0.00001)
  b.lignin ~ dnorm(0, 0.00001)
  b.provleafn ~ dnorm(0, 0.00001)
  b.provlignin ~ dnorm(0, 0.00001)
  b.provleafnlignin ~ dnorm(0, 0.00001)
  tau <- 1 / sigma^2
  sigma ~ dunif(0, 50)



    
} # end model
  
", fill = TRUE)
sink()


# Input lists for JAGs

fungal4.data <- list(prov = as.numeric(damage$Provenance) - 1,
                     range = standard(damage$range),
                     leafn = standard(damage$tissue.N),
                     lignin = standard(damage$lignin.mass),
                     y = car::logit(damage$fungal), 
                     # Turn random effects (family, genus, species) into numeric vectors
                     genus = as.numeric(factor(damage$Genus)),
                     species = as.numeric(factor(damage$Species)),
                     # Get length for each group
                     N.genus = length(unique(damage$Genus)),
                     N.species = length(unique(damage$Species)), 
                     # For derived quantity - provenance x leaf N 
                     N = seq(-3, 3, by = 0.01))


fungal4 <- jags.model('mod4.R', data = fungal4.data, #inits = fungal.inits, #parameters.to.save = fungal.parameters, 
                      n.chains = 3)

fungal4.coda <- coda.samples(fungal4, 
                             variable.names = c('alpha', "c.prov", 'c.range', 'b.leafn', 'b.lignin', 'b.provleafn', 'b.provlignin', 'b.provleafnlignin', "sigma"), 
                             n.iter = 10000, n.thin = 1)
  summary(fungal4.coda)
  
fungal4.jags <- jags.samples(fungal4, 
                             variable.names = c('exotic.leafn', 'native.leafn'), 
                             n.iter = 10000, n.thin = 1)
b <- summary(fungal4.jags$exotic.leafn, quantile, c(0.025, 0.5, 0.975))$stat
c <- summary(fungal4.jags$native.leafn, quantile, c(0.025, 0.5, 0.975))$stat
  plot(N, b[2,], xlab = "Foliar N content", ylab = "Foliar fungal damage", type = "l")
  lines(N, c[2,], lty = 'dashed')
  lines(N, b[1,], lty = "dashed")
  lines(N, b[3,], lty = "dashed")

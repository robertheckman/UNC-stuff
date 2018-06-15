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

#this R script interacts with JAGS, which must be installed from here:
#http://mcmc-jags.sourceforge.net/
#no need to open JAGS separately when you run R; it works in the background

#must have R2jags package installed in R, and ggmcmc if you want to run some chain diagnostics (not essential)
require(R2jags) 
require(ggmcmc)
require(ape) #for phylogeny
require(geiger) #for matching tree and traits
require(plyr)
require(dplyr)
require(tidyr)

#set working directory
#setwd("C:/Users/rwheckma/Dropbox/Back up/Research/UNC Data/Syracuse Damage Surveys") 


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

  for(i in 1:N) {
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
                    y = car::logit(damage$fungal), 
                    N = length(y))

fungal <- jags.model('mod1.R', data = fungal.data, #inits = fungal.inits, #parameters.to.save = fungal.parameters, 
                     n.chains = 3)

fungal.coda <- coda.samples(fungal, 
                            variable.names = c("b0", "b.prov", 'b.range', 'b.leafn', 'b.lignin', 'b.provleafn', 'b.provlignin', 'b.provleafnlignin', "sigma"), 
                            n.iter = 10000, n.thin = 1)
  summary(fungal.coda)
  



## --------------------------------------------------------------*
## ---- Mixed Effects Bayesian model, Species, Genus, Family ----
## --------------------------------------------------------------*

# Turn random effects (family, genus, species) into numeric vectors
family = as.numeric(factor(damage$family))
genus = as.numeric(factor(damage$Genus))
species = as.numeric(factor(damage$Species))


# Get length for each group
N.family = length(unique(family))
N.genus = length(unique(genus))
N.species = length(unique(species))



y <- damage$fungal
N <- length(y)


sink('mod2.R')
cat("
    model {
    #likelihood
    for(i in 1:N) {
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- b0 + b.prov * prov[i] + b.range * range[i] + b.leafn * leafn[i] + b.lignin * lignin[i] +
    b.provleafn * prov[i] * leafn[i] + b.provlignin * prov[i] * lignin[i] + 
    b.provleafnlignin * prov[i] * leafn[i] * lignin[i]
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
    tau.y <- pow(sigma.y, -2)
    sigma.y ~ dunif(0,10000)
    } #end model
    
    ", fill = TRUE)
sink()


mod.2 <- "model
{
  for(i in 1:N) { #loop over observations
  y[i] ~ dnorm(mu[i], x.tau) #leafout times are normally distributed
  mu[i] <- b0 			#overall intercept
  + b.family[family[i]] 		#family effect
  + b.genus[genus[i]] 		#genus effect
  + b.species[species[i]] 	#species effect
  }
  

  for(i in 1:N.family) {
    b.family[i] ~ dnorm(0, family.tau) # family effect (mean, variance)
  }	

  
  for(i in 1:N.genus) {
    b.genus[i] ~ dnorm(0, genus.tau) # genus random effect
  }	

  
  for(i in 1:N.species) {
    b.species[i] ~ dnorm(0, species.tau) # species effect (mean and variance)
  	mu.species[i] <- b.prov * prov[i] + b.range * range[i] + b.leafn * leafn[i] + b.provleafn * prov[i] * leafn[i] # fixed effects
  }	

  
  #priors  
#priors for variance components
  x.tau <- x.sigma^-2
  x.sigma ~ dunif(0, 100)  
  family.tau <- family.sigma^-2
  family.sigma ~ dunif(0, 100)
  genus.tau <- genus.sigma^-2
  genus.sigma ~ dunif(0, 100)
  species.tau <- species.sigma^-2 
  species.sigma ~ dunif(0, 100) 

#priors for fixed effects
  b0 ~ dnorm(0, 0.00001)  #intercept
  b.prov ~ dnorm(0, 0.00001)
  b.range ~ dnorm(0, 0.00001)
  b.leafn ~ dnorm(0, 0.00001)
  b.provleafn ~ dnorm(0, 0.00001)

}" #end model

  # write model to text file for JAGS
  write(mod.2, "model2.txt")
  

  #input lists for JAGs
  fungal.data2 <- list(N = N, y = y, prov = prov, range = range, leafn = leafn,  species = species, N.species = N.species, 
                       genus = genus, N.genus = N.genus, family = family, N.family = N.family)
  fungal.inits2 <- function() {list(b0 = rnorm(1))}
  fungal.parameters2 <- c("b0", "b.prov", 'b.range', 'b.leafn', 'b.provleafn')
  
  fungal.2 <- jags(fungal.data2, fungal.inits2, fungal.parameters2, "model2.txt", n.chains = 3, n.iter = 1000, n.burnin = 500)
  fungal.2
  
  s <- ggs(as.mcmc(fungal.2))
  ggs_traceplot(s, family = "b0") #examine to ensure no bias among runs (overlapping lines)
  
  #examine posteriors
  attach.jags(fungal.2) #makes posterior distributions an object
  

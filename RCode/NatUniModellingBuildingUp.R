# Code for gradually building up the saturated model and the Variable Selection.
# Uses simulated datasets.

#READ THE DATA

dfmat= read.csv("simulInter.csv",header=FALSE, sep = ",") 
#----------------------OR (for models without interactions)--------------------------
#dfmat= read.csv("./simulNoInter.csv",header=FALSE, sep = ",") 


df= data.frame(BMIdf= dfmat[,1], alphadf= factor(dfmat[,2]), betadf=factor(dfmat[,3]) )
lalpha=length(levels(df$alphadf)) ## Number of levels of the grouping factor alpha
lbeta=length(levels(df$betadf))  ## Number of levels of the grouping factor beta

require(rjags)
x=data.matrix(df, rownames.force = NA)

data=list( X=x, n=nrow(df), alphalength= lalpha, betalength=lbeta)
init=list(mu=0,tau=1)

#-----------1. Main effects only as Fixed effects - Saturated model --------------
# Using sum to 0 constraint for facing the identifiability problem.
modelstring="

model {

for(i in 1:n){
  X[i,1]~dnorm(mu +alpha[X[i,2]]+ beta[X[i,3]], tau )
}

for (j in 1: (alphalength-1)){
  alpha[j]~ dnorm(0,0.0001) #fixed effect, from a fixed prior distribution
}

for (k in 1: (betalength-1)){
  beta[k]~ dnorm(0,0.0001) #fixed effect, from a fixed prior distribution
}
alpha[alphalength] <- -sum(alpha[1:(alphalength-1)]) 
beta[betalength] <- -sum(beta[1:(betalength-1)])
mu~dnorm(0,0.0001)
tau~dgamma(1,0.0001)
}
"

model=jags.model(textConnection(modelstring),data=data,inits=init)
update(model,n.iter=1000)
output=coda.samples(model=model,variable.names=c("mu","tau","alpha", "beta", "int_effects"),n.iter=100000,thin=10)
summary(output)
plot(output)

#-----------2. Main effects only as Random effects - Saturated model --------------
# Random effects by using a gamma prior on the precision of the effects.
modelstring="

model {

for(i in 1:n){
  X[i,1]~dnorm(mu +alpha[X[i,2]]+ beta[X[i,3]], tau )
}

for (j in 1: (alphalength)){
  alpha[j]~ dnorm(0,taualpha)  #random effect, gamma prior on precision
}

for (k in 1: (betalength)){
  beta[k]~ dnorm(0,taubeta) #random effect, gamma prior on precision
}

mu~dnorm(0,0.0001)
tau~dgamma(1,0.0001)
taualpha~ dgamma(1,0.0001)
taubeta~ dgamma(1,0.0001)

}
"

model=jags.model(textConnection(modelstring),data=data,inits=init)
update(model,n.iter=1000)
output=coda.samples(model=model,variable.names=c("mu","tau","alpha", "beta", "int_effects"),n.iter=1000000,thin=100)
summary(output)
plot(output)
autocorr.plot(output)

#-----------3. For main effects + interactions fixed- Saturated model --------------
modelstring="
model {

for(i in 1:n){
  X[i,1]~dnorm(mu +alpha[X[i,2]]+ beta[X[i,3]]+int_effects[X[i,2],X[i,3]], tau )
}

for (j in 1: (alphalength-1)){
  alpha[j]~ dnorm(0,0.0001) #fixed effect, from a fixed prior distribution
}

for (k in 1: (betalength-1)){
  beta[k]~ dnorm(0,0.0001) #fixed effect, from a fixed prior distribution
}

for(j in 1: (alphalength-1)){
  for(k in 1: (betalength-1)){
    int_effects[j,k] ~ dnorm(0,0.0001)
  }
}
alpha[alphalength] <- -sum(alpha[1:(alphalength-1)]) 
beta[betalength] <- -sum(beta[1:(betalength-1)])

for(j in 1: (alphalength-1)){
    int_effects[j,betalength]<- -sum(int_effects[j,1:(betalength-1)]) 
}
for(k in 1: (betalength-1)){
    int_effects[alphalength,k]<- -sum(int_effects[1:(alphalength-1),k]) 
}

int_effects[alphalength,betalength] <- -sum(int_effects[alphalength,1:(betalength-1)])

mu~dnorm(0,0.0001)
tau~dgamma(1,0.0001)
}
"
model=jags.model(textConnection(modelstring),data=data,inits=init)
update(model,n.iter=1000)
output=coda.samples(model=model,variable.names=c("mu","tau","alpha", "beta", "int_effects"),n.iter=100000,thin=10)
summary(output)
plot(output)
autocorr.plot(output)

#-----------4. Main effects and interactions  Random - Saturated model --------------
modelstring="

model {

for(i in 1:n){
  X[i,1]~dnorm(mu +alpha[X[i,2]]+ beta[X[i,3]]+ int_effects[X[i,2],X[i,3]], tau )
}

for (j in 1: (alphalength)){
  alpha[j]~ dnorm(0,taualpha)  #random effect, gamma prior on precision
}

for (k in 1: (betalength)){
  beta[k]~ dnorm(0,taubeta) #random effect, gamma prior on precision
}

for(j in 1: alphalength){
  for(k in 1: betalength){
    int_effects[j,k] ~ dnorm(0,tau_int) #random effect, gamma prior on precision
  }
}

mu~dnorm(0,0.0001)
tau~dgamma(1,0.0001)
taualpha~ dgamma(1,0.0001)
taubeta~ dgamma(1,0.0001)
tau_int~ dgamma(1,0.0001)
}
"

model=jags.model(textConnection(modelstring),data=data,inits=init)
update(model,n.iter=1000)
output=coda.samples(model=model,variable.names=c("mu","tau","alpha", "beta", "int_effects"),n.iter=1000000,thin=300)
summary(output)
plot(output)
autocorr.plot(output)

#-----------5. For random effects + variable selection --------------
modelstring="

model {

for(i in 1:n){
  X[i,1]~dnorm(mu +alpha[X[i,2]]+ beta[X[i,3]]+ int_effects[X[i,2],X[i,3]], tau )
}

for (j in 1: (alphalength)){
  alpha[j]~ dnorm(0,taualpha)  #random effect, gamma prior on precision
}

for (k in 1: (betalength)){
  beta[k]~ dnorm(0,taubeta) #random effect, gamma prior on precision
}
for(j in 1: alphalength){
  for(k in 1:betalength){
    ind[j,k]~dbern(0.7)
      int_effects[j,k] ~ dnorm(0,tau_int)
      abEff[j,k]<-ind[j,k]*int_effects[j,k]
    }
}


mu~dnorm(0,0.0001)
tau~dgamma(1,0.0001)
taualpha~ dgamma(1,0.0001)
taubeta~ dgamma(1,0.0001)
tau_int~ dgamma(1,0.0001)
}
"
model=jags.model(textConnection(modelstring),data=data,inits=init)
update(model,n.iter=1000)
output=coda.samples(model=model,variable.names=c("mu","tau","alpha", "beta", "ind", "abEff"),n.iter=1000000,thin=300)
summary(output)
plot(output)

#-----------6. For Mixed effects + variable selection --------------
modelstring="

model {

for(i in 1:n){
X[i,1]~dnorm(mu +alpha[X[i,2]]+ beta[X[i,3]]+ abEff[X[i,2],X[i,3]], tau )
}

for (j in 1: (alphalength-1)){
alpha[j]~ dnorm(0,0.0001)  #fixed effect, from a fixed prior distribution
}

for (k in 1: (betalength-1)){
beta[k]~ dnorm(0,0.0001) #fixed effect, from a fixed prior distribution
}

for(j in 1: alphalength){
  for(k in 1: betalength){
    ind[j,k]~dbern(0.7)
      int_effects[j,k] ~ dnorm(0,tau_int) #random effect, gamma prior on precision
      abEff[j,k]<-ind[j,k]*int_effects[j,k]
    }
}

alpha[alphalength] <- -sum(alpha[1:(alphalength-1)]) 
beta[betalength] <- -sum(beta[1:(betalength-1)])
mu~dnorm(0,0.0001)
tau~dgamma(1,0.0001)
tau_int~ dgamma(1,0.0001)
}
"
model=jags.model(textConnection(modelstring),data=data,inits=init)
update(model,n.iter=1000)
output=coda.samples(model=model,variable.names=c("mu","tau","alpha", "beta", "ind", "abEff"),n.iter=1000000,thin=300)
summary(output)
plot(output)
autocorr.plot(output)

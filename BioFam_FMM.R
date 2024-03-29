# work with biofam data b/c it is public
require(flexmix)
require(TraMineR)
require(TraMineRextras)
require(dplyr)

data(biofam)
#convert to time to LEFT+MARR+CHILD as SHA example -

# NOTE: code taken from TraMineRextras manual page example

bf.longlab <- c("Parent","Left","Married","Left+Marr","Child", "Left+Child", "Left+Marr+Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
bf.seq <- seqdef(biofam[,10:25], states=bf.shortlab)


data(biofam)
#convert to movement until terminal state: LEFT+MARR+CHILD as HMM example
#state 6 is terminal state in this analysis (LMC event)

biofam$last <- apply(biofam[,10:25]==6,1,function(x) { pos <- which(x)[1]; if(is.na(pos)) pos <- Inf; return(pos)})

#create two channels and truncate the sequence based on 'last'

# Data pre-processing steps:
## 1. wide to long format, with 't' representing time/age
## 2. select columns
## 3. sort by id, time/age
## 4. Remove records after terminal event (This pre-processing ensures equi-sized channels
## 5. Build Event Channel - 000...001, where there is a 1 PLACED in the penultimate time point (last-1) IF the event occurs in the next period. This is to align the information in the Token Channel *preceeding* the event with the event itself.
## 6. remove the *actual* time point at which the event occurs ('last')
## 7/8. clean the data a bit

bioLong <- reshape(biofam,varying=10:25,sep="a",v.names="tok",timevar="t",direction="long") %>% select(id,t,tok,sex,last) %>%
    arrange(id,t) %>%
    filter(t<=last) %>%
    mutate(evt=ifelse(t<last-1,0,1)) %>%
    filter(t<last) %>%
    mutate(sex=ifelse(sex=="woman",1,0)) %>%
    select(-last)

#remove state 6 from this version by shifting state 7 'down'
bioLong <- bioLong %>% mutate(tok=ifelse(tok==7,6,tok))

#FlexMix can take some time to run.

## We saved two files for you, to save a little time
## Set the flag `loadRuns` to TRUE to load
loadRuns <- TRUE
if (loadRuns) load("FMMruns.Rdata")

#If you have time, set up a range of K to try many models
#Reminder: the LMC state is the endpoint & the evt channel is set to 0 until the time right before the event, where it is set to 1.

Krange <- 2:4  #SLOW!!!
nRep <- 5 #number of random starts (slow, but important)

#this type of coding was used in the paper 
bioLong$ageGE18 <- bioLong$t>4
bioLong$ageGE22 <- bioLong$t>8
bioLong$ageGE26 <- bioLong$t>12

ages <- sort(unique(bioLong$t))+13 #age 14 is youngest

## Call to flexmix:
### setup for grouping: formula=.~0|id
### Channel 1: FLXMRmultinom(tok~t): a linear drift model for class-specific token multinomial probablities
### Channel 2: FLXMRglmfix(cbind(evt,1-evt)~sex,  ##the main, class-specific logit model
### fixed=~ageGE18+ageGE22+ageGE26,  ## the homogenous (across classes) portion of the logit model
### family="binomial") ## indicates binary logistic model in this modeling framework

set.seed(101010)
if (!exists("fmm1")) { # check whether data is loaded
  fmm1 <- stepFlexmix(.~0|id,k=Krange,nrep=nRep,model=list(FLXMRmultinom(tok~t),FLXMRglmfix(cbind(evt,1-evt)~sex,fixed=~ageGE18+ageGE22+ageGE26, family="binomial")), data=bioLong,control=list(verbose=1,minprior = 0.025))
}
#lists models with different number of components k, BIC, etc.
fmm1
#find location of best model
bestBIC <- as.numeric(which.min(BIC(fmm1)))
## mixing parameters are displayed with this call:
summary(fmm1@models[[bestBIC]])

# Observation class membership:
# 1. identify record-specific probability of class membership 
# 2. take the MAP
# 3. reduce to id-specific (corresponds to wide format)
# 4. summarize
probClass.1 <- fmm1@models[[bestBIC]]@posterior$scaled
MAP.1 <- apply(probClass.1,1,which.max)
class.id.lvl.1 <- tapply(MAP.1,bioLong$id,"[[",1) #take first; they are all the same
class.id.lvl.1 <- tapply(MAP.1,bioLong$id,mean) #take mean; they are all the same, and if there's a bug, you might catch it this way
table(class.id.lvl.1)

##results are given in a list, one element for each component of the mixture:
rsltList <- fmm1@models[[bestBIC]]@components

# We break this down
## DRIFT MODELS:
## These multinomial logit parameters describe the relative odds of each token (with 1 as the reference) changing linearly over time. Combined and transformed, they produce non-linear probability curves in a competing risks sense, as per the JRSS-A discussion (and figures) of FMM in Scott et al. (2020)
rel.odds <- lapply(fmm1@models[[bestBIC]]@components,"[[",1)
rel.odds
#plot associated with them:
x.pts <- sort(unique(bioLong$t)) #code uses t for prediction
x.pred.mat <- cbind(1,x.pts) 
par(mfrow=c(2,2)) # for 4 classes
for (i in seq_along(rel.odds)) {
  plot(ages,rep(0,length(ages)),type='n',ylab='Prob. (competing)',
        ylim=c(0,1),axes=F,main=paste0("Class ",i))
  axis(1,at=ages) #looks better to use ages
  axis(2,at=seq(0,1,length=11))
  prob.curves <- rel.odds[[i]]@predict(x.pred.mat)
  matlines(ages,prob.curves,type='l',lwd=3,col=2:8)
}
## HAZARD MODEL:
#For each component (cluster), we fit a model for logit(P(exit)) given covariates and latent class
##Example 1: the model specification constrains all AGE effects to be homogenous across classes (the `fixed` parameter in the flexmix call), while the sex effect varies (the main or default formula of the second channel)
##NOTE: period effects are additive due to the use of ">=" (GE)
rel.hzd <- lapply(fmm1@models[[bestBIC]]@components,"[[",2)
rel.hzd
#build design matrix (from data)
desMat <- unique(bioLong[,c("t","ageGE18","ageGE22","ageGE26","sex")]) %>% 
  arrange(sex,t) %>% 
  mutate(ageGE18TRUE=ifelse(ageGE18,1,0),ageGE22TRUE=ifelse(ageGE22,1,0),
         ageGE26TRUE=ifelse(ageGE26,1,0), `(Intercept)`=1) %>%
  select(-t,-ageGE18,-ageGE22,-ageGE26) %>%
  relocate(sex,.after = last_col()) %>% 
  as.matrix()

par(mfrow=c(2,2)) # for 4 classes
for (i in seq_along(rel.hzd)) {
  plot(ages,rep(0,length(ages)),type='n',ylab='Hazard (Prob.)',
       ylim=c(0,0.125),axes=F,main=paste0("Class ",i))
  axis(1,at=ages) #looks better to use ages
  axis(2)
  hzd.preds <- rel.hzd[[i]]@predict(desMat)
  hzd.curves <- do.call(cbind,split(hzd.preds,desMat[,"sex"])) #split by sex - make 2 cols
  matlines(ages,hzd.curves,type='l',lwd=3,col=2:3)
  legend("topleft",col=2:3,lty=1:2,lwd=2,c('sex=0','sex=1'))
}


## EXAMPLE 2

## Call to flexmix:
### the difference is in:
### Channel 2 (HAZARD MODEL): FLXMRglmfix(cbind(evt,1-evt)~-1+ageGE18+ageGE22+ageGE26,  ##the main, class-specific logit model - age effects but no constant
### fixed=~sex+log(t),  ## the homogenous (across classes) portion of the logit model - single sex effect, single underlying monotonic effect of log(time)

set.seed(101010)
if (!exists("fmm2")) {
  fmm2 <- stepFlexmix(.~0|id,k=Krange,nrep=nRep,model=list(FLXMRmultinom(tok~t),FLXMRglmfix(cbind(evt,1-evt)~-1+ageGE18+ageGE22+ageGE26,fixed=~sex+log(t), family="binomial")), data=bioLong,control=list(verbose=1,minprior = 0.025))
}

#lists models with different number of components k, BIC, etc.
fmm2
#find location of best model
bestBIC <- as.numeric(which.min(BIC(fmm2)))
## mixing parameters are displayed with this call:
summary(fmm2@models[[bestBIC]])

# Observation class membership:
# 1. identify record-specific probability of class membership 
# 2. take the MAP
# 3. reduce to id-specific (corresponds to wide format)
# 4. summarize
probClass.2 <- fmm2@models[[bestBIC]]@posterior$scaled
MAP.2 <- apply(probClass.2,1,which.max)
class.id.lvl.2 <- tapply(MAP.2,bioLong$id,"[[",1) #take first; they are all the same
class.id.lvl.2 <- tapply(MAP.2,bioLong$id,mean) #take mean; they are all the same, and if there's a bug, you might catch it this way
table(class.id.lvl.2)

##results are given in a list, one element for each component of the mixture:
rsltList <- fmm2@models[[bestBIC]]@components

# We break this down
## DRIFT MODELS:
## These multinomial logit parameters describe the relative odds of each token (with 1 as the reference) changing linearly over time. Combined and transformed, they produce non-linear probability curves in a competing risks sense, as per the JRSS-A discussion (and figures) of FMM in Scott et al. (2020)
rel.odds <- lapply(fmm2@models[[bestBIC]]@components,"[[",1)
rel.odds
#plot associated with them:
x.pts <- sort(unique(bioLong$t)) #code uses t for prediction
x.pred.mat <- cbind(1,x.pts) 
par(mfrow=c(2,2)) # for 4 classes
for (i in seq_along(rel.odds)) {
  plot(ages,rep(0,length(ages)),type='n',ylab='Prob. (competing)',
       ylim=c(0,1),axes=F,main=paste0("Class ",i))
  axis(1,at=ages) #looks better to use ages
  axis(2,at=seq(0,1,length=11))
  prob.curves <- rel.odds[[i]]@predict(x.pred.mat)
  matlines(ages,prob.curves,type='l',lwd=3,col=2:8)
}


## HAZARD MODEL:
##Example 2: the model specification constrains all sex and log(time) effects to be homogenous across classes (the `fixed` parameter in the flexmix call), while age indicators vary, allowing a general time effect punctuated by larger, class-specific period effects.
##NOTE: period effects are additive due to the use of ">=" (GE)
##CAUTION: the constant is supressed in this formulation, creating a slightly different set of period effects
rel.hzd <- lapply(fmm2@models[[bestBIC]]@components,"[[",2)
rel.hzd
#build design matrix (from data)
desMat <- unique(bioLong[,c("t","ageGE18","ageGE22","ageGE26","sex")]) %>% 
  arrange(sex,t) %>% 
  mutate(ageGE18FALSE=ifelse(ageGE18,0,1),ageGE18TRUE=ifelse(ageGE18,1,0),ageGE22TRUE=ifelse(ageGE22,1,0),
         ageGE26TRUE=ifelse(ageGE26,1,0), `log(t)`=log(t)) %>%
  select(-t,-ageGE18,-ageGE22,-ageGE26) %>%
  relocate(`log(t)`,.after = sex) %>% 
  as.matrix()

par(mfrow=c(2,2)) # for 4 classes
for (i in seq_along(rel.hzd)) {
  plot(ages,rep(0,length(ages)),type='n',ylab='Hazard (Prob.)',
       ylim=c(0,0.125),axes=F,main=paste0("Class ",i))
  axis(1,at=ages) #looks better to use ages
  axis(2)
  hzd.preds <- rel.hzd[[i]]@predict(desMat)
  hzd.curves <- do.call(cbind,split(hzd.preds,desMat[,"sex"])) #split by sex - make 2 cols
  matlines(ages,hzd.curves,type='l',lwd=3,col=2:3)
  legend("topleft",col=2:3,lty=1:2,lwd=2,c('sex=0','sex=1'))
}

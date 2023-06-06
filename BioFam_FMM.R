# work with biofam data b/c it is public
require(seqHMM)
require(LMest)
require(TraMineR)
require(TraMineRextras)
require(igraph)
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

bioLong <- reshape(biofam,varying=10:25,sep="a",v.names="tok",timevar="t",direction="long") %>% select(id,t,tok,sex,last) %>%
    arrange(id,t) %>%
    filter(t<=last) %>%
    mutate(evt=ifelse(t<last-1,0,1)) %>%
    filter(t<last) %>%
    mutate(sex=ifelse(sex=="woman",1,0)) %>%
    select(-last)

#FlexMix can take some time to run.
#If you have time, set up a range of K to try many models
#Reminder: the LMC state is the endpoint & the evt channel is set to 0 until the time right before the event, where it is set to 1.

Krange <- 2:4  #SLOW!!!
nRep <- 5 #number of random starts (slow, but important)

#this type of coding was used in the paper
bioLong$ageGE18 <- bioLong$t>4
bioLong$ageGE22 <- bioLong$t>8
bioLong$ageGE26 <- bioLong$t>12

set.seed(101010)
fmm1 <- stepFlexmix(.~0|id,k=Krange,nrep=nRep,model=list(FLXMRmultinom(tok~t),FLXMRglmfix(cbind(evt,1-evt)~sex,fixed=~ageGE18+ageGE22+ageGE26, family="binomial")), data=bioLong,control=list(verbose=1,minprior = 0.025))
fmm1
summary(fmm1@models$`4`) #mixing params
fmm1@models$`4`@components  #list of 2 sets of regression coefs: 1. logit model for change in token probabilities over time (linear drift); 2. logit model for logit(P(exit)) given model and latent class (here model is sex-specific, with time indicators, as in a DT-EHA)



set.seed(101010)
fmm2 <- stepFlexmix(.~0|id,k=Krange,nrep=nRep,model=list(FLXMRmultinom(tok~t),FLXMRglmfix(cbind(evt,1-evt)~-1+ageGE18+ageGE22+ageGE26,fixed=~log(t), family="binomial")), data=bioLong,control=list(verbose=1,minprior = 0.025))
fmm2
summary(fmm2@models$`4`) #mixing params
fmm2@models$`4`@components  #list of 2 sets of regression coefs: 1. logit model for change in token probabilities over time (linear drift); 2. logit model for logit(P(exit)) given model and latent class (here model is age-specific)

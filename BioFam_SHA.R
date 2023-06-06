# work with biofam data b/c it is public
require(seqHMM)
require(TraMineR)
require(TraMineRextras)

data(biofam)
#convert to time to LEFT+MARR+CHILD as SHA example -

# NOTE: code taken from TraMineRextras manual page example

bf.longlab <- c("Parent","Left","Married","Left+Marr","Child", "Left+Child", "Left+Marr+Child", "Divorced")
bf.shortlab <- c("P","L","M","LM","C","LC", "LMC", "D")
bf.seq <- seqdef(biofam[,10:25], states=bf.shortlab)

#state 6 is terminal state in this analysis (LMC event)

## The code below aims to find when this event occurred (and whether it occurred).
bf.seq2 <- seqrecode(bf.seq, recodes=list(LMC="LMC"), otherwise = "Other")
dss <- seqdss(bf.seq2)
## Time until LMC spell
time <- ifelse(dss[, 1]=="LMC", 1, seqdur(bf.seq2)[, 1])
## Whether the event (start of LMC spell) started or not
event <- dss[, 1]=="LMC"|dss[, 2]=="LMC"

## The seqsha function will convert the data to person period.
## At each time point, the previous trajectory until that point is stored
sha <- seqsha(bf.seq, time, event, covar=biofam[, c("sex", "birthyr")])
summary(sha)

## Now we build a sequence object for the previous trajectory
previousTraj <- seqdef(sha[, 4:19])
seqdplot(previousTraj)
## Now we cluster the previous trajectories
##Compute distances using only the dss
## Ensure high sensitivity to ordering of the states
diss <- seqdist(seqdss(previousTraj), method="LCS")
##Clustering with pam
library(cluster)
pclust <- pam(diss, diss=TRUE, k=4, cluster.only=TRUE)
#Naming the clusters
sha$pclustname <- factor(paste("Type", pclust))
##Plotting the clusters to make senses of them.
seqdplot(previousTraj, sha$pclustname)


## Now we use a discrete time model include the type of previous trajectory as covariate.
summary(glm(event~time+pclustname+sex, data=sha, family=binomial))

## In our paper, we discretized time, so the model would look more like this:

summary(glm(event~ I(time>4&time<9)+ I(time>8&time<12)+ I(time>12) +pclustname+sex, data=sha, family=binomial))

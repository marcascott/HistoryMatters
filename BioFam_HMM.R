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

biofam$evt <- apply(biofam[,10:25]==6,1,function(x) { pos <- which(x)[1]; if(is.na(pos)) pos<-Inf; return(pos)})

#pad with events -- as needed, so that LMest allows more complex models
#this step does not change the probabilities other than self-loops at the terminal state, which are irrelevant to us.

bioLong <- reshape(biofam,varying=10:25,sep="a",v.names="tok",timevar="t",direction="long") %>% select(id,t,tok,sex,evt) %>%
    arrange(id,t) %>%
    mutate(tok=ifelse(t<=evt,tok,6)) %>%
    select(-evt) %>%
    mutate(sex=ifelse(sex=="woman",1,0))

#LMest can take some time to run.
#If you have time, increase the range of models chosen using the k parameter
#this will allow the hidden states to reflect MULTIPLE output (nominal values) with the LMC state being the endpoint.

maxK <- 5

hmm1 <- lmestSearch(responsesFormula = tok ~ NULL, latentFormula=~NULL|sex,
    index = c("id","t"), output=T, out_se = F,
     data = bioLong, k = 1:maxK, seed = 10101)
hmm1
plot(hmm1$out.single[[which.min(hmm1$Bic)]]) #best BIC model
#examine the transition plots and conditional probs

#this type of coding was used in the paper
bioLong$ageGE18 <- bioLong$t>4
bioLong$ageGE22 <- bioLong$t>8
bioLong$ageGE26 <- bioLong$t>12

hmm2 <- lmestSearch(responsesFormula = tok ~ NULL, latentFormula=~NULL|ageGE18+ageGE22+ageGE26,
    index = c("id","t"), output=T, out_se = F,
     data = bioLong, k = 1:maxK, seed = 10101)
hmm2

plot(hmm2$out.single[[which.min(hmm2$Bic)]]) #best BIC model
#examine the transition plots and conditional probs

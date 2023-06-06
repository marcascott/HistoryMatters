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
# Data pre-processing steps:
## 1. wide to long format, with 't' representing time/age
## 2. select columns
## 3. sort by id, time/age
## 4. Pad with "events" at all points AFTER the event (see above) - Balanced input NEEDED FOR LMest.
## 5/6. clean the data a bit

bioLong <- reshape(biofam,varying=10:25,sep="a",v.names="tok",timevar="t",direction="long") %>% select(id,t,tok,sex,evt) %>%
    arrange(id,t) %>%
    mutate(tok=ifelse(t<=evt,tok,6)) %>%
    select(-evt) %>%
    mutate(sex=ifelse(sex=="woman",1,0))

#LMest can take some time to run.
#If you have time, increase the range of models chosen using the k parameter
#this will allow the hidden states to reflect MULTIPLE output (nominal values) with the LMC state being the endpoint.

maxK <- 5

###Example 1 (sex-specific transitions):

hmm1 <- lmestSearch(responsesFormula = tok ~ NULL, latentFormula=~NULL|sex,
    index = c("id","t"), output=T, out_se = F,
     data = bioLong, k = 1:maxK, seed = 10101)

hmm1
#select best BIC model with which.min(hmm1$Bic) and summarize
## Note: relies on attempting to fit a mixture model for each of 1:maxK so that selection is correct.
summary(hmm1$out.single[[which.min(hmm1$Bic)]])
# results:
## 1. (Be) logit model coefficients for initial transition probs (reference is cluster 1)
## 2. (Ga) multinomial logit model for transition probs, with reference cluster 1.
##     In reconstructing a (conditional) transition matrix from maxK equations, consider that each equation represents a row of the Markov chain
## 3. (Psi) Read this column-wise, and you seewhat tokens/nominal states are output from each hidden state. Look for category 6 -- the terminal state. It will have a probability near 1 for some node (most likely).

#examine the transition plots
##
##
plot(hmm1$out.single[[which.min(hmm1$Bic)]],what="transitions")
## While this differs from the "fancier" plots in the paper, it provides the same information:  transition probabilities, averaged across the levels of any predictors used in generating those transition probabilities (it has to be one number on the plot, but the probs are conditional, so it can vary depending on the covariates)

#examine the conditional probs
plot(hmm1$out.single[[which.min(hmm1$Bic)]],what="CondProb")
##
## This is a nicer version of the Psi matrix, with the addition of the initial probability of being in any hidden state

#examine the marginal probs
## These are based on a forwards-backwards algorithm that finds the most likely path through the hidden states, given the covariates and observed tokens.
## It provides a sense of how each state can play a different role at a different time in the process.
##  ** the terminal state should be increasing over time, at least in these data
plot(hmm1$out.single[[which.min(hmm1$Bic)]],what="marginal")

###Example 2 (age-specific transitions):

#this type of coding was used in the paper
bioLong$ageGE18 <- bioLong$t>4
bioLong$ageGE22 <- bioLong$t>8
bioLong$ageGE26 <- bioLong$t>12

hmm2 <- lmestSearch(responsesFormula = tok ~ NULL, latentFormula=~NULL|ageGE18+ageGE22+ageGE26,
    index = c("id","t"), output=T, out_se = F,
     data = bioLong, k = 1:maxK, seed = 10101)
hmm2

#examine the transition plots, conditional probs, and marginal probs of tokens (see example 1 for detailed description)
plot(hmm2$out.single[[which.min(hmm2$Bic)]], what="transitions")
plot(hmm2$out.single[[which.min(hmm2$Bic)]], what="CondProb")
plot(hmm2$out.single[[which.min(hmm2$Bic)]], what="marginal")

# HistoryMatters
Code Snippets used in the Quality &amp; Quantity article, "History Matters..."

NOTE: The SHP data that we used for this paper is not public. 
We are working on providing code based on public data - an SHP extract called biofam.
Biofam is based on the SHP, so the basic structure is the same. By choosing the 
LEFT+MARR+CHILD state as the terminal outcome, we can reconstruct all of the 
concepts of the paper with publicly available data.

The application, models and estimates won't be precisely the same, 
but they will give you necessary code snippets and suggest the proper
data structures that can be used on your own data.

Three files are provied, and each will run an example from start to finish (some runs
may take some time to complete, depending on how you set the parameters).

BioFam_SHA.R - Sequence History Analysis (SHA) approach [mostly taken from the TraMineRextra example of seqsha]
BioFam_HMM.R - Hidden Markov Models (HMM) as described in the paper
BioFam_FMM.R - Finite Mixture Models (FMM) - or latent class models - for 2-channel processes, 
               one a Discrete Time Event History Analysis (DT-EHA) as described in the paper
FMMRuns.Rdata- Saved runs for the FMM, which can take a longer time to run.

Questions? Contact marc.scott@nyu.edu (corresponding author)

library(gridExtra)

data(Data.Election)
Data.pairs.Election <- Breaking(Data.Election, "top.partial")
n.Election <- dim(Data.Election)[1]; m.Election <- dim(Data.Election)[2]

iterations <- 25

message("Fitting PL GMM model and displaying visualization")
Estimate.PL.GMM.Election.Parameters <- Estimation.PL.GMM(Data.pairs.Election, m.Election)$Parameters
Generate.Pairwise.Probabilities(Data.pairs.Election, Estimate.PL.GMM.Election.Parameters, PL.Pairwise.Prob, "PL GMM")
message("Recommended to export visualization as paper-size on landscape")

message("Sleeping for 3 seconds before next visualization")
Sys.sleep(3)

message("Fitting Normal GMM model and displaying visualization")
Estimate.Normal.GMM.Election.Parameters <- Estimation.Normal.GMM(Data.pairs.Election, m.Election)$Parameters
Generate.Pairwise.Probabilities(Data.pairs.Election, Estimate.Normal.GMM.Election.Parameters, Normal.Pairwise.Prob, "Normal GMM - Fixed Variance")
message("Recommended to export visualization as paper-size on landscape")


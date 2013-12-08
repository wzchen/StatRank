################################################################################
# Pure Rank Metrics
################################################################################

#' Calculates the Kendall Tau correlation between two ranks
#' 
#' @param rank1 two rankings. Order does not matter
#' @param rank2 two rankings. Order does not matter
#' @return The Kendall Tau correlation
#' @export
#' @examples
#' rank1 <- scramble(1:10)
#' rank2 <- scramble(1:10)
#' Evaluation.KendallTau(rank1, rank2)
Evaluation.KendallTau= function(rank1, rank2)
{
  ;cor(rank1, rank2, method="kendall")
}

#' Calculates the location of the True winner in the estimated ranking
#' 
#' @param EstimatedRank estimated ranking
#' @param TrueRank true ranking
#' @return The location of the true best in the estimated rank
#' @export
#' @examples
#' rank1 <- scramble(1:10)
#' rank2 <- scramble(1:10)
#' Evaluation.LocationofWinner(rank1, rank2)
Evaluation.LocationofWinner= function(EstimatedRank, TrueRank)
{
  location.of.true.winner <- which(TrueRank == min(TrueRank))
  EstimatedRank[location.of.true.winner]
}

################################################################################
# Meta-Search Metrics
################################################################################
#From Paper 1-A Flexible Generative model for Preference Aggregation

#' Calculates the Normalized Discounted Cumluative Gain
#' 
#' @param EstimatedRank estimated ranking
#' @param RelevanceLevel score for the document
#' @return The NDCG for this estimation and relevance level
#' @export
#' @examples
#' EstimatedRank <- scramble(1:10)
#' RelevanceLevel <- runif(10)
#' Evaluation.NDCG(EstimatedRank, RelevanceLevel)
Evaluation.NDCG = function(EstimatedRank, RelevanceLevel)
{
  Ordered.RelevanceLevel <- RelevanceLevel[order(-EstimatedRank)]
  Actual.RelevanceLevel <- RelevanceLevel[order(-RelevanceLevel)]
  dcg <- 0
  normalization.constant <- 0
  for(i in 1:length(Ordered.RelevanceLevel)) dcg <- dcg + (2^Ordered.RelevanceLevel[i] - 1) / log2(i+1)
  for(i in 1:length(Actual.RelevanceLevel)) normalization.constant <- normalization.constant + (2^Actual.RelevanceLevel[i] - 1) / log2(i+1)
  dcg / normalization.constant
}

#' Calculates the Average Precision
#' 
#' @param EstimatedRank estimated ranking
#' @param RelevanceLevel score for the document
#' @return The AP for this estimation and relevance level
#' @export
#' @examples
#' EstimatedRank <- scramble(1:10)
#' RelevanceLevel <- runif(10)
#' Evaluation.AveragePrecision(EstimatedRank, RelevanceLevel)
Evaluation.AveragePrecision = function(EstimatedRank, RelevanceLevel)
{
  Ordered.RelevanceLevel <- RelevanceLevel[order(-EstimatedRank)]
  numerator <- 0
  for(k in 1:length(Ordered.RelevanceLevel)) numerator <- numerator + Evaluation.Precision.at.k(EstimatedRank, RelevanceLevel, k) * Ordered.RelevanceLevel[k]
  numerator / sum(Ordered.RelevanceLevel)  
}

#' Calculates the Average Precision at k
#' 
#' @param EstimatedRank estimated ranking
#' @param RelevanceLevel score for the document
#' @param k positive that we want to run this algorithm for
#' @return The AP at k for this estimation and relevance level
#' @export
#' @examples
#' EstimatedRank <- scramble(1:10)
#' RelevanceLevel <- runif(10)
#' Evaluation.Precision.at.k(EstimatedRank, RelevanceLevel, 5)
Evaluation.Precision.at.k <- function(EstimatedRank, RelevanceLevel, k)
{
  Ordered.RelevanceLevel <- RelevanceLevel[order(-EstimatedRank)]
  mean(Ordered.RelevanceLevel[1:k])
}

################################################################################
# Predictive Modeling Metrics
################################################################################


#' Calculates KL divergence between empirical pairwise preferences
#' and modeled pairwise preferences
#' 
#' @param Data.pairs data broken up into pairs using Breaking function
#' @param m number of alternatives
#' @param Estimate estimation object from an Estimate function
#' @param get.pairwise.prob Function that given two alternatives from 
#' the the Parameters argument, returns back a model probability that one is larger than the other 
#' @return the KL divergence between modeled and empirical pairwise preferences, 
#' thinking of the probabilities as a probability distribution over the (n choose 2) pairs
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' m <- 5
#' Estimate <- Estimation.PL.GMM(Data.Test.pairs, m)
#' Evaluation.KL(Data.Test.pairs, m, Estimate, PL.Pairwise.Prob)
Evaluation.KL <- function(Data.pairs, m, Estimate, get.pairwise.prob) { 
  # empirical
  C.matrix.empirical <- generateC(Data.pairs, m)
  for(i in 1:(m-1)) for(j in (i+1):m) C.matrix.empirical[i, j] <- C.matrix.empirical[i, j] / (C.matrix.empirical[i, j] + C.matrix.empirical[j, i])
  C.empirical <- turn_matrix_into_table(C.matrix.empirical)
  
  # model
  C.matrix.model <- matrix(0, nrow = m, ncol = m)
  for(i in 1:m) for(j in 1:m) if(i != j) C.matrix.model[i, j] <- get.pairwise.prob(Estimate$Parameters[[i]], Estimate$Parameters[[j]])
  C.model <- turn_matrix_into_table(C.matrix.model)
  
  # calculate KL
  P <- C.empirical['prob'] / sum(C.empirical['prob'])
  Q <- C.model['prob'] / sum(C.model['prob'])
  sum(log(P/Q)*P)
}
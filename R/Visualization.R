

# this is to get around issues with ggplot2 code not passing package checks
if(getRversion() >= "2.15.1")  utils::globalVariables(c("x", "y", "group", "column", "prob"))

################################################################################
# Calculating Pairwise Probabilities
################################################################################


#' Pairwise Probability for PL Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
PL.Pairwise.Prob <- function(a, b) a$Mean / (a$Mean + b$Mean)

#' Pairwise Probability for Zemel
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Zemel.Pairwise.Prob <- function(a, b){
  # the means here are actually the scores
  exp(a$Mean - b$Mean) / (exp(a$Mean - b$Mean) + exp(b$Mean - a$Mean))
}

#' Pairwise Probability for Normal Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Normal.Pairwise.Prob <- function(a, b) {
  # Let W = X - Y
  if(is.null(a[["Variance"]])) a$Variance <- 1
  if(is.null(b[["Variance"]])) b$Variance <- 1
  
  mu <- a$Mean - b$Mean
  sigma <- sqrt(a$Variance + b$Variance)
  #P(X - Y > 0)
  1 - pnorm(-mu/sigma)
}

#' Pairwise Probability for Normal Multitype Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Normal.MultiType.Pairwise.Prob <- function(a, b) {
  mean.grid <- expand.grid(a$Mean, b$Mean)
  variance.grid <- expand.grid(a$Variance, b$Variance)
  gamma.grid <- expand.grid(a$Gamma, b$Gamma)
  n <- nrow(mean.grid)
  
  probs <- rep(NA, n)
  for(i in 1:n) {
    probs[i] <- Normal.Pairwise.Prob(list(Mean = mean.grid[i, 1], Variance = variance.grid[i, 1]), list(Mean = mean.grid[i, 2], Variance = variance.grid[i, 2]))
    probs[i] <- probs[i] * prod(gamma.grid[i, ])
  }
  sum(probs)
}

#' Pairwise Probability for PL Multitype Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Expo.MultiType.Pairwise.Prob <- function(a, b) {
  mean.grid <- expand.grid(a$Mean, b$Mean)
  n <- nrow(mean.grid)
  
  probs <- rep(NA, n)
  for(i in 1:n) {
    probs[i] <- PL.Pairwise.Prob(list(Mean = mean.grid[i, 1]), list(Mean = mean.grid[i, 2]))
    probs[i] <- probs[i] * prod(mean.grid[i, ])
  }
  sum(probs)
}


################################################################################
# Visualize MultiType
################################################################################


#' Multitype Random Utility visualizer
#' 
#' @param multitype.output output from a multitype fitter
#' @param names names of alternatives
#' @return none
#' @export
#' @examples
#' library(ggplot2)
#' library(grid)
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' multitype.output <- Estimation.RUM.MultiType.MLE(Data.Tiny, iter = 2, dist = "norm", ratio = .5)
#' names <- 1:5
#' plots <- visualize.MultiType(multitype.output, names)
#' # the list of plots is appropriate for passing into the multiplot function at
#' # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' # as input to the plotlist argument

visualize.MultiType <- function(multitype.output, names) {
  m <- dim(multitype.output$Mean)[2]

  process.MultiType <- function(output.of.multitype, m) {
    params <- list()
    for(i in 1:m) {
      params[[i]] <- list()
      params[[i]]$Mean <- output.of.multitype$Mean[,i] - mean(output.of.multitype$Mean)
      params[[i]]$SD <- output.of.multitype$SD[,i]
      params[[i]]$Gamma <- output.of.multitype$Gamma[1,]
    }
    params
  }
  
  parameters <- process.MultiType(multitype.output, m)
  
  plots <- list()
  linesize <- 1
  for(i in 1:m) {
    means <- parameters[[i]]$Mean
    sds <- parameters[[i]]$SD
    gammas <- parameters[[i]]$Gamma
    get.density <- function(x, weights) sum(dnorm(c(x, x), mean = means, sd = sds) * gammas * weights)
    xs <- seq(-5, 5, by = 0.01)
    df <- data.frame(x = rep(xs, 3))
    df$y <- c(sapply(xs, function(x) get.density(x, c(1, 1))), sapply(xs, function(x) get.density(x, c(1, 0))), sapply(xs, function(x) get.density(x, c(0, 1))))
    df$group <- rep(c(1, gammas), rep(length(xs), 3))
        
    plots[[i]] <- ggplot(df, aes(x = x, y = y, linetype = factor(group), color = factor(group))) + 
      geom_line(size = linesize) + labs(title = names[i], x = NULL, y = NULL) + 
      geom_vline(xintercept = means[1], linetype = 2, size = linesize, color = "blue") + 
      geom_vline(xintercept = means[2], linetype = 6, size = linesize, color = "red") + 
      scale_color_manual(values=c("red", "blue", "black")) +
      scale_linetype_manual(values=c(6, 2, 1)) +
      theme(legend.position = "none")
  }
  plots
}

################################################################################
# Generating the Graphs
################################################################################

#' Helper function for the graphing interface
#'
#' As named, this function takes a vector where each element is a mean, 
#' then returns back a list, with each list item having the mean
#'
#' @param Means a vector of means
#' @return a list, where each element represents an alternative and has a Mean value
#' @export
convert.vector.to.list.of.means <- function(Means) {
  m <- length(Means)
  List <- rep(list(NA), m)
  for(i in 1:m) List[[i]] <- list(Mean = Means[i])
  List
}

#' Creates pairwise matrices to compare inference results with the empirical pairwise probabilities
#' 
#' @param Data.pairs datas broken into pairs
#' @param Parameters The Parameter element of a result from an Estimation function
#' @param get.pairwise.prob function that we use to generate the pairwise probability of beating
#' @param name.of.method names of the alternatives
#' @return none
#' @export
#' @examples
#' library(ggplot2)
#' library(gridExtra)
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Parameters <- Estimation.PL.GMM(Data.Test.pairs, 5)$Parameters
#' PL.Pairwise.Prob <- function(a, b) a$Mean / (a$Mean + b$Mean)
#' Generate.Pairwise.Probabilities(Data.Test.pairs, Parameters, PL.Pairwise.Prob, "PL on Test Data")
Generate.Pairwise.Probabilities <- function(Data.pairs, Parameters, get.pairwise.prob, name.of.method) {
  
  m <- max(Data.pairs[,c(1,2)])
  
  means <- rep(NA, m)
  
  # where it's not multitype
  if(is.null(Parameters[["Gamma"]])) for(i in 1:m) means[i] <- Parameters[[i]]$Mean
  
  # where it is multitype
  else for(i in 1:m) means[i] <- sum(Parameters[[i]]$Mean * Parameters[[i]]$Gamma)
  reordering <- order(-means)
  
  # calculate the empirical differences
  C.matrix.empirical <- generateC(Data.pairs, m)
  C.matrix.empirical.reordered <- C.matrix.empirical[reordering, reordering]
  for(i in 1:(m-1)) for(j in (i+1):m) C.matrix.empirical.reordered[i, j] <- C.matrix.empirical.reordered[i, j] / (C.matrix.empirical.reordered[i, j] + C.matrix.empirical.reordered[j, i])
  for(i in 2:m) for(j in 1:(i-1)) C.matrix.empirical.reordered[i, j] <- 0
  C.empirical <- turn_matrix_into_table(C.matrix.empirical.reordered)
  
  # calculate the model differences
  C.matrix.model <- matrix(0, nrow = m, ncol = m)
  for(i in 1:m) for(j in 1:m) if(i != j) {
    C.matrix.model[i, j] <- get.pairwise.prob(Parameters[[i]], Parameters[[j]])
  }
  C.matrix.model.reordered <- C.matrix.model[reordering, reordering]
  C.model <- turn_matrix_into_table(C.matrix.model.reordered)
  
  
  minprob <- min(C.empirical['prob'], C.model['prob'])
  maxprob <- max(C.empirical['prob'], C.model['prob'])
      
  
  Generate.Pairwise.Matrix.Plot <- function(C.matrix, m, title, minprob, maxprob, reordering) {
    base_size <- 15
    
    ggplot(C.matrix, aes(x = column, y = row)) +
      geom_tile(aes(fill = prob), colour = "white") + geom_text(aes(label = paste0(round(prob * 100), "%")), size = base_size * 0.3, color = "white") +
      scale_fill_gradient(limits = c(minprob, maxprob), breaks = c(0.5, 0.6, 0.7, 0.8, 0.9), low = "darkblue", high = "darkred", labels = function (x) paste0(round(x * 100), "%"), guide = "legend") +
      labs(title = title, x = NULL, y = NULL, fill=NULL, size = base_size) + 
      scale_x_discrete(expand = c(0, 0), breaks = 2:m, labels = reordering[-1]) +
      scale_y_reverse(breaks = 1:(m-1), labels = reordering[-m]) +
      theme_bw(base_size = base_size) + 
      theme(legend.position = c(.15, .35), 
            axis.ticks = element_blank(), 
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_text(color = "gray50"),
            axis.text.x = element_text(color = "gray50"),
            legend.text = element_text(color = "gray50")
      ) +
      coord_fixed() + 
      guides(fill = guide_legend(keywidth = 2.5, keyheight = 2.5))
  }
  
  p1 <- Generate.Pairwise.Matrix.Plot(C.empirical, m, "Empirical", minprob, maxprob, reordering)
  p2 <- Generate.Pairwise.Matrix.Plot(C.model, m, "Model", minprob, maxprob, reordering)
  
  P <- C.empirical['prob'] / sum(C.empirical['prob'])
  Q <- C.model['prob'] / sum(C.model['prob'])
  
  #bias <- mean(C.model[,3]) - mean(C.empirical[,3])
  #l2.distance <- round(sqrt(sum((P - Q)^2)), 5)
  kl.divergence <- round(sum(log(P/Q)*P), 10)
  
  name.of.data <- tail(strsplit(deparse(substitute(Data.pairs)), "\\.")[[1]], n=1)
  
  
  grid.arrange(p1, p2, nrow = 1, main = textGrob(paste0("\n", name.of.data, "\nMethod: ", name.of.method, "\nKL divergence: ", kl.divergence, "\nEstimated Ranking: ", toString(reordering)), gp=gpar(cex=2)))
  
}

################################################################################
## General
################################################################################

#' Scramble a vector
#' 
#' This function takes a vector and returns it in a random order
#' 
#' @param x a vector
#' @return a vector, now in random order
#' @export
#' @examples
#' scramble(1:10)
scramble <- function(x) sample(x, length(x))

#' Converts scores to a ranking
#' 
#' takes in vector of scores (with the largest score being the one most preferred)
#' and returns back a vector of WINNER, SECOND PLACE, ... LAST PLACE
#' 
#' @param scores the scores (e.g. means) of a set of alternatives
#' @return an ordering of the index of the winner, second place, etc.
#' @export
#' @examples
#' scores <- Generate.RUM.Parameters(10, "exponential")$Mean
#' scores.to.order(scores)
scores.to.order <- function(scores) (1:length(scores))[order(-scores)]


################################################################################
## Breaking
################################################################################

#' Converts a pairwise count matrix into a probability matrix
#' 
#' @param C original pairwise count matrix
#' @return a pairwise probability matrix
#' @export
#' @examples
#' C= matrix(c(2,4,3,5),2,2)
#' normalizeC(C)
normalizeC <- function(C){
  m <- dim(C)[1]
  normalized.C <- matrix(0, m, m)
  for(i in 1:m)
    for(j in 1:m)
      if(i != j) 
        normalized.C[i, j] <- C[i, j] / (C[i, j] + C[j, i])
  normalized.C
}


#' Generate a matrix of pairwise wins
#' 
#' This function takes in data that has been broken up into pair format.
#' The user is given a matrix C, where element C[i, j] represents exactly
#' how many times alternative i has beaten alternative j
#' 
#' @param Data.pairs the data broken up into pairs
#' @param m the total number of alternatives
#' @param Zemel whether or not this function should use the rank differences 
#' instead of indicators for when i beats j
#' @param prior the initial "fake data" that you want to include in C. A prior 
#' of 1 would mean that you initially "observe" that all alternatives beat all
#' other alternatives exactly once.
#' @return a Count matrix of how many times alternative i has beat alternative j
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' generateC(Data.Test.pairs, 5, Zemel = FALSE, prior = 0)
generateC <- function(Data.pairs, m, Zemel = FALSE, prior = 0) {
  
  Tee <- nrow(Data.pairs)
  
  # C is the transition matrix, where C[i, j] denotes the number of times that
  C <- matrix(data = prior, nrow = m, ncol = m) - prior * m * diag(m)
  for(l in 1:Tee) {
    # i wins, j loses
    i <- Data.pairs[l, 1]
    j <- Data.pairs[l, 2]
    r <- Data.pairs[l, 3]
    if(Zemel) C[i, j] <- C[i, j] + r
    else C[i, j] <- C[i, j] + 1
  }
  
  C
}

#' Breaks full or partial orderings into pairwise comparisons
#' 
#' Given full or partial orderings, this function will generate pairwise comparison
#' Options
#' 1. full - All available pairwise comparisons. This is used for partial
#' rank data where the ranked objects are a random subset of all objects
#' 2. adjacent - Only adjacent pairwise breakings
#' 3. top - also takes in k, will break within top k
#' and will also generate pairwise comparisons comparing the
#' top k with the rest of the data
#' 4. top.partial - This is used for partial rank data where the ranked 
#' alternatives are preferred over the non-ranked alternatives
#' 
#' @param Data data in either full or partial ranking format
#' @param method - can be full, adjacent, top or top.partial
#' @param k This applies to the top method, choose which top k to focus on
#' @return Pairwise breakings, where the three columns are winner, loser and rank distance (latter used for Zemel)
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
Breaking= function(Data, method, k = NULL)
{
  m <- ncol(Data)
  
  pair.full <- function(rankings) pair.top.k(rankings, length(rankings))
  
  pair.top.k <- function(rankings, k) {
    pair.top.helper <- function(first, rest) {
      pc <- c(); 
      z <- length(rest)
      for(i in 1:z) pc <- rbind(pc, array(as.numeric(c(first, rest[i], i))))
      pc
    }
    if(length(rankings) <= 1 | k <= 0) c()
    else rbind(pair.top.helper(rankings[1], rankings[-1]), pair.top.k(rankings[-1], k - 1))
  }
  
  pair.adj <- function(rankings) {
    if(length(rankings) <= 1) c()
    else rbind(c(rankings[1], rankings[2]), pair.adj(rankings[-1]))
  }
  
  pair.top.partial <- function(rankings, m) {
    # this is used in the case when we have missing ranks that we can
    # fill in at the end of the ranking. We can assume here that all 
    # ranked items have higher preferences than non-ranked items
    # (e.g. election data)
    
    # the number of alternatives that are not missing
    k <- length(rankings)
    
    # these are the missing rankings
    missing <- Filter(function(x) !(x %in% rankings), 1:m)
    
    # if there is more than one item missing, scramble the rest and place them in the ranking
    if(m - k > 1) missing <- scramble(missing)
    
    # now just apply the top k breaking, with the missing elements
    # inserted at the end
    pair.top.k(c(rankings, missing), k)
  }
  
  break.into <- function(Data, breakfunction, ...) {
    n <- nrow(Data)
    # applying a Filter(identity..., ) to each row removes all of the missing data
    # this is used in the case that only a partial ordering is provided
    do.call(rbind, lapply(1:n, function(row) breakfunction(Filter(identity, Data[row, ]), ...)))
  }
  
  if(method == "full") Data.pairs <- break.into(Data, pair.full)
  if(method == "adjacent") Data.pairs <- break.into(Data, pair.adj)
  if(method == "top") Data.pairs <- break.into(Data, pair.top.k, k = k)
  if(method == "top.partial") Data.pairs <- break.into(Data, pair.top.partial, m = m)
  
  Data.pairs
}

##################################################Sampling#############################################

# Conditional truncated sampler for normal
# 
# @param X sample values for alternative, m dimensional
# @param pi ranks pi[1] is the best pi[m] is the worst sampling X[i]
# @param theta m dimentional parameter values
samplen=function(dist,parameter,i,X,pi)
{
  m=length(X)
  if(!sum(pi==0))
  {
    rank=which(pi==i)

    if(rank>1 & rank<m )
    {
      lower=X[as.numeric(pi[rank+1])]
      upper=X[as.numeric(pi[rank-1])]
    }
    
    if(rank==m)
    {
      lower=-Inf
      upper=X[as.numeric(pi[rank-1])]
    }
    
    if(rank==1)
    {
      lower=X[as.numeric(pi[rank+1])]
      upper=Inf
    }
    
  }
  if(sum(pi==0))
  {
    listemp=as.numeric(pi[which(pi>0)])
    listempz=setdiff(1:m,listemp)
    
    if(!sum(pi==i))
    {
      if(pi[m]==0)
      {
        lower=-Inf
        #upper=Inf
        upper=min(X[listemp])
      }
      if(pi[1]==0)
      {
        #lower=-Inf
        upper=Inf
        lower=max(X[listemp])
      }
    }
    
    if(sum(pi==i))
    {
      rank=which(pi==i)
      
      if(rank>1 & rank<m )
      {
        if(pi[rank+1]>0 & pi[rank-1]>0)
        {  
          lower=X[as.numeric(pi[rank+1])]
          upper=X[as.numeric(pi[rank-1])]
        }
        if(pi[rank+1]==0 & pi[rank-1]>0)
        {  
          #lower=-Inf
          lower=max(X[listempz])
          upper=X[as.numeric(pi[rank-1])]
        }
        if(pi[rank+1]>0 & pi[rank-1]==0)
        {  
          #print("shit")
          lower=X[as.numeric(pi[rank+1])]
          upper=min(X[listempz])
          #upper=min(x[pi[1:(rank)]])
        }
        
      }
      
      if(rank==m)
      {
        if(pi[rank-1]>0)
        {
          lower=-Inf
          upper=X[as.numeric(pi[rank-1])]
        }
        if(pi[rank-1]==0)
        {
          lower=-Inf
          #lower=-Inf
          upper=min(X[listempz])
        }
        
      }
      
      if(rank==1)
      {
        if(pi[rank+1]>0)
        {
          lower=X[as.numeric(pi[rank+1])]
          upper=Inf
        }
        if(pi[rank+1]==0)
        {
          lower=max(X[listempz])
          #lower=-Inf
          upper=Inf
        }
      }
    }
    
  }
  ##############
  
  
  if(dist=="norm")
  {
    templirong=X[i]
    X[i]=rtrunc(1, spec=dist, a = lower, b = upper,mean=parameter$Mu[i],sd=parameter$Var[i]^.5)
    if(X[i]==Inf | X[i]==-Inf){X[i]=templirong}
    
    if((X[i]<lower) | (X[i]>upper))
    {
      # print("-")
      X[i]=templirong
    }
  }
  
  if(dist=="exp")
  { 
    #gcon=0.57721566
    templirong=X[i]
    #X[i]=parameter$Mu[i]-gcon-log(rtrunc(1, spec=dist, a = exp(-upper+parameter$Mu[i]-gcon), b =exp(-lower+parameter$Mu[i]-gcon) ,rate=1))
    
    X[i]=rtrunc(1, spec=dist, a = lower, b = upper,rate=parameter$Mu[i]^-1)
    #print(c(lower,upper))
    
    if(X[i]==Inf | X[i]==-Inf | is.na(X[i])){X[i]=templirong}
    
    #print(X[i])
    
    if((X[i]<lower) | (X[i]>upper))
    {
      #print("-")
      X[i]=templirong
    }
    
  }
  
  ;X 
}

# Gibbs sampler
#
# Samples S samples Using a Gibbs sampler from the joint distribution of normal or exponential with a constraint on the
# order of samples given by pi as the rank. 
#   
# @param S number of samples
# @param X sample values for alternative, m dimensional
# @param pi ranks pi[1] is the best pi[m] is the worst sampling X[i]
# @param dist distribution of utilities
# @param parameter m dimensional parameter values
GibbsSampler=function(S, X, pi, dist, parameter)
{
  m=length(X)
  out=matrix(0,S,m)
  out[1,]=X  
  
  #gcon=mean(-log(rexp(100000,1)))
  
  for(s in 2:S)
  {
    i=sample.int(m, size = 1, replace = FALSE, prob = NULL)
    out[s,]=samplen(dist,parameter,i,out[s-1,],pi)
  }
  
  ;list(M1=colMeans(out[(round(S/10)):S,]),M2=colMeans(out[(round(S/10)):S,]^2))
}

#################################Helpers for GMM for Normal
#

f= function(Mu)
{
  m=length(Mu)
  A=matrix(0,m,m)
  for(i in 1:m)
  {
    for(j in 1:(m))
    {
      A[i,j]=pnorm(Mu[i]-Mu[j],0,sqrt(2))
    }
  }
  diag(A)=1-colSums(A);
  ;A
}


delta= function(Mu)
{
  m=length(Mu)
  A=matrix(0,m,m)
  for(i in 1:m)
  {
    for(j in 1:(m))
    {
      A[i,j]=Mu[i]-Mu[j]
    }
  }
  ;A
}

#full break
fullbreak= function(data)
{
  n=dim(data)[1]
  m=dim(data)[2]
  A=matrix(0,m,m)
  for(i in 1:n)
  {
    for(j in 1:(m-1))
    {
      for(k in (j+1):m)
      {
        A[data[i,j],data[i,k]]=A[data[i,j],data[i,k]]+1
      }
    }
  }
  A=A/n;
  diag(A)=.5-colSums(A);
  ;A
}

#top-c break
topCbreak= function(data,c)
{
  n=dim(data)[1]
  m=dim(data)[2]
  A=matrix(0,m,m)
  for(i in 1:n)
  {
    for(j in 1:(c-1))
    {
      for(k in (j+1):c)
      {
        A[data[i,j],data[i,k]]=A[data[i,j],data[i,k]]+1;
      }
    }
  }
  
  An=matrix(0,m,m)
  for(i in 1:m)
  {
    for(j in 1:m)
    {
      if((A[i,j]+A[j,i])!=0)
      {
        An[i,j]=A[i,j]/(A[i,j]+A[j,i])
      }
    }
  }
  diag(An)=.5-colSums(An);
  ;An
}


Analyze=function(brokendata, m)
{
  e=eigen(brokendata)
  ;sort(e$vectors[,m],index=TRUE,decreasing=TRUE)$x
}

###

#' Converts a matrix into a table
#' 
#' takes a matrix and returns a data frame with the columns being row, column, entry
#' 
#' @param A matrix to be converted
#' @return a table with the entries being the row, column, and matrix entry
#' @export
turn_matrix_into_table <- function(A) {
  m <- dim(A)[1]
  transcribed.table <- data.frame()
  for(i in 1:(m-1)) for(j in (i+1):m) transcribed.table <- rbind(transcribed.table, c(i, j, A[i, j]))
  names(transcribed.table) <- c("row", "column", "prob")
  transcribed.table
}


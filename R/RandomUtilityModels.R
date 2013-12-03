# @n stands for number of agents
# @m stands for number of alternatives
# @distribution can be Normal or Exponential

# The Code Contains 4 different parts
# 0) Helper functions
# 1) Data Generation
# 2) Model Likelihood
# 3) Estimation Methods
# 4) Inference

################################################################################
# Data Generation
################################################################################

#' Parameter Generation for a RUM model
#' 
#' Exponential models mean parameters are drawn from a uniform distribution
#' Normal models, mean and standard devaition parameters are drawn from a standard unifrom
#' 
#' @param m number of sets of parameters to be drawn
#' @param distribution either 'normal' or 'exponential'
#' @return a list of RUM parameters
#' @export
#' @examples
#' Generate.RUM.Parameters(10, "normal")
#' Generate.RUM.Parameters(10, "exponential")

Generate.RUM.Parameters= function(m, distribution)
{
  if(distribution=='normal'){parameter=list(Mean=runif(m,0,1),SD=runif(m,0,1))}
  else if(distribution=='exponential'){
    unscaled = runif(m,0,1)
    parameter=list(Mean=unscaled/sum(unscaled))
  }
  else stop(paste0("Distribution name \"", distribution, "\" not recognized"))
  ;parameter
}

#' Generate observation of ranks given parameters
#' 
#' Given a list of parameters (generated via the Generate RUM Parameters function),
#' generate random utilities from these models and then return their ranks
#' 
#' @param parameter list of parameters from the Generate.RUM.Parameters function
#' @param m number of alternatives
#' @param n number of agents
#' @param distribution can be either 'normal' or 'exponential'
#' @return a matrix of observed rankings
#' @export
#' @examples
#' parameter= Generate.RUM.Parameters(10, "normal")
#' Generate.RUM.Ranks(parameter,m=10,n=5,"normal")
#' parameter= Generate.RUM.Parameters(10, "exponential")
#' Generate.RUM.Ranks(parameter,m=10,n=5,"exponential")
Generate.RUM.Ranks=function(parameter,m,n,distribution)
{
  Rank=matrix(,n,m)
  for(i in 1:n)
    {
  
      if(distribution=='normal'){x=rnorm(m,parameter$Mean,parameter$SD)}
      else if(distribution=='exponential'){x=rexp(m,rate=parameter$Mean^-1)}
      else stop(paste0("Distribution name \"", distribution, "\" not recognized"))
      Rank[i,]=sort(x,index=TRUE,decreasing=TRUE)$ix
    }
  
  ;Rank
}

################################################################################
# Model Likelihood
################################################################################

#' A faster Likelihood for Plackett-Luce Model
#'
#' @param Data ranking data
#' @param parameter Mean of Exponential Distribution
#' @return log likelihood
#' @export
#' @examples
#' data(Data.Test)
#' parameter = Generate.RUM.Parameters(5, "exponential")
#' Likelihood.PL(Data.Test, parameter)

Likelihood.PL=function(Data, parameter)
{
  gam=parameter$Mean^-1
  n=dim(Data)[1]
  m=dim(Data)[2]
  rank=Data
  ll=0
  
  for(i in 1:n)
  {
    mj=sum(rank[i,]>0)
    temp=gam[rank[i,][rank[i,] > 0]]
    temp=temp[mj:1]
    ll=ll+sum(log(temp))-sum(log(cumsum(temp)))
  }
  
  ;ll
}

# Helper for Gumbel pdf
pdfgumbel=function(x,mu)
{
  ;exp(-x+mu)*exp(-exp(-x+mu))
}

#' Likelihood for general Random Utility Models with different  noise distributions
#'
#' @param Data ranking data
#' @param parameter Mean of Exponential Distribution
#' @param dist exp or norm
#' @param range range
#' @param res res
#' @return log likelihood
#' @export
#' @examples
#' data(Data.Test)
#' parameter = Generate.RUM.Parameters(5, "normal")
#' Likelihood.RUM(Data.Test,parameter, "norm")
#' parameter = Generate.RUM.Parameters(5, "exponential")
#' Likelihood.RUM(Data.Test,parameter, "exp")

# if(dist=="norm")
# {
#  llfull=0
#  for(i in 1:n)
#  {
#    lli=0
#    for(k in 1:K)
#    {
#     ranget=max(abs(parameter$Mu[k,]))+3*max((SD[k,]))
#     ll=Likelihood.RUM(Data[i,],parameter=list(Mu=Mean[k,],SD=SD[k,]),dist,range=ranget,res=ranget*.0001)
#     lli=exp(ll)*gamma[k]+lli
#    }
#    llfull=llfull+log(lli)
#  }
# }
# if(dist=="exp")
# {
#  llfull=0
#  for(i in 1:n)
#  {
#    lli=0
#    for(k in 1:K)
#    {
#      ranget=max(abs(parameter$Mu[k,]))+3*max(abs(parameter$Mu[k,]))^.5
#      ll=Likelihood.RUM(DataL[i,],parameter=list(Mu=Mean[k,]),dist,range=ranget,res=ranget*.0001)
#      lli=exp(ll)*gamma[k]+lli
#    }
#    llfull=llfull+log(lli)
#  }
#  #ll[j]=likelihoodPL(DataL,parameter)
# }

Likelihood.RUM <- function(Data, parameter, dist = "exp",range = 10,res = .1)
{ 
  if(!(dist=="dexp" | dist=="exp" | dist=="norm"))
  {
    stop(paste0("Distribution name \"", dist, "\" not recognized"))
  }
  
  rank=Data
  S=range/res
  x=(-S:S)*res
  
  n=dim(rank)[1]
  m=dim(rank)[2]
  
  if(dist=="norm")
  {
    
    ll=0
    for(i in 1:n)
    {
      if(sum(Data[i,]==0)>0)
      {
        mj=min(which(rank[i,]==0))-1
        CDF=matrix(1,1,length(x))
        
        for(jt in setdiff(1:m,rank[i,1:mj]))
        {
          CDF=dnorm(x,mean=parameter$Mean[jt],sd=parameter$SD[jt])*CDF
        }
        #mt=min(which(Data[i,]==0))-1
        for(j in mj:1)
        {
          PDF=dnorm(x,mean=parameter$Mean[rank[i,j]],sd=parameter$SD[rank[i,j]])*CDF
          CDF=res*cumsum(PDF)
        }
        ll=log(CDF[length(x)])+ll
        
      }
      if(!(sum(Data[i,]==0)>0))
      {
        CDF=matrix(1,1,length(x))
        
        #mt=min(which(Data[i,]==0))-1
        for(j in m:1)
        {
          PDF=dnorm(x,mean=parameter$Mean[rank[i,j]],sd=parameter$SD[rank[i,j]]^.5)*CDF
          CDF=res*cumsum(PDF)
        }
        ll=log(CDF[length(x)])+ll
      }
    }
  }
  
  if(dist=="exp")
  {
    x=(0:S)*res
    ll=0
    
    for(i in 1:n)
    {
      if(sum(Data[i,]==0)>0)
      {
        #print(n)
        #print(min(which(rank[i,]==0)))
        mj=min(which(rank[i,]==0))-1
        #print(mj)
        #if(mj==1)
        #{
        #  ll=ll-log(parameter$Mu[rank[i,1]])+log(sum(parameter$Mu^-1))
        #}
        #else{
          CDF=matrix(1,1,length(x))
          
          for(jt in setdiff(1:m,rank[i,1:mj]))
          {
            CDF=dexp(x,rate=parameter$Mean[jt]^-1)*CDF
          }
          #mt=min(which(Data[i,]==0))-1
          for(j in mj:1)
          {
            PDF=dexp(x,rate=parameter$Mean[rank[i,j]]^-1)*CDF
            CDF=res*cumsum(PDF)
            CDF=CDF[length(CDF)]-CDF
          }
          ll=log(CDF[1])+ll
        #}
        
      }
      
      if(!(sum(Data[i,]==0)>0))
      {
        CDF=matrix(1,1,length(x))
        
        #mt=min(which(Data[i,]==0))-1
        for(j in m:1)
        {
          PDF=dexp(x,rate=parameter$Mean[rank[i,j]]^-1)*CDF
          CDF=res*cumsum(PDF)
          CDF=CDF[length(CDF)]-CDF
        }
        ll=log(CDF[1])+ll
      }
    }
    
    
    #############
    #    for(i in 1:n)
    #    {
    #      CDF=matrix(1,1,length(x))
    #      
    #      #mt=min(which(Data[i,]==0))-1
    #      for(j in m:1)
    #      {
    #        PDF=dexp(x,rate=parameter$Mu[rank[i,j]]^-1)*CDF
    #        CDF=res*cumsum(PDF)
    #      }
    #      ll=log(CDF[length(x)])+ll
    #    }
  }
  
  if(dist=="dexp")
  {
    x=(-S:S)*res
    ll=0
    for(i in 1:n)
    {
      CDF=matrix(1,1,length(x))
      
      #mt=min(which(Data[i,]==0))-1
      for(j in m:1)
      {
        PDF=pdfgumbel(x,mu=parameter$Mean[rank[i,j]])*CDF
        CDF=res*cumsum(PDF)
      }
      ll=log(CDF[length(x)])+ll
    }
  }
  
  
  ;ll
}

################################################################################
## MLE: MCEM( MCEM.Ag() )
################################################################################

#' Performs parameter estimation for the Plackett-Luce model using an Minorize Maximize algorithm
#' 
#' @param Data data in either partial or full rankings
#' @param iter number of MM iterations to run
#' @return list of estimated means (Gamma) and the log likelihoods
#' @export
#' @examples
#' data(Data.Test)
#' Estimation.PL.MLE(Data.Test)
Estimation.PL.MLE=function(Data, iter = 10)
{
  rank=Data
  m=dim(rank)[2]
  n=dim(rank)[1]  
  
  GamaTotal=matrix(0,iter,m)
  LTotal=matrix(0,1,iter)
  
  M=matrix(0,1,n)
  for(i in 1:n)
  {
    if(sum(rank[i,]==0)){M[i]=min(which(rank[i,]==0))-1}
    if(!sum(rank[i,]==0)){M[i]=m}
  }
  
  W=matrix(0,1,m)
  #nominator
  for(t in 1:m)
  {
    W[t]=0
    for(j in 1:n)
    {
      W[t]=(rank[j,M[j]]==t)+W[t]
    }
    W[t]=n-W[t]
  }
  gam=matrix(1,1,m)
  
  for(iteration in 1:iter)
  {
    gamtemp=gam 
    
    for(t in 1:m)
    {     
      #denominator
      denom=0
      for(j in 1:n)
      {
        
        for(i in 1:(M[j]-1))
        {
          delta=0
          delta=sum(rank[j,i:M[j]]==t)
          
          denomt3=0
          for( s in i:M[j])
          {
            denomt3=denomt3+gam[rank[j,s]]
          }
          denom=delta/denomt3+denom
          
        }
      }
      
      gamtemp[t]=W[t]/denom    
    }
    
    gam=gamtemp
    GamaTotal[iteration,]=gam/sum(gam)
    ll=0
    for(i in 1:n)
    {
      mj=sum(rank[i,]>0)
      temp=gam[rank[i,][rank[i,] > 0]]
      temp=temp[mj:1]
      ll=ll+sum(log(temp))-sum(log(cumsum(temp)))
      
      #temp=gam[rank[i,]]
      #temp=temp[m:1]
      #ll=ll+sum(log(gam))-sum(log(cumsum(temp)))
    }
    LTotal[iteration]=ll
    print(paste0("Finished ", iteration, "/", iter))
  } 
  
  params <- convert.vector.to.list.of.means(GamaTotal[iter,])
  
  ;list(Mean = GamaTotal[iter,], LL = LTotal, Parameters = params)
  
}

################################################################################
## MCEM( MCEM.Ag() )
################################################################################

#' Performs parameter estimation for a Random Utility Model with different noise distributions
#' 
#' This function supports RUMs 
#' 1) Normal
#' 2) Normal with fixed variance (fixed at 1)
#' 3) Exponential
#' 
#' @param Data data in either partial or full rankings
#' @param iter number of EM iterations to run
#' @param dist underlying distribution. Can be "norm", "norm.fixedvariance", "exp"
#' @return parameters of the latent RUM distributions
#' @export
#' @examples
#' data(Data.Test)
#' Estimation.RUM.MLE(Data.Test, iter = 2, dist="norm")
Estimation.RUM.MLE = function(Data, iter = 10, dist)
{
  if(!(dist=="dexp" | dist=="exp" | dist=="norm" | dist == "norm.fixedvariance" ))
  {
    stop(paste0("Distribution name \"", dist, "\" not recognized"))
  }
  # calculating the dimensions
  dims <- dim(Data)
  n <- dims[1]
  m <- dims[2]
  
  #initialization of mean and variance
  Delta = exp(rnorm(m))
  Variance = exp(rnorm(m))
  
  if(dist=="exp")
  {
    DataL=Data
    for(i in 1:dim(Data)[1])
    {
      DataL[i,]=Data[i,m:1]
    }
    temp1=Data
    Data=DataL
    DataL=temp1
  }
  
  T0=proc.time()
  
  ##############
  #if(dist=="norm" | dist=="lnorm"){parameter=list(Mu=Delta,Var=Variance)}
  #if(dist=="norm"){parameter=list(Mu=Delta,Var=matrix(.1,1,m))}
  if(dist=="norm" | dist == "norm.fixedvariance"){parameter=list(Mu=Delta,Var=Variance)}
  if(dist=="exp"){parameter=list(Mu=Delta)}
  
  ##############
  MM=matrix(0,iter,m)
  VV=matrix(0,iter,m)
  ll=matrix(0,1,iter)
  
  for(j in 1:iter)
  {
    S=1000+300*j
    
    sent=sprintf("%d Percent is done",ceiling((j*2000+150*j^2)/(iter*2000+150*iter^2)*100))
    if(!ceiling((j*2000+150*j^2)/(iter*2000+150*iter^2)*100)==ceiling(((j-1)*2000+150*(j-1)^2)/(iter*2000+150*iter^2)*100))
    {
      print(sent)
    }
    
    U=matrix(0,n,m)
    U2=matrix(0,n,m)
    
    #E-Step-Parallelizable:
    for(k in 1:n)
    {
      initial=matrix(0,1,m)
      print(k)
      initial[Data[k,][Data[k,]>0]]=sort(runif(sum(Data[k,]>0)),decreasing=TRUE)
      Temp=GibbsSampler(S,initial,pi=Data[k,],dist,parameter)
      U[k,]=Temp$M1
      U2[k,]=Temp$M2
    }
    # Mstep
    
    ##############
    ##############
    
    
    Delta=1/n*colSums(U)
    if(dist != "exp") Delta[1]=1
    
    if(dist == "norm.fixedvariance") Variance=matrix(1,1,m)
    else if(dist == "exp") Variance = NA
    else{
      Variance=abs(1/n*colSums(U2)-Delta^2)
      Variance[1]=1
    }
        
    if(dist=="norm"){parameter=list(Mu=Delta,Var=Variance)}
    if(dist=="exp"){parameter=list(Mu=Delta/sum(Delta))}
    
    MM[j,]=Delta
    VV[j,]=Variance
    print(Delta)
    if(dist=="norm")(print(Variance))
#     
#     if(dist=="norm")
#     {
#       ranget=max(abs(parameter$Mu))+3*max((Variance)^.5)
#       ll[j]=likelihood(Data,parameter,dist,range=ranget,res=ranget*.0001)
#       print(ll)
#     }
#     if(dist=="exp")
#     {
#       ranget=max(abs(parameter$Mu))+3*max(abs(parameter$Mu))^.5
#       ll[j]=likelihood(DataL,parameter,dist,range=ranget,res=ranget*.0001)
#       print(ll)
#       #ll[j]=likelihoodPL(DataL,parameter)
#     }
#     
  }  
  
  #if(dist=="exp")(Variance=Delta^-1)
  
  RT=sort(MM[iter,],decreasing=TRUE,index=TRUE)
  DT=proc.time()-T0
  
  ### William's addition
  params <- rep(list(list()), m)
  for(i in 1:m) {
    if(dist == "exp") params[[i]]$Mean <- 1/Delta[i]
    else params[[i]]$Mean <- Delta[i]
    params[[i]]$SD <- Variance[i]^.5
  }
  ###
  
  #;list(Aggregated.Rank=RT$ix, Mean=MM,Variance=VV,LogLikelihood=ll,Time=DT) 
  ;list(Aggregated.Rank=RT$ix, Mean=MM, SD=VV^.5, Time=DT, Parameters = params) 
  
}

#' Performs parameter estimation for a Multitype Random Utility Model
#' 
#' This function supports RUMs 
#' 1) Normal
#' 2) Normal with fixed variance (fixed at 1)
#' 3) Exponential
#' 
#' @param Data data in either partial or full rankings
#' @param K number of components in mixture distribution
#' @param iter number of EM iterations to run
#' @param dist underlying distribution. Can be "norm", "norm.fixedvariance", "exp"
#' @param ratio parameter in the algorithm that controls the difference of the starting points, the bigger the ratio the more the distance
#' @return results from the inference
#' @export
#' @examples
#' data(Data.Test)
#' Estimation.RUM.MultiType.MLE(Data.Test, K=2, iter = 3, dist= "norm.fixedvariance", ratio=.2)
Estimation.RUM.MultiType.MLE = function(Data, K=2, iter = 10, dist, ratio)
{
  if(!(dist=="dexp" | dist=="exp" | dist=="norm" | dist == "norm.fixedvariance" ))
  {
    stop(paste0("Distribution name \"", dist, "\" not recognized"))
  }
  # calculating the dimensions
  dims <- dim(Data)
  n <- dims[1]
  m <- dims[2]
  
  
  #initialization of mean and mixture probabilities
  
  initiateMean = matrix(exp(rnorm(m*K)),K,m)
  for(k in 1:K)
  {
    initiateMean[k,] =  initiateMean[1,] + runif(m)*ratio
  }
  #initiateVar = matrix(exp(rnorm(m*K)),K,m)
  #initiateGamma = runif(K)
  #initiateGamma=initiateGamma/sum(initiateGamma)
  initiateGamma=matrix(1/K,1,K)
  initiateZ = matrix(sample(1:K, size=n, replace = TRUE, prob = NULL),1,n)
  
  if(dist=="exp")
  {
    DataL=Data
    for(i in 1:dim(Data)[1])
    {
      DataL[i,]=Data[i,m:1]
    }
    temp1=Data
    Data=DataL
    DataL=temp1
  }
  
  T0=proc.time()
  Delta=initiateMean
  #Variance=initiateVar
  Variance=matrix(1,K,m)
  Z=initiateZ
  gamma=initiateGamma
  
  ##############
  if(dist=="norm" | dist == "norm.fixedvariance"){parameter=list(Mu=Delta,Var=Variance)}
  if(dist=="exp"){parameter=list(Mu=Delta)}
  
  ##############
  MM=matrix(0,iter,m)
  VV=matrix(0,iter,m)
  ll=matrix(0,1,iter)
  
  for(j in 1:iter)
  {
    #S=1000+300*j
    S=100+100*j
    sent=sprintf("%d%% completed",ceiling((j*2000+150*j^2)/(iter*2000+150*iter^2)*100))
    if(!ceiling((j*2000+150*j^2)/(iter*2000+150*iter^2)*100)==ceiling(((j-1)*2000+150*(j-1)^2)/(iter*2000+150*iter^2)*100))
    {
      print(sent)
    }
    
    #X=matrix(0,n,m)
    z=matrix(0,1,n)
    ESTEP2=matrix(0,n,K)
    ESTEP1=matrix(0,K,m)
    #E-Step-Parallelizable:
    #############
    #############
    
    for(i in 1:n)
    {
      initial=matrix(0,1,m)
      #print(i)
      initial[Data[i,][Data[i,]>0]]=sort(runif(sum(Data[i,]>0)),decreasing=TRUE)
      X=initial
      for(t1 in 1:S)
      { 
        #print(parameter$Mu)
        ZProb=gamma*PdfModel(X,parameter)
        ZProb=ZProb/sum(ZProb)
        #print(ZProb)
        zi= sample(1:K, size=1, replace = TRUE, prob = ZProb)
        XTemp=GibbsSampler(3,X,pi=Data[i,],dist,parameter=list(Mu=Delta[zi,],Var=Variance))
        #X[i,]=XTemp$M1
        X=XTemp$M1
        ESTEP1[zi,]=ESTEP1[zi,]+XTemp$M1
        ESTEP2[i,zi]=ESTEP2[i,zi]+1
      }
    }
    
    # Mstep
    ##############
    ##############
    
    for(k in 1:K)
    {
      Delta[k,]=1/(sum(ESTEP2[,k]))*ESTEP1[k,]
    }
    Delta[,1]=matrix(1,1,K)
    #Variance=matrix(1,K,m)
    
    for(k in 1:K)
    {
      gamma[k]=sum(ESTEP2[,k])
    }
    gamma = gamma/sum(gamma)
    
    #gamma=gamma+matrix(1/K,1,K)
    
    gamma = gamma/sum(gamma)
    
    if(dist=="norm"){parameter=list(Mu=Delta,Var=Variance)}
    if(dist=="exp")
    { 
      for(k in 1:K)
      {
        Delta[k,]=Delta[k,]/sum(Delta[k,])
      }
      parameter=list(Mu=Delta)
    }
    
    #############
    #############
    
    #MM[j,]=Delta
    #VV[j,]=Variance
    #print("Delta")
    #print(Delta)
    #print("gamma")
    #print(gamma)
    #print("ESTEP2")
    #print(ESTEP2)
    
    #     if(dist=="norm")
    #     {
    #       ranget=max(abs(parameter$Mu))+3*max((Variance)^.5)
    #       ll[j]=likelihood(Data,parameter,dist,range=ranget,res=ranget*.0001)
    #       print(ll)
    #     }
    #     if(dist=="exp")
    #     {
    #       ranget=max(abs(parameter$Mu))+3*max(abs(parameter$Mu))^.5
    #       ll[j]=likelihood(DataL,parameter,dist,range=ranget,res=ranget*.0001)
    #       print(ll)
    #       #ll[j]=likelihoodPL(DataL,parameter)
    #     }
    #     
  }  
  
  #RT=sort(MM[iter,],decreasing=TRUE,index=TRUE)
  DT=proc.time()-T0
  
  #      if(dist=="norm")
  #      {
  #        ranget=max(abs(parameter$Mu))+3*max((Variance)^.5)
  #        ll[j]=likelihood(Data,parameter,dist,range=ranget,res=ranget*.0001)
  #        print(ll)
  #      }
  #      if(dist=="exp")
  #      {
  #        ranget=max(abs(parameter$Mu))+3*max(abs(parameter$Mu))^.5
  #        ll[j]=likelihood(DataL,parameter,dist,range=ranget,res=ranget*.0001)
  #        print(ll)
  #        #ll[j]=likelihoodPL(DataL,parameter)
  #      }
  
  params <- rep(list(list()), m)
  for(i in 1:m) {
    params[[i]]$Mean <- Delta[,i]
    params[[i]]$SD <- Variance[,i]^.5
    params[[i]]$Gamma <- gamma[1,]
  }
  ###
  
  
  ;list(Mean=Delta, SD=Variance^.5, Gamma=gamma, Personal=ESTEP2, Time=DT, Parameters = params) 
  
}



# The PDF of a set of independent normal random variables
# 
# 
# @param X vector of scalar values
# @param parameter list containing vetor of means for normal distributions, variances are set to one
# @return value of the PDF for vector X given the parameter
# @export
PdfModel = function(X,parameter)
{
  m=length(X)
  K=dim(parameter$Mu)[1]
  out=matrix(0,1,K)
  
  for(k in 1:K)
  {
    out[k]=prod(dnorm(X,mean=parameter$Mu[k,],sd=1))
  }
  ;out
}


################################################################################
## GMM
################################################################################

#' GMM Method for estimating Plackett-Luce model parameters
#' 
#' @param Data.pairs data broken up into pairs
#' @param m number of alternatives
#' @param prior magnitude of fake observations input into the model
#' @return Estimated mean parameters for distribution of underlying exponential
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Estimation.PL.GMM(Data.Test.pairs, 5)
Estimation.PL.GMM = function(Data.pairs, m, prior = 0)
{
  transition.matrix <- matrix(data = prior, nrow = m, ncol = m) - prior * m * diag(m)
  for(l in 1:nrow(Data.pairs)) {
    # i loses, j wins
    i <- Data.pairs[l, 2]    
    j <- Data.pairs[l, 1]
    transition.matrix[i, j] <- transition.matrix[i, j] + 1
  }
  
  diag(transition.matrix) <- diag(transition.matrix) - min(diag(transition.matrix))
  transition.matrix <- transition.matrix/rowSums(transition.matrix)
  
  first.eigenvector <- Re(eigen(t(transition.matrix))$vectors[,1])
  stationary.probability <- first.eigenvector / sum(first.eigenvector)
  
  # this is the estimated means
  list(Mean = stationary.probability, Parameters = convert.vector.to.list.of.means(stationary.probability))
}

#' GMM Method for Estimating Random Utility Model wih Normal dsitributions
#' 
#' @param Data.pairs data broken up into pairs
#' @param m number of alternatives
#' @param prior magnitude of fake observations input into the model
#' @return Estimated mean parameters for distribution of underlying normal (variance is fixed at 1)
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Estimation.Normal.GMM(Data.Test.pairs, 5)
Estimation.Normal.GMM = function(Data.pairs, m, prior = 0)
{
  
  muhat <- matrix(1,1,m)
  C <- normalizeC(generateC(Data.pairs, m, Zemel = FALSE, prior))
  for(iter in 1:1000)
  {
    alpha <- iter^-1
    muhat <- muhat + alpha * rowSums(exp(-delta(muhat)^2/4)*(C - f(muhat)))
    muhat <- muhat - min(muhat)
    #muhat <- muhat / sum(muhat)
  }
  
  list(Mean = muhat, SD = rep(1, m), Parameters = convert.vector.to.list.of.means(muhat))
}

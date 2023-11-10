#'@title Sample size calculation for Binomial data
#'@description
#'Find the sample size needed to have a desired false discovery rate and power for a large number of two-group comparisons under Binomial distribution
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param p1  probability in one group (vector)
#'@param p2  probability in other group (vector)
#'@param alternative one- or two-sided test
#'@param pi0.hat Approximation method for null proportion
#'@export

n.fdr.binomial=function(fdr,pwr,p1,p2,alternative="two.sided",pi0.hat="BH")
{
  pi0=mean(p1==p2)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.binomial, p1=p1,p2=p2,alternative=alternative)
  return(res)

}


#'@title Sample size calculation for Fisher's Exact tests
#'@description
#'Find the sample size needed to have a desired false discovery rate and power for a large number of fisher tests
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param p1  probability in one group (vector)
#'@param p2  probability in other group (vector)
#'@param alternative one- or two-sided test
#'@param pi0.hat Approximation method for null proportion
#'@export
n.fdr.fisher=function(fdr,pwr,p1,p2,alternative="two.sided",pi0.hat="BH")
{
  pi0=mean(p1==p2)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.fisher, p1=p1,p2=p2,alternative=alternative)
  return(res)

}

#'@title Sample size calculation for Poisson data
#'@description
#' Find the sample size needed to have a desired false discovery rate and power for a large number of two-group comparisons under Poisson distribution
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param rho Fold-change, usual null hypothesis is that rho=1 (vector)
#'@param mu0 Average count in control group
#'@param w Ratio of total number of
#'@param type Type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score.
#'@param pi0.hat Approximation method for null proportion
#'@references C-I Li, P-F Su, Y Guo, and Y Shyr (2013). Sample size calculation for differential expression analysis of RNA-seq data under Poisson distribution. Int J Comput Biol Drug Des 6(4). doi:10.1504/IJCBDD.2013.056830
#'@export
n.fdr.poisson=function(fdr,pwr,rho,mu0,w=1,type="w",pi0.hat="BH")
{
  pi0=mean(rho==1)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.li, rho=rho,mu0=mu0)
  return(res)

}

#'@title Sample size calculation for signed rank tests
#'@description
#' Find the sample size needed to have a desired false discovery rate and power for a large number of signed rank tests
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@param pi0.hat Approximation method for null proportion
#'@references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#'@export
n.fdr.signrank=function(fdr,pwr,p,pi0.hat="BH"){
  pi0=mean(p==0.25)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.signrank, p=p)
  return(res)

}

#'@title Sample size calculation for sign tests
#'@description
#' Find the sample size needed to have a desired false discovery rate and power for a large number of sign tests
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@param pi0.hat Approximation method for null proportion
#'@references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#'@examples
#'n.fdr.signtest(fdr=0.1,pwr=0.8,p=rep(c(0.8,0.5), c(100,9900)),pi0.hat="BH")
#'@export
n.fdr.signtest=function(fdr,pwr,p,pi0.hat="BH"){
  pi0=mean(p==0.5)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.signtest, p=p)
  return(res)

}

#'@title Sample size calculation for rank-sum tests
#'@description
#' Find the sample size needed to have a desired false discovery rate and average power for a large number of rank-sum tests
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@param pi0.hat Approximation method for null proportion
#'@references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#'@export
n.fdr.ranksum=function(fdr,pwr,p,pi0.hat="BH"){
  pi0=mean(p==0.5)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.ranksum, p=p)
  return(res)

}

#'@title Sample size calculation for one-way ANOVA
#'@description
#' Find the sample size needed to have a desired false discovery rate and power for a large number of one-way ANOVA tests
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param theta sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test (vector)
#'@param pi0.hat Approximation method for null proportion
#'@param k the number of groups to be compared
#'@export
n.fdr.oneway=function(fdr,pwr,theta,k,pi0.hat="BH"){
  pi0=mean(theta==0)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(alpha=a,pwr=pwr,avepow.func=average.power.oneway, theta=theta)
  return(res)

}


#'@title Sample size calculation for the Cox proportional hazards regression model
#'@description
#' Find number of events needed to have a desired false discovery rate and power for a large number of Cox regression models with nonbinary covariates
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param logHR log hazard ratio (vector)
#'@param v variance of predictor variable (vector)
#'@param pi0.hat Approximation method for null proportion
#'@references Hsieh, FY and Lavori, Philip W (2000) Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled Clinical Trials 21(6):552-560.
#'@export
n.fdr.coxph=function(fdr,pwr,logHR,v,pi0.hat="BH")

{
  pi0=mean(logHR==0)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(a,pwr,average.power.coxph,logHR=logHR,v=v)
  return(res)
}


#'@title Sample size calculation for negative binomial data
#'@description
#' Find the sample size needed to have a desired false discovery rate and power for a large number of negative binomial comparisons
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param log.fc log fold-change (vector), usual null hypothesis is log.fc=0
#'@param mu read depth per gene (vector, same length as log.fc)
#'@param sig coefficient of variation (CV) per gene (vector, same length as log.fc)
#'@param pi0.hat Approximation method for null proportion
#'@references SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013). Calculating Sample Size Estimates for RNA Sequencing Data. Journal of Computational Biology 20: 970-978.
#'@export
n.fdr.negbin=function(fdr,pwr,log.fc,mu,sig,pi0.hat="BH")

{
  pi0=mean(log.fc==0)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(a,pwr,average.power.hart,
                       log.fc=log.fc,mu=mu,sig=sig)
  return(res)
}


#'@title Sample size calculation for t-tests
#'@description
#' Find the sample size needed to have a desired false discovery rate and power for a large number of t-tests
#'@param fdr Desired FDR (scalar numeric)
#'@param pwr Power of single α-level test (scalar numeric)
#'@param delta difference of population means (vector)
#'@param sigma standard deviation (vector or scalar)
#'@param type type of t-test
#'@param pi0.hat Approximation method for null proportion
#'@export
n.fdr.ttest=function(fdr,pwr,delta,sigma=1,type="two.sample",pi0.hat="BH")
{

  pi0=mean(delta==0)
  a=alpha.power.fdr(fdr,pwr,pi0,pi0.hat)
  res=find.sample.size(a,pwr,average.power.t.test,
                       delta=delta,sigma=sigma,type=type)
  return(res)

}


#' @title Determines the sample size needed to achieve a desired power while controlling the FDR at a specified level.
#' @description
#' Determines the sample size needed to achieve a desired power while controlling the FDR at a specified level.
#' @param alpha Fixed p-value threshold  (scalar numeric)
#' @param pwr Power of single α-level test (scalar numeric)
#' @param avepow.func an R function to compute average power
#' @param n0 Lower limit for initial sample size range
#' @param n1 Upper limit for initial sample size range
#' @param max.its Number of iterations
#' @param ... Additional arguments to average power function
#' @export
find.sample.size=function(alpha,pwr,avepow.func,n0=3,n1=6,max.its=50,...)
{
  # compute average power for sample size limits n0 and n1
  pwr0=avepow.func(n0,alpha,...)
  pwr1=avepow.func(n1,alpha,...)

  # double sample size limit until desired average power is reached
  n.its=0
  while((pwr1<pwr)&&(n.its<max.its))
  {
    pwr0=pwr1
    n0=n1
    n1=2*n1
    pwr1=avepow.func(n1,alpha,...)
    n.its=n.its+1
  }

  # use bisection on sample size limits to get sample size
  while(((n1-n0)>1)&&(n.its<max.its))
  {
    n.mid=floor((n1+n0)/2)
    pwr.mid=avepow.func(n.mid,alpha,...)
    n.its=n.its+1
    if (pwr.mid>pwr) {n=n.mid; pwr=pwr.mid}
    else {n0=n.mid; pwr0=pwr.mid}
  }

  res=list(n=n1,
           computed.avepow=pwr1,
           desired.avepow=pwr,
           alpha=alpha,
           n.its=n.its,
           max.its=max.its,
           n0=n0,n1=n1)
           #other.input=list(...))

  return(res)



}

#' @title Compute p-value threshold for  given the proportion pi0 of tests with a true null
#' @description
#' Given the proportion pi0 of tests with a true null, find the p-value threshold that results in a desired FDR control and a desired power
#' @param fdr Desired FDR (scalar numeric)
#' @param pwr Power of single α-level test (scalar numeric)
#' @param pi0 The proportion of tests with a true null hypothesis.
#' @param method Approximation method includes "HH" (3-part histogram height) , "PC" (3-part histogram expected value),"BH"(Benjamini & Hochberg 1995),"Jung"(Jung 2005)
#' @references Pounds, Stan, and Cheng Cheng. "Sample size determination for the false discovery rate." Bioinformatics 21.23 (2005): 4263-4271. Gadbury GL, et al. (2004) Power and sample size estimation in high dimensional biology. Statistical Methods in Medical Research 13(4):325-38.
#'Jung, Sin-Ho."Sample size for FDR-control in microarray data analysis." Bioinformatics 21.14 (2005): 3097-3104.
#' @export

alpha.power.fdr=function(fdr,pwr,pi0,method="HH")

{
  res=NA

  # 3-part histogram height
  if (method=="HH")
  {
    a=-(1-fdr)*pi0
    b=(fdr*(1-pi0)*pwr+(1-fdr)*pi0+(1-pi0)*(1-pwr))
    c0= -(1-pi0)*pwr*fdr

    res=(-b+sqrt(b^2-4*a*c0))/(2*a)
    return(res)
  }

  # 3-part histogram expected value
  if (method=="PC")
  {
    a=(1-pi0)
    b=(1-pwr)*(1-pi0)+pi0*(1-fdr)
    c0= -fdr*(1-pi0)*pwr
    res=(-b+sqrt(b^2-4*a*c0))/(2*a)
    return(res)
  }

  # Benjamini & Hochberg 1995
  if (method=="BH")
  {
    res=(1-pi0)*pwr*fdr/(1-fdr*pi0)
  }

  # Jung 2005
  if (method=="Jung")
  {
    res=(1-pi0)*pwr*fdr/(pi0*(1-fdr))
  }

  if (is.na(res))
    stop('No calculation performed; method must be "HH", "PC", "BH", or "Jung".')

  return(res)
}

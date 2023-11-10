#'@title Compute the average power of many Cox regression models for a given number of events, p-value threshold, vector of effect sizes (log hazard ratio),and variance of predictor variables
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param logHR log hazard ratio (vector)
#'@param v variance of predictor variable (vector)
#'@examples
#'average.power.coxph(n=30,alpha=0.05,logHR=log(rep(c(1,2),c(9900,100))),v=rep(1,10000))
#'@export
average.power.coxph=function(n,alpha,logHR,v)

{
  m=length(logHR)
  if (length(v)==1) v=rep(v,m)
  alt=(logHR!=0)
  logHR=logHR[alt]
  v=v[alt]
  res=mean(power.cox(n,alpha,logHR,v))
  return(res)
}


#'@title Compute the power of a single-predictor Cox regression model
#'@description Use the formula of Hseih and Lavori (2000) to compute the power of a single-predictor Cox model
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param logHR log hazard ratio (vector)
#'@param v variance of predictor variable (vector)
#'@return vector of power estimates for two-sided test
#'@references Hsieh, FY and Lavori, Philip W (2000) Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled Clinical Trials 21(6):552-560.
#'@importFrom stats pnorm qnorm
#'@export
power.cox=function(n,
                   alpha,
                   logHR,
                   v)
{
  pnorm(qnorm(alpha/2),sqrt(n*v)*logHR)+
    1-pnorm(qnorm(1-alpha/2),sqrt(n*v)*logHR)
}



#'@title Compute average power for RNA-seq experiments assuming negative binomial distribution
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param log.fc log fold-change (vector), usual null hypothesis is log.fc=0
#'@param mu read depth per gene (vector, same length as log.fc)
#'@param sig coefficient of variation (CV) per gene (vector, same length as log.fc)
#'@details
#'This function is based on equation (1) of Hart et al (2013). It assumes a negative binomial model for RNA-seq read counts and equal sample size per group.
#'@export
average.power.hart=function(n,alpha,log.fc,mu,sig){
  alt=(log.fc!=0)
  log.fc=log.fc[alt]
  pwr=power.hart(n=n,alpha=alpha,log.fc=log.fc,mu=mu,sig=sig)
  res=mean(pwr)
  return(res)

}

#'@title Compute power for RNA-seq experiments assuming negative binomial distribution
#'@description
#'Use the formula of Hart et al (2013) to compute power for comparing RNA-seq expression across two groups assuming a negative binomial distribution
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param log.fc log fold-change (vector), usual null hypothesis is log.fc=0
#'@param mu read depth per gene (vector, same length as log.fc)
#'@param sig coefficient of variation (CV) per gene (vector, same length as log.fc)
#'@details
#'This function is based on equation (1) of Hart et al (2013). It assumes a negative binomial model for RNA-seq read counts and equal sample size per group.
#'@return vector of power estimates for the set of two-sided tests
#'@references SN Hart, TM Therneau, Y Zhang, GA Poland, and J-P Kocher (2013). Calculating Sample Size Estimates for RNA Sequencing Data. Journal of Computational Biology 20: 970-978.
#'@importFrom stats pnorm qnorm
#'@export
power.hart=function(n,
                    alpha,
                    log.fc,
                    mu,
                    sig)
{
  z.alpha=qnorm(alpha/2)
  res=pnorm(z.alpha,sqrt(n*log.fc^2/(2*(1/mu+sig^2))))+
    pnorm(-z.alpha,sqrt(n*log.fc^2/(2*(1/mu+sig^2))),lower.tail=F)
  return(res)
}


#'@title Compute average power of many t-tests
#'@description
#'Compute average power of many t-tests;Uses classical power formula for t-test;Assumes equal variance and sample size
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param delta difference of population means (vector)
#'@param sigma standard deviation (vector or scalar)
#'@param type type of t-test: two.sample, one.sample
#'@importFrom stats power.t.test
#'@export
average.power.t.test=function(n,alpha,delta,sigma,type="two.sample")
{
  alt=(delta!=0)
  delta=delta[alt]
  pwr=power.t.test(n,delta,sigma,alpha)
  res=mean(pwr$power)
  return(res)
}



#'@title Compute average power of Rank Sum tests
#'@param n sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@export
average.power.ranksum=function(n,alpha,p)

{
  alt=(p!=0.5)
  p=p[alt]
  pwr=power.ranksum(n,alpha,p)
  res=mean(pwr)
  return(res)

}

#'@title Compute power of Rank Sum tests
#'@description
#'Compute power of rank-sum test;Uses formula of Noether (JASA 1987)
#'@param n sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@details
#'In most applications, the null effect size will be designated by p = 0.5, which indicates that Thus, in the example below, the argument null.effect=0.5 is specified in the call to fdr.sampsize.
#'@return vector of power estimates for two-sided tests
#'@references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#'@importFrom stats qnorm pnorm
#'@export
power.ranksum=function(n,alpha,p)

{
  mu0=0.5*n*n
  mu1=p*n*n
  delta=(mu1-mu0)
  sig2=n*n*(2*n+1)/12
  sig=sqrt(sig2)
  z.reject=-abs(qnorm(alpha/2,0,sig))
  pow=pnorm(z.reject,delta,sig)+pnorm(-z.reject,delta,sig,lower.tail=F)
  return(pow)
}


#'@title Compute average power of many one-way ANOVA tests
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param theta sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test(vector)
#'@param k the number of groups to be compared
#'@export
average.power.oneway=function(n,alpha,theta,k)

{
  alt=(theta!=0)
  theta=theta[alt]
  pwr=power.oneway(n,alpha,theta,k)
  res=mean(pwr)
  return(res)

}
#'@title Compute power of one-way ANOVA
#'@description
#'Compute power of one-way ANOVA;Uses classical power formula for ANOVA;Assumes equal variance and sample size
#'@param n per-group sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param theta sum of ((group mean - overall mean)/stdev)^2 across all groups for each hypothesis test(vector)
#'@param k the number of groups to be compared, default k=2
#'@details
#'For many applications, the null effect is zero for the parameter theta described above
#'@returns vector of power estimates for test of equal means
#'@importFrom stats pf qf
#'@export
power.oneway=function(n,alpha,theta,k=2)

{
  1-pf(qf(1-alpha,k-1,k*(n-1)),k-1,k*(n-1),ncp=n*theta)
}



#'@title Compute average power of many sign tests
#'@param n sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@export
#'@importFrom stats qnorm pnorm
average.power.signtest=function(n,alpha,p)

{
  alt=(p!=0.5)
  p=p[alt]
  pwr=power.signtest(n,alpha,p)
  res=mean(pwr)
  return(res)

}


#'@title Compute power of the sign test
#'@description
#'Use the Noether (1987) formula to compute the power of the sign test
#'@param n sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@details
#' In most applications, the null effect size will be designated by p = 0.5 instead of p = 0. Thus, in the call to fdr.sampsize, we specify null.effect=0.5 in the example below.
#'@return vector of power estimates for two-sided tests
#'@references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#'@export
power.signtest=function(n,alpha,p)
{
  mu0=0.5*n
  mu1=p*n
  sig0=sqrt(n*0.25)
  sig1=sqrt(n*p*(1-p))
  pow=pnorm(qnorm(alpha/2,mu0,sig0),mu1,sig1)+
    pnorm(qnorm(1-alpha/2,mu0,sig0),mu1,sig1,lower.tail=F)
  return(pow)
}



#'@title Compute average power of many signrank tests
#'@param n sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@export
average.power.signrank=function(n,alpha,p)

{
  alt=(p!=0.5)
  p=p[alt]
  pwr=power.signrank(n,alpha,p)
  res=mean(pwr)
  return(res)

}


#'@title Compute power of signrank test
#'@description
#'Use the Noether (1987) formula to compute the power of the signrank test
#'@param n sample size (scalar)
#'@param alpha p-value threshold (scalar)
#'@param p Pr(Y>X), as in Noether (JASA 1987)
#'@references Noether, Gottfried E (1987) Sample size determination for some common nonparametric tests. Journal of the American Statistical Association, 82:645-647.
#'@importFrom stats pnorm qnorm
#'@export
power.signrank=function(n,
                        alpha,
                        p)
{
  mu0=0.25*n*(n+1)
  mu1=0.5*p*n*(n+1)
  delta=(mu1-mu0)
  sig2=n*(n+1)*(2*n+1)/24
  sig=sqrt(sig2)
  z.reject=-abs(qnorm(alpha/2,0,sig))
  pow=pnorm(z.reject,delta,sig)+pnorm(-z.reject,delta,sig,lower.tail=F)
  return(pow)
}



#'@title Compute Average Power for RNA-Seq Experiments Assuming Poisson Distribution
#'@description
#'Use the formula of Li et al (2013) to compute power for comparing RNA-seq expression across two groups assuming the Poisson distribution.
#'@param n per-group sample size
#'@param alpha p-value threshold (scalar)
#'@param rho Fold-change, usual null hypothesis is that rho=1 (vector)
#'@param mu0 Average count in control group
#'@param w Ratio of total number of
#'@param type Type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score
#'@details
#'This function computes the average power for a series of two-sided tests defined by the input parameters. The power is based on the sample size formulas in equations 10-13 of Li et al (2013). Also, note that the null.effect is set to 1 in the examples because the usual null hypothesis is that the fold-change = 1.
#'@return Average power estimates for a series of two-sided tests
#'@references C-I Li, P-F Su, Y Guo, and Y Shyr (2013). Sample size calculation for differential expression analysis of RNA-seq data under Poisson distribution. Int J Comput Biol Drug Des 6(4). doi:10.1504/IJCBDD.2013.056830
#'@export
average.power.li=function(n,alpha,rho,mu0,w=1,type="w")

{
  alt=(rho!=1)
  rho=rho[alt]
  pwr=power.li(n,alpha=alpha,rho=rho,mu0=mu0,w=w,type=type)
  res=mean(pwr)
  return(res)

}

#'@title Compute Power for RNA-Seq Experiments Assuming Poisson Distribution
#'@description
#'Use the formula of Li et al (2013) to compute power for comparing RNA-seq expression across two groups assuming the Poisson distribution
#'@param n per-group sample size
#'@param alpha p-value threshold (scalar)
#'@param rho Fold-change, usual null hypothesis is that rho=1 (vector)
#'@param mu0 Average count in control group
#'@param w Ratio of total number of
#'@param type Type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for log-transformed score
#'@details
#'This function computes the power for each of a series of two-sided tests defined by the input parameters. The power is based on the sample size formulas in equations 10-13 of Li et al (2013). Also, note that the null.effect is set to 1 in the examples because the usual null hypothesis is that the fold-change = 1.
#'@return vector of power estimates for two-sided tests
#'@references C-I Li, P-F Su, Y Guo, and Y Shyr (2013). Sample size calculation for differential expression analysis of RNA-seq data under Poisson distribution. Int J Comput Biol Drug Des 6(4). doi:10.1504/IJCBDD.2013.056830
#'@importFrom stats pnorm qnorm
#'@export
power.li=function(n, #> \item{n}{per-group sample size}
                  alpha, #> \item{alpha}{p-value threshold}
                  rho, #> \item{rho}{fold-change, usual null hypothesis is that rho=1}
                  mu0, #> \item{mu0}{average count in control group}
                  w=1, #> \item{w}{ratio of total number of }
                  type="w") #> \item{type}{type of test: "w" for Wald, "s" for score, "lw" for log-transformed Wald, "ls" for logtransformed score.}#> }
{
  z.alpha=qnorm(alpha/2)
  if (!is.element(type,c("w","lw","s","ls")))
    stop("type must be 'w', 's', 'ls', or 'lw'.")
  if (any(w<=0)) stop("w must be >0.")
  if (any(rho<=0)) stop("rho must be >0.")
  if (type=="w")
  {
    mu.z=sqrt(n*mu0*(rho-1)^2/(1+rho/w))
    res=pnorm(z.alpha,mu.z)+pnorm(-z.alpha,mu.z,lower.tail=F)
    return(res)
  }
  if (type=="lw")
  {
    mu.z=sqrt(n*mu0*(log(rho)^2)/(1+1/(rho*w)))
    res=pnorm(z.alpha,mu.z)+pnorm(-z.alpha,mu.z,lower.tail=F)
    return(res)
  }
  if (type=="s")
  {
    z.cut=z.alpha*sqrt((1+w*rho)/(w+rho))
    mu.z=sqrt(n*mu0*(rho-1)^2/(1+rho/w))
    res=pnorm(z.cut,mu.z)+pnorm(-z.cut,mu.z,lower.tail=F)
    return(res)
  }
  if (type=="ls")
  {
    z.cut=z.alpha*sqrt((2+w+1/w)/((1+w*rho)*(1+1/(w*rho))))
    mu.z=sqrt(n*mu0*(log(rho)^2)/(1+1/(rho*w)))
    res=pnorm(z.cut,mu.z)+pnorm(-z.cut,mu.z,lower.tail=F)
    return(res)
  }
}



#'@title Compute average power of many Fisher's exact tests
#'@param p1  probability of success in group1 (vector)
#'@param p2  probability of success in group1 (vector)
#'@param n per-group sample size
#'@param alpha p-value threshold
#'@param alternative one- or two-sided test
#'@importFrom stats dbinom fisher.test
#'@export
average.power.fisher=function(p1,p2,n,alpha,alternative) #n:per group sample size

{
  alt=(p1!=p2)
  p1=p1[alt]
  p2=p2[alt]
  n=rep(n,length(p1))

  pwr=0
  for (i in 1:length(p1))
  {
    temp=power.fisher(p1=p1[i],p2=p2[i],n=n[i],alpha=alpha,alternative=alternative)

    pwr=pwr+temp
  }
  pwr=pwr/length(p1)
  return(pwr)

}

#'@title Compute Power for Fisher's exact test
#'@param p1  probability in one group (vector)
#'@param p2  probability in other group (vector)
#'@param n per-group sample size
#'@param alpha p-value threshold
#'@param alternative one- or two-sided test
#'@export
#'@importFrom stats dbinom fisher.test
power.fisher=function(p1,p2,n,alpha,alternative)  #Calculation of power for Fisher's exact test for comparing two proportions by simulated dataset, n:per group sample size

{

  pr1=dbinom(0:n,n,p1,log=T)
  pr2=dbinom(0:n,n,p2,log=T)
  y1=rep(0:n,each=n+1)
  y2=rep(0:n,times=n+1)
  y=cbind(y1,n-y1,y2,n-y2)
  pr=exp(rep(pr1,each=n+1)+rep(pr2,times=n+1))
  #p.value <- rep(0,length(pr))
  res=0
  for (i in 1:length(pr))
  {
    pval <- fisher.test(matrix(y[i,],2,2),alternative=alternative)$p.value
    res=res+(pval<=alpha)*pr[i]
  }

  return(res)
}

#'@title Computer average power of many binomial comparisons
#'@param n  per-group sample size
#'@param p1  probability in one group
#'@param p2  probability in other group
#'@param alpha p-value threshold
#'@param alternative one- or two-sided test
#'@importFrom stats power.prop.test
#'@export
average.power.binomial=function(n,alpha,p1,p2,alternative)

{
  alt=(p1!=p2)
  p1=p1[alt]
  p2=p2[alt]
  pwr=power.prop.test(n=n, p1 = p1, p2 = p2, sig.level = alpha,alternative = alternative)$power
  res=mean(pwr)
  return(res)
}

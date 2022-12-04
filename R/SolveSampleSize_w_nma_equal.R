#' Calculate the optimal sample sizes for a new two-arm trial when analyze it with the existing network
#'
#' This function calculates the optimal sample size for each treatment group
#' when planning a new two-arm trial with 1) two arms exist in the existing network;
#' 2) the response variable is binary; 3) we analyze it with the existing network by NMA.
#' 4) the sample size for each treatment group is required to be equal.
#'
#' @param p1 Risk of treatment 1
#' @param p2 Risk of treatment 2
#' @param sigma Standard error of the estimated effect size (log odds ratio) between treatment 1 and treatment 2
#' @param power_level Power of test we want to obtain
#' @param sig.level Significance level, the default value is 0.05  
#' @return the total sample size for the new trial
#' @export
#' @example 
#' SolveSampleSize_Withprev_equal(p1 = 0.2, p2 = 0.3, sigma = 0.4, power_level = 0.8)

SolveSampleSize_Withprev_equal <- function(p1,p2,sigma,power_level,sig.level = 0.05){
  # cal the total sample size when we added one more condition: sample sizes are equal for both treatment groups
  beta1 = log(p1/(1-p1))
  beta2 = log(p2/(1-p2)) - log(p1/(1-p1))
  mu_1 <- exp(beta1)/(1+exp(beta1))^2
  mu_2 <- exp(beta1 + beta2)/(1+exp(beta1 + beta2))^2
  sigma_prev <- sigma
  n0=c(1,1)
  
  confun_withprev <- function(n){
    f = power_level-power_withprev(n)
    f = rbind(f,-n)
    return(list(ceq=NULL,c=f))
  }
  
  power_withprev <- function(n){
    size <- n/2
    var_inv <- 1/(1/(mu_1 * size) + 1/(mu_2 * size))+1/sigma_prev^2
    var <- 1/var_inv
    se <- sqrt(var)
    z <- beta2/se
    power <- pnorm(z-qnorm(1-sig.level/2))+pnorm(-z-qnorm(1-sig.level/2))
    return(power)
  }
  
  objfun=function(n){
    n[1]+n[2]
  }
  
  solution_temp <- NlcOptim::solnl(n0,objfun=objfun,confun=confun_withprev)$par
  solution_temp_int <- round(solution_temp[2],0)
  # get the even integer solution around this
  if(solution_temp_int %% 2 ==0){
    res <- solution_temp_int
  }else{
    n <- c(solution_temp_int+1,solution_temp_int-1)
    n <- n[n>0]
    power <- power_single(n)
    dat_para <- data.frame(n,power)
    dat_para <- dat_para[dat_para$power>=power_level,]
    res <- min(dat_para$n)
    #res <- solution_temp_int+1
  }
  return(res)
}


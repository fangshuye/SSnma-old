#' Calculate the optimal sample sizes for a new two-arm trial
#'
#' This function calculates the optimal sample size for each treatment group
#' when 1) the response variable is binary; 2) we analyze it in isolation (w/o the existing network);
#' 3) the sample size for each treatment group is required to be equal.
#'
#' @param p1 Risk of treatment 1
#' @param p2 Risk of treatment 2
#' @param power_level Power of test we want to obtain
#' @param sig.level Significance level, the default value is 0.05  
#' @return the total sample size for the new trial
#' @export
#' @example 
#' SolveSampleSize_Single_equal(p1 = 0.2, p2 = 0.3,power_level = 0.8)

SolveSampleSize_Single_equal <- function(p1,p2,power_level, sig.level = 0.05){
  beta1 = log(p1/(1-p1))
  beta2 = log(p2/(1-p2)) - log(p1/(1-p1))
  mu_1 <- exp(beta1)/(1+exp(beta1))^2
  mu_2 <- exp(beta1 + beta2)/(1+exp(beta1 + beta2))^2
  n0=c(1,1)
  
  power_single <- function(n){
    size <- n/2
    var <- 1/(mu_1 * size) + 1/(mu_2 * size)
    se <- sqrt(var)
    z <- beta2/se
    power <- pnorm(z-qnorm(1-sig.level/2))+pnorm(-z-qnorm(1-sig.level/2))
    return(power)
  }
  
  confun_single <- function(n){
    f = power_level-power_single(n)
    f = rbind(f,-n)
    return(list(ceq=NULL,c=f))
  }
  
  objfun=function(n){
    n[1]+n[2]
  }
  
  solution_temp <- NlcOptim::solnl(n0,objfun=objfun,confun=confun_single)$par
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

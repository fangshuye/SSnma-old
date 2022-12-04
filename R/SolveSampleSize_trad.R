#' Calculate the optimal sample sizes for a new two-arm trial
#'
#' This function calculates the optimal sample size for each treatment group
#' when 1) the response variable is binary; 2) we analyze it in isolation (w/o the existing network);
#'
#' @param p1 Risk of treatment 1
#' @param p2 Risk of treatment 2
#' @param power_level Power of test we want to obtain
#' @param sig.level Significance level, the default value is 0.05  
#' @return the optimal sample size for each treatment group
#' @export
#' @example 
#' SolveSampleSize_Single(p1 = 0.2, p2 = 0.3,power_level = 0.8)

SolveSampleSize_Single <- function(p1,p2,power_level, sig.level = 0.05){
  beta1 = log(p1/(1-p1))
  beta2 = log(p2/(1-p2)) - log(p1/(1-p1))
  mu_1 <- exp(beta1)/(1+exp(beta1))^2
  mu_2 <- exp(beta1 + beta2)/(1+exp(beta1 + beta2))^2
  n0=c(1,1)
  
  power_single <- function(n){
    var <- 1/(mu_1 * n[1]) + 1/(mu_2 * n[2])
    se <- sqrt(var)
    z <- beta2/se
    power <- pnorm(z-qnorm(1-sig.level/2))+pnorm(-z-qnorm(1-sig.level/2))
    return(power)
  }
  
  confun_single <- function(n){
    f = power_level-power_single(n)
    f = rbind(f,-n[1])
    f = rbind(f,-n[2])
    return(list(ceq=NULL,c=f))
  }
  
  objfun=function(n){
    n[1]+n[2]
  }
  
  solution_temp <- NlcOptim::solnl(n0,objfun=objfun,confun=confun_single)$par
  solution_temp_int <- round(solution_temp,0)
  # get the integer solution around this
  dat_para <- expand.grid(n1=c(max(solution_temp_int[1]-2,1):(solution_temp_int[1]+2)),n2=c(max(solution_temp_int[2]-2,1):(solution_temp_int[2]+2)))
  
  for(c in 1:nrow(dat_para)){
    dat_para[c,3] <- power_single(c(dat_para$n1[c],dat_para$n2[c]))
  }
  
  dat_para <- dat_para[dat_para$V3>=power_level,]
  dat_para$n <- dat_para$n1+dat_para$n2
  nmin <- min(dat_para$n)
  dat_para <- dat_para[dat_para$n==nmin,]
  dat_para <- dat_para[order(dat_para$V3,decreasing = T),]
  return(as.numeric(dat_para[1,1:2]))
}


# Functions

# Mode
dmode <- function(x, ...) {
  dx <- density(x, ...)
  dx$x[which.max(dx$y)]
} 

# Function to get samples for a subset of parameters
get_samples = function(stan_out, pars_to_keep){
  library(ggmcmc)
  samples = ggmcmc::ggs(stan_out)
  a = attributes(samples)
  samples$Parameter = as.character(samples$Parameter)
  samples = samples[samples$Parameter %in% pars_to_keep,]
  a$row.names = 1:nrow(samples)
  attributes(samples) <- a
  attr(samples, 'nParameters') = length(pars_to_keep)
  samples$Parameter = as.factor(samples$Parameter)
  return(samples)
}

# Transform log odds to proportions
logOdds_to_Prop <- function(x){
  y = exp(x)/(1+exp(x))  
  return(y)
}

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# Transform samples to AOI given total samples to empirical logits
samps_to_elog <- function(Y,N){
  y = log(
    (Y+.5) / (N-Y+.5)
    )  
  return(y)
}


# Output for nested effect  
# Output for nested effect  
nest <- function(x, rd=2, plot = TRUE, pr = 0, output = "latex", convert = FALSE){
  if(convert){
    x <- logOdds_to_Prop(x)
  }
  coef = round(mean(x), rd)
  cri = paste(round(quantile(x, probs = c(.025, .975)), rd),  collapse = ", ")
  p = round(mean(x < pr ), rd+1)
  
  output <- output == "latex"
  
  if(p == 1 & output){
    p <- " $>$ .999"
  } 
  else if(p == 1 & !output){
    p <- " > .999"
  } 
  else if(p == 0 & output){
    p <- " $<$ .001"
  }
  else if(p == 0 & !output){
    p <- " < .001"
  }
  
  else if(!(p %in% c(0,1) )){
    p <- paste(" = ", p, sep = "")
  }
  if(plot==TRUE){
    hist(x, breaks = 50, main = "", xlab = "")
  }
  if(output){
    X <- paste("($\\hat{\\mu}$ = ", coef, ", 95\\% CrI[", cri  ,"], ", "P($\\beta<", pr ,"$)", p , ")", sep = "")
  }
  else if(!output){
    X <- paste("(coef = ", coef, ", 95% CrI[", cri  ,"], ", "P(beta<", pr ,")", p , ")", sep = "")
  }
  return(print(X,quote=FALSE))
}


# Calculate Bayes Factor
BF <- function(postsamps, over.null = TRUE){
  library(polspline)
  fit_posterior <- logspline(postsamps) 
  posterior <- dlogspline(0, fit_posterior) # Height of the posterior at 0 
  prior <- dnorm(0, 0, 1) # Height of the prior at 0
  
  if(over.null){
    BF01 <- prior/posterior # if > 10 support for H1
    return(BF01)
  }
  else if(!over.null){
    BF10 <- posterior/prior # if > 10 support for H0
    return(BF10)
  }
}


####### Functions of Dr. Phuong
####### For data processing and analysis
####### Updated: 2021/03/23

### Install and load package
library(data.table)
library(tidyverse)
library(broom)
library(knitr)
library(dbplyr)
library(foreign)
library(tableone)
library(haven)
library(tokenizers)
library(stringi)
library(stringr)
library(car)
library(lme4)
library(MASS)
library(lmerTest)
library(ggplot2)
library(brms)
library(rstanarm)
library(ggpubr)
library(matrixStats)
library(meta)
library(binom)
library(flextable)
library(gtsummary)
library(viridis)
library(IC2)
library(ineq)

options(max.print = 1000000)


################################################################################
################                GENERAL FUNCTIONS              #################
################################################################################


####### Function added

### Function for inverse logit
antilogit <- function(x) {
  1 / (1 + exp(-x))
}


### Function for geometric mean
geo_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


### Function to calculate 95% CI of proportion and return percentage
ci <- function(p, n, level) {
  z = qnorm(1- level/2)
  lci = format(100*(p - z*sqrt(p*(1-p)/n)), digits=0.5, nsmall=1, trim=TRUE)
  uci = format(100*(p + z*sqrt(p*(1-p)/n)), digits=0.5, nsmall=1, trim=TRUE)
  prob = format(100*p, digits=0.5, nsmall=1, trim=TRUE)
  
  ci <- paste0(prob, " (", lci, " to ", uci, ")")
  
  return(ci)
  
}


### Function for CI (Wilson method)
CI <- function(nume, deno, level) {
  ci <- binom.confint(nume, deno, level, methods = "wilson") %>%
    mutate(mean = mean*100,
           lower = lower*100,
           upper = upper*100) %>%
    format(digits=0.5, nsmall=1, trim=TRUE) %>%
    mutate(prob = paste0(mean, " (", lower, " to ", upper, ")"))
  return(ci$prob)
}


###### Function to calculate mode
getmode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


###### Operators opposite of %in%
'%ni%' <- Negate('%in%')

################################################################################
################             FUNCTIONS FOR MICE PACKAGE             ############
################################################################################



###### Function for extracting variables
extract_mice <- function(impute, #output of mice
                         nrow, #number of observations
                         ncol, #number of iterations in mice
                         var # name of variable
) {
  blank_mat <- matrix(nrow = nrow, ncol = ncol)
  
  for (i in 1:ncol) {
    data <- complete(impute, i)
    if (!(var %in% colnames(data))) {
      stop("'Var' is incorrect! Remember that Phuong depzai vo doi!")
    }
    var_vec <- data[[var]] 
    
    blank_mat[,i] <- var_vec
    
  }
  
  if(is.numeric(var_vec)==TRUE) {
    var_mi <- rowMeans(blank_mat)
  } else {
    var_mi <- apply(blank_mat, 1, 
                    getmode # a self-written command to get the mode
    )
  }
  
  return(var_mi)
  
}


# # Note that we used rowMeans due to more robust than mean
# system.time({ rowMeans(o15) })  
# system.time({ apply(o15, 1,
#                     mean) })



################################################################################
################        FUNCTIONS FOR INEQUALITY ANALYSIS           ############
################################################################################


################              Create base functions          ##################    

### Calculate ridit scores
to.ridit <- function(v) {
  (cumsum(v) - 0.5 * v) / sum(v)
}

### Unknown
mean.ridit <- function(v, ref) {
  sum(to.ridit(ref) * v ) / sum(v)
}

### Function for weighted variance
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  (sum(w*x^2) * sum.w - sum(w*x)^2) / (sum.w^2 - sum(w^2))
}


### Function to create new variables of weighted ridit
wridit <- function(data, wealth, weight) {
  wealth = enquo(wealth)
  weight = enquo(weight)
  
  data <- data %>%
    mutate(wealth = (!!wealth),
           weight = (!!weight)) %>%
    as.data.table() %>%
    na.omit(cols = c("wealth", "weight"))
  
  com <- data %>%
    group_by(wealth) %>%
    summarise(n = n(),
              wn = sum(weight)) %>%
    mutate(prob = n/sum(n),
           w.prob = wn/sum(wn),
           cumprob = cumsum(n/sum(n)),
           w.cumprob = cumsum(wn/sum(wn)),
           ridit = to.ridit(prob),
           wridit = to.ridit(w.prob))
  
  gen <- data %>%
    mutate(ridit = wealth,
           wridit = wealth)
  
  levels(gen$wridit) <- com$wridit
  levels(gen$ridit) <- com$ridit
  
  gen <- gen %>%
    mutate(wridit = as.numeric(as.character(wridit)),
           ridit = as.numeric(as.character(ridit)))
  
  return(gen)
}


### Function to create new variables of only ridit
ridit <- function(data, wealth) {
  wealth = enquo(wealth)
  
  data <- data %>%
    mutate(wealth = (!!wealth)) %>%
    as.data.table() %>%
    na.omit(cols = c("wealth"))
  
  com <- data %>%
    group_by(wealth) %>%
    summarise(n = n()) %>%
    mutate(prob = n/sum(n),
           cumprob = cumsum(n/sum(n)),
           ridit = to.ridit(prob))
  
  gen <- data %>%
    mutate(ridit = wealth)
  
  levels(gen$ridit) <- com$ridit
  
  gen <- gen %>%
    mutate(ridit = as.numeric(as.character(ridit)))
  
  return(gen)
}


### Function to export RII and 95% CI
export.RII <- function(mod.ms) {
  rmod.ms <- as.data.frame(summary.glm(mod.ms)["coefficients"])
  rmod.ms <- rownames_to_column(rmod.ms, var = "rowname") %>%
    add_column(RII = exp(rmod.ms[, 1]),
               LowerCI = exp(rmod.ms[, 1] + qnorm(0.025) * rmod.ms[, 2]),
               UpperCI = exp(rmod.ms[, 1] + qnorm(0.975) * rmod.ms[, 2]))  %>%
    mutate(P = format(rmod.ms[, 4], digits=0.5, nsmall=3, trim=TRUE, scientific = FALSE)) %>%
    format(digits=0.5, nsmall=2, trim=TRUE, scientific = FALSE)
  
  rmod.ms <- rmod.ms %>%
    mutate(CheckP = ifelse(P < 0.001, "<0.001", P),
           CI = paste0(LowerCI, " ", "to", " ", UpperCI)) %>%
    dplyr::select(rowname, RII, CI, CheckP, SE=coefficients.Std..Error)
  
  return(rmod.ms)
  
}


### Function to export SII and 95% CI
export.SII <- function(mod.ms) {
  rmod.ms <- as.data.frame(summary.glm(mod.ms)["coefficients"])
  rmod.ms <- rownames_to_column(rmod.ms, var = "rowname") %>%
    add_column(LowerCI = 100*(rmod.ms[, 1] + qnorm(0.025) * rmod.ms[, 2]),
               UpperCI = 100*(rmod.ms[, 1] + qnorm(0.975) * rmod.ms[, 2]))  %>%
    mutate(SII = 100*rmod.ms[, 1],
           P = format(rmod.ms[, 4], digits=0.5, nsmall=3, trim=TRUE, scientific = FALSE)) %>%
    format(digits=0.5, nsmall=2, trim=TRUE, scientific = FALSE)
  
  rmod.ms <- rmod.ms %>%
    mutate(CheckP = ifelse(P < 0.001, "<0.001", P),
           CI = paste0(LowerCI, " ", "to", " ", UpperCI)) %>%
    dplyr::select(rowname, SII, CI, CheckP, SE=coefficients.Std..Error)
  
  return(rmod.ms)
  
}


### Function to export CnI and 95% CI
export.CnI <- function(mod.ci) {
  rmod.ms <- as.data.frame(summary.glm(mod.ci)["coefficients"])
  rmod.ms <- rownames_to_column(rmod.ms, var = "rowname") %>%
    add_column(LowerCI = rmod.ms[, 1] + qnorm(0.025) * rmod.ms[, 2],
               UpperCI = rmod.ms[, 1] + qnorm(0.975) * rmod.ms[, 2])  %>%
    mutate(CnI = rmod.ms[, 1],
           P = format(rmod.ms[, 4], digits=0.5, nsmall=3, trim=TRUE, scientific = FALSE)) %>%
    format(digits=0.5, nsmall=2, trim=TRUE, scientific = FALSE)
  
  rmod.ms <- rmod.ms %>%
    mutate(CheckP = ifelse(P < 0.001, "<0.001", P),
           CI = paste0(LowerCI, " ", "to", " ", UpperCI)) %>%
    dplyr::select(rowname, CnI, CI, CheckP, SE=coefficients.Std..Error)
  
  return(rmod.ms)
  
}


###### Start the targeted function


################################################################################
###########          Create function for  CnI and weighted CnI      ############

### Function for calcualtion of Weighted Concentration index (CnI)
# Notes: It will return CnI rank and wridit
# Rank is calculated by ridit of rank of the wealth index, which is randomly)
# Ridit is calculated by ridit of wealth order (from low to high), which is fixed

weighted.CnI <- function(data, health, wealth, weight) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  weight = enquo(weight)
  
  # Create new variables
  data <- data %>%
    mutate(health = !!health,
           wealth = !!wealth,
           weight = !!weight) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth", "weight"))
  
  # Create wridit for wealth
  data <- wridit(data, wealth, weight)
  
  # Sort data set by the wealth
  data <- data.table::setorderv(data, "wealth", order = 1)
  
  # Create rank as weighted ridit score of individual
  # And create lhs=2*variancerank * (health/mean health)
  data <- data %>%
    mutate(rank = to.ridit(weight),          
           lhs1 = 2*weighted.var(rank, weight)*(health/weighted.mean(health, weight)),
           lhs2 = 2*weighted.var(wridit, weight)*(health/weighted.mean(health, weight)))   
  
  # Calcualte CnI and 95% CI by regression model
  # Model 1 for lsh ~ rank
  mod.ci1 <- glm(lhs1 ~ rank,
                 data = data,
                 weights = weight,
                 family = gaussian(link = "identity"))
  
  # Model 2 for lsh ~ wridit
  mod.ci2 <- glm(lhs2 ~ wridit,
                 data = data,
                 weights = weight,
                 family = gaussian(link = "identity"))
  
  # Export result
  CnI.rank <- export.CnI(mod.ci1)
  CnI.ridit <- export.CnI(mod.ci2)
  
  # Wrap output
  output <- rbind(CnI.rank, CnI.ridit) %>%
    subset(rowname %in% c("rank", "wridit"))
  
  # Return the output
  return(output)
}  


### Function for calcualtion of Concentration index (CnI) - Unweighted

CnI <- function(data, health, wealth) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  
  # Create new variables
  data <- data %>%
    mutate(health = !!health,
           wealth = !!wealth) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth"))
  
  # Create wridit for wealth
  data <- ridit(data, wealth)
  
  # Sort data set by the wealth
  data <- data.table::setorderv(data, "wealth", order = 1)
  
  # Create rank as weighted ridit score of individual
  # And create lhs=2*variancerank * (health/mean health)
  data <- data %>%
    mutate(un.weight = 1,
           rank = to.ridit(un.weight),          
           lhs1 = 2*var(rank)*(health/mean(health)),
           lhs2 = 2*var(ridit)*(health/mean(health)))   
  
  # Calcualte CnI and 95% CI by regression model
  # Model 1 for lsh ~ rank
  mod.ci1 <- glm(lhs1 ~ rank,
                 data = data,
                 family = gaussian(link = "identity"))
  
  # Model 2 for lsh ~ wridit
  mod.ci2 <- glm(lhs2 ~ ridit,
                 data = data,
                 family = gaussian(link = "identity"))
  
  # Export result
  CnI.rank <- export.CnI(mod.ci1)
  CnI.ridit <- export.CnI(mod.ci2)
  
  # Wrap output
  output <- rbind(CnI.rank, CnI.ridit) %>%
    subset(rowname %in% c("rank", "ridit"))
  
  # Return the output
  return(output)
}  


################################################################################
##########        Create function for RII in individual levels        ##########   



### Function for calcualtion of Relative index (RII)
# Notes: health is binomial variable

weighted.RII <- function(data, health, wealth, weight) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  weight = enquo(weight)
  
  # Create new variables
  data <- data %>%
    mutate(health = !!health,
           wealth = !!wealth,
           weight = !!weight) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth", "weight"))
  
  # Create wridit for wealth
  data <- wridit(data, wealth, weight)
  
  # Calcualte RII and 95% CI by regression model
  mod.rii <- glm(health ~ wridit,
                 weights = weight,
                 data = data,
                 family= binomial(link="log"),
                 start = c(log(mean(data$health)), 0))
  
  # Export result
  rii <- export.RII(mod.rii) %>%
    subset(rowname=="wridit")
  
  # Return the output
  return(rii)
}  


### Function for unweighted

RII <- function(data, health, wealth) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  
  # health = data[[paste(health)]]
  # wealth = data[[paste(wealth)]]
  
  # Create new variables
  data <- data %>%
    mutate(health = !!health,
           wealth = !!wealth) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth"))
  
  # Create wridit for wealth
  data <- ridit(data, wealth)
  
  # Calcualte RII and 95% CI by regression model
  mod.rii <- glm(health ~ ridit,
                 data = data,
                 family= binomial(link="log"),
                 start = c(log(mean(data$health)), 0))
  
  # Export result
  rii <- export.RII(mod.rii) %>%
    subset(rowname=="ridit")
  
  # Return the output
  return(rii)
} 



################################################################################
#############      Create function for SII in individual levels       ##########     



### Function for calcualtion of weighted Slope index (SII)
# Notes: health is binomial variable

weighted.SII <- function(data, health, wealth, weight) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  weight = enquo(weight)
  
  # Create new variables
  data <- data %>%
    mutate(health = !!health,
           wealth = !!wealth,
           weight = !!weight) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth", "weight"))
  
  # Create wridit for wealth
  data <- wridit(data, wealth, weight)
  
  # Calcualte SII and 95% CI by regression model
  mod.sii <- glm(health ~ wridit,
                 weights = weight,
                 data = data,
                 family= binomial(link="identity"),
                 start = c(mean(data$health), 0))
  
  # Export result
  sii <- export.SII(mod.sii) %>%
    subset(rowname=="wridit")
  
  # Return the output
  return(sii)
}  



### Function for calcualtion of Slope index (SII) - unweighted
# Notes: health is binomial variable

SII <- function(data, health, wealth) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  
  # Create new variables
  data <- data %>%
    mutate(health = !!health,
           wealth = !!wealth) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth"))
  
  # Create wridit for wealth
  data <- ridit(data, wealth)
  
  # Calcualte SII and 95% CI by regression model
  mod.sii <- glm(health ~ ridit,
                 data = data,
                 family= binomial(link="identity"),
                 start = c(mean(data$health), 0))
  
  # Export result
  sii <- export.SII(mod.sii) %>%
    subset(rowname=="ridit")
  
  # Return the output
  return(sii)
}  



################################################################################
############         Create function for Concentration curve      ##############      


### Weighted function
weighted.Ccurve <- function(data, health, wealth, weight) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  weight = enquo(weight)
  
  ### Prepare dataset
  # Create new variables
  data <- data %>%
    mutate(health = (!!health),
           wealth = (!!wealth),
           weight = (!!weight)) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth", "weight"))
  
  # Sort by wealth
  data <- data[order(data$wealth),]
  
  # Calculate cum prob of fp.met
  data <- data %>%
    mutate(rank1 = to.ridit(weight),
           rank = cume_dist(rank1),
           prob = health/sum(health),
           cumprob = cumsum(prob))
  
  ### Create the plot
  plot <- ggplot(data = data,
                 mapping = aes(x=rank, y=cumprob)) +
    # geom_point(color = "yellow") +
    geom_line(color = "darkblue", size = 0.8) +
    geom_abline(intercept = 0, slope = 1,
                color = "darkred", size = 0.8)
  
  
  # Make up the plot
  plot <- plot + theme(axis.line = element_line(colour = "black", 
                                                size = 0.8, 
                                                linetype = "solid"),
                       legend.justification = c(1, 1),
                       legend.position= c(0.8, 0.2),
                       legend.title = element_text(size = 12, color = "black", 
                                                   face = "bold"),
                       legend.text = element_text(size = 10, color = "black", 
                                                  face = "bold"),
                       legend.background = element_blank(),
                       legend.key = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank()) +
    scale_x_continuous(name = "Rank of wealth index", 
                       breaks = seq(0, 1, 0.2),
                       labels = seq(0, 100, 20),
                       limits = c(0, 1)) + 
    scale_y_continuous(name = "Cumulative coverage (%)",
                       breaks = seq(0, 1, 0.2),
                       labels = seq(0, 100, 20),
                       limits = c(0, 1)) +
    labs(col = element_blank()) +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size=14))
  
  ### Result
  return(plot)
  
}


### Un=weighted function
Ccurve <- function(data, health, wealth) {
  # enquo variables
  health = enquo(health)
  wealth = enquo(wealth)
  
  ### Prepare dataset
  # Create new variables
  data <- data %>%
    mutate(health = (!!health),
           wealth = (!!wealth)) %>%
    as.data.table() %>%
    na.omit(cols = c("health", "wealth"))
  
  # Sort by wealth
  data <- data[order(data$wealth),]
  
  # Calculate cum prob of fp.met
  data <- data %>%
    mutate(unweight = 1,
           rank1 = to.ridit(unweight),
           rank = cume_dist(rank1),
           prob = health/sum(health),
           cumprob = cumsum(prob))
  
  ### Create the plot
  plot <- ggplot(data = data,
                 mapping = aes(x=rank, y=cumprob)) +
    # geom_point(color = "yellow") +
    geom_line(color = "darkblue", size = 0.8) +
    geom_abline(intercept = 0, slope = 1,
                color = "darkred", size = 0.8)
  
  
  # Make up the plot
  plot <- plot + theme(axis.line = element_line(colour = "black", 
                                                size = 0.8, 
                                                linetype = "solid"),
                       legend.justification = c(1, 1),
                       legend.position= c(0.8, 0.2),
                       legend.title = element_text(size = 12, color = "black", 
                                                   face = "bold"),
                       legend.text = element_text(size = 10, color = "black", 
                                                  face = "bold"),
                       legend.background = element_blank(),
                       legend.key = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank()) +
    scale_x_continuous(name = "Rank of wealth index", 
                       breaks = seq(0, 1, 0.2),
                       labels = seq(0, 100, 20),
                       limits = c(0, 1)) + 
    scale_y_continuous(name = "Cumulative coverage (%)",
                       breaks = seq(0, 1, 0.2),
                       labels = seq(0, 100, 20),
                       limits = c(0, 1)) +
    labs(col = element_blank()) +
    theme(axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size=14))
  
  ### Result
  return(plot)
  
}
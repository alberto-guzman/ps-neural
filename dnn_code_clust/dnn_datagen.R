# Data generation --------------------------------------------------------------------------------------------------------------------------------------------------------------------

#########################
# helper functions
#########################

# quick negation 

`%notin%` <- Negate(`%in%`)


#########################################
############## funcov(size, num_variables,nrep, scenarioT, scenarioY)
#########################################

#inputs: sample size n, scenario, num_variables
#output: a data set of size n in the chosen scenario, containing:
#covariates v1....vp
#treatment T
#outcome Y


funcov <- function(size, num_variables,nrep, scenarioT, scenarioY){


# n = 500
# p = 20
# nrep = 1

# induce correlations on non-normal variables

n = size
p = num_variables


  p_nor <- round((1/2*(p)),digits = 0)
  p_unif <- round((1/4*(p)),digits = 0)
  p_bino <- round((1/4*(p)),digits = 0)
  covars <- c()
  master_covar <- c()
  beta <- c()
  alpha <- c()
  i <- 1
  t <- 0

  # generates normal vars
  while (i <= p_nor) {
    x <- rnorm(n)
    t = t + 1
    assign(paste0('v',t), x)
    covar <- paste0('v',t)
    covars[[length(covars) + 1]] <- covar
    master_covar[[length(master_covar) + 1]] <- covar
    rm(x)
    i = i + 1
    t <- t
  }

  # generates unif vars
  i <- 1
  while (i <= p_unif) {
    x <- runif(n)
    t = t + 1
    assign(paste0('v',t), x)
    covar <- paste0('v',t)
    covars[[length(covars) + 1]] <- covar
    master_covar[[length(master_covar) + 1]] <- covar
    rm(x)
    i = i + 1
    t <- t
  }

  # generates binom vars
  i <- 1
  while (i <= p_bino) {
    g <- rbeta(1,1,1)
    x <- rbinom(n, 1, g)
    t = t + 1
    assign(paste0('v',t), x)
    covar <- paste0('v',t)
    master_covar[[length(master_covar) + 1]] <- covar
    rm(x)
    i = i + 1
    t <- t
  }


  # induce correlations
  covars_cor <- sample(covars,((.5)*length(covars)))
  
  for (i in covars_cor) {
    t = t + 1
    g <- rbeta(1,1,1)
    x <- F.sample.cor(get(i), g)
    assign(paste0('v',t), x)
    covar <- paste0('v',t)
    master_covar[[length(master_covar) + 1]] <- covar
    rm(x)
    t <- t
  }
  

  
  # ~~~~~~~~~~~~~~Global Variables~~~~~~~~~~~~~~~~~~~~#
  
  # ~~ coefficients for the treatment and potential outcomes equations
  
  b0 <- 0
  
  for (i in seq_len(length(master_covar))) {
    g <- round(rbeta(1,1,1),2)
    assign(paste0('b', i), g)
    b <- paste0('b',i)
    beta[[length(beta) + 1]] <- b
  }
    
  # ~~ coefficients for potential outcomes equation Y(0) and Y(1)
  # ~~ the intercept in Y(1) has been fixed to obtain Y(1)-Y(0)=-0.4 in each scenario
  
  a0 <- -3.85
  a1 <- 0.3
  
  for (i in 2:length(master_covar)) {
    g <- round(rbeta(1,1,1),2)
    assign(paste0('a', i), g)
    a <- paste0('a',i)
    alpha[[length(alpha) + 1]] <- a
  }
  
  
  #########################################
  # ~~ scenarios for data generation models
  #########################################
  
  # A: model with additivity and linearity
  # [1] "(1 + exp(-(0 + b1 * v1 + b2 * v2 + b16 * v16 + b14 * v14 + b4 * v4 + b7 * v7 + b8 * v8 + b13 * v13 + b5 * v5 + b10 * v10 + b9 * v9 + b3 * v3 + b15 * v15 )))^-1"
  
  # B: model with non-linearity
  # [1] "(1 + exp(-(0 + b1 * v1 + b2 * v2 + b16 * v16 + b14 * v14 + b4 * v4 + b7 * v7 + b8 * v8 + b13 * v13 + b5 * v5 + b10 * v10 + b9 * v9 + b3 * v3 + b15 * v15  +  b13*v13*v13 + b4*v4*v4 )))^-1"

  # C: model with non-additivity
  # [1] "(1 + exp(-(0 + b1 * v1 + b2 * v2 + b16 * v16 + b14 * v14 + b4 * v4 + b7 * v7 + b8 * v8 + b13 * v13 + b5 * v5 + b10 * v10 + b9 * v9 + b3 * v3 + b15 * v15  +  b15*0.5*v8*v15 + b8*0.5*v3*v8 )))^-1"

  # D: model with non-additivity and non-linearity
  # [1] "(1 + exp(-(0 + b1 * v1 + b2 * v2 + b16 * v16 + b14 * v14 + b4 * v4 + b7 * v7 + b8 * v8 + b13 * v13 + b5 * v5 + b10 * v10 + b9 * v9 + b3 * v3 + b15 * v15  +  v13*v13  +  b15*0.5*v8*v15 + v4*v4  +  b7*0.5*v13*v7 )))^-1"

  #########################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~ scenarios for treatment assignment
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #########################################
  
  toSample_master_covar <- master_covar[-1]
  
  
  #take sample of covars
  cov_80 <- sample(toSample_master_covar,((.8)*length(toSample_master_covar)))
  cov_20 <- sample(cov_80,((.2)*length(cov_80)))
  b <- sub(".*v","",cov_80)
  element <- paste0('b',b,' * ', cov_80)
  form <- paste("(1 + exp(-(0 + b1 * v1 +", paste(element, collapse = " + ")) 
  
  #cov_20 <- sample(element,((.2)*length(element)))

  if (scenarioT == "A") {
    
  #additivity and linearity 

  f_form <- paste("(1 + exp(-(0 + b1 * v1 +", paste(element, collapse = " + "),")))^-1") 
  trueps <- eval(parse(text = f_form))

    } else
  
      if (scenarioT == "B") {
        
      #non-lineary function
      add_lin <- paste0(cov_20, '*',cov_20)
      b_2 <- sub(".*v","",cov_20)
      f_form2 <- paste(form," + ", paste0('b',b_2,'*',add_lin, collapse = " + "),")))^-1")
      trueps <- eval(parse(text = f_form2))
      
        } else

          if (scenarioT == "C") {
            
          #non-addivity function
          n1 <- 0.5
          add_cov <- sample(cov_80,((.2)*length(cov_80)))
          b_3 <- sub(".*v","",add_cov)
          b_4 <- sample(b,length(add_cov))
          add <- paste0('b',b_3,'*',n1,'*','v',b_4,'*',add_cov)
          f_form3 <- paste(form," + ", paste(add, collapse = " + "),")))^-1") 
          trueps <- eval(parse(text = f_form3))
          
            } else
              
              if (scenarioT == "D") {
                
              #non-lineary + aditivity function
              n1 <- 0.5
              add_lin <- paste0(cov_20, '*',cov_20)
              add_cov <- sample(cov_80,((.2)*length(cov_80)))
              b_5 <- sub(".*v","",add_cov)
              b_6 <- sample(b,length(add_cov))
              add <- paste0('b',b_5,'*',n1,'*','v',b_6,'*',add_cov)
              f_form4 <- paste(form," + ", paste(add_lin, ' + ',add, collapse = " + "),")))^-1") 
              trueps <- eval(parse(text = f_form4))
              }
  
  #########################################
  # ~~ binary treatment T
  #########################################
  
  unif1 <- runif(n, 0, 1)
  T <- ifelse(trueps > unif1, 1, 0) # there is a probability of unif1 that T=1

  #########################################
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ~~ scenarios for outcome
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #########################################
  
  cov_60 <- sample(cov_80,((.6)*length(cov_80)))
  cov_rem <- master_covar[which(toSample_master_covar %notin% cov_80)]
  
  a_1 <- sub(".*v","",cov_60)
  element <- paste0('a',a_1,' * ', cov_60)
  
  a_2 <- sub(".*v","",cov_rem)
  element_rem <- paste0('a',a_2,' * ', cov_rem)

  form_y0 <- paste("Y0 <- a0 + a1 * v1 +", paste0(element, collapse = " + "),"+", paste0(element_rem, collapse = " + ")) 
  

  
  if (scenarioY == "a") {

    # linearity and aditivity
    Y0 <- eval(parse(text = form_y0))
    
    form_y1 <- paste("Y1 <- (a0 - 0.4) + (a1 + 0) * v1 + ", paste0('(a', a_1, '+0)','*','v',a_1, collapse = " + "),'+',paste0('(a', a_2, '+0)','*','v',a_2, collapse = " + ")) 
    Y1 <- eval(parse(text = form_y1))

    } else
  
      if (scenarioY == "b") {

        # non-linearity
        Y0 <- eval(parse(text = form_y0))
        
        cov_30 <- sample(cov_60,((.2)*length(cov_60)))
        a_3 <- sub(".*v","",cov_30)
        a_4 <- a_1[which(a_1 %notin% a_3)]
        form_y1_2 <- paste("Y1 <- (a0 - 0.4) + (a1 + 0) * v1 + ", paste0('(a', a_4, '+0)','*','v',a_4, collapse = " + "),'+',paste0('(a', a_3, '+0.25',"*",'a',a_3,')*','v',a_3, collapse = " + ")) 
        inter <- paste(form_y1_2,'+', paste0('(a', a_2, '+0)','*','v',a_2, collapse = " + "))
        Y1 <- eval(parse(text = inter))
                         
      } else

         if (scenarioY == "c") {

         # non-additivity
           cov_30 <- sample(cov_60,((.2)*length(cov_60)))
           a_3 <- sub(".*v","",cov_30)
           a_4 <- a_1[which(a_1 %notin% a_3)]
           form_y0_2 <- paste0(form_y0, '+',paste0('0.5', '*', 'a', a_3, '*','v',a_3,'^2', collapse = ' + '))
           Y0 <- eval(parse(text = form_y0_2))  
           
           form_y1_3 <- paste("Y1 <- (a0 - 0.2) + (a1 + 0) * v1 + ", paste0('(a', a_4, '+0)','*','v',a_4, collapse = " + "),'+',paste0('a', a_3, '*','v',a_3,'^2',"+",'(a', a_2, '+0)','*','v',a_2, collapse = " + "))
           Y1 <- eval(parse(text = form_y1_3))
  
         } else
         
          if (scenarioY == "d") {
            
            # non-linearity and non-adivitiy
            cov_30 <- sample(cov_60,((.2)*length(cov_60)))
            a_3 <- sub(".*v","",cov_30)
            a_4 <- a_1[which(a_1 %notin% a_3)]
            form_y0_3 <- paste0(form_y0,"+", paste0('0.5', '*', 'a', a_3, '*','v',a_3,'^2','+','a', a_2,'*','v',a_2, collapse = ' + '))
            Y0 <- eval(parse(text = form_y0_3))
            
            cov_30 <- sample(cov_60,((.2)*length(cov_60)))
            a_3 <- sub(".*v","",cov_30)
            a_4 <- a_1[which(a_1 %notin% a_3)]
            form_y1_4 <- paste("Y1 <- (a0 - 0.4) + (a1 + 0) * v1 + ", paste0('(a', a_4, '+0)','*','v',a_4, collapse = " + "),'+',paste0('(a', a_3, '+0.25',"*",'a',a_3,')*','v',a_3, collapse = " + ")) 
            form_y1_5 <- paste(form_y1_4,'+', paste0('(a', a_4, '+0)','*','v',a_4, collapse = " + "),'+',paste0('a', a_3, '*','v',a_3,'^2', collapse = " + ")) 
            inter <- paste(form_y1_5,'+', paste0('(a', a_2, '+0)','*','v',a_2, collapse = " + "))
            Y1 <- eval(parse(text = inter))
            
            } 
            

  # continuous outcome Y
  
  Y <- T * Y1 + (1 - T) * Y0
  
  # true individual effect of T on Y

  indeff <- Y1 - Y0
  
  # form simulation tibble

  v_list <- mget(paste0("v", 1:length(master_covar)))
  sim <- as_tibble(v_list)
  sim$T <- T
  sim$Y <- Y
  sim$indeff <- indeff
  return(sim)
  }
  



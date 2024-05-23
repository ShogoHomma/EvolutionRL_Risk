# Run simulation of model with two learning rates
# mainly used for the simulation in the tasks where the expected values are the same

RunQL_dual <- function(params, task_struct, trial_N) {
  
  # params
  ap <- params %>% dplyr::pull(ap)
  an <- params %>% dplyr::pull(an)
  beta <- params %>% dplyr::pull(bt)
  print(paste0("N of agent: ", length(ap)))
  
  # task struct
  m1 <- task_struct[1] # risky option
  sd1 <- task_struct[2]
  m2 <- task_struct[3] # non-risky option
  sd2 <- task_struct[4]
  
  print(paste0("task: N(", m1, ", ", sd1, "), N(", m2, ", ", sd2, ")"))
  
  sim_N <- nrow(params) # number of agents
  
  Choice <- matrix(0, nrow = sim_N, ncol = trial_N)
  Reward <- matrix(0, nrow = sim_N, ncol = trial_N)
  deltas <- matrix(0, nrow = sim_N, ncol = trial_N)
  
  p <- matrix(0, nrow = sim_N, ncol = trial_N) # p[sim_i, trial]
  Q <- array(0, dim = c(sim_N, trial_N + 1, 2)) 
  # Q[, t, 1] : risky option、Q[, t, 2] : non-risky option
  
  for (t in 1:trial_N) {
    
    ## choice
    p[, t] <- 1 / (1 + exp(-beta*( Q[, t, 1] - Q[, t, 2] )))
    
    sim_choice <- as.integer(rbernoulli(sim_N, p = p[, t]))
    sim_choice[sim_choice == 0] <- 2 # convert choice 0 -> choice 2
    Choice[, t] <- sim_choice
    
    ## get payoff
    reward_list <- rep(0, times = sim_N)
    
    Choices_1 <- Choice[, t] == 1 # agents who chose option 1
    if (length(which(Choices_1)) > 0) {
      
      sim_reward <- rnorm(n = length(which(Choices_1)), mean = m1, sd = sd1)
      reward_list[Choices_1] <- sim_reward
      
    }
    Choices_2 <- Choice[, t] == 2 # agents who chose option 2
    if (length(which(Choices_2)) > 0) {
      
      sim_reward <- rnorm(n = length(which(Choices_2)), mean = m2, sd = sd2)
      reward_list[Choices_2] <- sim_reward
      
    }
    
    Reward[, t] <- reward_list
    
    #-- agents who chose option 1
    if (length(which(Choices_1)) > 0) {
      
      deltas[Choices_1, t] <- Reward[Choices_1, t]-Q[Choices_1, t, 1]
      
      deltas_one_plus <- (deltas[, t] >= 0) & Choices_1 
      deltas_one_minus <- (deltas[, t] < 0) & Choices_1
      
      ## update Q value (Q1: high risk, Q2: low risk)
      # if RPE >= 0
      if (length(which(deltas_one_plus)) > 0) {
        
        Q[deltas_one_plus, t+1, 1] <- Q[deltas_one_plus, t, 1] + ap[deltas_one_plus] * deltas[deltas_one_plus, t] # for chosen option
        Q[deltas_one_plus, t+1, 2] <- Q[deltas_one_plus, t, 2] # for non-chosen option
        
      }
      
      # if RPE < 0
      if (length(which(deltas_one_minus)) > 0) {
        
        Q[deltas_one_minus, t+1, 1] <- Q[deltas_one_minus, t, 1] + an[deltas_one_minus] * deltas[deltas_one_minus, t] # for chosen option
        Q[deltas_one_minus, t+1, 2] <- Q[deltas_one_minus, t, 2] # for non-chosen option
        
      }
      
    }
    
    #-- agents who chose option 2
    if (length(which(Choices_2)) > 0) {
      
      deltas[Choices_2, t] <- Reward[Choices_2, t]-Q[Choices_2, t, 2]
      
      deltas_two_plus <- (deltas[, t] >= 0) & Choices_2
      deltas_two_minus <- (deltas[, t] < 0) & Choices_2
      
      ## update Q value (Q1: high risk, Q2: low risk)
      # if RPE >= 0
      if (length(which(deltas_two_plus)) > 0) {
        
        Q[deltas_two_plus, t+1, 2] <- Q[deltas_two_plus, t, 2] + ap[deltas_two_plus] * deltas[deltas_two_plus, t] # for chosen option
        Q[deltas_two_plus, t+1, 1] <- Q[deltas_two_plus, t, 1] # for non-chosen option
        
      }
      
      # if RPE < 0
      if (length(which(deltas_two_minus)) > 0) {
        
        Q[deltas_two_minus, t+1, 2] <- Q[deltas_two_minus, t, 2] + an[deltas_two_minus] * deltas[deltas_two_minus, t] # for chosen option
        Q[deltas_two_minus, t+1, 1] <- Q[deltas_two_minus, t, 1] # for non-chosen option
        
      }
      
    }
    
  }
  
  return(
    list(Q = Q, p = p, Reward = Reward, Choice = Choice)
  )
  
}

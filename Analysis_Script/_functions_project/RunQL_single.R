# 1つの学習率の強化学習シミュレーションを行う

RunQL_single <- function(params, task_struct, trial_N) {
  
  # params
  alpha <- params %>% dplyr::pull(alpha)
  beta <- params %>% dplyr::pull(bt)
  
  # task struct
  m1 <- task_struct[1] # リスキー
  sd1 <- task_struct[2]
  m2 <- task_struct[3] # 非リスキー
  sd2 <- task_struct[4]
  
  sim_N <- nrow(params) # number of agents
  
  Choice <- matrix(0, nrow = sim_N, ncol = trial_N)
  Reward <- matrix(0, nrow = sim_N, ncol = trial_N)
  deltas <- matrix(0, nrow = sim_N, ncol = trial_N)
  
  p <- matrix(0, nrow = sim_N, ncol = trial_N) # p[sim_i, trial]
  Q <- array(0, dim = c(sim_N, trial_N + 1, 2)) 
  # Q[, t, 1] がリスキー選択肢、Q[, t, 2] が非リスキー選択肢
  
  for (t in 1:trial_N) {
    
    ## choice
    p[, t] <- 1 / (1 + exp(-beta*( Q[, t, 1] - Q[, t, 2] )))
    
    sim_choice <- as.integer(rbernoulli(sim_N, p = p[, t]))
    sim_choice[sim_choice == 0] <- 2 # convert choice 0 -> choice 2
    Choice[, t] <- sim_choice
    
    ## get payoff
    reward_list <- rep(0, times = sim_N)
    
    Choices_1 <- Choice[, t] == 1 # 1を選んだ個体
    if (length(which(Choices_1)) > 0) {
      
      sim_reward <- rnorm(n = length(which(Choices_1)), mean = m1, sd = sd1)
      reward_list[Choices_1] <- sim_reward
      
    }
    Choices_2 <- Choice[, t] == 2 # 2を選んだ個体
    if (length(which(Choices_2)) > 0) {
      
      sim_reward <- rnorm(n = length(which(Choices_2)), mean = m2, sd = sd2)
      reward_list[Choices_2] <- sim_reward
      
    }
    
    Reward[, t] <- reward_list
    
    #-- agents who chose option 1
    if (length(which(Choices_1)) > 0) {
      
      deltas[Choices_1, t] <- Reward[Choices_1, t]-Q[Choices_1, t, 1]
      
      deltas_one <- deltas[, t] & Choices_1 
      
      ## update Q value (Q1: high risk, Q2: low risk)
      if (length(which(deltas_one)) > 0) {
        
        Q[deltas_one, t+1, 1] <- Q[deltas_one, t, 1] + alpha[deltas_one] * deltas[deltas_one, t] # for chosen option
        Q[deltas_one, t+1, 2] <- Q[deltas_one, t, 2] # for non-chosen option
        
      }
      
    }
    
    #-- agents who chose option 2
    if (length(which(Choices_2)) > 0) {
      
      deltas[Choices_2, t] <- Reward[Choices_2, t]-Q[Choices_2, t, 2]
      
      deltas_two <- deltas[, t] & Choices_2
      
      ## update Q value (Q1: high risk, Q2: low risk)
      if (length(which(deltas_two)) > 0) {
        
        Q[deltas_two, t+1, 2] <- Q[deltas_two, t, 2] + alpha[deltas_two] * deltas[deltas_two, t] # for chosen option
        Q[deltas_two, t+1, 1] <- Q[deltas_two, t, 1] # for non-chosen option
        
      }
      
    }
    
  }
  
  return(
    list(Q = Q, p = p, Reward = Reward, Choice = Choice)
  )
  
}

# Random課題の_nat.csvで、最初と最後の世代の情報だけ取り出す

Read_FirstLast_Rand <- function(path_i, data_path, GeneN, AgentN, TaskN) {
  
  read_df <- 
    data.table::fread(data_path[path_i]) %>% 
    dplyr::mutate(
      sim_i = path_i,
      generation = rep(0:(GeneN - 1), each = AgentN * TaskN),
      agent = rep(rep(0:(AgentN - 1), each = TaskN), times = GeneN)
    ) %>% 
    dplyr::filter(generation %in% c(0, GeneN-1)) %>% 
    dplyr::select(sim_i, generation, agent, task, everything())
  
  print(paste0("--- finish:", data_path[path_i]))
  
  return(read_df)
  
}

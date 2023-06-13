# read_trial data of random multiple-task simulation

MakeDfTrial_Rand <- function(dir_pattern, sim_name, DataDir, DataDir_sub, csv_pattern, Extract_SimN) {
  
  data_path <- make_DataPath(DataDir, DataDir_sub, dir_pattern, csv_pattern)
  data_path <- data_path[1:Extract_SimN]
  
  # to save memory, select specidfic columns and make summarised df for each simulation 
  df <- purrr::map_dfr(data_path, 
                       ~{print(.x); 
                         data.table::fread(.x) %>% 
                          dplyr::rename(
                            sim_i = sim
                          ) %>% 
                          dplyr::select(sim_i, generation, task, trial, p) %>% 
                          dplyr::group_by(sim_i, generation, task, trial) %>% 
                          dplyr::summarise(
                            mean_p = mean(p),
                            mean_risk_aversion = 1 - mean_p,
                            sd_p = sqrt(variance(p))
                          )
                         })　
  
  df <- df %>% 
    dplyr::rename(task_i = task) %>% 
    dplyr::mutate(sim = sim_name)
  
  return(df)
  
}

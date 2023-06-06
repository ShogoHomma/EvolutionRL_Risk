# Random Multiple Taskのリスク回避率を網羅的に載せる
# df_agent_randmulti_summary_bysimを渡し、横軸：シミュレーション・課題、縦軸：リスク回避率
# sim_range <- c(first, last)

VisFirstLastGene_RandMulti_RiskAversion_bySim <- function(df, sim_filt, tasktype, sim_range, main_title, sub_title, legend = TRUE) {
  
  all_sim <- df$sim %>% unique(.) 
  
  if (sim_filt %in% all_sim == FALSE) {
    
    stop("Error in sim_filt! Only 'sim_X/X can be accepted. See VisFirstLastGene_RandMulti_RiskAversion_bySim & sim_filt.")
    
  } 
  
  df_rev <-
    df %>% 
    dplyr::filter(sim == sim_filt) %>% 
    dplyr::mutate(generation_rev = generation + 1) %>% 
    dplyr::filter((sim_i >= sim_range[1]) & (sim_i <= sim_range[2]))

  
  if (tasktype == "risk-aversion") {
    
    df_rev <- 
      df_rev %>% 
      dplyr::filter(task_type == "risk_aversion")
    
  } else if (tasktype == "risk-seeking") {
    
    df_rev <- 
      df_rev %>% 
      dplyr::filter(task_type == "risk_seeking")
    
  } else {
    
    stop("Error in tasktype! Only 'risk-aversion' or 'risk-seeking' can be accepted. See VisFirstLastGene_RandMulti_RiskAversion_bySim(). & tasktype")
    
  }
  
  g <-
    ggplot(data = df_rev, aes(x = sim_n_sequence, y = mean_risk_aversion, 
                              fill = as.factor(generation_rev), group = generation_rev)) + 
    geom_linerange(aes(ymin = mean_risk_aversion - SD_rate_choice1, ymax = mean_risk_aversion + SD_rate_choice1), 
                   size = 0.9, position = position_dodge(width = 0.6)) +
    geom_point(size = 2.0, stroke = 1.3, shape = 21, position = position_dodge(width = 0.6)) + 
    my_theme2 + 
    theme(axis.text.x = element_text(size = 11, angle = 55, hjust = 1), 
          axis.text.y = element_text(size = 17), 
          axis.title = element_text(size = 28, face = "bold"), 
          strip.text.x = element_text(size = 15),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.background = element_rect(fill = "transparent"), 
          legend.key = element_rect(fill = "white", size = 15), 
          legend.position = c(0.98, 0.38), 
          legend.justification = c(1, 1),
    ) + 
    coord_cartesian(ylim = c(0, 1)) + 
    scale_fill_manual(values = c("white", "black")) + 
    guides(fill = guide_legend(override.aes = list(size = 4))) + 
    labs(x = "simulation number - task number", y = "Rate of risk aversion", 
         title = main_title, subtitle = sub_title, fill = "Generation") 
  
  if (legend == FALSE) {
    
    g <- g + guides(fill = "none")
    
  }
  
  return(g)
  
}

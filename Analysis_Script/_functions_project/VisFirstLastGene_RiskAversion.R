# horizontal axis : task, vertical axis : risk aversion rate of the first and last generation

VisFirstLastGene_RiskAversion <- function(df, target_simexp, sub_title, tag, legend = TRUE) {
  
  g <-
    df %>% 
    dplyr::filter(generation %in% c(0, 4999)) %>% 
    dplyr::filter(task_exp %in% target_simexp) %>% 
    dplyr::mutate(
      task_sd = paste0("Variance of Risky Option = ", sd1),
      generation = generation + 1
      ) %>% 
    
    ggplot(aes(x = task_exp, y = mean_risk_aversion, fill = as.factor(generation), group = as.factor(generation))) + 
    geom_linerange(aes(ymin = mean_risk_aversion - sd_choice1, ymax = mean_risk_aversion + sd_choice1), 
                   position = position_dodge(width = 0.8), size = 1.0) + 
    geom_point(size = 2.5, stroke = 1.3, position = position_dodge(width = 0.8), shape = 21) + 
    my_theme2 + 
    # old theme 
    # theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
    #       axis.text.y = element_text(size = 15), 
    #       axis.title = element_text(size = 23, face = "bold"), 
    #       strip.text.x = element_text(size = 17),
    #       legend.title = element_text(size = 16), #change legend title font size
    #       legend.text = element_text(size = 16),
    #       legend.background = element_rect(fill = "transparent"), 
    #       legend.key = element_rect(fill = "white"), 
    #       legend.position = c(1, 0.36), 
    #       legend.justification = c(1, 1),
    #       plot.tag = element_text(size = 30)
    # ) + 
    theme(axis.text.x = element_text(size = 13.5, angle = 30, hjust = 1), 
          axis.text.y = element_text(size = 20), 
          axis.title = element_text(size = 25, face = "bold"), 
          strip.text.x = element_text(size = 18),
          legend.title = element_text(size = 19), #change legend title font size
          legend.text = element_text(size = 19),
          legend.background = element_rect(fill = "transparent"), 
          legend.key = element_rect(fill = "white"), 
          legend.position = c(1, 0.40), 
          legend.justification = c(1, 1),
          plot.tag = element_text(size = 30)
    ) + 
    scale_fill_manual(values=c("white", "black")) + 
    coord_cartesian(ylim = c(0, 1)) + 
    scale_x_discrete(limits = target_simexp) + 
    facet_grid(.~task_sd) + 
    labs(x = "Location of two distributions (risky vs safe)", 
         y = "Rate of risk aversion",
         fill = "Generation",
         subtitle = sub_title,
         tag = tag)
  
  if (legend == FALSE) {
    
    g <- g + guides(fill = "none")
  }
  
  return(g)
  
}

# horizontal axis : task, vertical axis : parameter value
# colored area : SD averaged across replications

VisLastGene_Para_RunMean <- function(df, target_simexp, params, sub_title, tag, legend = TRUE) {
  
  if (params == "all") {  # ap, an, & bt
    
    df_rev <-
      df %>% 
      dplyr::filter(generation == 4999) %>% # 最終世代のみ取り出す
      dplyr::filter(task_exp %in% target_simexp) %>% 
      dplyr::mutate(task_sd = paste0("Variance of Risky Option = ", sd1))
    
  } else if (params == "alpha") {
    
    df_rev <-
      df %>% 
      dplyr::filter(generation == 4999) %>% # 最終世代のみ取り出す
      dplyr::filter(task_exp %in% target_simexp) %>% 
      dplyr::mutate(task_sd = paste0("Variance of Risky Option = ", sd1)) %>% 
      dplyr::filter(para == "an" | para == "ap")
    
  } else if (params == "beta") {
    
    df_rev <-
      df %>% 
      dplyr::filter(generation == 4999) %>% # 最終世代のみ取り出す
      dplyr::filter(task_exp %in% target_simexp) %>% 
      dplyr::mutate(task_sd = paste0("Variance of Risky Option = ", sd1)) %>% 
      dplyr::filter(para == "bt")
    
    
  } else {
    
    stop("Error! See VisLastGene_Para_RunMean(). & params")
    
  }
  
  # 共通部分のグラフを作成
  
  g <- 
    ggplot(data = df_rev, aes(x = task_exp, y = mean_rate, color = para, group = para)) + 
    geom_linerange(aes(ymin = mean_rate - sd_rate, ymax = mean_rate + sd_rate), size = 1.4, position = position_dodge(width = 0.5)) + 
    geom_line(position = position_dodge(width=0.7), size = 1.2, linetype = "twodash") + 
    geom_point(size = 2.5, position = position_dodge(width = 0.5), stroke = 2, shape = 21, fill = "white") + 
    
    my_theme2 + 
    # old theme
    # theme(
    #   axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
    #   axis.text.y = element_text(size = 15), 
    #   axis.title = element_text(size = 23, face = "bold"), 
    #   strip.text.x = element_text(size = 17), 
    #   legend.text = element_text(family = "serif", face = "italic"), 
    #   legend.background = element_rect(fill = "transparent"), 
    #   legend.key = element_rect(fill = "white"), 
    #   legend.position = c(1, 0.38), 
    #   legend.justification = c(1, 1),
    #   plot.tag = element_text(size = 30)
    #   ) + 
    theme(
      axis.text.x = element_text(size = 13.5, angle = 30, hjust = 1), 
      axis.text.y = element_text(size = 20), 
      axis.title = element_text(size = 25, face = "bold"), 
      strip.text.x = element_text(size = 18),
      legend.text = element_text(family = "serif", face = "italic"), 
      legend.background = element_rect(fill = "transparent"), 
      legend.key = element_rect(fill = "white"), 
      legend.position = c(1, 0.40), 
      legend.justification = c(1, 1),
      plot.tag = element_text(size = 30)
    ) + 
    scale_x_discrete(limits = target_simexp) + 
    facet_grid(.~task_sd) + 
    labs(x = "Location of two distributions (risky vs safe)", 
         y = "Parameter value",
         subtitle = sub_title,
         tag = tag) 
  
  if (params == "all") {  # ap, an, & bt
    
    g_rev <- 
      g + 
      scale_color_manual(name = "", values = c("#1F78B4", "#FF7F00", "#33A02C"), labels = c(expression(paste(alpha[n])), expression(paste(alpha[p])), expression(beta))) +
      coord_cartesian(ylim = c(0, 1.0))
    
  } else if (params == "alpha") {
    
    g_rev <- 
      g + 
      scale_color_manual(name = "", values = c("#1F78B4", "#FF7F00"), labels = c(expression(paste(alpha[n])), expression(paste(alpha[p])))) +
      coord_cartesian(ylim = c(0, 1.0))
    
  } else if (params == "beta") {
    
    g_rev <- 
      g + 
      scale_color_manual(name = "", values = c("#33A02C"), labels = expression(beta)) +
      coord_cartesian(ylim = c(0, 0.5))
    
  } else {
    
    stop("Error! See VisLastGene_Para_RunMean(). & colors")
    
  }
  
  if (legend == FALSE) {
    
    g_rev <- g_rev + guides(color = "none")
  }
  
  return(g_rev)
  
  
}

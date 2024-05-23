# comprehensively display the evolutionary dynamics of parameters in the single-task simulation

# horizontal axis：generation, vertical axis：parameter value, colored area: SD averaged across replications
# row：non-risky option (seven)、column：risk of risky option (five)

VisEvol_Para_RunMean_facet_sd1_nonrisky <- function(df, params, tasktype, abs_Dvalue, maintext_sim, legend = TRUE) {
  
  df_rev <-
    df %>% 
    dplyr::mutate(
      risky = paste0('N(', m1, ',', sd1, ')'),
      safe = factor(paste0('N(', m2, ',', sd2, ')')) %>% forcats::fct_reorder(., m2),
      sim_sort_tmp = factor(paste0('N(', m1, ',', sd1, ') vs N(', m2, ',', sd2, ')')),
      D = m1 - m2,
      abs_D = abs(D)
    ) %>% 
    dplyr::filter(abs_D == abs_Dvalue) %>% 
    dplyr::mutate(
      maintext_sim_sign = if_else(sim %in% maintext_sim, " *", ""),
      sim_sort = paste0(sim_sort_tmp, maintext_sim_sign)
    )
  
  # --- filter by params
  
  if (params == "all") {  # ap, an, & bt
    
    print("We set params to all.")
    
  } else if (params == "alpha") {
    
    df_rev <-
      df_rev %>% 
      dplyr::filter(para == "an" | para == "ap")
    
  } else if (params == "beta") {
    
    df_rev <-
      df_rev %>% 
      dplyr::filter(para == "bt")
    
  } else {
    
    stop("Error! See VisEvol_RunMean_sd1_nonrisky(). & params")
    
  }
  
  # --- tasktypeに応じてfilter
  
  if (tasktype == "risk-aversion") {
    
    df_rev <-
      df_rev %>% 
      dplyr::filter(m2 > m1) 
    
    #main_title <- "Risk-aversion task"
    sub_title <- paste0("Risk-aversion task (D = -", abs_Dvalue, ")")
    
  } else if (tasktype == "risk-seeking") {
    
    df_rev <-
      df_rev %>% 
      dplyr::filter(m2 < m1) 
    
    #main_title <- "Risk-seeking task"
    sub_title <- paste0("Risk-seeking task (D = +", abs_Dvalue, ")")
    
  } else {
    
    stop("Error in tasktype! Only 'risk-aversion' or 'risk-seeking' can be accepted. See VisEvol_RunMean_sd1_nonrisky(). & tasktype")
    
  }

  # label部分
  
  label_df <-
    df_rev %>% 
    dplyr::distinct(sd1, m2, safe, sim_sort)
  
  
  # 共通部分のグラフを作成
  
  g <- 
    ggplot(data = df_rev, aes(x = generation, y = mean_rate)) + 
    geom_line(aes(x = generation, y = mean_rate, color = para), size = 1.0) + 
    geom_ribbon(aes(x = generation, ymin = mean_rate - sd_rate, ymax = mean_rate + sd_rate, fill = para), alpha = 0.5) + 
    my_theme2 + 
    
    theme(
      axis.text.x = element_text(size = 15), 
      axis.text.y = element_text(size = 15), 
      axis.title = element_text(size = 30, face = "bold"), 
      strip.text.x = element_text(size = 19),
      strip.text.y = element_text(size = 17),
      legend.text = element_text(family = "serif", face = "italic"), 
      legend.background = element_rect(fill = "transparent"), 
      legend.key = element_rect(color = "transparent", fill = "white"), 
      legend.position = c(1, 0.35), 
      legend.justification = c(1, 1)
    ) + 
    #labs(x= 'Generation', y= 'Parameter value', title = main_title, subtitle = sub_title) + 
    labs(x= 'Generation', y= 'Parameter value', subtitle = sub_title) + 
    scale_x_continuous(breaks = c(0, 2500, 5000)) + 
    
    facet_grid(safe~sd1, labeller = labeller(sd1 = label_both, safe = label_both))
    
  
  if (params == "all") {  # ap, an, & bt
    
    g_rev <- 
      g + 
      scale_color_manual(name = "", values = c("#1F78B4", "#FF7F00", "#33A02C"), labels = c(expression(paste(alpha[n])), expression(paste(alpha[p])), expression(beta))) +
      scale_fill_manual(name = "", values = c("#1F78B4", "#FF7F00", "#33A02C"), labels = c(expression(paste(alpha[n])), expression(paste(alpha[p])), expression(beta))) +
      coord_cartesian(ylim = c(0, 1.0))
    
  } else if (params == "alpha") {
    
    g_rev <- 
      g + 
      scale_color_manual(name = "", values = c("#1F78B4", "#FF7F00"), labels = c(expression(paste(alpha[n])), expression(paste(alpha[p])))) +
      scale_fill_manual(name = "", values = c("#1F78B4", "#FF7F00"), labels = c(expression(paste(alpha[n])), expression(paste(alpha[p])))) +
      coord_cartesian(ylim = c(0, 1.0))
    
  } else if (params == "beta") {
    
    g_rev <- 
      g + 
      scale_color_manual(name = "", values = c("#33A02C"), labels = expression(beta)) +
      scale_fill_manual(name = "", values = c("#33A02C"), labels = expression(beta)) + 
      coord_cartesian(ylim = c(0, 0.5))
    
  } else {
    
    stop("Error! See VisEvol_Para_RunMean_facetSim(). & colors")
    
  }
  
  if (tasktype == "risk-aversion") { # リスク回避課題なら
    
    g_rev <- 
      g_rev + geom_text(x = 2500, y = 0.05, aes(label = sim_sort), data = label_df, size = 5.5, fontface = "bold")
    
  } else if (tasktype == "risk-seeking") { # リスク追求課題なら
    
    if (params == "all") {
      y_coord <- 0.98
    } else if (params == "alpha") {
      y_coord <- 0.98
    } else if (params == "beta") {
      y_coord <- 0.04
    }
    
    g_rev <- 
      g_rev + geom_text(x = 2500, y = y_coord, aes(label = sim_sort), data = label_df, size = 5.5, fontface = "bold")
    
  } else {
    
    stop("Error in tasktype! Only 'risk-aversion' or 'risk-seeking' can be accepted. See VisEvol_RunMean_sd1_nonrisky(). & tasktype")
    
  }
  
  
  if (legend == FALSE) {
    
    g_rev <- g_rev + guides(color = "none", fill = "none")
    
  }
  
  return(g_rev)
  
}

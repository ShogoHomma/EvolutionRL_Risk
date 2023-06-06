# 世代内シミュレーションのdfを使い、ヒートマップ：横軸 ap, 縦軸 an

VisHeatmap_RiskAversion_facet_sd1_nonrisky <- function(df, tasktype, abs_Dvalue, maintext_sim, bt_filt) {
  
  
  if (bt_filt %in% c(0.10, 0.25, 0.40) == FALSE) {
   
    stop("Error in bt_filt! 0.1, 0.25, 0.4 can be accepted. See VisHeatmap_RiskAversion_facet_sd1_nonrisky(). & bt_filt")
    
  }
  
  df_rev <- 
    df %>% 
    dplyr::filter(bt == bt_filt) %>% 
    dplyr::mutate(
      m1 = as.integer(m1),
      sd1 = as.integer(sd1),
      m2 = as.integer(m2),
      sd2 = as.integer(sd2)
    ) %>%
    dplyr::mutate(
      risky = paste0('N(', m1, ',', sd1, ')'),
      nonrisky = paste0('N(', m2, ',', sd2, ')') %>% forcats::fct_reorder(., m2),
      sim_sort_tmp = factor(paste0('N(', m1, ',', sd1, ') vs N(', m2, ',', sd2, ')')),
      D = m1 - m2,
      abs_D = abs(D)
    ) %>% 
    dplyr::filter(abs_D == abs_Dvalue) %>% 
    dplyr::mutate(
      maintext_sim_sign = if_else(sim %in% maintext_sim, " *", ""),
      sim_sort = paste0(sim_sort_tmp, maintext_sim_sign)
    )
    
  
  # --- tasktypeに応じてfilter
  
  if (tasktype == "risk-aversion") {
    
    df_rev <-
      df_rev %>% 
      dplyr::filter(m2 > m1) # リスク回避課題なのでm2の方が大きい課題を取り出す
    
    main_title <- paste0("Risk-aversion task (D = -", abs_Dvalue, ")")
    sub_title <- paste0("beta = ", bt_filt)
    
  } else if (tasktype == "risk-seeking") {
    
    df_rev <-
      df_rev %>% 
      dplyr::filter(m2 < m1) # リスク追求課題なのでm1の方が大きい課題を取り出す
    
    main_title <- paste0("Risk-seeking task (D = +", abs_Dvalue, ")")
    sub_title <- paste0("beta = ", bt_filt)
    
  } else {
    
    stop("Error in tasktype! Only 'risk-aversion' or 'risk-seeking' can be accepted. See VisHeatmap_RiskAversion_facet_sd1_nonrisky(). & tasktype")
    
  }
  
  
  # label部分
  label_df <-
    df_rev %>% 
    dplyr::distinct(sd1, m2, nonrisky, sim_sort)
  
  
  # Heatmapを描画
  g <-
    ggplot(data = df_rev, aes(x = ap, y = an)) + 
    geom_raster(aes(fill = mean_risk_aversion)) + 
    geom_text(x = 0.50, y = 0.95, aes(label = sim_sort), data = label_df, size = 4.7, fontface = "bold") + 
    
    my_theme2 + 
    theme(
      axis.text.x = element_text(size = 15, angle = 30, hjust = 1), 
      axis.text.y = element_text(size = 15), 
      axis.title = element_text(size = 25, face = "bold"), 
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 16),
      legend.text = element_text(size = 16), 
      legend.title = element_text(size = 16), 
      panel.grid = element_blank()
    ) + 
    scale_fill_gradient2(midpoint = 0.5, low = "#FF7F00", mid = "white", high = "#1F78B4", limits = c(0, 1)) + 
    labs(x = expression(alpha[p]), y = expression(alpha[n]), fill = "Risk aversion",
         title = main_title, subtitle = sub_title) + 
    facet_grid(nonrisky~sd1, labeller = labeller(sd1 = label_both, nonrisky = label_both)) 
  
  return(g)
  
}

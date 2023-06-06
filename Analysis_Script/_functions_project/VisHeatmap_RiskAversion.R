# 世代内シミュレーションのdfを使い、ヒートマップ：横軸 ap, 縦軸 an

VisHeatmap_RiskAversion <- function(df, target_sim, sub_title) {
  
  df_rev <- 
    df %>%
      dplyr::filter(sim %in% target_sim) %>% 
      dplyr::mutate(
        ap = round(ap, 4),
        an = round(an, 4),
        bt = round(bt, 4),
        rate_risk_aversion = 1 - choice_rate1,
        sim_fac = factor(sim, levels = target_sim) %>% as.numeric()
      ) %>% 
      dplyr::mutate(
        sim_sort = factor(paste0('N(', m1, ',', sd1, ') vs N(', m2, ',', sd2, ')')) %>% forcats::fct_reorder(., sim_fac, .desc = TRUE)
      ) %>%
      dplyr::group_by(sim_sort, ap, an, bt) %>%
      dplyr::summarise(
        mean_pvalue = mean(average_pvalue), 
        sd_pvalue = sqrt(variance(average_pvalue)),
        mean_risk_aversion = mean(rate_risk_aversion),
        sd_risk_aversion = sqrt(variance(rate_risk_aversion))
      ) %>%
      dplyr::filter(bt == 0.25)
    
  g <-
    ggplot(data = df_rev, aes(x = ap, y = an)) + 
    geom_raster(aes(fill = mean_risk_aversion)) + 
    my_theme2 + 
    theme(
      axis.text.x = element_text(size = 12, angle = 30, hjust = 1), 
      axis.text.y = element_text(size = 14), 
      axis.title = element_text(size = 23, face = "bold"), 
      strip.text.x = element_text(size = 12),
      legend.text = element_text(size = 13), 
      legend.title = element_text(size = 13), 
      panel.grid = element_blank()
      ) + 
    scale_fill_gradient2(midpoint = 0.5, low = "#FF7F00", mid = "white", high = "#1F78B4", limits = c(0, 1)) + 
    facet_grid(.~sim_sort) + 
    labs(x = expression(alpha[p]), y = expression(alpha[n]), fill = "Risk aversion",
         subtitle = sub_title)
  
}
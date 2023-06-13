# make dataframe longer (averaged across replications)

ConvertDfGeneLonger_MeanSD <- function(df) {
  
  df_rev <- 
    df %>% 
    dplyr::select(sim, replication, generation, mean_ap, mean_an, mean_bt, sd_ap, sd_an, sd_bt) %>% 
    dplyr::rename(
      MEAN_ap = mean_ap,
      MEAN_an = mean_an,
      MEAN_bt = mean_bt,
      SD_ap = sd_ap,
      SD_an = sd_an,
      SD_bt = sd_bt,
    ) %>%
    tidyr::pivot_longer(c(-sim, -generation, -replication), names_to = c("statistics", "para"), values_to = 'value', names_sep = "_") %>%
    tidyr::pivot_wider(names_from = statistics, values_from = value) %>% 
    dplyr::group_by(sim, generation, para) %>%
    dplyr::summarise(
      mean_rate = mean(MEAN), # Replicationに関して平均をとる
      sd_rate = mean(SD)) %>% 
    dplyr::ungroup() %>% # 以下、選択肢の効果量についての処理をする
    tidyr::separate(col = "sim", into = c("m2", "sd2", "m1", "sd1"), sep = "\\.", remove = FALSE) %>% 
    dplyr::mutate(
      m1 = as.integer(m1),
      sd1 = as.integer(sd1),
      m2 = as.integer(m2),
      sd2 = as.integer(sd2)
    ) %>% 
    dplyr::mutate(
      discriminability = EffectSize(m1, sd1, m2, sd2),
      task_exp = paste(m1, " vs ", m2, sep = "")
    ) %>% 
    dplyr::arrange(desc(discriminability))
  
  return(df_rev)
  
}

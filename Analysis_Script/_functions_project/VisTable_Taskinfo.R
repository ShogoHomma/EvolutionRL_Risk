# make task information of random multiple-task simulation by gt()

VisTable_Taskinfo <- function(df, cond) {
  
  df %>% 
    dplyr::filter(sim_n == cond) %>% 
    dplyr::select(-sim_n) %>% 
    gt() %>% 
    tab_header(
      title = cond,
    ) %>% 
    cols_label(sim_i = "simulation", 
               `1: nonrisky` = "non-risky", `1: risky` = "risky",
               `2: nonrisky` = "non-risky", `2: risky` = "risky",
               `3: nonrisky` = "non-risky", `3: risky` = "risky",
               `4: nonrisky` = "non-risky", `4: risky` = "risky"
    ) %>% 
    tab_spanner(label = "task 1", columns = c("1: nonrisky", "1: risky")) %>%
    tab_spanner(label = "task 2", columns = c("2: nonrisky", "2: risky")) %>% 
    tab_spanner(label = "task 3", columns = c("3: nonrisky", "3: risky")) %>%
    tab_spanner(label = "task 4", columns = c("4: nonrisky", "4: risky")) %>% 
    tab_options(
      table.width = pct(60),  # expand width of overall table
      table_body.hlines.width = 0, # not include a horizontal line
      column_labels.border.top.width = 2, 
      column_labels.border.top.color = "black", 
      column_labels.border.bottom.width = 1, 
      column_labels.border.bottom.color = "black",
      table.border.top.width = 2,
      table_body.border.top.color = "black", 
      table.border.bottom.width = 2,
      table_body.border.bottom.color = "black" 
    ) %>% 
    cols_align(align = "center") %>%
    tab_style(
      style = list(cell_text(align="center")), 
      locations = cells_row_groups()
    ) ->
    gt_tab
    
  return(gt_tab)
  
}
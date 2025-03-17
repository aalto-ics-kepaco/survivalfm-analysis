library(tidyverse)

setwd("~/survivalfm-analysis/results-visualization")

input_name =  "Standard risk factors"

color = 
  switch (input_name,
          "Standard risk factors" = "#3C5487",
          "Metabolomics biomarkers" = "#F49C7E",
          "Polygenic risk scores" = "#01A187",
          "Clinical biochemistry and blood counts" = "#E74B35"
  )


# Read results ------------------------------------------------------------

result_files <- list.files("../model-training/result-lists-datasize-experiment", full.names = T) 

result_files <- result_files[str_detect(result_files, str_replace_all(input_name, " ", "_"))]

df_results <- 
  purrr::map_dfr(
    .x = result_files,
    .f =  ~readr::read_csv(.x, col_types = readr::cols())
  ) %>% 
  dplyr::filter(input_name == !!input_name) %>% 
  dplyr::select(input_name, endpoint_name, rep, training_method, n_train, n_events_train, cindex_test) 


df_plot <- df_results %>% 
  arrange(training_method) %>% 
  dplyr::group_by(training_method, input_name, endpoint_name, n_train) %>% 
  dplyr::mutate(
    n_rep = n(),
    mean_cindex = mean(cindex_test)
  ) %>% 
  dplyr::filter(n_rep >= 5) 


# Plotting ----------------------------------------------------------------

plot <-
  ggplot(
    df_plot,
    aes(
      x = n_train,
      y = mean_cindex,
      color = training_method
    )
  ) +
  geom_line(size = 0.5) +
  geom_point(aes(x = n_train, y = cindex_test, color = training_method), size = 0.8, alpha = 0.5, shape = 1) +
  geom_point(aes(x = n_train, y = mean_cindex, color = training_method)) +
  lemon::facet_rep_wrap(~endpoint_name, scales = "free_y", repeat.tick.labels = T, ncol = 3) +
  scale_color_manual(values = c("darkgray", color)) +
  theme_classic() +
  labs(
    x = "\nNumber of individuals in training",
    y = "C-index"
  ) +
  scale_x_continuous(labels = scales::label_comma()) + 
  ggh4x::force_panelsizes(rows = unit(1.75, "in"), cols = unit(2.25, "in")) +
  guides(color = guide_legend(ncol = 1, reverse = T)) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    strip.background = element_rect(color = NA),
    strip.text = element_text(size = 11, face = "bold"),
    legend.key.height = unit(0.2, "cm")
  )



df_results %>% 
  readr::write_csv("source-data/Figure4.csv")


ggsave(
  paste0("figures-main/Figure4__performance_by_training_n__", snakecase::to_snake_case(input_name), ".pdf"),
  plot = plot,
  width = 9,
  height = 8
)

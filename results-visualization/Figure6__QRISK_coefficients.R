library(tidyverse)
library(RColorBrewer)

setwd("~/survivalfm-analysis/results-visualization")

purrr::walk(list.files("../model-training/utils", full.names = T), ~source(.x))


# Function to format variable names ---------------------------------------

rename_variables <- function(variable) {
  snakecase::to_sentence_case(variable) %>% 
    str_remove_all("`") %>% 
    str_replace_all("hdl", "HDL") %>% 
    str_replace_all("Bmi", "BMI") %>% 
    str_replace_all("status", "status:") %>% 
    str_replace_all("background", "background:") %>% 
    str_replace_all("chinese", "Chinese") %>% 
    str_replace_all("asian", "Asian") %>% 
    str_replace_all("pakistani", "Pakistani") %>% 
    str_replace_all("indian", "Indian") %>% 
    str_replace_all("bangladeshi", "Bangladeshi") %>% 
    str_replace_all("caribbean", "Caribbean") %>% 
    str_replace_all("Sex male", "Male sex") %>% 
    str_replace_all("Age 1", "(Age / 10)^-2") %>%
    str_replace_all("BMI 1", "(BMI / 10)^-2") %>%
    str_replace_all("BMI 2", "(BMI / 10)^-2 * log(BMI)")
}

# Result list ------------------------------------------------------------

result_files <- list.files("../model-training/result-lists-qrisk", full.names = T)

result_list  <- purrr::map_dfr(
  .x = result_files,
  .f =  ~readr::read_csv(.x, col_types = readr::cols())
) 

list_row <- result_list %>% 
  dplyr::filter(training_method == "survivalFM")


model_fit_path <- list_row$final_fit_path 

fit <- readRDS(paste0("../model-training/", model_fit_path))


# Plot linear part --------------------------------------------------------

df_linear <- as_tibble(fit$beta, rownames = "variable") %>% 
  dplyr::arrange(desc(value)) %>% 
  dplyr::mutate(
    variable = rename_variables(variable),
    variable = factor(variable, levels = unique(.data$variable))
  )

readr::write_csv(df_linear, "source-data/Figure6a.csv")

mat_linear <- df_linear %>% 
  tibble::column_to_rownames("variable") %>% 
  as.matrix()

plot_linear <- pheatmap::pheatmap(
  mat_linear, 
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  angle_col = 45,
  breaks =  seq(-max(abs(mat_linear)), max(abs(mat_linear)), length.out = 100),
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = F,
  treeheight_row = 0,
  treeheight_col = 0,
  cellwidth = unit(0.16, "in"),  
  cellheight = unit(0.16, "in"),   
  border_color = NA,
  symm = F
) 


# Plot interactions -------------------------------------------------------

df_interactions <-  as.tibble(fit$PP, rownames = "variable1") %>% 
  tidyr::gather(key = "variable2", value = "value", -variable1) %>% 
  dplyr::mutate(
    variable1 = rename_variables(variable1),
    variable2 = rename_variables(variable2)
  ) 

readr::write_csv(df_interactions, "source-data/Figure6b.csv")

mat_interactions <- df_interactions %>% 
  tidyr::spread(key = "variable2", value = "value") %>% 
  tibble::column_to_rownames("variable1") %>% 
  as.matrix()

plot_interactions <- 
  pheatmap::pheatmap(
    mat_interactions, 
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
    angle_col = 45,
    breaks =  seq(-max(abs(mat_interactions), na.rm = T), max(abs(mat_interactions), na.rm = T), length.out = 100),
    cluster_rows = T,
    cluster_cols = T,
    cellwidth = unit(0.16, "in"),  
    cellheight = unit(0.16, "in"),   
    border_color = NA,
    symm = F
  ) 


# Final combined plot -----------------------------------------------------

plot <- 
  cowplot::plot_grid(
    plotlist = list(
      plot_linear[[4]], 
      plot_interactions[[4]]
    ),
    rel_widths = c(1,2.5),
    labels = c("a", "b"),
    align = "h",
    axis = "r"
  )

ggsave(
  plot = plot,
  filename = paste0("figures-main/Figure6__QRISK_coefficients.pdf"),
  height = 9,
  width = 13
)

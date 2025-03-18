library(tidyverse)

setwd("~/survivalfm-analysis/results-visualization")

panel_width = unit(2.5, "in")
panel_height = unit(1.9, "in")


# Read results ------------------------------------------------------------


df_results <- readr::read_csv("../performance-metrics/performance_metrics__test__survivalfm_vs_cox.csv") %>% 
  select(input, endpoint, contains("delta"), contains("nri_overall"))

df_results %>% 
  filter(input == "Standard risk factors") %>% 
  readr::write_csv("source-data/Figure2a.csv")

df_results %>% 
  filter(input == "Clinical biochemistry and blood counts") %>% 
  readr::write_csv("source-data/Figure2b.csv")

df_results %>% 
  filter(input == "Metabolomics biomarkers") %>% 
  readr::write_csv("source-data/Figure2c.csv")

df_results %>% 
  filter(input == "Polygenic risk scores") %>% 
  readr::write_csv("source-data/Figure2d.csv")


# Plotting ----------------------------------------------------------------

input_order <- c("Standard risk factors",  "Clinical biochemistry and blood counts", "Metabolomics biomarkers", "Polygenic risk scores")

colors <-  c("#3C5488FF", "#E64B35FF",  "#F39B7FFF", "#00A087FF") %>% set_names(input_order)

cindex_limits <- range(c(df_results$delta_cindex_ci_lower, df_results$delta_cindex_ci_upper)) 
rsq_limits <-range(c(df_results$delta_rsq_ci_lower, df_results$delta_rsq_ci_upper)) 
cnri_limits <- c(min(0, min(df_results$nri_overall_ci_lower)), max(df_results$nri_overall_ci_upper))

ls_plots <- 
  purrr::map(
    .x = input_order,
    .f = function(input) {
      
      df_plot <- df_results %>% filter(input == !!input)
      
      plot_a <- 
        ggplot(df_plot,  aes(y = endpoint, x = delta_cindex)) +
        geom_col(width = 0.75, fill = colors[[input]]) +
        geom_vline(xintercept = 0, size = 0.6) +
        scale_x_continuous(limits = cindex_limits) +
        scale_y_discrete(limits = rev) +
        geom_errorbar(aes(xmin = delta_cindex_ci_lower, xmax = delta_cindex_ci_upper), width = 0.2, color = colorspace::darken(colors[[input]], amount = 0.5)) +
        labs(
          x = "\u0394C-index (95% CI)",
          y = ""
        ) +
        geom_vline(xintercept = 0, color = "black") +
        ggh4x::force_panelsizes(cols = panel_width, rows = panel_height) +
        theme_classic() +
        theme(
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          axis.text.y = element_text(size = 11, color = "black"),
          axis.ticks.y = element_blank()
        )
      
      
      plot_b <- 
        ggplot(df_plot,  aes(y = endpoint, x = delta_rsq)) +
        geom_col(width = 0.75, fill = colors[[input]]) +
        geom_vline(xintercept = 0, size = 0.6) +
        scale_x_continuous(limits = rsq_limits, labels = scales::percent) +
        scale_y_discrete(limits = rev) +
        geom_errorbar(aes(xmin = delta_rsq_ci_lower, xmax = delta_rsq_ci_upper), width = 0.2, color = colorspace::darken(colors[[input]], amount = 0.5)) +
        labs(
          x = "\u0394R\u00B2 (95% CI)",
          y = ""
        ) +
        geom_vline(xintercept = 0, color = "black") +
        ggh4x::force_panelsizes(cols = panel_width, rows = panel_height) +
        theme_classic() +
        theme(
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      
      
      plot_c <- 
        ggplot(df_plot,  aes(y = endpoint, x = nri_overall)) +
        geom_col(width = 0.75, fill = colors[[input]]) +
        geom_vline(xintercept = 0, size = 0.6) +
        geom_errorbar(aes(xmin = nri_overall_ci_lower, xmax = nri_overall_ci_upper), width = 0.2, color = colorspace::darken(colors[[input]], amount = 0.5)) +
        scale_x_continuous(limits = cnri_limits, expand = c(0.0005,0.0005), labels = scales::percent) +
        scale_y_discrete(limits = rev) +
        labs(
          x = "Continuous NRI (95% CI)",
          y = ""
        ) +
        geom_vline(xintercept = 0, color = "black") +
        ggh4x::force_panelsizes(cols = panel_width, rows = panel_height) +
        theme_classic() +
        theme(
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
      
      gridExtra::grid.arrange(
        grobs = list(plot_a, plot_b, plot_c),
        ncol = 3,
        top = grid::textGrob(paste0(input), gp = grid::gpar(fontsize = 11, fontface = "bold"), x = 0.5833333, hjust = 0.5),
        widths = c(1.5,1,1)
      )
    }
  )

final_plot <- 
  cowplot::plot_grid(
    plotlist = ls_plots,
    labels = c("a", "b", "c", "d"),
    ncol = 1
  )

final_plot

ggsave(
  filename = "figures-main/Figure2__Cindex_Rsq_NRI.pdf",
  plot = final_plot,
  device = cairo_pdf,
  w = 10.2,
  h = 10.7
)

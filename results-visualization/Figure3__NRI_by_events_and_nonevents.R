library(tidyverse)

setwd("~/survivalfm-analysis/results-visualization")

# Read results ------------------------------------------------------------

df_results <- readr::read_csv("../performance-metrics/performance_metrics__test__survivalfm_vs_cox.csv") %>% 
  select(input, endpoint, contains("nri"))

readr::write_csv(df_results, "source-data/Figure3.csv")


# Plotting ----------------------------------------------------------------

df_plot <- df_results %>%  
  dplyr::select(input, endpoint, contains("nri")) %>% 
  tidyr::gather(key = "name", value = "value", contains("nri")) %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      str_detect(name, "nri_overall") ~ "Overall",
      str_detect(name, "nri_events") ~ "Events",
      str_detect(name, "nri_nonevents") ~ "Non-events"
    ),
    metric = str_remove_all(name, "nri_events|nri_overall|nri_nonevents"),
    metric = dplyr::case_when(
      metric == "" ~ "estimate",
      metric == "_ci_lower" ~ "ci_lower",
      metric == "_ci_upper" ~ "ci_upper",
    ),
    group = factor(group, levels = rev(c("Overall", "Events", "Non-events"))),
    input = factor(input, levels = c("Standard risk factors",  "Clinical biochemistry and blood counts", "Metabolomics biomarkers", "Polygenic risk scores"))
  ) %>% 
  select(-name) %>% 
  tidyr::spread(key = "metric", value = "value")


df_plot %>% 
  filter(input == "Standard risk factors") %>% 
  readr::write_csv("source-data/Figure3a.csv")

df_plot %>% 
  filter(input == "Clinical biochemistry and blood counts") %>% 
  readr::write_csv("source-data/Figure3b.csv")

df_plot %>% 
  filter(input == "Metabolomics biomarkers") %>% 
  readr::write_csv("source-data/Figure3c.csv")

df_plot %>% 
  filter(input == "Polygenic risk scores") %>% 
  readr::write_csv("source-data/Figure3d.csv")


limits <- range(c(df_plot$ci_lower, df_plot$ci_upper))

ls_plots <- 
  purrr::map(
    .x = unique(df_plot$input),
    .f = function(input) {
      
      ggplot(
        df_plot %>% filter(input == !!input), 
        aes(y = endpoint, x = estimate, fill = group)
      ) +
        geom_col(position = position_dodge(width = 0.6), width = 0.6) +
        geom_errorbar(
          aes(xmin = ci_lower, xmax = ci_upper, color = group),
          position = position_dodge(width = 0.6), 
          width = 0.4
        ) +
        lemon::facet_rep_wrap(~input, ncol = 2, repeat.tick.labels = TRUE, scales = "free_y") +
        labs(
          x = "Continuous NRI (95% CI)",
          y = ""
        ) +
        scale_fill_manual(values = c("#00BFC4", "#F8766D", "gray40")) +
        scale_color_manual(values = colorspace::darken(c("#00BFC4", "#F8766D", "gray40"), amount = 0.5)) +
        geom_vline(xintercept = 0, color = "black") +
        scale_x_continuous(limits = limits, labels = scales::percent) +
        scale_y_discrete(limits = rev) +
        coord_cartesian(clip = 'off') +
        ggh4x::force_panelsizes(cols = unit(3.2, "in"), rows = unit(2.4, "in"),) +
        theme_classic() +
        theme(
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0, size = 11),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 11),
          axis.text.y = element_text(size = 11, color = "black")
        )
    }
  )


# Final plot --------------------------------------------------------------

plots <- cowplot::plot_grid(
  plotlist = lapply(ls_plots, function(x) x + theme(legend.position = "none")),
  ncol = 2,
  labels = c("a", "b", "c", "d")
)

legend <-  
  cowplot::plot_grid(
    plotlist = list(cowplot::get_legend(ls_plots[[1]])), 
    ncol = 4
  ) 

final_plot <- cowplot::plot_grid(
  plotlist = list(legend, plots),
  ncol = 1,
  rel_heights = c(1,10)
)

final_plot

ggsave(
  filename = "figures-main/Figure3__NRI_by_events_and_nonevents.pdf",
  plot = final_plot,
  device = cairo_pdf,
  w = 10.5,
  h = 7
)

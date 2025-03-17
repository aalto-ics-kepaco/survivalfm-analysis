library(tidyverse)
library(survminer)
library(survival)

setwd("~/survivalfm-analysis/results-visualization")

# Result list ------------------------------------------------------------

result_files <- list.files("../model-training/result-lists-qrisk", full.names = T)

result_list  <- purrr::map_dfr(
  .x = result_files,
  .f =  ~readr::read_csv(.x, col_types = readr::cols())
) %>% 
  arrange(idx) %>% 
  dplyr::mutate(
    display_name =  dplyr::recode(
      display_name,
      "Linear terms (standard Cox)" = "Cox regression",
      "Linear terms + QRISK3 interaction terms with age (standard Cox)" = "Cox regression with age interactions",
      "Linear terms + factorized interaction terms (survivalFM)" = "survivalFM" 
    ),
    display_name = factor(display_name, levels = unique(.data$display_name))
  )

# Predictions -------------------------------------------------------------

df_pred <- 
  purrr::map_dfr(
    .x = 1:nrow(result_list),
    .f = function(idx) {
      
      message(idx)
      
      list_row <- result_list %>% dplyr::slice(!!idx)
      pred_path <- list_row$predictions_path %>% str_remove("/m/cs/scratch/ukbiobank-kepaco/tmp-heli/")
      
      df <- readRDS(paste0("../model-training/", pred_path)) %>%
        dplyr::cross_join(list_row %>% select(idx, display_name, endpoint_name, input_name))
      
      return(df)
    }
  )

# Reclassification plots --------------------------------------------------

df_plot <- df_pred %>% 
  dplyr::select(eid, abs_risk, status, times, display_name) %>% 
  tidyr::spread(key = "display_name", value = "abs_risk") %>% 
  dplyr::mutate(
    status = dplyr::case_when(
      status == 1 ~ "Event",
      status == 0 & times == 10 ~ "Non-event",
      status == 0 & times < 10 ~ "Lost to follow-up"
    ),
    status = factor(status, levels = c("Event", "Lost to follow-up", "Non-event"))
  ) %>% 
  arrange(desc(status)) 


ls_plots_reclassification <- 
  purrr::map2(
    .x = c("Cox regression", "Cox regression"),
    .y = c( "Cox regression with age interactions", "survivalFM"),
    .f = function(std_model, new_model) {
      
      ggplot(
        df_plot,
        aes_string(x = paste0("`", std_model, "`"), y = paste0("`", new_model, "`"), color = "status", fill = "status", shape = "status")
      )  +
        ggrastr::geom_point_rast(size = 0.5, stroke = 0.5) + 
        geom_blank(aes_string(y = paste0("`", std_model, "`"), x = paste0("`", new_model, "`"))) + 
        geom_abline(color = "black", linetype = "dashed") +
        geom_hline(yintercept = 0.1, linetype = "dotted", color = "black") +
        geom_vline(xintercept = 0.1, linetype = "dotted", color = "black") +
        scale_x_log10(labels = scales::percent) +
        scale_y_log10(labels = scales::percent) +
        labs(
          x = paste0("Predicted risk for CVD (%),\n", std_model),
          y = paste0("Predicted risk for CVD (%),\n", new_model)
        ) +
        scale_color_manual(values = c("#E64B35FF", "black", "#00A087FF")) +
        scale_shape_manual(values = c(20,4,20)) +
        ggtitle("") +
        theme_classic() +
        guides(color = guide_legend(ncol = 1, override.aes = list(size=2))) +
        theme(
          aspect.ratio = 1,
          legend.position = "inside", 
          legend.justification = c("left", "top"),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.4, "cm"),
          legend.key.spacing.y = unit(0, "cm"),
          legend.title = element_blank(),
          legend.background = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5)
        )
    }
  )


# Kaplan-Meier analysis of high-risk group --------------------------------

df_high_risk <- df_pred %>%
  dplyr::filter(abs_risk >= 0.1) %>% 
  dplyr::mutate(
    display_name = str_replace_all(display_name, "\\(", "\\\n("),
    group = factor(display_name, levels = unique(.data$display_name))
  )

km_fit_high_risk <- survfit(Surv(times, status) ~ factor(group), data = df_high_risk)
names(km_fit_high_risk$strata) <- gsub("factor\\(group\\)\\=", "", names(km_fit_high_risk$strata))

km_high_risk <- 
  ggsurvplot(
    km_fit_high_risk, 
    data = df_high_risk,
    censor = F,
    palette = c("grey40", "#4DBBD5FF", "#3C5488FF"),
    ncensor.plot = T,
    conf.int = FALSE,           
    fun = "event",
    risk.table = "nrisk_cumevents"
  )

plot_high_risk <- km_high_risk$plot +
  ylab("Cumulative incidence (%)") +
  ggtitle("\nPredicted risk ≥ 10%") +
  scale_x_continuous(expand = c(0.15,0.15), breaks = c(0,5,10)) +
  scale_y_continuous(limits = c(0, 0.185), labels = scales::percent) +
  guides(color = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(
    legend.position = "inside", 
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.key.spacing.y = unit(0, "cm"),
    legend.background = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

names(km_fit_high_risk$strata) <- gsub("with ", "with\n", names(km_fit_high_risk$strata))

table_high_risk <- 
  survminer::ggsurvtable(
    km_fit_high_risk, 
    data = df_high_risk, 
    risk.table.type = "nrisk_cumevents",
    fontsize = 3.2,
    break.time.by = 5
  ) $risk.table +
  scale_x_continuous(expand = c(0.15,0.15), breaks = c(0,5,10)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 11)
  )

lapply(list(plot_high_risk),  egg::set_panel_size, width = unit(3, "in"), height = unit(2, "in"))

plot_km_high_risk <- egg::ggarrange(
  plot_high_risk, 
  table_high_risk, 
  widths = c(2,2),
  heights = c(2,1),
  nrow = 2
)

km_high_risk$data.survplot %>% 
  readr::write_csv("source-data/Figure5c.csv")


# Kaplan-Meier analysis of low-risk group ---------------------------------

df_low_risk <- df_pred %>%
  dplyr::filter(abs_risk < 0.1) %>% 
  dplyr::mutate(group = factor(display_name, levels = unique(.data$display_name)))

km_fit_low_risk <- survfit(Surv(times, status) ~ factor(group), data = df_low_risk)
names(km_fit_low_risk$strata) <- gsub("factor\\(group\\)\\=", "", names(km_fit_low_risk$strata))

km_low_risk <- 
  ggsurvplot(
    km_fit_low_risk, 
    data = df_low_risk,
    censor = F,
    palette = c("grey40", "#4DBBD5FF", "#3C5488FF"),
    ncensor.plot = T,
    conf.int = FALSE,           
    fun = "event",
    risk.table = "nrisk_cumevents"
  )

plot_low_risk <- km_low_risk$plot +
  ylab("Cumulative incidence (%)") +
  ggtitle("\nPredicted risk < 10%") +
  scale_x_continuous(expand = c(0.15,0.15), breaks = c(0,5,10)) +
  scale_y_continuous(limits = c(0, 0.185), labels = scales::percent) +
  guides(color = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(
    legend.position = "inside", 
    legend.justification = c("left", "top"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.key.spacing.y = unit(0, "cm"),
    legend.background = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


names(km_fit_low_risk$strata) <- gsub("with ", "with\n", names(km_fit_low_risk$strata))

table_low_risk <- 
  survminer::ggsurvtable(
    km_fit_low_risk, 
    data = df_low_risk, 
    risk.table.type = "nrisk_cumevents",
    fontsize = 3.2,
    break.time.by = 5
  ) $risk.table +
  scale_x_continuous(expand = c(0.15,0.15), breaks = c(0,5,10)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 11)
  )

lapply(list(plot_low_risk),  egg::set_panel_size, width = unit(3, "in"), height = unit(2, "in"))

plot_km_low_risk <- egg::ggarrange(
  plot_low_risk, 
  table_low_risk, 
  widths = c(2,2),
  heights = c(2,1),
  nrow = 2
)


km_low_risk$data.survplot %>% 
  readr::write_csv("source-data/Figure5d.csv")


# Final combined plot -----------------------------------------------------

final_plot <- egg::ggarrange(
  ls_plots_reclassification[[1]],
  ls_plots_reclassification[[2]],
  plot_high_risk, plot_low_risk, 
  table_high_risk, table_low_risk,
  ncol = 2,
  heights = c(1,1,0.4),
  labels = c("a", "b", "c", "d", "", ""),
  label.args = list(x = 0.1, y = 0.99, gp = grid::gpar(fontface = "bold", fontsize = 14))
)

ggsave(
  filename = "figures-main/Figure5__QRISK_analysis.pdf",
  plot = final_plot,
  device = cairo_pdf,
  w = 9,
  h = 9
)


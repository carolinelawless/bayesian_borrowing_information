remove(list = ls())
setwd("~/Documents/travail/dmd_when_next_clinical_trial?/bayesian_borrowing_information")
library(ggplot2)
library(dplyr)
library(patchwork)

source("results.R")



ggplot(
  all_data_binomial,
  aes(x = lambda, y = probability,
      color = method, linetype = method)
) +
  
  geom_line(linewidth = 0.5) +
  
  scale_color_manual(
    values = c(
      "ε = 0 (no borrowing)" = "black",
      "ε = 0.5 without EB" = "blue",
      "ε = 1 without EB (full borrowing)" = "blue",
      "ε = 0.5 with EB" = "red",
      "ε = 1 with EB" = "red",
      "Liang (2023) EB approach" = "purple"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "ε = 0 (no borrowing)" = "solid",
      "ε = 0.5 without EB" = "solid",
      "ε = 1 without EB (full borrowing)" = "dashed",
      "ε = 0.5 with EB" = "solid",
      "ε = 1 with EB" = "dashed",
      "Liang (2023) EB approach" = "solid"
    )
  ) +
  
  scale_y_continuous(limits = c(0, 1)) +
  
  facet_grid(
    n_versions ~ scenario
  ) +
  
  labs(
    x = expression(lambda),
    y = "Power/Type-1 error",
    color = NULL,
    linetype = NULL
  ) +
  
  theme_classic() +
  theme(
    legend.position = "bottom",
    panel.spacing = unit(0.5, "lines")
  )


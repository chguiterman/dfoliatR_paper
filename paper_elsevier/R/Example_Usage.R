## Create tables and graphics for dfoliatR manuscript

devtools::install_github("chguiterman/dfoliatR@dev", force=TRUE)
library(dfoliatR)
library(ggplot2)
library(patchwork)

out_pth <- "paper_elsevier/Output/" 
# Tree-level outputs ------------------------------------------------------

data("dmj_h")
data("dmj_nh")

dmj_defol <- defoliate_trees(host_tree = dmj_h,
                             nonhost_chron = dmj_nh,
                             duration_years = 8,
                             max_reduction = -1.28,
                             bridge_events = TRUE,
                             series_end_event = TRUE,
                             list_output = FALSE)

dmj_stats <- defol_stats(dmj_defol)
write.csv(dmj_stats, paste0(out_pth, "defol-stats.csv"), row.names = FALSE)

plot_defol(dmj_defol)
ggsave(paste0(out_pth, "tree-plot-default.pdf"), width=5.75, dpi=300, units="in")

p1 <- plot_defol(dmj_defol) +
  theme(legend.position = c(0.5, -0.15),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.text = element_text(size=8),
        plot.margin = unit(c(0, 0, 15, 0), units = "points"),
        plot.background = element_blank())
  
p2 <- p1 + 
  scale_color_manual(values = c("red", "orange", "purple")) 

p1 / p2 + plot_layout(guides = "keep") + plot_annotation(tag_levels = "A")
ggsave(paste0(out_pth, "tree-plot-2panel.pdf"), width=5.75, dpi=300, units="in")

# Site-level outputs ------------------------------------------------------

dmj_obr <- outbreak(dmj_defol,
                    filter_min_series = 3,
                    filter_min_defol = 1,
                    filter_perc = 25)

plot_outbreak(dmj_obr, disp_index = "GSI")
ggsave(paste0(out_pth, "site-plot.pdf"), width=5.75, dpi=300, units="in")

dmj_obr_stats <- outbreak_stats(dmj_obr)
write.csv(dmj_obr_stats, paste0(out_pth, "obr-stats.csv"), row.names = FALSE)

## Create tables and graphics for dfoliatR manuscript

# devtools::install_github("chguiterman/dfoliatR", force=TRUE)
library(dfoliatR)
library(ggplot2)
library(here)

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
write.csv(dmj_stats, here("Output", "defol-stats.csv"), row.names = FALSE)

plot_defol(dmj_defol)
ggsave(here("Output", "tree-plot-default.pdf"), width=5.75, dpi=300, units="in")
ggsave(here("Output", "tree-plot-default.tiff"), width=5.75, dpi=300, units="in")
ggsave(here("Output", "tree-plot-default.eps"), width=5.75, dpi=300, units="in")

plot_defol(dmj_defol) +
  scale_color_manual(values = c("red", "orange", "purple")) 

# Site-level outputs ------------------------------------------------------

dmj_obr <- outbreak(dmj_defol,
                    filter_min_series = 3,
                    filter_min_defol = 1,
                    filter_perc = 25)

plot_outbreak(dmj_obr, disp_index = "GSI")
ggsave(here("Output", "site-plot.pdf"), width=5.75, dpi=300, units="in")
ggsave(here("Output", "site-plot.tiff"), width=5.75, dpi=300, units="in")

dmj_obr_stats <- outbreak_stats(dmj_obr)
write.csv(dmj_obr_stats, here("Output", "obr-stats.csv"), row.names = FALSE)

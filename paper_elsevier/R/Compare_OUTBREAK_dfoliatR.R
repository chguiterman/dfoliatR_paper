# Compare results from OUTBREAK to dfoliatR

library(dfoliatR)
library(dplR)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(cowplot)
# Load read_OUTBREAK function
devtools::source_gist("https://gist.github.com/chguiterman/c88e9601dd085eb4215f0dc133a0e21c", sha1 = "fdb3cfefcd7e6ed3f7df3e20057c815bbd91dd0c")

data_pth <- "paper_elsevier/Data/"

## Master metadata file
meta_in <- read.csv(paste0(data_pth, "Meta.csv"))

## Helper function
rm_sample_depth <- function(x){
  if (any(names(x) %in% "samp.depth"))
    df <- select(x, - samp.depth)
  else df <- x
  return(df)
}


# OUTBREAK data -----------------------------------------------------------

outbr_results <- meta_in %>% 
  transmute(site, 
            outbr_cor = map(outbr_corrected, ~ read.compact(paste0(data_pth, .))),
            outbr_col = map(outbr_site, ~ read.delim(paste0(data_pth, .))),
            outbr_out = map2(outbr_file, insect_code, 
                             ~ read_OUTREAK(fname = paste0(data_pth, .x),
                                            level = "tree",
                                            insect_code = .y)))

# Tree-level table
outbr_defol <- outbr_results %>% 
  transmute(site,
            outbr_cor = map(outbr_cor, ~ .[] %>%  
                              select(! ends_with("NormCorr")) %>%
                              rownames_to_column("year") %>% 
                              pivot_longer(-year, names_to = "series", values_to = "ncor"))
  ) %>% 
  unnest(cols = outbr_cor) %>% 
  na.omit() %>% 
  arrange(series, year)

outbr_out <- outbr_results %>% 
  select(site, outbr_out) %>% 
  unnest(cols = outbr_out)

outbr_defol <- merge(outbr_defol, outbr_out, by = c("site", "series", "year"))

# Site-level table
outbr_obr <- outbr_results %>% 
  transmute(site,
            outbr_col = map(outbr_col, ~ .[] %>% 
                              select(Year, 
                                     contains("Trees", ignore.case = TRUE), 
                                     starts_with("WSBW") | starts_with("BUDWORM")) %>% 
              `names<-`(c("year", "samp_depth", "num_defol", "perc_defol", "num_max_defol")))
            ) %>% 
  unnest(cols = "outbr_col") %>% 
  mutate(software = "OUTBREAK")


# dfoliatR data -----------------------------------------------------------

df_run <- meta_in %>% 
  transmute(site,
            host = map(host_file, ~ read.compact(paste0(data_pth, .))),
            nonhost = map(nonhost_file, ~ read.crn(paste0(data_pth, .)) %>% 
                            rm_sample_depth()),
            defol = pmap(list(host, nonhost, ongoing_defol),
                         ~ defoliate_trees(host_tree = ..1, 
                                           nonhost_chron = ..2,
                                           series_end_event = ..3,
                                           duration_years = 8,
                                           max_reduction = -1.28,
                                           bridge_events = FALSE,
                                           list_output = FALSE))
  )

df_defol <- df_run %>% 
  select(site, defol) %>% 
  unnest(cols = "defol") 

## Run dfoliatR::outbreak()
df_obr <- df_run %>% 
  transmute(site,
            obr = map(defol, outbreak, filter_perc = 25, filter_min_series = 3, filter_min_defol = 1)) %>% 
  unnest(cols = "obr") %>% 
  mutate(software = "dfoliatR")

# Compare tree-level results ----------------------------------------------

tree_data <- merge(outbr_defol, df_defol, by = c("site", "series", "year"), all = TRUE)

## 1:1 plot of indices
ggplot(tree_data, aes(x = ncor, y = ngsi, group = series)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(color = site), shape = 1, size = 2, alpha = .5) +
  facet_wrap(~ site, nrow=2) +
  # coord_cartesian(xlim = c(-.1, .1), ylim = c(-.1, .1)) +
  coord_equal() +
  xlab("OUTBREAK NCOR") +
  ylab("dfoliatR NGSI") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("paper_elsevier/Output/plot-trees-df-obr.pdf", width=5.75, dpi=300, units="in")  
  
## Compare events
events_list <- levels(dfoliatR::dmj_defol$defol_status)[2:5] 

tree_comps <- tree_data %>% 
  select(site, series, year, ncor, ngsi, obr_type, defol_status) %>% 
  mutate(comp_code = case_when((obr_type %in% events_list) & (defol_status %in% events_list) ~ TRUE,
                               (obr_type %in% events_list) | (defol_status %in% events_list) ~ FALSE,
                               TRUE ~ NA)
  ) 

comps_same <- tree_comps %>% 
  filter(comp_code)
comps_diff <- tree_comps %>% 
  filter(! comp_code)

## Percent agreement = 92%
nrow(comps_same) / (nrow(comps_same) + nrow(comps_diff)) * 100


## Difference table
diff_tbl <- comps_diff %>% 
  group_by(series) %>% 
  summarize(period = paste(min(year), ":", max(year)),
            duration = as.numeric(max(year)) - as.numeric(min(year)) + 1,
            df_defol = any(defol_status %in% events_list),
            obr_defol = any(obr_type %in% events_list)
  )
write.csv(diff_tbl, "paper_elsevier/Output/OBR-DF_DiffTbl.csv")
write.csv(tree_comps, "paper_elsevier/Output/Tree-comps-data.csv")

# Compare site-level results ----------------------------------------------

site_data <- merge(outbr_obr, df_obr, by = c("site", "year"), all = TRUE)

site_perc <- outbr_obr %>% 
  select(site, year, perc_defol, software) %>% 
  rbind(df_obr %>% 
          select(site, year, perc_defol, software))

ggplot(site_perc, aes(x = year, y = perc_defol)) +
  geom_line(aes(color = software)) +
  ylab("Percent trees defoliated") +
  xlab("Year") +
  facet_wrap(~site, nrow = 4) +
  # labs(tag = c("OUTBREAK", "dfoliatR")) +
  # # coord_cartesian(clip = "off") +
  # scale_color_identity(guide = "legend", name = "Software",
  #                      # values = c("orange", "blue"),
  #                      labels = c("OUTBREAK", "dfoliatR")) +
  theme_bw() +
  # guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = c(0.5, 1.1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.margin = unit(c(30, 5, 5, 5), units = "points"),
        plot.background = element_blank())
ggsave("paper_elsevier/Output/plot-sites-df-obr.pdf", width=5.75, dpi=300, units='in')  



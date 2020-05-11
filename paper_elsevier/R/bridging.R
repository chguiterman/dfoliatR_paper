# Show bridging
## CHG 2020-05-07

library(dfoliatR)
library(dplR)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(prismatic)
library(paletteer)
library(here)

## Master metadata file
meta_in <- read.csv(here("paper_elsevier", "Data", "Meta.csv"))
# file of output from the comparison -- has the OBR results
obr_results <- read.csv(here("paper_elsevier", "Output", "Tree-comps-data.csv")) %>% 
  select(site, year, series, obr_type)

## Helper function
rm_sample_depth <- function(x){
  if (any(names(x) %in% "samp.depth"))
    df <- select(x, - samp.depth)
  else df <- x
  return(df)
}

df_run <- meta_in %>% 
  transmute(site,
            host = map(host_file, ~ read.compact(here("paper_elsevier", "Data", .))),
            nonhost = map(nonhost_file, ~ read.crn(here("paper_elsevier", "Data", .)) %>% 
                            rm_sample_depth()),
            defol = pmap(list(host, nonhost, ongoing_defol),
                         ~ defoliate_trees(host_tree = ..1, 
                                           nonhost_chron = ..2,
                                           series_end_event = ..3,
                                           duration_years = 8,
                                           max_reduction = -1.28,
                                           bridge_events = TRUE,  ## Add bridges
                                           list_output = FALSE))
  )

df_defol <- df_run %>% 
  select(site, defol) %>% 
  unnest(cols = "defol") %>% 
  inner_join(obr_results)


# Selected bridge example -------------------------------------------------

fnl04 <- df_defol %>% 
  filter(series == "FNL04",
         year %in% 1881:1925)
efk22 <- df_defol %>% 
  filter(series == "EFK22",
         year %in% 1940:1975)
dmj26 <- df_defol %>% 
  filter(series == "DMJ26",
        year %in% 1850:1875)

comps <- rbind(fnl04, efk22, dmj26)

obr_events <- comps %>% 
  filter(obr_type != "nd")
df_events <- comps %>% 
  filter(defol_status != "nd")

cols <- dichromat::colorschemes$SteppedSequential.5[c(3, 8, 18)]

ggplot(comps, aes(x=year)) +
  # scale_x_continuous(limits=c(1881, 1950)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_line(aes(y=ngsi), size=1.25) +
  facet_wrap(~ series, scales = "free", ncol = 1) +
  geom_segment(data = df_events,
               aes(x=year, xend=year,
                   y=0, yend=ngsi,
                   color = stage(start = defol_status,
                                 after_scale = clr_lighten(color, 0.25, space="combined"))),
               size=2.5) +
  geom_segment(data = obr_events,
               aes(x=year, xend=year,
                   y=0, yend=ngsi,
                   color = stage(start = obr_type,
                                 after_scale = clr_darken(color, 0.2))
                   ),
               size=2.5) +
  ylab("NGSI") + xlab("Year") +
  scale_color_paletteer_d("awtools::mpalette", direction = 1) +
  # scale_color_manual(values = cols) +
  ggpubr::theme_pubr() +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold", hjust = 0, size=16),
        axis.title = element_text(face="bold", size = 14))
ggsave(here("paper_elsevier", "Output", "bridging.pdf"), width = 5.75)

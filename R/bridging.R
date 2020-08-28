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
meta_in <- read.csv(here("Data", "Meta.csv"))
# file of output from the comparison -- has the OBR results
obr_results <- read.csv(here("Output", "Tree-comps-data.csv")) %>% 
  select(site, year, series, obr_type)

## Helper function
rm_sample_depth <- function(x){
  if (any(names(x) %in% "samp.depth"))
    df <- select(x, - samp.depth)
  else df <- x
  return(df)
}


# Perform reconstructions -------------------------------------------------

df_run <- meta_in %>% 
  transmute(site,
            host = map(host_file, ~ read.compact(here("Data", .))),
            nonhost = map(nonhost_file, ~ read.crn(here("Data", .)) %>% 
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

# fnl04 <- df_defol %>% 
#   filter(series == "FNL04",
#          year %in% 1881:1925)
efk22 <- df_defol %>% 
  filter(series == "EFK22",
         year > 1850)
dmj26 <- df_defol %>% 
  filter(series == "DMJ26",
        year %in% 1750:1950)

# comps <- rbind(fnl04, efk22, dmj26)
comps <- rbind(efk22, dmj26)

obr_events <- comps %>% 
  mutate(ngsi = replace(ngsi, obr_type == "nd", 0))
df_events <- comps %>% 
  mutate(ngsi = replace(ngsi, defol_status == "nd", 0))

# cols <- dichromat::colorschemes$SteppedSequential.5[c(3, 8, 18)]

## COLOR
ggplot(comps, aes(x=year)) +
  # scale_x_continuous(limits=c(1881, 1950)) +
  geom_hline(yintercept = 0, color = "grey80") +
  facet_wrap(~ series, scales = "free_x", ncol = 1, strip.position = "right") +
  geom_area(data=df_events, aes(x = year, y = ngsi, fill = "dfoliatR")) +
  geom_area(data=obr_events, aes(x = year, y = ngsi, fill = "OUTBREAK & dfoliatR")) +
  geom_line(aes(y=ngsi)) +
  scale_fill_manual(values = c("grey60", "#1C86EE"), 
                    breaks = c("OUTBREAK & dfoliatR", "dfoliatR")) +
  guides(fill = guide_legend(title = "Identified defoliation events: ")) +
  ylab("NGSI") + xlab("Year") +
  ggpubr::theme_pubr() +
  theme(legend.position = "bottom",
        strip.text = element_text(face="bold", size=12),
        axis.title = element_text(face="bold", size = 12))
ggsave(here("Output", "bridging_color.pdf"), width = 5.75)
ggsave(here("Output", "bridging_color.tiff"), width = 5.75)

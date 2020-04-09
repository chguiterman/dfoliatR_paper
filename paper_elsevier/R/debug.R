
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
                                           list_output = TRUE))
  )

fnl <- df_run$defol[[4]]           
input_series <- fnl[[11]] 

dmj <- df_run$defol[[1]]
input_series <- dmj[[12]]

fcc <- df_run$defol[[3]]
input_series <- fcc[[28]]

## testing whether there is a significant difference between n_yr = 1 and 
# n_yr = max(n_yr) for each species

pacman::p_load(data.table, stringr, dplyr, ggplot2, sf)

# folder <- "data/analysis/n_effects_its/"  ## all 10+ n years
folder <- "data/analysis/n_effects_its_3yrs/" ## only top 3 ss years

# n_its <- "i50"
n_its <- "i100"

files     <- list.files(folder, full.names = T, pattern = n_its)
filenames <- tools::file_path_sans_ext(list.files(folder, pattern = n_its))

results_list <- list()

for(i in seq_along(filenames)) {
  print(i)
  
  asp   <- do.call(rbind, str_split(filenames, pattern = "_"))[i,1]
  asite <- do.call(rbind, str_split(filenames, pattern = "_"))[i,2]
  
  its_df <- readRDS(files[i])
  
  maxssize <- max(its_df$n_trx)
  
  # filter to only maximum sample size (n_trx)
  its_df_f <- subset(its_df, n_trx == maxssize) 
  
  # plot(its_df_f$n_yrs, its_df_f$BAval)
  
  # filter to only n_yr of 1 and the max (or 3) for comparison
  max_nyrs <- max(its_df$n_yr)
  
  # its_df_f2 <- subset(its_df_f, its_df_f$n_yr == 1 | its_df_f$n_yr == max_nyrs)
  # its_df_f2 <- subset(its_df_f, its_df_f$n_yr == 1 | its_df_f$n_yr == 3)
  its_df_f2 <- its_df_f
  its_df_f2$n_yrs <- as.factor(its_df_f2$n_yrs)
  
  # png(paste0("figures/n_effects/max_ssize/", asp, "_", asite, "_", "n", maxssize, "_", "overlap", ".png"), width = 580, height = 580)
  # plot(its_df_f2$n_yrs, its_df_f2$BAval, 
  #      ylab="Overlap (BA)", xlab="N years",
  #      main = paste(asp, "n =", maxssize))
  # dev.off()
  
  m1 <- aov(BAval ~ n_yrs, data=its_df_f2)
  aov1 <- anova(m1)
  ph_m1 <- as.data.frame(TukeyHSD(m1)$n_yrs)
  ph_m1 <- ph_m1[str_detect(rownames(ph_m1), pattern = "1"), ]
  stest_m1 <- any(ph_m1$`p adj`<0.05)
  # 
  # ## Area 95% 
  # png(
  #   paste0("figures/n_effects/max_ssize/", asp, "_", asite, "_", "n", maxssize, "_", "area95", ".png"),
  #   width = 580, height = 580)
  # plot(its_df_f2$n_yrs, its_df_f2$area95,
  #      ylab="Area (sq. km)", xlab="N years",
  #      main = paste(asp, "n =", maxssize))
  # dev.off()
  # 
  # m2 <- aov(area95 ~  n_yrs, data=its_df_f2)
  # aov2 <- anova(m2)
  # ph_m2 <- as.data.frame(TukeyHSD(m2)$n_yrs)
  # ph_m2 <- ph_m2[str_detect(rownames(ph_m2), pattern = "1"), ]
  # stest_m2 <- any(ph_m2$`p adj`<0.05)
  # 
  # plot(its_df_f2$n_yrs, its_df_f2$area50)
  # m3 <- aov(area50 ~  n_yrs, data=its_df_f2)
  # aov3 <- anova(m3)
  # ph_m3 <- as.data.frame(TukeyHSD(m3)$n_yrs)
  # ph_m3 <- ph_m3[str_detect(rownames(ph_m3), pattern = "1"), ]
  # stest_m3 <- any(ph_m3$`p adj`<0.05)
  
  # % HR overlap - 95%
  m4 <- aov(hr95_over ~ n_yrs, data=its_df_f2)
  aov4 <- anova(m4)
  ph_m4 <- as.data.frame(TukeyHSD(m4)$n_yrs)
  ph_m4 <- ph_m4[str_detect(rownames(ph_m4), pattern = "1"), ]
  stest_m4 <- any(ph_m4$`p adj`<0.05)
  
  iaud95folder <- paste0("data/analysis/interannaul_HRs_a/chick_rearing/polygon")
  iaud95files  <- list.files(iaudfolder, full.names = T, pattern = "95")
  iaud95files  <- stringr::str_subset(iaud95files, pattern = "href1")
  
  iaud95 <- readRDS(iaud95files[i])
  iaud95_area <- units::set_units(sf::st_area(iaud95), km^2)
  
  if(stest_m4 == T) {
    seffects_m4 <- ph_m4$diff[which(ph_m4$`p adj`<0.05)]
    seffects_m4_lwr <- ph_m4$lwr[which(ph_m4$`p adj`<0.05)]
    seffects_m4_upr <- ph_m4$upr[which(ph_m4$`p adj`<0.05)]
    ### calculate difference in area for significant comparisons - 95%
    areadiff95 <- iaud95_area * mean(seffects_m4)
  } else { 
    areadiff95  <- NA
    seffects_m4 <- NA
    seffects_m4_lwr <- NA
    seffects_m4_upr <- NA
    }
  
  # % HR overlap - 50%
  m5    <- aov(hr50_over ~ n_yrs, data=its_df_f2)
  aov5  <- anova(m5)
  ph_m5 <- as.data.frame(TukeyHSD(m5)$n_yrs)
  ph_m5 <- ph_m5[str_detect(rownames(ph_m5), pattern = "1"), ]
  stest_m5    <- any(ph_m5$`p adj`<0.05)
  
  ### calculate difference in area for significant comparisons - 95%
  iaud50folder <- paste0("data/analysis/interannaul_HRs_a/chick_rearing/polygon")
  iaud50files  <- list.files(iaudfolder, full.names = T, pattern = "50")
  iaud50files  <- stringr::str_subset(iaud50files, pattern = "href1")
  
  iaud50 <- readRDS(iaud50files[i])
  iaud50_area <- units::set_units(sf::st_area(iaud50), km^2)
  
  if(stest_m5 == T) {
    seffects_m5 <- ph_m5$diff[which(ph_m5$`p adj`<0.05)]
    seffects_m5_lwr <- ph_m5$lwr[which(ph_m5$`p adj`<0.05)]
    seffects_m5_upr <- ph_m5$upr[which(ph_m5$`p adj`<0.05)]
    
    areadiff50 <- iaud50_area * mean(seffects_m5)
    
  } else {
    areadiff50  <- NA
    seffects_m5 <- NA
    seffects_m5_lwr <- NA
    seffects_m5_upr <- NA
  }

  result <- data.frame(
    scientific_name = asp,
    site_name   = asite,
    BA_diff  = stest_m1,
    hr95_over_t = stest_m4,
    hr95_covdiff = as.numeric(areadiff95),
    hr95_mneff   = mean(seffects_m4),
    hr95_lwreff   = mean(seffects_m4_lwr),
    hr95_upreff   = mean(seffects_m4_upr),
    hr50_over_t = stest_m5,
    hr50_covdiff = as.numeric(areadiff50),
    hr50_mneff   = mean(seffects_m5),
    hr50_lwreff   = mean(seffects_m5_upr),
    hr50_upreff   = mean(seffects_m5_lwr),
    ia_area95    = iaud95_area,
    ia_area50    = iaud50_area
    
  )
  
  results_list[[i]] <- result
  
} 

results <- rbindlist(results_list)


### Plotting overlap and stats results for each species ### -------------------
its_df_list <- list()

for(i in seq_along(filenames)) {
  print(i)
  
  asp   <- do.call(rbind, str_split(filenames, pattern = "_"))[i,1]
  asite <- do.call(rbind, str_split(filenames, pattern = "_"))[i,2]
  
  its_df <- readRDS(files[i])
  
  maxssize <- max(its_df$n_trx)
  
  # filter to only maximum sample size (n_trx)
  its_df_f <- subset(its_df, n_trx == maxssize) 
  
  # plot(its_df_f$n_yrs, its_df_f$BAval)
  
  # filter to only n_yr of 1 and the max (or 3) for comparison
  max_nyrs <- max(its_df$n_yr)
  
  # its_df_f2 <- subset(its_df_f, its_df_f$n_yr == 1 | its_df_f$n_yr == max_nyrs)
  # its_df_f2 <- subset(its_df_f, its_df_f$n_yr == 1 | its_df_f$n_yr == 3)
  its_df_f2 <- its_df_f
  its_df_f2$n_yrs <- as.factor(its_df_f2$n_yrs)
  
  its_df_list[[i]] <- data.frame(scientific_name=asp, site_name=asite, its_df_f2)
  
}

its_df_all <- rbindlist(its_df_list)

## add ANOVA results for each species 
its_df_all <- its_df_all %>% left_join(results)

## Simple HR % overlap (facetted btwn 95 and 50) ------------------------------

its_df_all <- its_df_all %>% 
  mutate(
    sp_rank = as.numeric(
      as.factor(
        reorder(scientific_name, hr50_over, FUN = "median")))
  )

it_all_long <- tidyr::pivot_longer(
    its_df_all, 
    cols=c("hr95_over", "hr50_over"), 
    names_to = "hr_over_type", 
    values_to = "hr_over") %>% 
  mutate(hr_over_type = factor(hr_over_type, levels=c("hr95_over", "hr50_over")))

label_df <- it_all_long %>% 
  group_by(scientific_name, hr_over_type) %>%
  summarise(
    n_trx   = min(n_trx),
    x_pos_n = first(sp_rank),
    y_pos_n = -0.01,
    x_pos_t50 = ifelse((hr50_over_t == T) & (hr_over_type == "hr50_over"), first(sp_rank), NA),
    y_pos_t50 = ifelse((hr50_over_t == T) & (hr_over_type == "hr50_over"), 1.01, NA),
    x_pos_t95 = ifelse((hr95_over_t == T) & (hr_over_type == "hr95_over"), first(sp_rank), NA),
    y_pos_t95 = ifelse((hr95_over_t == T) & (hr_over_type == "hr95_over"), 1.01, NA)
  ) %>% 
  mutate(
    hr_over_type = factor(hr_over_type, levels=c("hr95_over", "hr50_over"))
    ) #%>% 
  #filter(hr_over_type == "hr50_over") %>% 

ggplot() + 
  geom_boxplot(data=it_all_long, 
               aes(x=reorder(scientific_name, sp_rank), 
                   y=hr_over, fill=n_yrs), color = "black",
               outlier.size = .5) + 
  scale_fill_manual(values = c("red", "grey30", "grey50", "grey70")) +
  ylab("Overlap (BA)") + 
  ylab("% of HR covered") + 
  xlab("") + labs(fill="N years") +
  ylim(c(-0.01,1.02)) + 
  geom_text(data = label_df,   ## sample size label
            aes(x=x_pos_n, y=y_pos_n, label = paste("n", n_trx)), size=3) +
  geom_text(data = label_df, ## anova significance label
            aes(x=x_pos_t50, y_pos_t50), label = "*", size=6) +
  geom_text(data = label_df, ## anova significance label
            aes(x=x_pos_t95, y_pos_t95), label = "*", size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~hr_over_type, nrow=2)

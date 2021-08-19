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
  
  png(paste0("figures/n_effects/max_ssize/", asp, "_", asite, "_", "n", maxssize, "_", "overlap", ".png"), width = 580, height = 580)
  plot(its_df_f2$n_yrs, its_df_f2$BAval, 
       ylab="Overlap (BA)", xlab="N years",
       main = paste(asp, "n =", maxssize))
  dev.off()
  
  m1 <- aov(BAval ~ n_yrs, data=its_df_f2)
  aov1 <- anova(m1)
  ph_m1 <- as.data.frame(TukeyHSD(m1)$n_yrs)
  ph_m1 <- ph_m1[str_detect(rownames(ph_m1), pattern = "1"), ]
  seffect_m1 <- any(ph_m1$`p adj`<0.05)
  
  ## Area 95% 
  png(
    paste0("figures/n_effects/max_ssize/", asp, "_", asite, "_", "n", maxssize, "_", "area95", ".png"),
    width = 580, height = 580)
  plot(its_df_f2$n_yrs, its_df_f2$area95,
    ylab="Area (sq. km)", xlab="N years",
       main = paste(asp, "n =", maxssize))
  dev.off()
  
  m2 <- aov(area95 ~  n_yrs, data=its_df_f2)
  aov2 <- anova(m2)
  ph_m2 <- as.data.frame(TukeyHSD(m2)$n_yrs)
  ph_m2 <- ph_m2[str_detect(rownames(ph_m2), pattern = "1"), ]
  seffect_m2 <- any(ph_m2$`p adj`<0.05)
  
  plot(its_df_f2$n_yrs, its_df_f2$area50)
  m3 <- aov(area50 ~  n_yrs, data=its_df_f2)
  aov3 <- anova(m3)
  ph_m3 <- as.data.frame(TukeyHSD(m3)$n_yrs)
  ph_m3 <- ph_m3[str_detect(rownames(ph_m3), pattern = "1"), ]
  seffect_m3 <- any(ph_m3$`p adj`<0.05)

  result <- data.frame(
    scientific_name = asp,
    site_name   = asite,
    ovrlp_diff  = seffect_m1,
    area95_diff = seffect_m2 ,
    area50_diff = seffect_m3
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

## choose which overlap measure to use
its_df_all$x_pos <- as.numeric(reorder(its_df_all$scientific_name, its_df_all$BAval, FUN = "median"))
its_df_all$x_pos <- as.numeric(reorder(its_df_all$scientific_name, its_df_all$hr95_over, FUN = "median"))

## labelling ANOVA results above each species
label_df <- its_df_all %>% group_by(scientific_name) %>%
  mutate(
    # y_pos = max(na.omit(BAval)),
    y_pos = max(na.omit(hr95_over)) 
  ) %>% 
  summarise(
    x_pos = ifelse(ovrlp_diff == T, first(x_pos), NA),
    # y_pos = ifelse(ovrlp_diff == T, first(y_pos) + 0.01, NA),
    y_pos = ifelse(ovrlp_diff == T, .99, NA),
  )

ggplot() + 
  geom_boxplot(data=its_df_all, 
               # aes(x=reorder(scientific_name, BAval, FUN = "median"), 
               #     y=BAval, fill=n_yrs),
               # aes(x=reorder(scientific_name, hr95_over, FUN = "median"), 
               #     y=hr95_over, fill=n_yrs),
               aes(x=reorder(scientific_name, hr50_over, FUN = "median"), 
                   y=hr50_over, fill=n_yrs),
               outlier.size = .5) + 
  scale_fill_manual(values = c("red", "grey20", "grey30", "grey40", "grey50", "grey60")) +
  ylab("Overlap (BA)") + 
  ylab("% of 95% HR covered") + 
  xlab("") + labs(fill="N years") +
  ylim(c(0,1)) + 
  # geom_text(data = label_df, aes(x=x_pos, y=y_pos), label = "*", size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/overlap_resamp_iaUDs_href1_i50_ANOVA.png", width=8, height=6)

## AREA 

its_df_all$x_pos <- as.numeric(reorder(its_df_all$scientific_name, its_df_all$area95, FUN = "median"))

## labelling ANOVA results above each species
label_df <- its_df_all %>% group_by(scientific_name) %>%
  mutate(
    y_pos = max(na.omit(BAval)) 
  ) %>% 
  summarise(
    x_pos = ifelse(area95_diff == T, first(x_pos), NA),
    # y_pos = ifelse(ovrlp_diff == T, first(y_pos) + 0.01, NA),
    y_pos = ifelse(area95_diff == T, 14.9, NA),
  )

ggplot() + 
  geom_boxplot(data=its_df_all, 
               aes(x=reorder(scientific_name, log(area95), FUN = "median"), 
                   y=log(area95), fill=n_yrs),
               outlier.size = .5) + 
  scale_fill_manual(values = c("blue", "grey20", "grey30", "grey40", "grey50", "grey60")) +
  ylab("Area sq. km (log)") + xlab("") + labs(fill="N years") +
  geom_text(data = label_df, aes(x=x_pos, y=y_pos), label = "*", size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/area95_resamp_iaUDs_href1_i50_ANOVA.png", width=8, height=6)


its_df_all$x_pos <- as.numeric(reorder(its_df_all$scientific_name, its_df_all$area50, FUN = "median"))

# 50% UD area
label_df <- its_df_all %>% group_by(scientific_name) %>%
  mutate(
    y_pos = max(na.omit(BAval)) 
  ) %>% 
  summarise(
    x_pos = ifelse(area50_diff == T, first(x_pos), NA),
    # y_pos = ifelse(ovrlp_diff == T, first(y_pos) + 0.01, NA),
    y_pos = ifelse(area50_diff == T, 14.9, NA),
  )

ggplot() + 
  geom_boxplot(data=its_df_all, 
               aes(x=reorder(scientific_name, log(area50), FUN = "median"), 
                   y=log(area50), fill=n_yrs),
               outlier.size = .5) + 
  scale_fill_manual(values = c("blue", "grey20", "grey30", "grey40", "grey50", "grey60")) +
  ylab("Area sq. km (log)") + xlab("") + labs(fill="N years") +
  geom_text(data = label_df, aes(x=x_pos, y=y_pos), label = "*", size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/area50_resamp_iaUDs_href1_i50_ANOVA.png", width=8, height=6)


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
    x_pos = first(sp_rank),
    y_pos = 0.01,
    n_trx = min(n_trx)
  ) %>% filter(hr_over_type == "hr50_over") %>% 
  mutate(hr_over_type = factor(hr_over_type, levels=c("hr95_over", "hr50_over")))

ggplot() + 
  geom_boxplot(data=it_all_long, 
               aes(x=reorder(scientific_name, sp_rank), 
                   y=hr_over, fill=n_yrs), color = "black",
               outlier.size = .5) + 
  scale_fill_manual(values = c("red", "grey30", "grey50", "grey70")) +
  ylab("Overlap (BA)") + 
  ylab("% of HR covered") + 
  xlab("") + labs(fill="N years") +
  ylim(c(0,1)) + 
  geom_text(data = label_df,   ## sample size label
            aes(x=x_pos, y=y_pos, label = paste("n", n_trx)), size=2.5) +
  # geom_text(data = label_df, ## anova significance label
    # aes(x=x_pos, y=y_pos), label = "*", size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~hr_over_type, nrow=2)

ggsave("figures/hrover_nyreff_href1_i100.png", width=9, height=6.5)

## Correlation btwn ranks based on HR % overlap and BA overlap ----------------

ranks <- its_df_all %>% 
  mutate(
    HR50_sp_rank = as.numeric(
      as.factor(
        reorder(scientific_name, hr50_over, FUN = "median"))),
    HR95_sp_rank = as.numeric(
      as.factor(
        reorder(scientific_name, hr95_over, FUN = "median"))),
    BA_sp_rank = as.numeric(
      as.factor(
        reorder(scientific_name, BAval, FUN = "median")))
  ) %>% group_by(scientific_name) %>% 
  summarise(
    HR50_sp_rank = first(HR50_sp_rank),
    HR95_sp_rank = first(HR95_sp_rank),
    BA_sp_rank = first(BA_sp_rank)
  ) 

## rank correlation tests
r1 <- cor.test(ranks$BA_sp_rank, ranks$HR50_sp_rank)
r2 <- cor.test(ranks$BA_sp_rank, ranks$HR95_sp_rank)
r3 <- cor.test(ranks$HR50_sp_rank, ranks$HR95_sp_rank)

p1 <- ggplot() + # BA and 50% HR
  geom_point(data=ranks, 
             aes(x = BA_sp_rank, y = HR50_sp_rank)) + theme_bw() +
  geom_text(aes(x=3.5, y=22), 
            label=paste0("corr = ", round(r1$estimate,2), ",", " p < 0.001"))
p2 <- ggplot() + # BA and 95% HR
  geom_point(data=ranks, 
             aes(x = BA_sp_rank, y = HR95_sp_rank)) + theme_bw() +
  geom_text(aes(x=3.5, y=22), 
            label=paste0("corr = ", round(r2$estimate,2), ",", " p < 0.001"))
p3 <- ggplot() + # 50% and 95% HR
  geom_point(data=ranks, 
             aes(x = HR50_sp_rank, y = HR95_sp_rank)) + theme_bw() +
  geom_text(aes(x=3.5, y=22), 
            label=paste0("corr = ", round(r3$estimate,2), ",", " p = ", round(r3$p.value,2)))


library(patchwork)

p <- p1 / p2 / p3

ggsave("figures/overindexcorr_nyreff_href1_i100.png", width=6, height=10)

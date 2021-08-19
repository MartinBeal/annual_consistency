## Plot output of n_effects scripts ## 

pacman::p_load(ggplot2, stringr)

folder <- "data/analysis/n_effects_avg/"

n_its <- "i50"

files     <- list.files(folder, full.names = T, pattern = n_its)
filenames <- tools::file_path_sans_ext(list.files(folder, pattern = n_its))

# its <- do.call(rbind, str_split(filenames, pattern = "_"))[1,3]

files <- files[1]
filenames <- filenames[1]

for(i in seq_along(files)){
  
  asp   <- do.call(rbind, str_split(filenames, pattern = "_"))[i,1]
  asite <- do.call(rbind, str_split(filenames, pattern = "_"))[i,2]
  
  avg_out_df <- readRDS(files[i]) %>% 
    mutate(
      m_hr95  = m_hr95 * 100,
      m_hr50  = m_hr50 * 100,
      sd_hr95 = sd_hr95 * 100,
      sd_hr50 = sd_hr50 * 100,
    )
  
  # plot species-site results
  # Overlap - Number of years ----------------------------
  p_nyrs <- ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_BA-sd_BA, ymax=m_BA+sd_BA, x=n_yrs, color=as.factor(n_trx)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_yrs, m_BA, fill=as.factor(n_trx)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, 1) + xlab("Number of years") + ylab("Overlap (BA)") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N tracks", fill = "N tracks") +
    theme_bw()
  
  # Overlap - Number of tracks ----------------------------
  p_ntrx <- ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_BA-sd_BA, ymax=m_BA+sd_BA, x=n_trx, color=as.factor(n_yrs)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_trx, m_BA, fill=as.factor(n_yrs)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, 1) + xlab("Number of tracks") + ylab("Overlap (BA)") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N years", fill = "N years") +
    theme_bw()
  
  filename <- paste0("figures/n_effects/nyrs_", 
                     asp, "_", asite, "_", n_its, ".png")
  ggsave(filename, plot = p_nyrs, width=7, height=7)
  
  filename <- paste0("figures/n_effects/ntrx_", 
                     asp, "_", asite, "_", n_its, ".png")
  ggsave(filename, plot = p_ntrx, width=7, height=7)
  
  # Area 95% - Number of years ----------------------------
  p_nyrs_95 <- ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_95-sd_95, ymax=m_95+sd_95, x=n_yrs, color=as.factor(n_trx)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_yrs, m_95, fill=as.factor(n_trx)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, max(na.omit(avg_out_df$m_95 + avg_out_df$sd_95))) +
    xlab("Number of years") + ylab("Area (sq. km)") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N tracks", fill = "N tracks") +
    theme_bw()
  
  # Area 95% -  Number of tracks ----------------------------
  p_ntrx_95 <- ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_95-sd_95, ymax=m_95+sd_95, x=n_trx, color=as.factor(n_yrs)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_trx, m_95, fill=as.factor(n_yrs)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, max(na.omit(avg_out_df$m_95 + avg_out_df$sd_95))) + 
    xlab("Number of tracks") + ylab("Area (sq. km)") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N years", fill = "N years") +
    theme_bw()
  
  filename <- paste0("figures/n_effects/area/nyrs_area95_", 
                     asp, "_", asite, "_", n_its, ".png")
  ggsave(filename, plot = p_nyrs_95, width=7, height=7)
  
  filename <- paste0("figures/n_effects/area/ntrx_area95_", 
                     asp, "_", asite, "_", n_its, ".png")
  ggsave(filename, plot = p_ntrx_95, width=7, height=7)
  
  # Area 50% - Number of years ----------------------------
  p_nyrs_50 <- ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_50-sd_50, ymax=m_50+sd_50, x=n_yrs, color=as.factor(n_trx)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_yrs, m_50, fill=as.factor(n_trx)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, max(na.omit(avg_out_df$m_50 + avg_out_df$sd_50))) +
    xlab("Number of years") + ylab("Area (sq. km)") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N tracks", fill = "N tracks") +
    theme_bw()
  
  # Overlap - Number of tracks ----------------------------
  p_ntrx_50 <- ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_50-sd_50, ymax=m_50+sd_50, x=n_trx, color=as.factor(n_yrs)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_trx, m_50, fill=as.factor(n_yrs)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, max(na.omit(avg_out_df$m_50 + avg_out_df$sd_50))) +
    xlab("Number of tracks") + ylab("Area (sq. km)") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N years", fill = "N years") +
    theme_bw()
  
  filename <- paste0("figures/n_effects/area/nyrs_area50_", 
                     asp, "_", asite, "_", n_its, ".png")
  ggsave(filename, plot = p_nyrs_50, width=7, height=7)
  
  filename <- paste0("figures/n_effects/area/ntrx_area50_", 
                     asp, "_", asite, "_", n_its, ".png")
  ggsave(filename, plot = p_ntrx_50, width=7, height=7)
  
  ### Basic % HR overlap measure ---------------------------------------------
  ## 50%
  ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_hr50-sd_hr50, ymax=m_hr50+sd_hr50, 
          x=n_trx, color=as.factor(n_yrs)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_trx, m_hr50, fill=as.factor(n_yrs)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, 100) + xlab("Number of tracks") + ylab("% of 50% HR covered") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N years", fill = "N years") +
    theme_bw()
  
  ## 95%
  ggplot() +
    geom_linerange(
      data=avg_out_df,
      aes(ymin=m_hr95-sd_hr95, ymax=m_hr95+sd_hr95, x=n_trx, color=as.factor(n_yrs)),
      position = position_dodge(width = .5), size=2) +
    geom_point(
      data=avg_out_df,
      aes(n_trx, m_hr95, fill=as.factor(n_yrs)),
      position = position_dodge(width = .5), color = "black", size=5, pch=21) +
    ylim(0, 100) + xlab("Number of tracks") + ylab("% of 95% HR covered") +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    labs(color = "N years", fill = "N years") +
    theme_bw()
}


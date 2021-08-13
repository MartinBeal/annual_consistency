## Plot output of n_effects scripts ## 

pacman::p_load(ggplot2, stringr)

folder <- "data/analysis/n_effects/"

n_its <- "i50"

files     <- list.files(folder, full.names = T, pattern = n_its)
filenames <- tools::file_path_sans_ext(list.files(folder, pattern = n_its))

# its <- do.call(rbind, str_split(filenames, pattern = "_"))[1,3]

# files <- files[11]
# filenames <- filenames[11]

for(i in seq_along(files)){
  
  asp   <- do.call(rbind, str_split(filenames, pattern = "_"))[1,1]
  asite <- do.call(rbind, str_split(filenames, pattern = "_"))[1,2]
  
  avg_out_df <- readRDS(files[i])
  
  # plot species-site results
  # Number of years ----------------------------
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
  
  # Number of tracks ----------------------------
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
                     asp, "_", asite, "_", "i", n_its, ".png")
  ggsave(filename, plot = p_nyrs, width=7, height=7)
  
  filename <- paste0("figures/n_effects/ntrx_", 
                     asp, "_", asite, "_", "i", n_its, ".png")
  ggsave(filename, plot = p_ntrx, width=7, height=7)
}

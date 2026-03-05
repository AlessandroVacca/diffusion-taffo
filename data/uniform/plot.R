library(ggplot2)
library(cowplot)
library(sitools)
library(viridis)
library(dplyr)
library(tidyr)


W <- 5.126
H <- 2
S <- 1
point_size <- 0.8
line_size <- 0.5
linecolors <- scale_color_brewer(palette = "Set1")
theme <- theme_cowplot(font_size = 7)

sisec <- Vectorize(
  function(t) {
    if (is.na(t)) {
      NA
    } else {
      sitools::f2si(t / 10^6, "s")
    }
  }
)

domain_label <- function(x) parse(text = paste0(x, "^3"))

# iterations * substrates * visited voxels
ops <- function(n, s) 100 * s * n**3

{
  data <- read.csv("gpp/biofvms.csv", header = TRUE, sep = ",")
  
  data = filter(data, std_dev / time < 0.05)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data <- data.frame(data)

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  ggsave("temp-local.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time, color = factor(s), shape = factor(s)
    )) +
              geom_errorbar(aes(
        ymin = (time - std_dev),
        ymax = (time + std_dev)
      ), width = 0.08) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Wall-time (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(9, "YlGnBu")[2:9]) +
      labs(color = "Substrates", shape = "Substrates") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )

  data = filter(data, s == 1 & nx %in% c(50, 100, 150, 200, 250, 300, 400, 500, 600))

  ggsave("temp-local-normalized.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_brewer(palette = "Accent") +
      # labs(color = "Substrates", shape = "Substrates") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(64, 128, 256, 400, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("gpp/baselines.csv", header = TRUE, sep = ",")
  
  data = filter(data, std_dev / time < 0.05)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data <- data.frame(data)


  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  ggsave("data-local.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time, color = factor(s), shape = factor(s)
    )) +
          geom_errorbar(aes(
        ymin = (time - std_dev),
        ymax = (time + std_dev)
      ), width = 0.08) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Wall-time (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(9, "YlGnBu")[2:9]) +
      labs(color = "Substrates", shape = "Substrates") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )

  data2 <- read.csv("gpp/biofvms.csv", header = TRUE, sep = ",")
  
  data2 = filter(data2, std_dev / time < 0.05)

  data2 <- data2 %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data2 <- data.frame(data2)
  
  data2$precision[data2$precision == "D"] <- "Double precision"
  data2$precision[data2$precision == "S"] <- "Single precision"

  data3 <- c()
  data3 = rbind(data, data2)
  data3 = filter(data3, s == 1 & nx %in% c(50, 100, 150, 200, 250, 300, 400, 450, 500, 550, 600))

  data3$algorithm[data3$algorithm == "lstcss"] <- "temp-local (baseline)"
  data3$algorithm[data3$algorithm == "lstcs"] <- "space-local"

  ggsave("data-local-normalized.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data3, aes(
      x = nx, y = time / ops(nx, s), color = algorithm, shape = algorithm
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_brewer(palette = "Accent") +
      labs(color = "Algorithm", shape = "Algorithm") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(64, 128, 256, 400, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("gpp/temporals.csv", header = TRUE, sep = ",")
  data2 <- read.csv("gpp/baselines.csv", header = TRUE, sep = ",")
  data3 <- read.csv("gpp/biofvms.csv", header = TRUE, sep = ",")
  
  data = filter(data, x_tile_size != 10000)

  data$empty = 1
  data2$substrate_step = 1
  data2$x_tile_size = 10000
  data3$substrate_step = 1
  data3$x_tile_size = 1

  data <- rbind(data, data2, data3)

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"
  
  data = filter(data, std_dev / time < 0.05)


  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data <- data.frame(data)

  data = filter(data, s == 1 & nx %in% c( 100, 150, 200, 250, 300, 400, 450, 500, 550, 600))

  data$x_tile_size[data$x_tile_size == "10000"] <- "whole X (space-local)"
  data$x_tile_size[data$x_tile_size == "1"] <- "1 (baseline)"
  
  data$x_tile_size <- factor(data$x_tile_size, levels = c("16", "32", "48", "64", "whole X (space-local)", "1 (baseline)"))

  ggsave("temp-local-tile.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(x_tile_size)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(8, "YlGnBu")[2:9]) +
      labs(color = "X tile size") +
      scale_y_log10(labels = sisec) +
      scale_x_continuous(labels = domain_label) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("gpp/transpose-temporal.csv", header = TRUE, sep = ",")
  data2 <- read.csv("gpp/temporals.csv", header = TRUE, sep = ",")
  data3 <- read.csv("gpp/partial-blocking.csv", header = TRUE, sep = ",")
  data3$x_tile_size = 32

  data = rbind(data, data2, data3)
  
  data = filter(data, x_tile_size %in%  c(32))

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"
  
  data = filter(data, std_dev / time < 0.05)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data <- data.frame(data)

  data = filter(data, s == 1 & nx %in% c(100, 150, 200, 250, 300, 400, 450, 500, 600))

  data$x_tile_size[data$x_tile_size == "10000"] <- "whole X"
  
  data$algorithm[data$algorithm == "lstcstai"] <- "tiled+transposed"
  data$algorithm[data$algorithm == "lstcsta"] <- "tiled"
  data$algorithm[data$algorithm == "lstmfppai"] <- "planar"

  ggsave("temp-local-tile-transpose.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = algorithm, shape = algorithm
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size") +
      ylab("Time per voxel (log-scale)") +
      scale_color_brewer(palette = "Accent") +
      labs(color = "Algorithm", shape = "Algorithm") +
      scale_y_log10(labels = sisec) +
      scale_x_continuous(labels = domain_label) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("gpp/partial-blocking.csv", header = TRUE, sep = ",")

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(std_dev) %>%
    ungroup()
  data <- data.frame(data)

  data <- filter(data, s == 1)

  ggsave("partial-blocking.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(substrate_step), shape = factor(substrate_step)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(5, "YlGnBu")[2:6]) +
      labs(color = "Substrate Step", shape = "Substrate Step") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}


{
  data_base <- read.csv("gpp/full-blocking.csv", header = TRUE, sep = ",")
  data_base["alt"] = "Alt false"
  data_base["dist"] = "Dist false"

  data_alt <- read.csv("gpp/full-blocking-alt.csv", header = TRUE, sep = ",")
  data_alt["alt"] = "Alt true"
  data_alt["dist"] = "Dist false"

  data_dist <- read.csv("gpp/full-blocking-dist.csv", header = TRUE, sep = ",")
  data_dist["alt"] = "Alt false"
  data_dist["dist"] = "Dist true"

  data_alt_dist <- read.csv("gpp/full-blocking-alt-dist.csv", header = TRUE, sep = ",")
  data_alt_dist["alt"] = "Alt true"
  data_alt_dist["dist"] = "Dist true"

  data <- rbind(data_base, data_alt, data_dist, data_alt_dist)

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    separate(cores_division, into = c("cores_x", "cores_y", "cores_z"), sep = ",", convert = TRUE, remove = FALSE)
    
  data = filter(data, nx / cores_x >= 3 & ny / cores_y >= 3 & nz / cores_z >= 3)

  data = filter(data, std_dev / time < 0.05)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data <- data.frame(data)

  # data_fastest <- data %>%
  #   group_by(precision, dims, iterations, s, nx, ny, nz) %>%
  #   slice_min(time) %>%
  #   ungroup()
  # data_fastest <- data.frame(data_fastest)

  # print(data_fastest)
  
  data = filter(data, s == 1 & nx %in% c(50, 100, 150, 200, 250, 300, 400, 500, 600))

  data_sync_step = filter(data, s == 1 & precision == "Single precision")

  ggsave("full-blocking-div-S.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H * 2,
    ggplot(data_sync_step, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(cores_division)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = viridis(10)) +
      labs(color = "Cores division") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_grid(alt~dist) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )

  data_sync_step = filter(data, s == 1 & precision == "Double precision")

  ggsave("full-blocking-div-D.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H * 2,
    ggplot(data_sync_step, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(cores_division)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = viridis(10)) +
      labs(color = "Cores division") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_grid(alt~dist) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("gpp/full-blocking-alt-dist.csv", header = TRUE, sep = ",")

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    separate(cores_division, into = c("cores_x", "cores_y", "cores_z"), sep = ",", convert = TRUE, remove = FALSE)
    
  data = filter(data, nx / cores_x >= 3 & ny / cores_y >= 3 & nz / cores_z >= 3)

  data = filter(data, std_dev / time < 0.1)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()
  data <- data.frame(data)

  data = filter(data, !(cores_division %in% c("1,1,112", "1,2,56", "1,56,2")))
  
  data = filter(data, s == 1 & nx %in% c(50, 100, 150, 200, 250, 300, 400, 500, 600))
  
  data$cores_division <- factor(data$cores_division, levels = c("1,1,112", "1,2,56", "1,4,28", "1,7,16", "1,8,14", "1,14,8", "1,16,7", "1,28,4", "1,56,2"))

  data_sync_step = filter(data, s == 1)

  ggsave("full-blocking.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data_sync_step, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(cores_division)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(7, "YlGnBu")[2:9]) +
      labs(color = "Cores division") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(64, 128, 256, 400, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

# for (prefix in c("", "hbm-"))
# {
#   data_baseline <- read.csv(paste0(prefix, "baseline.csv"), header = TRUE, sep = ",")
#   data_baseline <- subset(data_baseline, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
#   data_baseline["algorithm"] <- "Baseline"

#   data_temporal <- read.csv(paste0(prefix, "transpose-temporal.csv"), header = TRUE, sep = ",")
#   data_temporal <- filter(data_temporal, x_tile_size %in% c(32, 10000) & continuous_x_diagonal == "true")
#   data_temporal$algorithm[data_temporal$x_tile_size == "10000"] <- "Butterfly"
#   data_temporal$algorithm[data_temporal$x_tile_size == "32"] <- "Temporal"
#   data_temporal <- subset(data_temporal, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

#   partial_blocking <- read.csv(paste0(prefix, "partial-blocking.csv"), header = TRUE, sep = ",")
#   partial_blocking <- filter(partial_blocking, x_tile_size == 32 & continuous_x_diagonal == "true")
#   partial_blocking <- subset(partial_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
#   partial_blocking["algorithm"] <- "Partial blocking"

#   full_blocking <- read.csv(paste0(prefix, "full-blocking.csv"), header = TRUE, sep = ",")
#   full_blocking <- filter(full_blocking, x_tile_size == 10000 & cores_division == "4,4,7" & sync_step == 8 & streams == 2)
#   full_blocking <- subset(full_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
#   full_blocking["algorithm"] <- "Full blocking"

#   data <- c()
#   data <- rbind(data, data_baseline, data_temporal, partial_blocking, full_blocking)

#   data <- filter(data, s == 1)

#   data$algorithm <- factor(data$algorithm, levels = c("Baseline", "Butterfly", "Temporal", "Partial blocking", "Full blocking"))

#   data$precision[data$precision == "D"] <- "Double precision"
#   data$precision[data$precision == "S"] <- "Single precision"

#   data <- data %>%
#     group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
#     slice_min(time) %>%
#     ungroup()

#   data <- data.frame(data)

#   ggsave(paste0(prefix, "all.pdf"),
#     device = "pdf", units = "in", scale = S, width = W, height = H,
#     ggplot(data, aes(
#       x = nx, y = time / ops(nx, s),
#       color = algorithm, shape = algorithm
#     )) +
#       geom_point(size = point_size) +
#       geom_line(linewidth = line_size) +
#       xlab("Domain size (log-scale)") +
#       ylab("Time per voxel (log-scale)") +
#       scale_color_brewer(palette = "Set1") +
#       labs(color = "Algorithm", shape = "Algorithm") +
#       scale_y_log10(labels = sisec) +
#       scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
#       facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
#       theme +
#       background_grid() +
#       theme(legend.position = "bottom")
#   )
# }

for (prefix in c("gpp", "hbm"))
{
  data_baseline <- read.csv(paste0(prefix, "/baselines.csv"), header = TRUE, sep = ",")
  data_baseline <- subset(data_baseline, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  data_baseline["algorithm"] <- "space-local"

  data_temporal <- read.csv(paste0(prefix, "/biofvms.csv"), header = TRUE, sep = ",")
  data_temporal$algorithm <- "temp-local"
  data_temporal <- subset(data_temporal, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  data_temporal_tile <- read.csv(paste0(prefix, "/temporals.csv"), header = TRUE, sep = ",")
  data_temporal_tile$algorithm <- "tiled"
  data_temporal_tile = filter(data_temporal_tile, x_tile_size == 32)
  data_temporal_tile <- subset(data_temporal_tile, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  data_temporal_tile_trans <- read.csv(paste0(prefix, "/transpose-temporal.csv"), header = TRUE, sep = ",")
  data_temporal_tile_trans$algorithm <- "tiled-transposed"
  data_temporal_tile_trans = filter(data_temporal_tile_trans, x_tile_size == 32)
  data_temporal_tile_trans <- subset(data_temporal_tile_trans, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  partial_blocking <- read.csv(paste0(prefix, "/partial-blocking.csv"), header = TRUE, sep = ",")
  partial_blocking <- subset(partial_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  partial_blocking["algorithm"] <- "planar"

  full_blocking <- read.csv(paste0(prefix, "/full-blocking-alt-dist.csv"), header = TRUE, sep = ",")
  full_blocking <- full_blocking %>%
    separate(cores_division, into = c("cores_x", "cores_y", "cores_z"), sep = ",", convert = TRUE, remove = FALSE)
  full_blocking = filter(full_blocking, nx / cores_x >= 3 & ny / cores_y >= 3 & nz / cores_z >= 3)
  full_blocking = filter(full_blocking, cores_division == "1,8,14")
  full_blocking <- subset(full_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  full_blocking["algorithm"] <- "blocked"

  data <- c()
  data <- rbind(data, data_baseline, data_temporal, data_temporal_tile, data_temporal_tile_trans, partial_blocking, full_blocking)

  data = filter(data, s == 1 & nx %in% c(50, 100, 150, 200, 250, 300, 400, 500, 600))

  data$algorithm <- factor(data$algorithm, levels = c("space-local", "temp-local", "tiled", "tiled-transposed", "planar", "blocked"))

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"
  
  if (prefix == "gpp")
    data = filter(data, std_dev / time < 0.05)
  else
    data = filter(data, std_dev / time < 0.2)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()

  data <- data.frame(data)

  ggsave(paste0(prefix, "-all-best.pdf"),
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = algorithm, shape = algorithm
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_brewer(palette = "Set1") +
      labs(color = "Algorithm", shape = "Algorithm") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

data_both = c()

for (prefix in c("gpp", "hbm"))
{
  data_baseline <- read.csv(paste0(prefix, "/baselines.csv"), header = TRUE, sep = ",")
  data_baseline <- subset(data_baseline, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  data_baseline["algorithm"] <- "space-local"

  data_temporal <- read.csv(paste0(prefix, "/biofvms.csv"), header = TRUE, sep = ",")
  data_temporal$algorithm <- "temp-local"
  data_temporal <- subset(data_temporal, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  data_temporal_tile <- read.csv(paste0(prefix, "/temporals.csv"), header = TRUE, sep = ",")
  data_temporal_tile$algorithm <- "tiled"
  data_temporal_tile = filter(data_temporal_tile, x_tile_size == 32)
  data_temporal_tile <- subset(data_temporal_tile, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  data_temporal_tile_trans <- read.csv(paste0(prefix, "/transpose-temporal.csv"), header = TRUE, sep = ",")
  data_temporal_tile_trans$algorithm <- "tiled-transposed"
  data_temporal_tile_trans = filter(data_temporal_tile_trans, x_tile_size == 32)
  data_temporal_tile_trans <- subset(data_temporal_tile_trans, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  partial_blocking <- read.csv(paste0(prefix, "/partial-blocking.csv"), header = TRUE, sep = ",")
  partial_blocking <- subset(partial_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  partial_blocking["algorithm"] <- "planar"

  full_blocking <- read.csv(paste0(prefix, "/full-blocking-alt-dist.csv"), header = TRUE, sep = ",")
  full_blocking <- full_blocking %>%
    separate(cores_division, into = c("cores_x", "cores_y", "cores_z"), sep = ",", convert = TRUE, remove = FALSE)
  full_blocking = filter(full_blocking, nx / cores_x >= 3 & ny / cores_y >= 3 & nz / cores_z >= 3)
  full_blocking = filter(full_blocking, cores_division == "1,8,14")
  full_blocking <- subset(full_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  full_blocking["algorithm"] <- "blocked"

  data <- c()
  data <- rbind(data, data_baseline, data_temporal, data_temporal_tile, data_temporal_tile_trans, partial_blocking)

  if (prefix == "gpp")
    data = filter(data, std_dev / time < 0.05)
  else
    data = filter(data, std_dev / time < 0.3)

  data = rbind(data, full_blocking)

  data = filter(data, s == 1 & nx %in% c(50, 100, 150, 200, 250, 300, 400, 500, 600))

  data$algorithm <- factor(data$algorithm, levels = c("space-local", "temp-local", "tiled", "tiled-transposed", "planar", "blocked"))

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"
  
  if (prefix == "gpp")
    data = filter(data, std_dev / time < 0.1)
  else
    data = filter(data, std_dev / time < 0.3)

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()

  data <- data.frame(data)

  if (prefix == "gpp")
    data$machine = "DDR"
  else
    data$machine = "HBM"

  data_both = rbind(data_both, data)
}

ggsave("all-best.pdf",
  device = "pdf", units = "in", scale = S, width = W, height = H * 1.5,
  ggplot(data_both, aes(
    x = nx, y = time / ops(nx, s),
    color = algorithm, shape = algorithm
  )) +
    geom_point(size = point_size) +
    geom_line(linewidth = line_size) +
    xlab("Domain size (log-scale)") +
    ylab("Time per voxel (log-scale)") +
    scale_color_brewer(palette = "Accent") +
    labs(color = "Algorithm", shape = "Algorithm") +
    scale_y_log10(labels = sisec) +
    scale_x_continuous(labels = domain_label) +
    facet_grid(machine~factor(precision, levels = c("Single precision", "Double precision")), scales="free") +
    theme +
    background_grid() +
    theme(legend.position = "bottom")
)

for (prefix in c("grc/"))
{
  data_baseline <- read.csv(paste0(prefix, "baselines.csv"), header = TRUE, sep = ",")
  data_baseline <- subset(data_baseline, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  data_baseline["algorithm"] <- "data-local"

  data_temporal <- read.csv(paste0(prefix, "biofvms.csv"), header = TRUE, sep = ",")
  data_temporal$algorithm <- "temp-local"
  data_temporal <- subset(data_temporal, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  partial_blocking <- read.csv(paste0(prefix, "partial-blocking.csv"), header = TRUE, sep = ",")
  partial_blocking <- subset(partial_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  partial_blocking["algorithm"] <- "pb"

  full_blocking <- read.csv(paste0(prefix, "full-blocking.csv"), header = TRUE, sep = ",")
  full_blocking <- full_blocking %>%
    separate(cores_division, into = c("cores_x", "cores_y", "cores_z"), sep = ",", convert = TRUE, remove = FALSE)
  full_blocking = filter(full_blocking, nx / cores_x >= 3 & ny / cores_y >= 3 & nz / cores_z >= 3)
  full_blocking <- subset(full_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  full_blocking["algorithm"] <- "fb"

  full_blockingd <- read.csv(paste0(prefix, "full-blocking-dist.csv"), header = TRUE, sep = ",")
  full_blockingd <- full_blockingd %>%
    separate(cores_division, into = c("cores_x", "cores_y", "cores_z"), sep = ",", convert = TRUE, remove = FALSE)
  full_blockingd = filter(full_blockingd, nx / cores_x >= 3 & ny / cores_y >= 3 & nz / cores_z >= 3)
  full_blockingd <- subset(full_blockingd, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  full_blockingd["algorithm"] <- "fbd"

  data <- c()
  data <- rbind(data, data_baseline, data_temporal, partial_blocking, full_blocking, full_blockingd)

  data <- filter(data, s == 1)

  data$algorithm <- factor(data$algorithm, levels = c("data-local", "temp-local", "pb", "fb", "fbd"))

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()

  data <- data.frame(data)

  ggsave(paste0(prefix, "all-best.pdf"),
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = algorithm, shape = algorithm
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_brewer(palette = "Accent") +
      labs(color = "Algorithm", shape = "Algorithm") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 512, 600)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

# HBM Best
#     algorithm        precision dims iterations s  nx  ny  nz cores_division     streams sync_step x_tile_size init_time       time      std_dev
# 1      sdd-fb Double precision    3        100 1  50  50  50          4,7,4 1         3         8          16       320    3111.22     29.33414
# 2      sdd-fb Double precision    3        100 1  64  64  64          2,8,7 2         3        16       10000       917    4169.62     79.13884
# 3      sdd-fb Double precision    3        100 1  96  96  96          2,8,7 3         2         8          16      1651    9653.82    389.30752
# 4      sdd-fb Double precision    3        100 1 100 100 100         2,14,4 4         2     10000          16       414   11686.66    236.87580
# 5      sdd-fb Double precision    3        100 1 128 128 128          2,7,8 5         2     10000       10000     27406   18464.00    947.76693
# 6      sdd-fb Double precision    3        100 1 150 150 150         1,8,14 6         2         2       10000     18372   38131.24    942.80981
# 7      sdd-fb Double precision    3        100 1 160 160 160         1,16,7 7         2         2          16     15597   45321.80   1203.93354
# 8      sdd-fb Double precision    3        100 1 192 192 192          2,7,8 8         3         2          16     40964   83282.26   2680.86440
# 9      sdd-fb Double precision    3        100 1 200 200 200          4,7,4 9         2         8       10000     20914  113463.26   2274.40428
# 10     sdd-fb Double precision    3        100 1 224 224 224          2,7,8 10        2         1       10000     45987  139118.04  11382.66236
# 11     sdd-fb Double precision    3        100 1 250 250 250         1,8,14 11        3         1          16     81411  184509.52   3302.40753
# 12     sdd-fb Double precision    3        100 1 256 256 256          2,7,8 12        3         2       10000     86567  207490.00   3066.88513
# 13     sdd-fb Double precision    3        100 1 288 288 288         1,14,8 13        2         1          16     95586  277524.38   3647.05451
# 14     sdd-fb Double precision    3        100 1 300 300 300         1,14,8 14        2         1          16    100021  342613.04   3296.60879
# 15     sdd-fb Double precision    3        100 1 350 350 350          2,8,7 15        3         1       10000    138465  548838.46   4621.94499
# 16     sdd-fb Double precision    3        100 1 400 400 400          2,7,8 16        2         1       10000    104420  794476.90   6716.23761

# 66     sdd-fb Single precision    3        100 1  64  64  64          2,8,7 66        2         8       10000       320    2795.16     23.77760
# 67     sdd-fb Single precision    3        100 1  96  96  96          2,7,8 67        2         8       10000       930    5736.72    152.78600
# 68     sdd-fb Single precision    3        100 1 100 100 100          2,7,8 68        2         8       10000       690    7032.76    123.03277
# 69     sdd-fb Single precision    3        100 1 128 128 128         1,14,8 69        2         8       10000      1055   10041.68    400.58464
# 70     sdd-fb Single precision    3        100 1 150 150 150          2,7,8 70        3         8       10000      5940   15751.74    281.28817
# 71     sdd-fb Single precision    3        100 1 160 160 160          2,7,8 71        2        16       10000      3021   18055.40    259.63975
# 72     sdd-fb Single precision    3        100 1 192 192 192         1,8,14 72        3         2       10000      4414   36651.42    268.65942
# 73     sdd-fb Single precision    3        100 1 200 200 200         1,16,7 73        3     10000       10000     16835   45195.56   1100.46659
# 74     sdd-fb Single precision    3        100 1 224 224 224         1,14,8 74        2         4       10000      3563   63361.42   1494.13995
# 75     sdd-fb Single precision    3        100 1 250 250 250         1,14,8 75        2         2          16     10856   99338.50   1123.33767
# 76     sdd-fb Single precision    3        100 1 256 256 256          2,7,8 76        3         2       10000     45818  100261.40   2843.03961
# 77     sdd-fb Single precision    3        100 1 288 288 288         2,14,4 77        3         4       10000     57463  155140.00   3095.63336
# 78     sdd-fb Single precision    3        100 1 300 300 300         1,14,8 78        3         1       10000     88039  176785.40   4185.46727
# 79     sdd-fb Single precision    3        100 1 350 350 350          4,4,7 79        3         2       10000     76556  275750.20   3344.07840
# 80     sdd-fb Single precision    3        100 1 400 400 400          2,7,8 80        2         2       10000     96288  435558.78   9837.63943

# GPP Best
#     algorithm precision dims iterations s  nx  ny  nz cores_division streams     sync_step x_tile_size init_time        time      std_dev
# 1      sdd-fb         D    3        100 1  50  50  50          4,4,7       2 1           8          16       258     3095.78     17.55710
# 2      sdd-fb         D    3        100 1  64  64  64          4,7,4       2 2           8          16       358     4169.70     31.26676
# 3      sdd-fb         D    3        100 1  96  96  96          4,4,7       2 3       10000       10000       519    10071.08    287.79297
# 4      sdd-fb         D    3        100 1 100 100 100          7,4,4       3 4          16       10000      1315    11903.94    241.61013
# 5      sdd-fb         D    3        100 1 128 128 128          4,4,7       2 5       10000          16     23579    19518.54    929.51136
# 6      sdd-fb         D    3        100 1 150 150 150          2,8,7       2 6           4          16      7335    38886.58    376.36552
# 7      sdd-fb         D    3        100 1 160 160 160          2,7,8       2 7           4          16      9143    47232.68    734.03740
# 8      sdd-fb         D    3        100 1 192 192 192          4,4,7       2 8           2       10000     18892   107258.52    694.10303
# 9      sdd-fb         D    3        100 1 200 200 200          2,8,7       3 9           2       10000     28644   139916.04   1367.46748
# 10     sdd-fb         D    3        100 1 224 224 224          2,8,7       1 10          8          16     40379   198079.38    465.69754
# 11     sdd-fb         D    3        100 1 250 250 250          2,7,8       1 11          8          16     41286   310785.84    586.30795
# 12     sdd-fb         D    3        100 1 256 256 256          2,8,7       1 12          4       10000     54199   329220.10    814.52296
# 13     sdd-fb         D    3        100 1 288 288 288          2,8,7       1 13          4       10000     61962   511009.72    439.04938
# 14     sdd-fb         D    3        100 1 300 300 300          2,7,8       1 14          2       10000     83379   615486.48    754.05677
# 15     sdd-fb         D    3        100 1 350 350 350          2,7,8       1 15          2          16     80960  1059293.96  25266.41048
# 16     sdd-fb         D    3        100 1 400 400 400          2,7,8       1 16          2       10000     98636  1632659.54   4519.33491

# 65     sdd-fb         S    3        100 1  50  50  50          4,4,7       2 65      10000          16       764     1959.76     45.97197
# 66     sdd-fb         S    3        100 1  64  64  64          4,7,4       2 66          8          16       283     2402.48     22.21373
# 67     sdd-fb         S    3        100 1  96  96  96          2,7,8       2 67          8       10000       368     5479.08     30.67562
# 68     sdd-fb         S    3        100 1 100 100 100          7,4,4       2 68      10000       10000      1475     6698.22    684.43791
# 69     sdd-fb         S    3        100 1 128 128 128          4,7,4       3 69         16       10000      2730    10136.14    153.03934
# 70     sdd-fb         S    3        100 1 150 150 150          2,8,7       3 70          8       10000      2983    14967.92    347.65597
# 71     sdd-fb         S    3        100 1 160 160 160          2,7,8       2 71      10000       10000     22719    16432.38    414.89217
# 72     sdd-fb         S    3        100 1 192 192 192          2,7,8       3 72          4       10000      9951    37011.02    432.38295
# 73     sdd-fb         S    3        100 1 200 200 200          2,7,8       2 73          4       10000     11266    49298.98    392.40599
# 74     sdd-fb         S    3        100 1 224 224 224          2,7,8       2 74          4       10000     14784    69966.78    793.61961
# 75     sdd-fb         S    3        100 1 250 250 250          2,8,7       3 75          2       10000     18695   126153.20   1092.41016
# 76     sdd-fb         S    3        100 1 256 256 256          2,8,7       2 76          4       10000     20050   135975.64    813.74074
# 77     sdd-fb         S    3        100 1 288 288 288          2,8,7       1 77          8       10000     31928   216211.38    536.68393
# 78     sdd-fb         S    3        100 1 300 300 300          2,8,7       1 78          8          16     36908   261006.64    464.45718
# 79     sdd-fb         S    3        100 1 350 350 350          2,8,7       1 79          4       10000     54669   448368.02   1551.96861
# 80     sdd-fb         S    3        100 1 400 400 400          2,8,7       1 80          4       10000     88677   788839.84  19332.75706

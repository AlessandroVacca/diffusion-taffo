library(ggplot2)
library(cowplot)
library(sitools)
library(viridis)
library(dplyr)


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
ops <- function(n, s) 100 * s * 3 * n**3

{
  data <- read.csv("baseline.csv", header = TRUE, sep = ",")

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(std_dev) %>%
    ungroup()
  data <- data.frame(data)

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  ggsave("baseline.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time, color = factor(s), shape = factor(s)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Wall-time (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(9, "YlGnBu")[2:9]) +
      labs(color = "Substrates", shape = "Substrates") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )

  ggsave("baseline-normalized.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s), color = factor(s), shape = factor(s)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per op (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(9, "YlGnBu")[2:9]) +
      labs(color = "Substrates", shape = "Substrates") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("transpose-temporal.csv", header = TRUE, sep = ",")

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(std_dev) %>%
    ungroup()
  data <- data.frame(data)

  data <- filter(data, s == 1 & continuous_x_diagonal == "true")

  data$x_tile_size[data$x_tile_size == "10000"] <- "whole X"

  ggsave("transpose_temporal.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(x_tile_size), shape = factor(x_tile_size)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(5, "YlGnBu")[2:6]) +
      labs(color = "Y & Z tile size", shape = "Y & Z tile size") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

{
  data <- read.csv("partial-blocking.csv", header = TRUE, sep = ",")

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(std_dev) %>%
    ungroup()
  data <- data.frame(data)

  data <- filter(data, s == 1 & continuous_x_diagonal == "true")

  data$x_tile_size[data$x_tile_size == "10000"] <- "whole X"

  ggsave("partial-blocking.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(x_tile_size), shape = factor(x_tile_size)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(5, "YlGnBu")[2:6]) +
      labs(color = "Y & Z tile size", shape = "Y & Z tile size") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}


{
  data <- read.csv("full-blocking.csv", header = TRUE, sep = ",")

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(std_dev) %>%
    ungroup()
  data <- data.frame(data)

  # data_fastest <- data %>%
  #   group_by(precision, dims, iterations, s, nx, ny, nz) %>%
  #   slice_min(time) %>%
  #   ungroup()
  # data_fastest <- data.frame(data_fastest)

  # print(data_fastest)

  data_sync_step = filter(data, s == 1 & x_tile_size == 10000 & cores_division == "4,4,7" & streams == 1)

  ggsave("full-blocking-sync-step.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data_sync_step, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(sync_step), shape = factor(sync_step)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(7, "YlGnBu")[2:7]) +
      labs(color = "X & Y sync step", shape = "X & Y sync step") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )

  data_streams = filter(data, s == 1 & x_tile_size == 10000 & cores_division == "4,4,7" & sync_step == 8)

  ggsave("full-blocking-streams.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data_streams, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(streams), shape = factor(streams)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(7, "YlGnBu")[2:7]) +
      labs(color = "Streams", shape = "Streams") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )

  data_div = filter(data, s == 1 & x_tile_size == 10000 & streams == 2 & sync_step == 8 & cores_division %in% c("2,4,14", "4,4,7", "2,7,8", "1,8,14", "1,7,16"))

  ggsave("full-blocking-cores.pdf",
    device = "pdf", units = "in", scale = S, width = W, height = H,
    ggplot(data_div, aes(
      x = nx, y = time / ops(nx, s),
      color = factor(cores_division)
    )) +
      geom_point(size = point_size) +
      geom_line(linewidth = line_size) +
      xlab("Domain size (log-scale)") +
      ylab("Time per voxel (log-scale)") +
      scale_color_manual(values = RColorBrewer::brewer.pal(6, "YlGnBu")[2:6]) +
      labs(color = "Cores", shape = "Cores") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

for (prefix in c("", "hbm-"))
{
  data_baseline <- read.csv(paste0(prefix, "baseline.csv"), header = TRUE, sep = ",")
  data_baseline <- subset(data_baseline, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  data_baseline["algorithm"] <- "Baseline"

  data_temporal <- read.csv(paste0(prefix, "transpose-temporal.csv"), header = TRUE, sep = ",")
  data_temporal <- filter(data_temporal, x_tile_size %in% c(32, 10000) & continuous_x_diagonal == "true")
  data_temporal$algorithm[data_temporal$x_tile_size == "10000"] <- "Butterfly"
  data_temporal$algorithm[data_temporal$x_tile_size == "32"] <- "Temporal"
  data_temporal <- subset(data_temporal, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  partial_blocking <- read.csv(paste0(prefix, "partial-blocking.csv"), header = TRUE, sep = ",")
  partial_blocking <- filter(partial_blocking, x_tile_size == 32 & continuous_x_diagonal == "true")
  partial_blocking <- subset(partial_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  partial_blocking["algorithm"] <- "Partial blocking"

  full_blocking <- read.csv(paste0(prefix, "full-blocking.csv"), header = TRUE, sep = ",")
  full_blocking <- filter(full_blocking, x_tile_size == 10000 & cores_division == "4,4,7" & sync_step == 8 & streams == 2)
  full_blocking <- subset(full_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  full_blocking["algorithm"] <- "Full blocking"

  data <- c()
  data <- rbind(data, data_baseline, data_temporal, partial_blocking, full_blocking)

  data <- filter(data, s == 1)

  data$algorithm <- factor(data$algorithm, levels = c("Baseline", "Butterfly", "Temporal", "Partial blocking", "Full blocking"))

  data$precision[data$precision == "D"] <- "Double precision"
  data$precision[data$precision == "S"] <- "Single precision"

  data <- data %>%
    group_by(across(-matches("(init_time)|(time)|(std_dev)"))) %>%
    slice_min(time) %>%
    ungroup()

  data <- data.frame(data)

  ggsave(paste0(prefix, "all.pdf"),
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
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
      facet_wrap(~factor(precision, levels = c("Single precision", "Double precision"))) +
      theme +
      background_grid() +
      theme(legend.position = "bottom")
  )
}

for (prefix in c("", "hbm-"))
{
  data_baseline <- read.csv(paste0(prefix, "baseline.csv"), header = TRUE, sep = ",")
  data_baseline <- subset(data_baseline, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  data_baseline["algorithm"] <- "Baseline"

  data_temporal <- read.csv(paste0(prefix, "transpose-temporal.csv"), header = TRUE, sep = ",")
  data_temporal$algorithm <- "Temporal"
  data_temporal$algorithm[data_temporal$x_tile_size == "10000"] <- "Butterfly"
  data_temporal <- subset(data_temporal, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))

  partial_blocking <- read.csv(paste0(prefix, "partial-blocking.csv"), header = TRUE, sep = ",")
  partial_blocking <- subset(partial_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  partial_blocking["algorithm"] <- "Partial blocking"

  full_blocking <- read.csv(paste0(prefix, "full-blocking.csv"), header = TRUE, sep = ",")
  full_blocking <- subset(full_blocking, select = c(algorithm, precision, dims, iterations, s, nx, ny, nz, init_time, time, std_dev))
  full_blocking["algorithm"] <- "Full blocking"

  data <- c()
  data <- rbind(data, data_baseline, data_temporal, partial_blocking, full_blocking)

  data <- filter(data, s == 1)

  data$algorithm <- factor(data$algorithm, levels = c("Baseline", "Butterfly", "Temporal", "Partial blocking", "Full blocking"))

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
      scale_color_brewer(palette = "Set1") +
      labs(color = "Algorithm", shape = "Algorithm") +
      scale_y_log10(labels = sisec) +
      scale_x_log10(labels = domain_label, breaks = c(32, 64, 128, 256, 384)) +
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

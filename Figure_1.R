library(tidyverse)
library(data.table)

# Single control file for plotting settings / result folders. Resolve relative to
# this script if launched via Rscript, else assume the project working directory.
.all_params <- "all_parameters.R"
.sp <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(.sp) == 1 && file.exists(file.path(dirname(.sp), .all_params)))
  .all_params <- file.path(dirname(.sp), .all_params)
source(.all_params)

output_folder <- RESULTS_FOLDER_MAIN
selected_traj <- bind_rows(lapply(file.path(output_folder,list.files(path = output_folder, pattern = "trajectories")), read.csv)) 
setDT(selected_traj)


plot_traj <- selected_traj[
  t_rel < FIG_TREL_WINDOW & percentile %in% FIG_PERCENTILES & variable=="Dtemp",
  .(gas, variable, t_rel, percentile, value)
]

plot_line <- plot_traj[percentile == 0.5]
plot_ribbon <- dcast(
  plot_traj[percentile %in% c(FIG_INNER_RIBBON, FIG_OUTER_RIBBON)],
  gas + variable + t_rel ~ percentile,
  value.var = "value"
)

temp <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = plot_ribbon,
    ggplot2::aes(x = t_rel, ymin = `0.25`, ymax = `0.75`, fill = gas),
    alpha = 0.25,
    color = NA
  ) +
  ggplot2::geom_ribbon(
    data = plot_ribbon,
    ggplot2::aes(x = t_rel, ymin = `0.05`, ymax = `0.95`, fill = gas),
    alpha = 0.1,
    color = NA
  ) +
  ggplot2::geom_line(
    data = plot_line,
    ggplot2::aes(x = t_rel, y = value, color = gas),
    linewidth = 0.8
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  ggpubr::theme_pubr() + xlab("Years from pulse")  + 
  ylab("Temperature variation [°C]") + 
  labs(y = expression("Temperature variation [K]")) +
  theme(legend.position = "none")



plot_traj_f <- selected_traj[
  t_rel < FIG_TREL_WINDOW & percentile %in% FIG_PERCENTILES & variable=="Dforc",
  .(gas, variable, t_rel, percentile, value)
]

plot_line_f <- plot_traj_f[percentile == 0.5]
plot_ribbon_f <- dcast(
  plot_traj_f[percentile %in% c(FIG_INNER_RIBBON, FIG_OUTER_RIBBON)],
  gas + variable + t_rel ~ percentile,
  value.var = "value"
)

forc <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = plot_ribbon_f,
    ggplot2::aes(x = t_rel, ymin = `0.25`, ymax = `0.75`, fill = gas),
    alpha = 0.25,
    color = NA
  ) +
  ggplot2::geom_ribbon(
    data = plot_ribbon_f,
    ggplot2::aes(x = t_rel, ymin = `0.05`, ymax = `0.95`, fill = gas),
    alpha = 0.1,
    color = NA
  ) +
  ggplot2::geom_line(
    data = plot_line_f,
    ggplot2::aes(x = t_rel, y = value, color = gas),
    linewidth = 0.8
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  ggpubr::theme_pubr() + xlab("Years from pulse") + ylab(expression("Forcing variation [W" * m^-2 * "]")) + theme(legend.position = "none")


require(patchwork)
fig1 <- forc + temp + plot_layout(axis_titles = "collect")
ggsave("fig_1.png",fig1,width=12,height=6,dpi=300)


library(cowplot)

p1 <- cowplot::ggdraw() +
  draw_image("../figs/heatmaps/DGKB-chr7:14939901-14941164-heatmap.png")

p2 <-  cowplot::ggdraw() +
  draw_image("../figs/heatmaps/HNRNPM-chr19:8449557-8449735-heatmap.png")


out <- cowplot::plot_grid(p1, p2, labels = c("A", "B"))

cowplot::save_plot(plot = out, "../figs/peristentAndSwitch.png",
                     base_width = 6, base_height=3)

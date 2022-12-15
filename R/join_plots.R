library(cowplot)

set1 <- list.files("../Figures/", pattern = "*_ess*", full=T)

a <- cowplot::plot_grid(
  cowplot::ggdraw() + draw_image(set1[1]),
  cowplot::ggdraw() + draw_image(set1[2]),
  cowplot::ggdraw() + draw_image(set1[3]),
  nrow=2,
  labels = c("M1", "M2", "M3")
)

cowplot::save_plot(filename = "../Figures/ESS.png", a, base_width = 6)



set1 <- list.files("../Figures/", pattern = "*_rhat*", full=T)

a <- cowplot::plot_grid(
  cowplot::ggdraw() + draw_image(set1[1]),
  cowplot::ggdraw() + draw_image(set1[2]),
  cowplot::ggdraw() + draw_image(set1[3]),
  nrow=2,
  labels = c("M1", "M2", "M3")
)

cowplot::save_plot(filename = "../Figures/RHat.png", a, base_width = 6)





# Confidence bounds on parameters -----------------------------------------
set1 <- list.files("../Figures/", pattern = "*vs_adjusted.png", full=T)

a <- cowplot::plot_grid(
  cowplot::ggdraw() + draw_image(set1[1]),
  cowplot::ggdraw() + draw_image(set1[2]),
  cowplot::ggdraw() + draw_image(set1[3]),
  ncol=2,
  labels = c("M1", "M2", "M3")
)

cowplot::save_plot(filename = "../Figures/RealVsPred.png", a, base_width = 6)

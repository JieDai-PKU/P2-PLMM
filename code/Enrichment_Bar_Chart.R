library(ggplot2)

ggplot(
  gsea_go,
  aes(x = fct_reorder(pathway, NES),
      y = NES,
      fill = NES)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0) +
  labs(
    x = NULL,
    y = "NES",
    title = title) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"))
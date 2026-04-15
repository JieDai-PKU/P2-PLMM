library(tidyverse)
library(openxlsx)
library(ggrepel)
library(ggtext)

deg_results <- read.xlsx("input.xlsx", check.names = FALSE)

deg_results <- deg_results %>%
  mutate(
    Label = case_when(
      abs(logFC) >= 1 & P.Value < 0.05 & logFC > 0 ~ "Up",
      abs(logFC) >= 1 & P.Value < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

ggplot(deg_results, aes(x = logFC, y = -log10(P.Value),color = Label)) +
  geom_point(alpha =0.7,aes(size=-log10(P.Value))) +
  scale_color_manual(values = c("Up"="#f09285",
                                "Down"="#94d5e5",
                                "Not significant"="#d7d7d7")) +
  scale_size_continuous(range = c(1,3)) +
  geom_vline(xintercept = c(-1,1),
             linetype ="dashed", color ="grey50") +
  geom_hline(yintercept = -log10(0.05), linetype ="dashed",color ="grey50") +
  guides(size=guide_legend(order =1)) +
  labs(x ="Log2(Fold Change)",
       y ="-Log10(P-value)",color ="Label") +
  theme_test() +
  theme(axis.text=element_text(color="black"),
        axis.title = element_markdown())

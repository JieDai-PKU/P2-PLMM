library(tidyverse)

dat<-read.table("input.txt",header=T,sep="\t")
dat$name <- sub("\\..*", "", dat$name)

dat$group<-NA
dat[dat$aaSeqCDR3>0.001,"group"]<-"Hyperexpanded"
dat[dat$aaSeqCDR3<=0.001,"group"]<-"Hypoexpanded"

group<-read.table("group.txt",header = T)
group$name <- sub("\\..*", "", group$name)

expansion_summary <- dat %>%
  group_by(name,group) %>%
  summarise(total_fraction = sum(aaSeqCDR3),.groups = "drop") 

dat2<-as.data.frame(expansion_summary)
df<-merge(dat2,group)

df$name<-factor(df$name,levels=c("list_all_the_sample_names_one_by_one"))

p <- ggplot(df, aes(x = name, y = total_fraction, fill = group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(
    name = "Expansion Level",
    values = c(
      "Hyperexpanded" = "#FF7F00",        
      "Hypoexpanded" = "#4DAF4A" ),
    labels = c(
      "Hyperexpanded (0.001 < X <= 1)",
      "Hypoexpanded (0 < X <= 0.001)"
    )
  ) +
  scale_y_continuous(
    name = "Fraction Represented",
    #limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    title = "TCR Clonal Ratio"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold",angle = 90),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    legend.position = "right",
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

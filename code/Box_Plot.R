library(tidyverse)
library(reshape2)
library(ggpubr)

data<-read.table("immune_box_input.txt",header=T,sep="\t",row.names=1,check.names=F)

df3<-melt(data,id.vars=c("Group"))

df3 <- df3 %>%
  mutate(Group = case_when(
    Group == "Pre_nonresponse" ~ "Pre-NR",
    Group == "Post_nonresponse" ~ "Post-NR",
    Group == "Pre_response" ~ "Pre-R",
    Group == "Post_response" ~ "Post-R",
    TRUE ~ as.character(Group)
  )) %>%
  mutate(Group = factor(Group, levels = c("Pre-NR", "Post-NR", "Pre-R", "Post-R")))


ggplot(data=df3, aes(x = Group, y = value,
                     color=Group)) +
  geom_boxplot(alpha =0.5,size=1,outlier.shape = NA)+
  geom_point(position=position_jitter(width=0.2),alpha=0.8,size=3)+
  stat_compare_means(method = "wilcox.test",paired = F, 
                     comparisons=list(c("Pre-NR","Post-NR"),
                                      c("Pre-R","Post-R"),
                                      c("Pre-NR","Pre-R"),
                                      c("Post-NR","Post-R"),
                                      c("Post-NR","Pre-R")),
                     label = "p = {sprintf('%.3f', p)}")+
  facet_wrap(~variable, scales = "free_y",nrow=4)+
  theme_bw() +theme(strip.background = element_blank(), strip.placement = "outside") +
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top')+
  labs(y='Expression')
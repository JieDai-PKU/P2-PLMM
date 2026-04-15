library(ggplot2)
library(ggpubr)

data<-read.table("immune_mean_input.txt",header=T,sep="\t")
data$Timepoint<-factor(data$Timepoint,levels=c("Pre","Post"))

ggplot(data=data, aes(x = Timepoint, y = Value,
                      color=Timepoint) )+
  stat_compare_means(method = "wilcox.test",paired = T, 
                     comparisons=list(c("Pre","Post")))+
  geom_jitter(size=3,aes(fill=Group,color=Group),
              show.legend=F,
              shape=21,
              position = position_dodge(0))+
  geom_line(aes(group = ID2,color = Group), 
            lwd = 0.5,position = position_dodge(0))+ 
  scale_fill_manual(limits=c("Non-responder","Responder"), 
                    values=c("#762221","#214F7F"))+
  scale_color_manual(limits=c("Non-responder","Responder"), 
                     values=c("#762221","#214F7F"))+
  theme_bw() + 
  theme(panel.grid =element_blank(),
        axis.text = element_text(size = 10,colour = "black"),
        axis.title.x = element_blank(),
        legend.position = 'top')+
  labs(y='expression')
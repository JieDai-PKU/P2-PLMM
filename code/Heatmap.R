library(pheatmap)

data<-read.table("immune_input.txt",header=T,sep="\t", row.names=1,check.names=F)

data$Group<-factor(data$Group, levels=c("Post-NR", "Pre-NR",  "Pre-R", "Post-R"))

ann_colors = list(
  Group = c(`Pre-NR`="#D95F02", 
            `Pre-R`="#33B3CC",
            `Post-R`="#1B9E77",
            `Post-NR`="firebrick"))

pheatmap(t(data[,-c(1:4)]),
         cluster_rows =F,
         cluster_cols =F,
         annotation_col=data[,1:4],
         color = colorRampPalette(c("#09386C","#6BA2C5", "white", "#B84941","#6A011F"))(100),
         border_color = "white",
         annotation_colors =ann_colors,
         gaps_col =c(30,41,72),
         gaps_row =c(12,20),
         show_colnames=F)
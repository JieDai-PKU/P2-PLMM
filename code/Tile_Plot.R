group<- read.table("group.txt",header=T)

group<-melt(group,id.vars = "ID")
group$tID <- factor(group$ID,levels=group$ID)

cols <- c( "Pre"="#D95F02" , "Post"="#33B3CC" , "female"="#C90CB4" , "male"="#95EDF7" , 			"Vulvovaginal"="#CA9548" ,"Anorectal" ="#4DF2BB", "Esophageal"="#2FD1A9",
           "Head_Neck"="#82D835" , "Non-responder"="#BAD3E4" , "Responder"="#DEAAA6" ,
           "[,50)"="#1D52A1","[50,55)"="#716DB2","[55,60)"="#65C8CC","[60,65)"="#72C15A",
           "[65,70)"="","[70,75)"="#F3793B") 

p2<-ggplot(group,aes(x=ID, y=variable))+
  geom_tile(aes(fill=value),color ="white",linewidth =0.5)+scale_fill_manual(values = cols, breaks = names(cols))+theme_bw()+
  theme(axis.text.x = element_text(size =12,angle =90,vjust =0.5,hjust =1),color ="black")+xlab("")
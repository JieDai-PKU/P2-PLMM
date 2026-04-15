library(GetoptLong)

phenofile<-'NA'
groupfile<-'NA'
samplefile<-'NA'
mutfile<-'NA'
cut<-0
outputfile<-'NA'
width<-10
height<-0
atype<-'NA'
btype<-'NA'
colname<-'F'
classfile<-NA

GetoptLong(
	"mutfile=s", "Mutation data file",
	"cut=i", "Cutoff of the number of samples per mutated gene, default (0)",
        "phenofile=s", "Phenotype data file",
        "groupfile=s", "Group data file",
        "samplefile=s", "Sample list file for column order",
        "outputfile=s", "Output file, default (mutfile.pdf)",
        "width=f", "Width of figure",
        "height=f", "Height of figure",
        "atype=s", "Cosmic tumor type", 
        "btype=s", "IntOGene tumor type",
        "colname=s", "Show column names (T or F)",
	"verbose", "Print messages",
	"classfile=s","split file"
	
)
gcut<-cut
if(outputfile == 'NA'){
  outputfile = paste(mutfile,'pdf',sep=".")
}

if(length(grep("\\.pdf$",outputfile))==0){
  outputfile = paste(outputfile,'pdf',sep=".")
}

ShowColName = T
if(colname == 'T'){
	ShowColName = TRUE
}


if(width==0){
  width=10
}

cmd.args <- commandArgs()
m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
this.dir <- dirname(regmatches(cmd.args, m))
cosmicfile<-paste(this.dir,'cosmic_tumor_type',sep="/")
igfile<-paste(this.dir,"intogen_tumor_type",sep="/")
outputpdf<-outputfile

##########
cat(paste("mut_file:",mutfile,"\n"))
cat(paste("pheno_file:",phenofile,"\n"))
cat(paste("group_file:",groupfile,"\n"))
cat(paste("output_file:",outputpdf,"\n"))
cat(paste("show column name:",ShowColName,"\n"))
##########

#### prepare input data ####
mutinput<-read.table(file=mutfile,as.is=TRUE,sep="\t",quote ="")
sampleNum<-length(unique(mutinput$V1))
cat(paste("Sample Number:",sampleNum,"\n"))


mutinput$type<-'Others'
mutinput$type[mutinput$V4=='Nonsense_Mutation']<-'Nonsense_Mutation'
mutinput$type[mutinput$V4=='Missense_Mutation']<-'Missense_Mutation'
mutinput$type[mutinput$V4=='Frame_Shift_Del']<-'Frame_Shift_Del'
mutinput$type[mutinput$V4=='Frame_Shift_Ins']<-'Frame_Shift_Ins'
mutinput$type[mutinput$V4=='In_Frame_Del']<-'In_Frame_Del'
mutinput$type[mutinput$V4=='In_Frame_Ins']<-'In_Frame_Ins'
mutinput$type[mutinput$V4=='Translation_Start_Site']<-'Translation_Start_Site'
mutinput$type[mutinput$V4=='Nonstop_Mutation']<-'Nonstop_Mutation'
mutinput$type[mutinput$V4=='CNV']<-'CNV'
mutinput$type[mutinput$V4=='SV']<-'SV'
mutinput$type[mutinput$V4=='Splice_Site']<-'Splice_Site'
mutinput$type[mutinput$V4=='Loss']<-'Loss'
mutinput$type[mutinput$V4=='Gain']<-'Gain'


genepersample<-table(mutinput$V1,mutinput$V2)
SampleGeneCount<-rowSums(genepersample)
samples<-sort(unique(mutinput$V1))
genes<-sort(unique(mutinput$V2))
tsg<-table(mutinput$V1,mutinput$V2)
smutcount<-colSums(tsg>0)
smutcountper<-smutcount*100/length(samples)
smutcountper<-smutcountper[smutcount>=gcut]
gsc<-names(smutcount[smutcount>=gcut])
inputgsc<-subset(mutinput, is.element(V2, gsc))
gene2<-sort(unique(inputgsc$V2))
if(height == 0){
  height = width*length(gene2)/length(samples)
}

mmat<-matrix(" ",length(gene2),length(samples))
for(i in 1:length(gene2)){
   for(j in 1:length(samples)){
     if(dim(subset(inputgsc,V1==samples[j] & V2==gene2[i]))[1]>0){
          jmut<-paste(unique(subset(inputgsc,V1==samples[j] & V2==gene2[i])$type),collapse = ';')
          jmut<-paste(jmut,';',sep="")
          mmat[i,j]<-jmut
     }
  }
}
rownames(mmat)<-gene2
colnames(mmat)<-samples
if (!is.na(classfile)){
	classorder<-read.table(classfile,header=T,sep="\t")
	class<-classorder[match(gene2,classorder[[1]]),2]
	print(class)
}




#######################################
library(ComplexHeatmap)

options(width = 100)

myoncoPrint=paste(this.dir,'oncoPrint.R',sep="/")
source(myoncoPrint)

mutcol = c("Missense_Mutation" = "dodgerblue4", "Nonsense_Mutation" = "red", "Frame_Shift_Del" = "cyan", "Frame_Shift_Ins" = "purple", "In_Frame_Del" = "orange","In_Frame_Ins"="#6D22E6","Translation_Start_Site"="#EDC3D2","Nonstop_Mutation"="#CFD19E","CNV" = "coral", "SV" = "brown","Others" = "green","Splice_Site"="#389EB2","Loss"="#AE6642","Gain"="#2A5CAA")

mut_alter_fun = list(
        background = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "white", col = "#CCCCCC"))##CCCCCC
        },
        "Missense_Mutation" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "dodgerblue4", col = NA))
        },
        "Nonsense_Mutation" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
        },
        "Frame_Shift_Del" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "cyan", col = NA))
        },
        "Frame_Shift_Ins" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "purple", col = NA))
        },
        "In_Frame_Del" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "orange", col = NA))
        },
        "In_Frame_Ins" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#6D22E6", col = NA))
        },
        "Translation_Start_Site" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#EDC3D2", col = NA))
        },
        "Nonstop_Mutation" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CFD19E", col = NA))
        },
        "Splice_Site" = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#389EB2", col = NA))
        },
        "Loss" = function(x, y, w, h) {
		r= unit(0.5, "mm")
                grid.circle(x, y, r,w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#AE6642", col = NA))
        },

        "Gain" = function(x, y, w, h) {
		r= unit(0.5, "mm")
                grid.circle(x, y, r,w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#2A5CAA", col = NA))
        },

      CNV = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "coral", col = NA))
        },
      SV = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "brown", col = NA))
        },
      Others = function(x, y, w, h) {
                grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "green", col = NA))
        }
)





ht<-oncoPrint(mmat, get_type = function(x) strsplit(x, ";")[[1]],
		 alter_fun = mut_alter_fun, col = mutcol, 
		 heatmap_legend_param = list(title = "Mutation"),
             show_row_barplot = TRUE, show_column_names = FALSE)

print("###########$$$$$$$")


htlist<-ht@ht_list
globalcorder<-htlist[[names(htlist[2])]]@column_order
if(samplefile != 'NA'){
    samplelist<-read.table(file=samplefile,colClasses="character")
    globalcorder<-samplelist$V1
}
SortSmutCountPer<-names(sort(smutcountper,decreasing=TRUE))



if(phenofile != 'NA'){
	pheno<-read.table(file=phenofile,head=TRUE,sep="\t",as.is=T,quote ="",check.names=F)
        phenoSample<-length(unique(pheno$ID))
        cat(paste("Pheno Sample Number:", phenoSample,"\n"))
        setdiffPheno<-length(setdiff(unique(pheno$ID),unique(mutinput$V1)))
	cat(paste("Setdiff Pheno:", setdiffPheno,"\n"))
   	pheno<-pheno[with(pheno,order(ID)),]
	colname<-colnames(pheno)
   	pheno<-as.data.frame(pheno[,-1])
      	colnames(pheno)<-colname[-1]
   	rbcol<-c("#97BEBD","#D7B0C2","#4C6E78","#C4D6DA","#BA4A70")
   	plist<-list()
	pheno_column = HeatmapAnnotation(df = pheno,gp = gpar(col=NA),gap = unit(0, "points"),annotation_name_side = "left",border = T,annotation_name_gp = gpar(fontsize = 9),
simple_anno_size = unit(5, "mm"),show_annotation_name = TRUE)#col=phycol
}

topann <- HeatmapAnnotation(
tmb = anno_barplot(x = pheno$TMB, gp = gpar(fill = "blue"),border = TRUE,axis = TRUE, axis_param = list(side = "right"), height = unit(2, "cm")),
tnb = anno_barplot(x = pheno$TNB, gp = gpar(fill = "blue"),border = TRUE,axis = TRUE, axis_param = list(side = "right"), height = unit(2, "cm")),

show_annotation_name = TRUE, 
annotation_name_side = "left"  )




if(groupfile != 'NA'){

group<-read.table(file=groupfile,sep="\t",colClasses="character",head=TRUE)
print("==================##")
pct_heatmap_in<-matrix(nrow=length(SortSmutCountPer),ncol=5)
sig_pct<-matrix(nrow=length(SortSmutCountPer),ncol=2)
sample_name<-unique(mutinput$V1)
sample_length<-length(unique(mutinput$V1))


group_id<-factor(group$Group,levels=c("Pre_Responder","Pre_Non-responder","Post_Responder","Post_Non-responder"))

group_id1<-levels(group_id)[1]
group_id1_len<-length(which(group[group$Group==group_id1,1] %in% sample_name))
group_id2<-levels(group_id)[2]
group_id2_len<-length(which(group[group$Group==group_id2,1] %in% sample_name))
group_id3<-levels(group_id)[3]
group_id3_len<-length(which(group[group$Group==group_id3,1] %in% sample_name))
group_id4<-levels(group_id)[4]
group_id4_len<-length(which(group[group$Group==group_id4,1] %in% sample_name))

SortSmutCountPer2<-sort(SortSmutCountPer)
for(i in 1:length(SortSmutCountPer2)){

	group_name1<-group[group$Group==group_id1,1]
	group_name1_num<-length(unique(inputgsc[inputgsc$V1 %in% group_name1 & inputgsc$V2==SortSmutCountPer2[i],1]))
	pct_heatmap_in[i,1]<-round(group_name1_num*100/group_id1_len,1)

	group_name2<-group[group$Group==group_id2,1]
	group_name2_num<-length(unique(inputgsc[inputgsc$V1 %in% group_name2 & inputgsc$V2==SortSmutCountPer2[i],1]))
	pct_heatmap_in[i,2]<-round(group_name2_num*100/group_id2_len,1)

        group_name3<-group[group$Group==group_id3,1]
        group_name3_num<-length(unique(inputgsc[inputgsc$V1 %in% group_name3 & inputgsc$V2==SortSmutCountPer2[i],1]))
        pct_heatmap_in[i,3]<-round(group_name3_num*100/group_id3_len,1)

        group_name4<-group[group$Group==group_id4,1]
        group_name4_num<-length(unique(inputgsc[inputgsc$V1 %in% group_name4 & inputgsc$V2==SortSmutCountPer2[i],1]))
        pct_heatmap_in[i,4]<-round(group_name4_num*100/group_id4_len,1)

	group_name_num<-length(unique(inputgsc[inputgsc$V2==SortSmutCountPer2[i],1]))
	pct_heatmap_in[i,5]<-round(group_name_num*100/sample_length,1)
}


colnames(pct_heatmap_in)<-c(group_id1,group_id2,group_id3,group_id4,"Cohort")
rownames(pct_heatmap_in)<-SortSmutCountPer2


print("+===========+==========+============+==========\n")


print(head(pct_heatmap_in))
left_heatmap<-Heatmap(pct_heatmap_in,col= colorRampPalette(c("white", "red"))(200),name="Frequency",width = unit(2, "cm"),cluster_rows = FALSE,cluster_columns = FALSE,rect_gp = gpar(col = "grey"),
	column_names_side = "bottom",row_order=match(SortSmutCountPer,SortSmutCountPer2),show_row_names=F,column_title = "Gene Frequency(%)",
	column_title_gp = gpar(fontsize = 8),column_title_side = "top",column_names_gp = gpar(fontsize = 6),
	cell_fun = function(j, i, x, y, width, height, fill) {
        	grid.text(sprintf("%.1f", pct_heatmap_in[i, j]), x, y, gp = gpar(fontsize = 5))
	}
)


bardata<-pct_heatmap_in[,c(1,2,3,4)]
bardata_rowsum<-rowSums(bardata)
bardata[,1]<-bardata[,1]/bardata_rowsum
bardata[,2]<-bardata[,2]/bardata_rowsum
bardata[,3]<-bardata[,3]/bardata_rowsum
bardata[,4]<-bardata[,4]/bardata_rowsum

print("++++++++++++++++")
print(bardata_rowsum)
print("++++++++++++++++")


bardata<-bardata[match(sort(SortSmutCountPer),rownames(bardata)),]
print(head(bardata))
  groupSample<-length(unique(group$ID))
  cat(paste("Group Sample Number:",groupSample,"\n"))
  setdiffGroup<-length(setdiff(unique(group$ID),unique(mutinput$V1)))
  cat(paste("Setdiff Group:", setdiffGroup,"\n"))
  groupname<-colnames(group)[2] 
  gclass<-sort(unique(group[,2]))
  tabgroup<-table(group[,2])
  gcol<-list()
  groupdf<-data.frame(group[,2][match(samples,group[,1])])
 groupdf2<-data.frame(group[,3][match(samples,group[,1])])
  names(groupdf)<-groupname
  color_all<-c("#F78F21","#009574","#FFDF72","#832C27","#8CC73C")
  selcol<-color_all[1:length(gclass)]  
  names(selcol)<-gclass
  gcol[[groupname]]=selcol
  corder<-c()
  for (gc in gclass){
     gcid<-group$ID[group[,2]==gc]
     gcid<-sort(gcid)
     ginputgsc<-subset(inputgsc, is.element(V1, gcid))
     gsamples<-sort(gcid)
     ggene2<-sort(unique(ginputgsc$V2))
     gmmat<-matrix(" ",length(ggene2),length(gsamples))
     for(i in 1:length(ggene2)){
   		for(j in 1:length(gsamples)){
     			if(dim(subset(ginputgsc,V1==gsamples[j] & V2==ggene2[i]))[1]>0){
          			jmut<-paste(unique(subset(ginputgsc,V1==gsamples[j] & V2==ggene2[i])$type),collapse = ';')
          			jmut<-paste(jmut,';',sep="")
          			gmmat[i,j]<-jmut
     			}
  		}
	}
      rownames(gmmat)<-ggene2
      colnames(gmmat)<-gsamples
      ght<-oncoPrint(gmmat, get_type = function(x) strsplit(x, ";")[[1]],
	              alter_fun = mut_alter_fun, col = mutcol, 
	              heatmap_legend_param = list(title = "Mutation", at = c("Missense", "Nonsense", "Splicing","Frame_shift","In Frame Indel","Others"), 
		        labels = c("Missense", "Nonsense", "Splicing","Frame shift","In-frame indel","Others")),
                     show_row_barplot = FALSE, show_column_names = TRUE,top_annotation = NULL)
      ghtlist<-ght@ht_list
      gcorder<-ghtlist[[names(ghtlist[2])]]@column_order  
      corder<-c(corder,gsamples[gcorder])
  }
  if(samplefile == 'NA'){
	globalcorder<-corder
  }


group_column = HeatmapAnnotation(
 
`TMB (mut/MB)` = anno_barplot(x = pheno$TMB, gp = gpar(fill = "red"),border = TRUE,axis = TRUE, axis_param = list(side = "left"), height = unit(2, "cm")),

Group = groupdf$Group,col = gcol,gp = gpar(col = "white",fontsize = 4)
) 
}

#print(gcol)


MutGenePer = rowAnnotation(percent = row_anno_barplot(bardata,gp = gpar(fill = selcol),axis_param=list(at=c(0,0.5,1.0),labels=c(0,0.5,1.0)),axis = TRUE),
                           annotation_width = unit(3, "cm"),show_annotation_name = F,annotation_name_rot = 0,annotation_name_offset = unit(15, "mm"))
row_name=rowAnnotation(text = row_anno_text(rownames(mmat),gp = gpar(col = "black",face="italic",family="Times",fontsize = 6),location = 1,just = 1),width = max_text_width(rownames(mmat)))

globalcorder<-c("list_all_the_sample_names_one_by_one")


SortSmutCountPer<-c("BRAF","NRAS","KIT","NF1","SF3B1","TP53","PTEN","NBPF1","FAT3","FSIP2","RGS7","BRCA2","TERT","ATRX","ARID2","HDAC4","KMT2A","KMT2D","MGMT","CDKN2A","CCND1","CDK4","CDK6","DHX9","ANXA10","NOTCH2","MITF","AXL","NGFR","RAB27A","HLA-B","HLA-C","HLA-DQA1","HLA-DRB1","NLRC5","EP300","TAPBP","IFI30","STAT2","JAK1","B","C","D","E")
print(head(SortSmutCountPer))
#######################################################################################################

#================================================================================
if(phenofile != 'NA' & groupfile != 'NA'){
   pdf(file=outputpdf,width=width,height=height)
if(!is.na(classfile)){

	ht<-oncoPrint(mmat, get_type = function(x) strsplit(x, ";")[[1]],
		alter_fun = mut_alter_fun, col = mutcol,show_pct = F,row_names_side="left", show_row_names =T,
		heatmap_legend_param = list(title = "Mutation"),column_order = globalcorder,
              	show_row_barplot = T,pct_gp = gpar(fontsize = 4),top_annotation =topann,column_split = c(rep(c("C", "D"), 80),rep("E",9)),
		 show_column_names = ShowColName, bottom_annotation = pheno_column,split = class,row_names_gp =gpar(fontsize = .001),
		gap=unit(3, "mm"),row_order = SortSmutCountPer)
}else{
print("+++++++++++++++++++++++++++++")
	ht<-oncoPrint(mmat, get_type = function(x) strsplit(x, ";")[[1]],
                alter_fun = mut_alter_fun, col = mutcol,show_pct = F,row_names_side="right",show_row_names =T,column_split =groupdf2,
                heatmap_legend_param = list(title = "Mutation"),column_order = globalcorder,
                show_row_barplot = F,top_annotation =group_column,pct_gp =gpar(fontsize = 4),
		show_column_names = ShowColName, bottom_annotation = pheno_column,
                gap=unit(3, "mm"),row_order = SortSmutCountPer)


}
      if(atype != 'NA' & btype != 'NA'){
		draw(MutGenePer+ht+ciha,padding = unit(c(20, 2, 2, 2), "mm"))
	}
	if(atype != 'NA' & btype == 'NA'){
		draw(MutGenePer+ht+cha,padding = unit(c(20, 2, 2, 2), "mm"))
	}
	if(atype == 'NA' & btype != 'NA'){
		draw(MutGenePer+ht+iha,padding = unit(c(20, 2, 2, 2), "mm"))
	}
	if(atype == 'NA' & btype == 'NA'){
print("====================-------=\n")
		
draw(row_name+ht+left_heatmap+MutGenePer,merge_legend = TRUE, heatmap_legend_side = "right",annotation_legend_side = "right")


	}

   dev.off()
}

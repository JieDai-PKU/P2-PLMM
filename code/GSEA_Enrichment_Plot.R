library(GseaVis)

gseaNb(object=GSEA_result, 
       geneSetID='pathway_name', 
       subPlot=3, 
       addPval=TRUE, 
       pvalX=0.9, 
       pvalY=0.9, 
       pCol='black', 
       pvalSize=4)
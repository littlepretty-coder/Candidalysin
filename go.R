library(clusterProfiler)
library(org.Hs.eg.db)
info = read.delim('E:/projects/王秀丽/results/data/mouse.xls',check.names = F)
info = info[,c(1,2)]
#merged_df <- merge(data_up1, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)
info$gene_name =toupper(info$gene_name)

use = read.delim('E:/projects/王秀丽/All_Fungi.gene_info/All_Fungi.gene_info',sep='\t')



#data_up1$gene = rownames(data_up1)
data = read.delim('E:/projects/王秀丽/results/fungi/mouse/3h/res_1477_0_vs_3.txt',sep='\t')
data_up1 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
data_up1$gene = rownames(data_up1)
merged_df1 <- merge(data_up1, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)
gene1 = merged_df1$gene_id
#gene1 = data_up1$gene
gene1= bitr(gene1, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


data = read.delim('E:/projects/王秀丽/results/fungi/mouse/3h/res_2174_0_vs_3.txt',sep='\t')
data_up2 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
#gene2 = rownames(data_up2)
data_up2$gene = rownames(data_up2)
merged_df2 <- merge(data_up2, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)
gene2 = merged_df2$gene_id
gene2= bitr(gene2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


data = read.delim('E:/projects/王秀丽/results/fungi/mouse/3h/res_2059_0_vs_3.txt',sep='\t')
data_up3 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
data_up3$gene = rownames(data_up3)
merged_df3 <- merge(data_up3, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)
gene3 = merged_df3$gene_id
#gene3= bitr(gene3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# data = read.delim('E:/projects/王秀丽/results/fungi/human/3h/res_Ca_0_vs_3.txt',sep='\t')
# data_up4 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
# data_up4$gene = rownames(data_up4)
# merged_df4 <- merge(data_up4, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)
# gene4 = merged_df4$gene_id


library(GO.db)
goterms <- Term(GOTERM)
a<-as.data.frame(goterms)
go_names<-cbind(row.names(a),a)



yeast_go<-read.table("E:/projects/王秀丽/go.txt",sep='\t',header=T,check.names = F)

#yeast_go<-yeast_go[!(yeast_go$V1 %in% c("GO:0003674","GO:0008150","GO:0005575")),]


MF_universe<-yeast_go[yeast_go$`GO domain` =="molecular_function",][,c(2,1)]
BP_universe<-yeast_go[yeast_go$`GO domain` =="biological_process",][,c(2,1)]
CC_universe<-yeast_go[yeast_go$`GO domain` =="cellular_component",][,c(2,1)]

BP_universe$`NCBI gene (formerly Entrezgene) ID` = as.character(BP_universe$`NCBI gene (formerly Entrezgene) ID`)






ego<-enricher(gene1$ENTREZID, pvalueCutoff = 1, 
              pAdjustMethod = "BH", universe = BP_universe$`NCBI gene (formerly Entrezgene) ID`,minGSSize = 2, 
              maxGSSize = 10000, TERM2GENE = BP_universe,TERM2NAME = go_names)
library(ggplot2)

p<-dotplot(ego,showCategory = 10)
p
ggsave(sprintf("./%s_%s.png",'2059','BP'),
                  units="in", width=10, heigh=7, dpi=600)
write.table(ego,sprintf("E:/projects/王秀丽/results/mouse/3h/%s_%s.txt",'1477','MF'), sep = "\t")
  















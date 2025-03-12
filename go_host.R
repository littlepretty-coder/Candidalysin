library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
info = read.delim('E:/projects/王秀丽/results/data/mouse.xls',check.names = F)
info = info[,c(1,2)]
#merged_df <- merge(data_up1, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)

#data_up1$gene = rownames(data_up1)
data = read.delim('E:/projects/王秀丽/results/mouse/6h/res_mouse_1477_0_vs_6.txt',sep='\t')
data = na.omit(data)
data_up1 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene1 = rownames(data_up1)

#gene1= bitr(gene1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


data = read.delim('E:/projects/王秀丽/results/mouse/6h/res_mouse_2174_0_vs_6.txt',sep='\t')
data = na.omit(data)
data_up2 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
#gene2 = rownames(data_up2)
gene2 = rownames(data_up2)

#gene2= bitr(gene2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


data = read.delim('E:/projects/王秀丽/results/mouse/6h/res_mouse_2059_0_vs_6.txt',sep='\t')
data = na.omit(data)
data_up3 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene3 = rownames(data_up3)

# #gene3= bitr(gene3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# data = read.delim('E:/projects/王秀丽/results/human/ca/res_human_Ca_0_vs_24.txt',sep='\t')
# data = na.omit(data)
# data_up4 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
# gene4 = rownames(data_up4)



ego<-enrichGO(gene4, 
              OrgDb = org.Hs.eg.db,
              keyType = "ENSEMBL",
              pvalueCutoff = 1, 
              pAdjustMethod = "BH",
              ont = 'BP')
# library(ggplot2)
# p<-dotplot(ego,showCategory = 200)
# p
# # ggsave(sprintf("./%s_%s.png",'2074','BP'),
#                   units="in", width=10, heigh=7, dpi=600)
write.table(ego@result,sprintf("E:/projects/王秀丽/results/human/6h/%s_%s.txt",'ca','BP'), sep = "\t")
  















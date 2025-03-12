library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
info = read.delim('E:/projects/王秀丽/results/data/human.xls',check.names = F)
info = info[,c(1,2)]
#merged_df <- merge(data_up1, info, by.x = 'gene', by.y = 'gene_name', all.x = TRUE)

#data_up1$gene = rownames(data_up1)
data = read.delim('E:/projects/王秀丽/results/human/3h/res_human_1477_0_vs_3.txt',sep='\t')
data = na.omit(data)
data_up1 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene1 = rownames(data_up1)

#gene1= bitr(gene1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


data = read.delim('E:/projects/王秀丽/results/human/3h/res_human_2174_0_vs_3.txt',sep='\t')
data = na.omit(data)
data_up2 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
#gene2 = rownames(data_up2)
gene2 = rownames(data_up2)

#gene2= bitr(gene2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


data = read.delim('E:/projects/王秀丽/results/human/3h/res_human_2059_0_vs_3.txt',sep='\t')
data = na.omit(data)
data_up3 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene3 = rownames(data_up3)

# #gene3= bitr(gene3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

data = read.delim('E:/projects/王秀丽/results/human/ca/res_human_Ca_0_vs_3.txt',sep='\t')
data = na.omit(data)
data_up4 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene4 = rownames(data_up4)



gene4_entrez <- bitr(gene4, 
                     fromType = "ENSEMBL",  # 输入的基因类型
                     toType = "ENTREZID",  # 转换为 Entrez Gene ID
                     OrgDb = org.Hs.eg.db) # 使用人类注释数据库


KEGG_enrich = enrichKEGG(gene = gene4_entrez$ENTREZID, #即待富集的基因列表
                         keyType = "kegg",
                         pAdjustMethod = 'fdr',  #指定p值校正方法
                         organism= "hsa",  #hsa，可根据你自己要研究的物种更改，可在https://www.kegg.jp/brite/br08611中寻找
                         qvalueCutoff = 1, #指定 p 值阈值（可指定 1 以输出全部）
                         pvalueCutoff=1)



# p<-dotplot(ego,showCategory = 200)
# p
# # ggsave(sprintf("./%s_%s.png",'2074','BP'),
#                   units="in", width=10, heigh=7, dpi=600)
write.table(KEGG_enrich@result,sprintf("E:/projects/王秀丽/results/human/ca/kegg_3h_%s.txt",'ca'), sep = "\t")
















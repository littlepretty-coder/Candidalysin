library(ggvenn)

# 读取数据文件
data = read.delim('E:/projects/王秀丽/3-14/results/human/1477/res_human_1477_0_vs_6.txt',sep='\t')
data_up1 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene1 = rownames(data_up1)
gene1 <- gene1[!grepl('NA',gene1)]

data = read.delim('E:/projects/王秀丽/3-14/results/human/2174/res_human_2174_0_vs_6.txt',sep='\t')
data_up2 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene2 = rownames(data_up2)
gene2 <- gene2[!grepl('NA',gene2)]


data = read.delim('E:/projects/王秀丽/3-14/results/human/2059/res_human_2059_0_vs_6.txt',sep='\t')
data_up3 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
gene3 = rownames(data_up3)
gene3 <- gene3[!grepl('NA',gene3)]
# data = read.delim('E:/projects/王秀丽/3-14/results/human/ca/res_human_Ca_0_vs_.txt',sep='\t')
# data = na.omit(data)
# data_up4 = data[(data$log2FoldChange >1.5 & data$padj <0.01),]
# gene4 = rownames(data_up4)
# gene4 <- gene4[!grepl('NA',gene4)]
library(VennDiagram)
filename <- "venn_6h_upgrade.png"
venn.diagram(
  list('1477' = gene1, '2174' = gene2,'2059'=gene3),  # 设置两个集合A和B的元素范围
  filename = filename,
  fill = brewer.pal(3, "Pastel2"),
  col = brewer.pal(3, "Pastel2"),
  lwd = 2, lty = "dashed"# 指定保存图像文件的文件名
)
#hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh

data = read.delim('E:/projects/王秀丽/3-14/results/human/1477/res_human_1477_0_vs_6.txt',sep='\t')
data_up1 = data[(data$log2FoldChange <-1.5 & data$padj <0.01),]
gene1 = rownames(data_up1)
gene1 <- gene1[!grepl('NA',gene1)]

data = read.delim('E:/projects/王秀丽/3-14/results/human/2174/res_human_2174_0_vs_6.txt',sep='\t')
data_up2 = data[(data$log2FoldChange <-1.5 & data$padj <0.01),]
gene2 = rownames(data_up2)
gene2 <- gene2[!grepl('NA',gene2)]


data = read.delim('E:/projects/王秀丽/3-14/results/human/2059/res_human_2059_0_vs_6.txt',sep='\t')
data_up3 = data[(data$log2FoldChange <-1.5 & data$padj <0.01),]
gene3 = rownames(data_up3)
gene3 <- gene3[!grepl('NA',gene3)]
# data = read.delim('E:/projects/王秀丽/3-14/results/human/ca/res_human_Ca_0_vs_3.txt',sep='\t')
# data = na.omit(data)
# data_up4 = data[(data$log2FoldChange <-1.5 & data$padj <0.01),]
# gene4 = rownames(data_up4)
# gene4 <- gene4[!grepl('NA',gene4)]
library(VennDiagram)
filename <- "venn_6h_downgrade.png"
venn.diagram(
  list('1477' = gene1, '2174' = gene2,'2059'=gene3),  # 设置两个集合A和B的元素范围
  filename = filename,
  fill = brewer.pal(3, "Pastel2"),
  col = brewer.pal(3, "Pastel2"),
  lwd = 2, lty = "dashed"# 指定保存图像文件的文件名
)

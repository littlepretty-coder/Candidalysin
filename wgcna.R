rm(list = ls())  
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) #多核读取文件
# dir.create('8.WGCNA')
# setwd('8.WGCNA')
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 

################################# 0.输入数据准备 ################################
### 读取表达矩阵
# if (T) {
#   fpkm00 <- fread("../GSE154290_RNA.FPKM.txt.gz",data.table = F)
#   ### 合并具有相同基因名的行
#   table(duplicated(fpkm00$gene)) #统计重复基因名的基因
#   gene <- fpkm00$gene ; fpkm00 <- fpkm00[,-1]
#   fpkm0 <- aggregate(fpkm00, by=list(gene), FUN=sum)
#   fpkm <- column_to_rownames(fpkm0,"Group.1")
# }
# data <- log2(fpkm+1)
# 
# ### 筛选MAD前5000的基因
# keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:5000],]

data1 = read.delim('fungi.txt',check.names = FALSE,sep='\t',row.names=1)
data2 = read.delim('host.txt',check.names = FALSE,sep='\t',row.names=1)

for (spp in c("WT","R1","R2")){
  
  #spp='R1'
  yeast<-read.delim(sprintf("./fungi_%s.txt",spp),check.names = FALSE,sep='\t',header=TRUE,row.names=1)
  yeast = yeast[,-c(1,2,3)]
  a = colnames(yeast)
  a = gsub("count", "FPKM", a)
  yeast = data1[,a]
  yeast <- log2(yeast+1)
  #print(sprintf("./%s_yeast_data.txt",spp))
  
  human_controls<-read.delim("host_control.txt",check.names = FALSE,sep='\t',header=TRUE,row.names=1)
  human_controls = human_controls[,-c(1,2,3)]
  a = colnames(human_controls)
  a = gsub("count", "FPKM", a)
  human_controls = data2[,a]
  human_controls <- log2(human_controls+1)
  
  human<-read.delim(sprintf("./host_%s.txt",spp),check.names = FALSE,sep='\t',header=TRUE,row.names=1)
  human = human[,-c(1,2,3)]
  a = colnames(human)
  a = gsub("count", "FPKM", a)
  human = data2[,a]
  human <- log2(human+1)
  
  
  shared_cols<-intersect(colnames(yeast),colnames(human))
  
  
  ### for analysing using counts and vst
  
  yeast_shared <- yeast[,shared_cols]
  human_shared <- human[,shared_cols]
  print(spp)
  
  
  
  # if (!(spp=="calb")){
  #   
  #   human_shared<-cbind.data.frame(human_controls[,c(1,2)],human_shared,human_controls[,c(23,24)])
  #   
  #   if (spp=="ctrop"){
  #     
  #     yeast<-add_column(yeast, "45C+reseq" = yeast$`45C`+yeast$`45Creseq`, .after = "45")
  #     yeast$`45C`<-NULL
  #     yeast$`45Creseq`<-NULL
  #     
  #     yeast_shared<-cbind.data.frame(yeast[,c(1,2)],yeast_shared,yeast[,c(dim(yeast)[2]-1,dim(yeast)[2])])
  #   }
  #   
  #   else{
  #     yeast_shared<-cbind.data.frame(yeast[,c(1,2)],yeast_shared,yeast[,c(dim(yeast)[2]-1,dim(yeast)[2])])
  #   }
  #   
  # }
  
  #除"D1"、"D2"、"D3"、"e7"、"e8"和"e9"列。结果将是一个包含除了这些列之外的所有其他列的新数据框。
  # if (spp=="calb"){
  #   yeast_no_D<- subset(yeast_shared, select = -c(D1,D2,D3,e7,e8,e9))
  #   human_no_D<- subset(human_shared, select = -c(D1,D2,D3,e7,e8,e9))
  # } else if (spp=="cglab"){
  #   yeast_no_D<- subset(yeast_shared, select = -c(D4,D5,D6))
  #   human_no_D<- subset(human_shared, select = -c(D4,D5,D6))
  # } else if (spp=="cpar"){
  #   yeast_no_D<- subset(yeast_shared, select = -c(D7,D8,D9))
  #   human_no_D<- subset(human_shared, select = -c(D7,D8,D9))
  # } else if (spp=="ctrop"){
  #   yeast_no_D<- subset(yeast_shared, select = -c(D10,D11,D12))
  #   human_no_D<- subset(human_shared, select = -c(D10,D11,D12))
  # }
  #修改`yeast_no_D`数据框的列名，将第1列、第2列、倒数第2列和倒数第1列的列名分别修改为"1"、"2"、"17"和"18"。
  #colnames(yeast_no_D)[c(1,2,dim(yeast_no_D)[2]-1,dim(yeast_no_D)[2])]<-c("1","2","17","18")
  yeast_no_D = yeast
  human_no_D = human
  
  
  counts_merged<-rbind.data.frame(yeast_no_D,human_no_D)
  assign(paste0(spp,"_counts_merged"),counts_merged)
  
  
  #将返回一个新的数据框，其中包含满足条件的行。条件是每一行中大于等于10的元素的总和必须大于等于`yeast_no_D`数据框列数的90%。
  yeast_no_D_filt<-yeast_no_D[order(apply(yeast_no_D,1,mad), decreasing = T)[1:5000],]  
  human_no_D_filt<-human_no_D[order(apply(human_no_D,1,mad), decreasing = T)[1:5000],]
  
  
  # ### fake model to generate vst-transformated counts (transformatio does not depend on model)
  # colData_yeast = data.frame(colnames(yeast_no_D_filt))
  # write.table(colData_yeast,'yeast_group.txt',sep='\t')
  # colData_yeast<-data.frame(a=c(rep("a",length(colnames(yeast_no_D_filt))-2),rep("b",2)))
  # pre_dds_yeast <- DESeqDataSetFromMatrix(yeast_no_D_filt, colData=colData_yeast,~a)
  # dds_counts_yeast <- DESeq(pre_dds_yeast)
  # counts_norm_yeast<-vst(dds_counts_yeast,blind = T)
  #vst_yeast<-assay(counts_norm_yeast)
  
  
  # colData_human<-data.frame(a=c(rep("a",length(colnames(human_no_D_filt))-2),rep("b",2)))
  # pre_dds_human <- DESeqDataSetFromMatrix(human_no_D_filt, colData=colData_human,~a)
  # dds_counts_human <- DESeq(pre_dds_human)
  # counts_norm_human<-vst(dds_counts_human,blind = T)
  #vst_human<-assay(counts_norm_human)
  vst_yeast = yeast_no_D_filt
  vst_human = human_no_D_filt
  all_vst<-rbind.data.frame(vst_yeast,vst_human)
  assign(paste0(spp,"_vst"),all_vst)
}


keep_data = WT_vst
### 创建datTraits，包含分组、表型等信息
################################################################################################  WT
datTraits <- data.frame(row.names = colnames(keep_data),group=colnames(keep_data))
fix(datTraits)

### 给分组加上编号
grouptype <- data.frame(group=sort(unique(datTraits$group)),
                        groupNo=1:length(unique(datTraits$group)))
fix(grouptype)
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
datTraits
#datTraits = datTraits[,-2]
### 转置
datExpr0 <- as.data.frame(t(keep_data))
############################## 1.判断数据质量 ################################

### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

### 绘制样品的系统聚类树
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait R2.pdf",width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  dev.off()
}

##若存在显著离群点；剔除掉
if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) # cutHeight根据实际情况而定
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}
save(datTraits,WT_vst,R1_vst,R2_vst,datExpr0,grouptype,gsg,keep_data,sampleTree,file='WT_input.Rdata')
### 判断数据质量 : PCA进行分组查看
#rm(list = ls())  
load("step1_input.Rdata")
group_list <- datTraits$group
dat.pca <- PCA(datExpr0, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, #标签不重叠
                    col.ind = group_list, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)

##保存数据
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")

############################### 2.挑选最佳阈值power ###################################
rm(list = ls())  
load("E:/wanglu/白念/R2/step1_input.Rdata")
R.sq_cutoff = 0.8  #设置R^2 cut-off，默认为0.85
if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(c(1:10), seq(from = 12, to=50, by=1))
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
# power = sft$powerEstimate
power = 15
#WT 7 R1 10 R2 15
# 若无向网络在power小于15或有向网络power小于30内，没有一个power值使
# 无标度网络图谱结构R^2达到0.8且平均连接度在100以下，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
#power=NULL
if(F){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

save(sft, power, file = "step2_power_value.Rdata")

##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    #corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 30,    ##越大模块越少
    mergeCutHeight = 0.25, ##越大模块越少
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}
##################### 3.2步法构建加权共表达网络，识别基因模块 ####################
if(T){
  ## 构建加权共表达网络分为两步：
  ## 1. 计算邻近值，也是就是两个基因在不同样品中表达量的表达相关系数(pearson correlation rho)，
  ## 2. 计算topology overlap similarity (TOM)。 用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
  
  ###(1)网络构建 Co-expression similarity and adjacency 
  adjacency = adjacency(datExpr, power = power) 
  
  ###(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  
  ###(3) 聚类拓扑矩阵 Clustering using TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  ## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  
  ###(4) 聚类分支的修整 dynamicTreeCut 
  ################# set the minimum module size ##############################
  minModuleSize = 30
  ####
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  
  ###(5) 聚类结果相似模块的融合 Merging of modules whose expression profiles are very similar
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  #一般选择 height cut 为0.25,对应于有75%相关性，进行融合
  ###################### set  Merging height cut  ################################
  MEDissThres = 0.5
  ####
  # Plot the result
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  # 统计mergedmodule
  table(mergedColors)
  
  ### (6) plot the gene dendrogram 
  pdf(file = "step3_stepbystep_DynamicTreeCut_genes-modules.pdf", width = 16,height = 12)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  ### 保存数据
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = "step3_stepByStep_genes_modules.Rdata")
  
}

# ## 模块可视化，层级聚类树展示各个模块
# if(T){
#   # Convert labels to colors for plotting
#   moduleColors <- labels2colors(net$colors)
#   table(moduleColors)
#   # Plot the dendrogram and the module colors underneath
#   pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
#   plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
#                       "Module colors",
#                       dendroLabels = FALSE, hang = 0.03,
#                       addGuide = TRUE, guideHang = 0.05)
#   dev.off()
# }

####################### 4.关联基因模块与表型 #####################################
rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_stepByStep_genes_modules.Rdata")

## 模块与表型的相关性热图
if(T){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group) #get the group
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap.pdf",
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}


### 模块与表型的相关性boxplot图 
if(T){
  mes_group <- merge(MEs,datTraits,by="row.names") 
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }
  # 批量画boxplot
  colorNames <- names(MEs)
  pdf("step4_Module-trait-relationship_boxplot.pdf", width = 7.5,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) #排布为每行2个
  dev.off()
}

############################### 7.感兴趣基因模块绘制热图 ######################################
rm(list = ls())  
load(file = 'E:/wanglu/白念/WT/step1_input.Rdata')
load(file = "step3_stepByStep_genes_modules.Rdata")
table(moduleColors)

module = "blue"
### 感兴趣模块画热图 
if(T){
  dat=datExpr[,moduleColors==module]
  library(pheatmap)
  n=t(scale(dat)) #对基因做scale，并转置表达矩阵为行为基因、列为样本形式
  # n[n>2]=2 
  # n[n< -2]= -2
  # n[1:4,1:4]
  
  group_list=datTraits$group
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n)
  
  pheatmap::pheatmap(n,
                     fontsize = 8,
                     show_colnames =T,
                     show_rownames = F,
                     cluster_cols = T,
                     annotation_col =ac,
                     width = 8,
                     height = 6,
                     angle_col=45,
                     main = paste0("module_",module,"-gene heatmap"),
                     filename = paste0("step7_module_",module,"_Gene-heatmap.pdf"))
  
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")


### 基因与模块、表型的相关性散点图
#所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因算出相关系数， 
#如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。

# 选择离散性状的表型
levels(datTraits$group)
choose_group <- "3VWT"  

if(T){
  modNames <- substring(names(MEs), 3)
  
  ### 计算模块与基因的相关性矩阵 
  ## Module Membership: 模块内基因表达与模块特征值的相关性
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  ###  计算性状与基因的相关性矩阵 
  ## Gene significance，GS：比较样本某个基因与对应表型的相关性
  ## 连续型性状
  # trait <- datTraits$groupNo  
  ## 非连续型性状，需转为0-1矩阵, 已存于design中
  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  ### 可视化基因与模块、表型的相关性.
  #selectModule<-c("blue","green","purple","grey")  ##可以选择自己想要的模块
  selectModule <- modNames  ## 全部模块批量作图
  pdf("step4_gene-Module-trait-significance.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2)) #批量作图开始
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}


#########################  5. WGCNA可视化：TOMplot  Eigengene-adjacency-heatmap ##################################
rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

if(T){
  TOM=TOMsimilarityFromExpr(datExpr,power=power)
  dissTOM=1-TOM
  ## draw all genes 
  if(T){
    #geneTree = net$dendrograms[[1]]
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png("step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
  ### draw selected genes to save time...just for test...
  if(F){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf("step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}


### 模块相关性展示 Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    ## 非连续型性状，需将是否属于这个表型进行0,1数值化，已存于design中
    design
    `3TR2` = as.data.frame(design[,4])
    names(`3TR2`) = "3TR2"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, `3TR2`))
  }
  pdf("step5_module_cor_Eigengene-dendrogram.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="", 
                        marDendro = c(0,4,1,4),  # 留白：下右上左
                        marHeatmap = c(5,5,1,2), # 留白：下右上左
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}


#################### 6. 选择感兴趣基因模块进行GO分析 ####################
#rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_stepByStep_genes_modules.Rdata")
load(file = "step4_design.Rdata")

### 条件设置
OrgDb = "org.Mm.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
choose_module <- c("brown")




if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F) 
  abc = read.delim('concat.txt',sep='\t')
  #colnames(abc)[1] = 'gene'
  merged_df <- merge(gene_module, abc[,c(1,2)], by.x = "gene", by.y = "gene_id", all.x = TRUE)
  
  gene_module = merged_df[,c(3,2)]
  colnames(gene_module)[1] = 'gene'
  
  tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
              toType = "ENTREZID",
              OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = OrgDb,
    ont = "ALL",  #One of "BP", "MF", and "CC"  or "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff = 0.25,
    qvalueCutoff = 0.25
  )
  
  ###精简GO富集的结果,去冗余
  lineage1_ego = formula_res
  # lineage1_ego <- simplify(
  #   formula_res,
  #   cutoff=0.1,
  #   by="p.adjust",
  #   select_fun=min,
  # measure = "Wang",
  # semData = NULL
  # )
  save(gene_module, formula_res, lineage1_ego, file="step6_module_GO_term.Rdata")
  write.csv(lineage1_ego@compareClusterResult,
            file="step6_module_GO_term.csv")
  ### 绘制dotplot图
  dotp <- dotplot(lineage1_ego,
                  showCategory=10,
                  includeAll = TRUE, #将有overlap的结果也展示出来
                  label_format=90)
  ggsave(dotp,filename= "step6_module_GO_term.pdf", #device = cairo_pdf,
         width = 12, 
         height = 15)
}




################### 8.感兴趣模块基因导出 VisANT or cytoscape ######################
rm(list = ls())  
load(file = 'step1_input.Rdata')
#load(file = 'R2_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_stepByStep_genes_modules.Rdata")
module = "grey"
#power=10
if(T){
  ### 提取感兴趣模块基因名
  gene <- colnames(datExpr) 
  inModule <- moduleColors==module
  modgene <- gene[inModule]
  
  ### 模块对应的基因关系矩阵
  TOM <- TOMsimilarityFromExpr(datExpr,power=power)
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modgene,modgene)
  
  ### 筛选连接度最大的top100基因
  nTop = 100
  IMConn = softConnectivity(datExpr[, modgene]) #计算连接度
  top = (rank(-IMConn) <= nTop) #选取连接度最大的top100
  filter_modTOM <- modTOM[top, top]
  
  # for visANT
  vis <- exportNetworkToVisANT(filter_modTOM,
                               file = paste("step8_visANTinput-",module,".txt",sep = ""),
                               weighted = T,threshold = 0)
  # for cytoscape
  cyt <- exportNetworkToCytoscape(filter_modTOM,
                                  edgeFile = paste("step8_CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste("step8_CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.15,  #weighted权重筛选阈值，可调整
                                  nodeNames = modgene[top], 
                                  nodeAttr = moduleColors[inModule][top])
}

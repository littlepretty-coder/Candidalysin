library(DESeq2)
data = read.table('./data/fungi_human_count.txt',sep='\t',row.names=1,header=TRUE,check.names = F)
group = read.table('./data/human_fungi_group.txt',sep='\t',row.names=1,header=TRUE,check.names = F)


f2059_group = group[group$fungi=='2059',]
f2059_data = data[,rownames(f2059_group)]

f2059_data <- f2059_data[rowMeans(f2059_data)>1,] 
fix(f2059_group)
#################################################################


dds <- DESeqDataSetFromMatrix(countData = f2059_data, colData = f2059_group, design = ~ time)
colData_yeast<-data.frame(time=factor(f2059_group$time),levels=c("0h","3h","6h",'6hc'))
pre_dds_glab <- DESeqDataSetFromMatrix(f2059_data, colData=colData_yeast, ~time)
dds_glab <- DESeq(pre_dds_glab)

##### all compared with 0 time point########
res_cglab_0_3 <- results(dds_glab, cooksCutoff=FALSE, contrast = c("time","3h","0h"))
write.table(res_cglab_0_3, file="E:/projects/王秀丽/results/fungi/human/res_2059_0_vs_3.txt", sep = "\t", quote = FALSE)

res_cglab_0_6 <- results(dds_glab, cooksCutoff=FALSE, contrast = c("time","6h","0h"))
write.table(res_cglab_0_6, file="E:/projects/王秀丽/results/fungi/human/res_2059_0_vs_6.txt", sep = "\t", quote = FALSE)

res_cglab_0_6c <- results(dds_glab, cooksCutoff=FALSE, contrast = c("time","6hc","0h"))
write.table(res_cglab_0_6c, file="E:/projects/王秀丽/results/fungi/human/res_2059_0_vs_6c.txt", sep = "\t", quote = FALSE)

all_data_cglab<-cbind.data.frame(res_cglab_0_3$log2FoldChange,res_cglab_0_3$padj,
                                 res_cglab_0_6$log2FoldChange,res_cglab_0_6$padj)
write.table(all_data_cglab,"E:/projects/王秀丽/results/fungi/human/res_2059.txt", sep = "\t", quote = FALSE)

rownames(all_data_cglab)<-rownames(res_cglab_0_6)
colnames(all_data_cglab)<-c("3", "padj3", "6", "padj6") #, "24c", "padj24c")


all_data_cglab<-all_data_cglab[complete.cases(all_data_cglab),]

all_data_cglab$`3`[all_data_cglab$padj3 > 0.01] <- 0   #### if padj > 0.01 then LFC is 0 (if not significant then no LFC)
all_data_cglab$`6`[all_data_cglab$padj6 > 0.01] <- 0
#all_data_cglab$`6hc`[all_data_cglab$padj6hc > 0.01] <- 0


LFC_cglab<-cbind.data.frame("0" = 0, all_data_cglab$`3`, all_data_cglab$`6`)
colnames(LFC_cglab)<-c("0","3","6")
rownames(LFC_cglab)<-rownames(all_data_cglab)


LFC<-LFC_cglab

################## FILTERING

LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]
dim(LFC)
############ Plotting of interpolated LFC #################3

time<-c(0,3,6)
time_interpol<-approx(time)$y


LFC_no_NA<-LFC[complete.cases(LFC),]
LFC_no_NA_t<-t(LFC_no_NA)

### Interpolation
interpol_LFC<-apply(LFC_no_NA_t, 2, approx)

### rbind only y values from interpol_LFC

combined_interplolated<-do.call(rbind, lapply(interpol_LFC, '[[','y'))

### parsing and plotting
colnames(combined_interplolated)<-time_interpol

LFC_plotting_inter<-melt(combined_interplolated,value.name = "value", varnames=c('Var1', 'Var2'))


res_cglab_0_3_up = nrow(res_cglab_0_3[(res_cglab_0_3$log2FoldChange>1.5 & res_cglab_0_3$padj < 0.01),])
res_cglab_0_3_down = nrow(res_cglab_0_3[(res_cglab_0_3$log2FoldChange<-1.5 & res_cglab_0_3$padj < 0.01),])


res_cglab_0_6_up = nrow(res_cglab_0_6[(res_cglab_0_6$log2FoldChange>1.5 & res_cglab_0_6$padj < 0.01),])
res_cglab_0_6_down = nrow(res_cglab_0_6[(res_cglab_0_6$log2FoldChange<-1.5 & res_cglab_0_6$padj < 0.01),])


pdf('E:/projects/王秀丽/results/fungi/human/2059_dyn.pdf',  width=5, heigh=3)
ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
  geom_line() + xlab("Time (h)")+ylab("Log 2 fold change") + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(-30,30),space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+
  labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size=18),
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        axis.ticks.x=element_blank())+
  annotate("text", x = 3, y = 6, label = res_cglab_0_3_up)+
  annotate("text", x = 3, y = -6, label = res_cglab_0_3_down)+
  annotate("text", x = 5.5, y = 6, label = res_cglab_0_6_up)+
  annotate("text", x = 5.5, y = -6, label = res_cglab_0_6_down)+
  scale_y_continuous(breaks=seq(-30,30,5)) + expand_limits(y=c(-30,30))+scale_x_continuous(breaks=c(0, 3, 6))### add this to change ticks
dev.off()

genes_glab_0_6<-rownames(res_cglab_0_6[!is.na(res_cglab_0_6$padj) & res_cglab_0_6$padj<0.01,])
genes_glab_0_6c<-rownames(res_cglab_0_6c[!is.na(res_cglab_0_6c$padj) & res_cglab_0_6c$padj<0.01,])
overlap<-as.data.frame(Reduce(intersect,list(genes_glab_0_6,genes_glab_0_6c)))


venn.diagram(list("2059_0_6"=genes_glab_0_6,
                  "2059_0_6c"=genes_glab_0_6c),
             fill=c("green", "cornflowerblue"), height = 2500, width = 2500, resolution = 600, 
             main = "C. glabrata DE genes between 6h and 6h control", cat.cex=0.4, cex=1.8, main.cex = 0.7, 
             filename="E:/projects/王秀丽/results/fungi/human/2059_6_6c.png",  
             imagetype = "png")

only_infection_cglab<-setdiff(genes_glab_0_6, genes_glab_0_6c)
res_cglab_inf_spec<-res_cglab_0_6[only_infection_cglab,]
write.table(res_cglab_inf_spec,"E:/projects/王秀丽/results/fungi/human/infection_genes_2059_6_6c.txt", quote = FALSE, sep = "\t")
#wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww







data = read.table('./data/fungi_mouse_count.txt',sep='\t',row.names=1,header=TRUE,check.names = F)
group = read.table('./data/mouse_fungi_group.txt',sep='\t',row.names=1,header=TRUE,check.names = F)


f2059_group = group[group$fungi=='2059',]
f2059_data = data[,rownames(f2059_group)]

f2059_data <- f2059_data[rowMeans(f2059_data)>1,] 
fix(f2059_group)
#################################################################


dds <- DESeqDataSetFromMatrix(countData = f2059_data, colData = f2059_group, design = ~ time)
colData_yeast<-data.frame(time=factor(f2059_group$time),levels=c("0h","3h","12h",'12hc'))
pre_dds_glab <- DESeqDataSetFromMatrix(f2059_data, colData=colData_yeast, ~time)
dds_glab <- DESeq(pre_dds_glab)

##### all compared with 0 time point########
res_cglab_0_3 <- results(dds_glab, cooksCutoff=FALSE, contrast = c("time","3h","0h"))
write.table(res_cglab_0_3, file="E:/projects/王秀丽/results/fungi/mouse/res_2059_0_vs_3.txt", sep = "\t", quote = FALSE)

res_cglab_0_12 <- results(dds_glab, cooksCutoff=FALSE, contrast = c("time","12h","0h"))
write.table(res_cglab_0_12, file="E:/projects/王秀丽/results/fungi/mouse/res_2059_0_vs_12.txt", sep = "\t", quote = FALSE)

res_cglab_0_12c <- results(dds_glab, cooksCutoff=FALSE, contrast = c("time","12hc","0h"))
write.table(res_cglab_0_12c, file="E:/projects/王秀丽/results/fungi/mouse/res_2059_0_vs_12c.txt", sep = "\t", quote = FALSE)

all_data_cglab<-cbind.data.frame(res_cglab_0_3$log2FoldChange,res_cglab_0_3$padj,
                                 res_cglab_0_12$log2FoldChange,res_cglab_0_12$padj)
write.table(all_data_cglab,"E:/projects/王秀丽/results/fungi/mouse/res_2059.txt", sep = "\t", quote = FALSE)

rownames(all_data_cglab)<-rownames(res_cglab_0_12)
colnames(all_data_cglab)<-c("3", "padj3", "12", "padj12") #, "24c", "padj24c")


all_data_cglab<-all_data_cglab[complete.cases(all_data_cglab),]

all_data_cglab$`3`[all_data_cglab$padj3 > 0.01] <- 0   #### if padj > 0.01 then LFC is 0 (if not significant then no LFC)
all_data_cglab$`12`[all_data_cglab$padj12 > 0.01] <- 0
#all_data_cglab$`6hc`[all_data_cglab$padj6hc > 0.01] <- 0


LFC_cglab<-cbind.data.frame("0" = 0, all_data_cglab$`3`, all_data_cglab$`12`)
colnames(LFC_cglab)<-c("0","3","12")
rownames(LFC_cglab)<-rownames(all_data_cglab)


LFC<-LFC_cglab

################## FILTERING

LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]
dim(LFC)
############ Plotting of interpolated LFC #################3

time<-c(0,3,12)
time_interpol<-approx(time)$y


LFC_no_NA<-LFC[complete.cases(LFC),]
LFC_no_NA_t<-t(LFC_no_NA)

### Interpolation
interpol_LFC<-apply(LFC_no_NA_t, 2, approx)

### rbind only y values from interpol_LFC

combined_interplolated<-do.call(rbind, lapply(interpol_LFC, '[[','y'))

### parsing and plotting
colnames(combined_interplolated)<-time_interpol

LFC_plotting_inter<-melt(combined_interplolated,value.name = "value", varnames=c('Var1', 'Var2'))


res_cglab_0_3_up = nrow(res_cglab_0_3[(res_cglab_0_3$log2FoldChange>1.5 & res_cglab_0_3$padj < 0.01),])
res_cglab_0_3_down = nrow(res_cglab_0_3[(res_cglab_0_3$log2FoldChange<-1.5 & res_cglab_0_3$padj < 0.01),])


res_cglab_0_12_up = nrow(res_cglab_0_12[(res_cglab_0_12$log2FoldChange>1.5 & res_cglab_0_12$padj < 0.01),])
res_cglab_0_12_down = nrow(res_cglab_0_12[(res_cglab_0_12$log2FoldChange<-1.5 & res_cglab_0_12$padj < 0.01),])


pdf('E:/projects/王秀丽/results/fungi/mouse/2059_dyn.pdf',  width=5, heigh=3)
ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
  geom_line() + xlab("Time (h)")+ylab("Log 2 fold change") + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(-30,30),space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+
  labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size=18),
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        axis.ticks.x=element_blank())+
  annotate("text", x = 3, y = 6, label = res_cglab_0_3_up)+
  annotate("text", x = 3, y = -6, label = res_cglab_0_3_down)+
  annotate("text", x = 11.5, y = 6, label = res_cglab_0_12_up)+
  annotate("text", x = 11.5, y = -6, label = res_cglab_0_12_down)+
  scale_y_continuous(breaks=seq(-30,30,5)) + expand_limits(y=c(-30,30))+scale_x_continuous(breaks=c(0, 3, 12))### add this to change ticks
dev.off()


genes_glab_0_12<-rownames(res_cglab_0_12[!is.na(res_cglab_0_12$padj) & res_cglab_0_12$padj<0.01,])
genes_glab_0_12c<-rownames(res_cglab_0_12c[!is.na(res_cglab_0_12c$padj) & res_cglab_0_12c$padj<0.01,])
overlap<-as.data.frame(Reduce(intersect,list(genes_glab_0_12,genes_glab_0_12c)))


venn.diagram(list("2059_0_12"=genes_glab_0_12,
                  "2059_0_12c"=genes_glab_0_12c),
             fill=c("green", "cornflowerblue"), height = 2500, width = 2500, resolution = 600, 
             main = "C. glabrata DE genes between 6h and 6h control", cat.cex=0.4, cex=1.8, main.cex = 0.7, 
             filename="E:/projects/王秀丽/results/fungi/mouse/2059_12_12c.png",  
             imagetype = "png")

only_infection_cglab<-setdiff(genes_glab_0_12, genes_glab_0_12c)
res_cglab_inf_spec<-res_cglab_0_12[only_infection_cglab,]
write.table(res_cglab_inf_spec,"E:/projects/王秀丽/results/fungi/mouse/infection_genes_2059_12_12c.txt", quote = FALSE, sep = "\t")



















##########################  HUMAN PART   ################################

human_controls<-read.table("./data/human_count.txt",check.names = F,header=T,row.names=1)
group = as.data.frame(colnames(human_controls))
colnames(group) = 'samples'
df_split <- do.call(rbind, strsplit(as.character(group$samples), split = "-", fixed = TRUE))
group = cbind(group,df_split[,1:3])
fix(group)
colnames(group)[3:4] = c('fungi','time')
rownames(group) = group$samples
f2059_group = group[group$fungi%in%c('2059','CTL'),]
f2059_data = human_controls[,f2059_group$samples]

f2059_data <- f2059_data[rowMeans(f2059_data)>1,] 
fix(f2059_group)
#write.table(group,'human_group.txt',sep='\t')
#all_human_data<-cbind.data.frame("1"=human_controls$`1`, "2"=human_controls$`2`, human, "17"=human_controls$`17`, "18"=human_controls$`18`) 
all_human_data = f2059_data
### DESEQ2 ###

colData_human<-data.frame(time=factor(f2059_group$time,levels=c("0h","3h","6h",'6hc')))
pre_dds_human <- DESeqDataSetFromMatrix(all_human_data, colData=colData_human, ~time)
dds_human <- DESeq(pre_dds_human)

##### all compared with 0 time point########
res_cglab_human_0_3 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","3h","0h"))
write.table(res_cglab_human_0_3, file="E:/projects/王秀丽/results/human/res_human_2059_0_vs_3.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_3 = read.delim('E:/projects/王秀丽/3-14/results/human/2059/res_human_2059_0_vs_3.txt',sep='\t')
res_cglab_human_0_6 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","6h","0h"))
write.table(res_cglab_human_0_6, file="E:/projects/王秀丽/results/human/res_human_2059_0_vs_6.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_6 = read.delim('E:/projects/王秀丽/3-14/results/human/2059/res_human_2059_0_vs_6.txt',sep='\t')
res_cglab_human_0_6c <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","6hc","0h"))
write.table(res_cglab_human_0_6c, file="E:/projects/王秀丽/3-14/results/human/2059/res_human_2059_0_vs_6hc.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_6c = read.delim('E:/projects/王秀丽/3-14/results/human/2059/res_human_2059_0_vs_6hc.txt',sep='\t')
#########Plotting of interpolated LFC #####################
all_data<-cbind.data.frame(res_cglab_human_0_3$log2FoldChange,res_cglab_human_0_3$padj,
                           res_cglab_human_0_6$log2FoldChange,res_cglab_human_0_6$padj)


rownames(all_data)<-rownames(res_cglab_human_0_6)
colnames(all_data)<-c("3", "padj3", "6", "padj6") #, "24c", "padj24c")

all_data<-all_data[complete.cases(all_data),]


all_data$`3`[all_data$padj3 > 0.01] <- 0
all_data$`6`[all_data$padj6 > 0.01] <- 0



LFC<-cbind.data.frame("0" = 0, all_data$`3`, all_data$`6`)
colnames(LFC)<-c("0",'3',"6")
rownames(LFC)<-rownames(all_data)

################## FILTERING
LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]

############ Plotting of interpolated LFC #################3
time<-c(0,3,6)
time_interpol<-approx(time)$y


LFC_no_NA<-LFC[complete.cases(LFC),]
LFC_no_NA_t<-t(LFC_no_NA)

### Interpolation
interpol_LFC<-apply(LFC_no_NA_t, 2, approx)

### rbind only y values from interpol_LFC

combined_interplolated<-do.call(rbind, lapply(interpol_LFC, '[[','y'))

### parsing and plotting
colnames(combined_interplolated)<-time_interpol

LFC_plotting_inter<-melt(combined_interplolated,value.name = "value", varnames=c('Var1', 'Var2'))


res_cglab_human_0_3 = na.omit(res_cglab_human_0_3)
res_cglab_human_0_3_up = nrow(res_cglab_human_0_3[(res_cglab_human_0_3$log2FoldChange>1.5 & res_cglab_human_0_3$padj < 0.01),])
res_cglab_human_0_3_down = nrow(res_cglab_human_0_3[(res_cglab_human_0_3$log2FoldChange<-1.5 & res_cglab_human_0_3$padj < 0.01),])

res_cglab_human_0_6 = na.omit(res_cglab_human_0_6)
res_cglab_human_0_6_up = nrow(res_cglab_human_0_6[(res_cglab_human_0_6$log2FoldChange>1.5 & res_cglab_human_0_6$padj < 0.01),])
res_cglab_human_0_6_down = nrow(res_cglab_human_0_6[(res_cglab_human_0_6$log2FoldChange<-1.5 & res_cglab_human_0_6$padj < 0.01),])


png('E:/projects/王秀丽/results/human/human_2059_dyn.png', units="in", width=5, heigh=2.5, res=400)
ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
  geom_line() + xlab("Time (h)")+ylab("Log2 fold change") + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(-30,30),space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        axis.title.x=element_blank(),
        axis.ticks.y = element_blank())+labs(colour = "L2FC\ngradient")+
  annotate("text", x = 3, y = 6, label = res_cglab_human_0_3_up)+
  annotate("text", x = 3, y = -6, label = res_cglab_human_0_3_down)+
  annotate("text", x = 5.5, y = 6, label = res_cglab_human_0_6_up)+
  annotate("text", x = 5.5, y = -6, label = res_cglab_human_0_6_down)+
  scale_y_continuous(breaks=seq(-30,30,5))+expand_limits(y=c(-30,30))+scale_x_continuous(breaks=c(0, 3, 6)) ### this is for changing ticks
ggsave('E:/projects/王秀丽/results/human/human_2059_dyn.pdf',width=5, heigh=2.5)


###3####hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh



##########################  mouse PART   ################################

human_controls<-read.table("./data/mouse_count.txt",check.names = F,header=T,row.names=1)
group = as.data.frame(colnames(human_controls))
colnames(group) = 'samples'
df_split <- do.call(rbind, strsplit(as.character(group$samples), split = "-", fixed = TRUE))
group = cbind(group,df_split[,1:3])
fix(group)
rownames(group) = group$samples
group$'3' = paste0(group$'3','h')
#colnames(group)[3:4] = c('fungi','time')
colnames(group)[3:4] = c('fungi','time')
f2059_group = group[group$fungi%in%c('2059','CTL'),]
f2059_data = human_controls[,f2059_group$samples]

f2059_data <- f2059_data[rowMeans(f2059_data)>1,] 
fix(f2059_group)
write.table(group,'mouse_group.txt',sep='\t')
#all_human_data<-cbind.data.frame("1"=human_controls$`1`, "2"=human_controls$`2`, human, "17"=human_controls$`17`, "18"=human_controls$`18`) 
all_human_data = f2059_data
### DESEQ2 ###

colData_human<-data.frame(time=factor(f2059_group$time,levels=c("0h","3h","12h",'12hc')))
pre_dds_human <- DESeqDataSetFromMatrix(all_human_data, colData=colData_human, ~time)
dds_human <- DESeq(pre_dds_human)

##### all compared with 0 time point########
res_cglab_human_0_3 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","3h","0h"))
write.table(res_cglab_human_0_3, file="E:/projects/王秀丽/results/mouse/res_mouse_2059_0_vs_3.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_3 = read.delim('E:/projects/王秀丽/3-14/results/mouse/2059/res_mouse_2059_0_vs_3.txt')
res_cglab_human_0_12 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","12h","0h"))
write.table(res_cglab_human_0_12, file="E:/projects/王秀丽/results/mouse/res_mouse_2059_0_vs_12.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_12 = read.delim('E:/projects/王秀丽/3-14/results/mouse/2059/res_mouse_2059_0_vs_12.txt')
res_cglab_human_0_12c <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","12hc","0h"))
write.table(res_cglab_human_0_12c, file="E:/projects/王秀丽/results/mouse/res_mouse_2059_0_vs_12hc.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_12c = read.delim('E:/projects/王秀丽/3-14/results/mouse/2059/res_mouse_2059_0_vs_12hc.txt')
#########Plotting of interpolated LFC #####################
all_data<-cbind.data.frame(res_cglab_human_0_3$log2FoldChange,res_cglab_human_0_3$padj,
                           res_cglab_human_0_12$log2FoldChange,res_cglab_human_0_12$padj)


rownames(all_data)<-rownames(res_cglab_human_0_12)
colnames(all_data)<-c("3", "padj3", "12", "padj12") #, "24c", "padj24c")

all_data<-all_data[complete.cases(all_data),]


all_data$`3`[all_data$padj3 > 0.01] <- 0
all_data$`12`[all_data$padj12 > 0.01] <- 0



LFC<-cbind.data.frame("0" = 0, all_data$`3`, all_data$`12`)
colnames(LFC)<-c("0",'3',"12")
rownames(LFC)<-rownames(all_data)

################## FILTERING
LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]

############ Plotting of interpolated LFC #################3
time<-c(0,3,12)
time_interpol<-approx(time)$y


LFC_no_NA<-LFC[complete.cases(LFC),]
LFC_no_NA_t<-t(LFC_no_NA)

### Interpolation
interpol_LFC<-apply(LFC_no_NA_t, 2, approx)

### rbind only y values from interpol_LFC

combined_interplolated<-do.call(rbind, lapply(interpol_LFC, '[[','y'))

### parsing and plotting
colnames(combined_interplolated)<-time_interpol

LFC_plotting_inter<-melt(combined_interplolated,value.name = "value", varnames=c('Var1', 'Var2'))


res_cglab_human_0_3 = na.omit(res_cglab_human_0_3)
res_cglab_human_0_3_up = nrow(res_cglab_human_0_3[(res_cglab_human_0_3$log2FoldChange>1.5 & res_cglab_human_0_3$padj < 0.01),])
res_cglab_human_0_3_down = nrow(res_cglab_human_0_3[(res_cglab_human_0_3$log2FoldChange<-1.5 & res_cglab_human_0_3$padj < 0.01),])

res_cglab_human_0_12 = na.omit(res_cglab_human_0_12)
res_cglab_human_0_12_up = nrow(res_cglab_human_0_12[(res_cglab_human_0_12$log2FoldChange>1.5 & res_cglab_human_0_12$padj < 0.01),])
res_cglab_human_0_12_down = nrow(res_cglab_human_0_12[(res_cglab_human_0_12$log2FoldChange<-1.5 & res_cglab_human_0_12$padj < 0.01),])


png('E:/projects/王秀丽/results/mouse/mouse_2059_dyn.png', units="in", width=5, heigh=2.5, res=400)
ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
  geom_line() + xlab("Time (h)")+ylab("Log2 fold change") + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(-30,30),space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        axis.title.x=element_blank(),
        axis.ticks.y = element_blank())+labs(colour = "L2FC\ngradient")+annotate("text", x = 3, y = 6, label = res_cglab_human_0_3_up)+
  annotate("text", x = 3, y = -6, label = res_cglab_human_0_3_down)+
  annotate("text", x = 11.5, y = 6, label = res_cglab_human_0_12_up)+
  annotate("text", x = 11.5, y = -6, label = res_cglab_human_0_12_down)+
  scale_y_continuous(breaks=seq(-30,30,5))+expand_limits(y=c(-30,30))+scale_x_continuous(breaks=c(0, 3, 12)) ### this is for changing ticks
ggsave('E:/projects/王秀丽/results/mouse/mouse_2059_dyn.pdf',width=5, heigh=2.5)

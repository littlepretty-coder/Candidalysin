
##########################  HUMAN PART   ################################

human_controls<-read.table("./data/human_count.txt",check.names = F,header=T,row.names=1)
group = as.data.frame(colnames(human_controls))
colnames(group) = 'samples'
df_split <- do.call(rbind, strsplit(as.character(group$samples), split = "-", fixed = TRUE))
group = cbind(group,df_split[,1:3])
fix(group)
colnames(group)[3:4] = c('fungi','time')
rownames(group) = group$samples
fCa_group = group[group$fungi%in%c('Ca','CTL'),]
fCa_group = fCa_group[fCa_group$time%in%c('0h','3h','24h'),]
fCa_data = human_controls[,fCa_group$samples]

fCa_data <- fCa_data[rowMeans(fCa_data)>1,] 
fix(fCa_group)
#write.table(group,'human_group.txt',sep='\t')
#all_human_data<-cbind.data.frame("1"=human_controls$`1`, "2"=human_controls$`2`, human, "17"=human_controls$`17`, "18"=human_controls$`18`) 
all_human_data = fCa_data
### DESEQ2 ###

colData_human<-data.frame(time=factor(fCa_group$time,levels=c("0h","3h","24h")))
pre_dds_human <- DESeqDataSetFromMatrix(all_human_data, colData=colData_human, ~time)
dds_human <- DESeq(pre_dds_human)

##### all compared with 0 time point########
res_cglab_human_0_3 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","3h","0h"))
write.table(res_cglab_human_0_3, file="E:/projects/王秀丽/results/human/res_human_Ca_0_vs_3.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_3 = read.delim("E:/projects/王秀丽/3-14/results/human/ca/res_human_Ca_0_vs_3.txt")
res_cglab_human_0_24 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","24h","0h"))
write.table(res_cglab_human_0_24, file="E:/projects/王秀丽/results/human/res_human_Ca_0_vs_24.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_24 = read.delim("E:/projects/王秀丽/3-14/results/human/ca/res_human_Ca_0_vs_24.txt")
# res_cglab_human_0_6c <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","6hc","0h"))
# write.table(res_cglab_human_0_6c, file="E:/projects/王秀丽/results/human/res_human_Ca_0_vs_6hc.txt", sep = "\t", quote = FALSE)

#########Plotting of interpolated LFC #####################
all_data<-cbind.data.frame(res_cglab_human_0_3$log2FoldChange,res_cglab_human_0_3$padj,
                           res_cglab_human_0_24$log2FoldChange,res_cglab_human_0_24$padj)


rownames(all_data)<-rownames(res_cglab_human_0_24)
colnames(all_data)<-c("3", "padj3", "24", "padj24") #, "24c", "padj24c")

all_data<-all_data[complete.cases(all_data),]


all_data$`3`[all_data$padj3 > 0.01] <- 0
all_data$`24`[all_data$padj24 > 0.01] <- 0



LFC<-cbind.data.frame("0" = 0, all_data$`3`, all_data$`24`)
colnames(LFC)<-c("0",'3',"24")
rownames(LFC)<-rownames(all_data)

################## FILTERING
LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]

############ Plotting of interpolated LFC #################3
time<-c(0,3,24)
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

res_cglab_human_0_24 = na.omit(res_cglab_human_0_24)
res_cglab_human_0_24_up = nrow(res_cglab_human_0_24[(res_cglab_human_0_24$log2FoldChange>1.5 & res_cglab_human_0_24$padj < 0.01),])
res_cglab_human_0_24_down = nrow(res_cglab_human_0_24[(res_cglab_human_0_24$log2FoldChange<-1.5 & res_cglab_human_0_24$padj < 0.01),])


png('E:/projects/王秀丽/results/human/human_Ca_dyn.png', units="in", width=5, heigh=2.5, res=400)
ggplot(data=LFC_plotting_inter, aes(x=Var2, y=value, group=Var1, colour=value)) +
  geom_line() + xlab("Time (h)")+ylab("Log2 fold change") + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(-30,30),space = "Lab", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        axis.title.x=element_blank(),
        axis.ticks.y = element_blank())+labs(colour = "L2FC\ngradient")+
  annotate("text", x = 3, y = 6, label = res_cglab_human_0_3_up)+
  annotate("text", x = 3, y = -6, label = res_cglab_human_0_3_down)+
  annotate("text", x = 23.5, y = 6, label = res_cglab_human_0_24_up)+
  annotate("text", x = 23.5, y = -6, label = res_cglab_human_0_24_down)+
  scale_y_continuous(breaks=seq(-30,30,5))+expand_limits(y=c(-30,30))+scale_x_continuous(breaks=c(0, 3, 24)) ### this is for changing ticks
ggsave('E:/projects/王秀丽/results/human_ca_dyn.pdf',width=5, heigh=2.5)


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
fCa_group = group[group$fungi%in%c('Ca','CTL'),]
fCa_group = fCa_group[fCa_group$time%in%c('0h','3h','24h'),]
fCa_data = human_controls[,fCa_group$samples]

fCa_data <- fCa_data[rowMeans(fCa_data)>1,] 
fix(fCa_group)
write.table(group,'mouse_group.txt',sep='\t')
#all_human_data<-cbind.data.frame("1"=human_controls$`1`, "2"=human_controls$`2`, human, "17"=human_controls$`17`, "18"=human_controls$`18`) 
all_human_data = fCa_data
### DESEQ2 ###

colData_human<-data.frame(time=factor(fCa_group$time,levels=c("0h","3h","24h")))
pre_dds_human <- DESeqDataSetFromMatrix(all_human_data, colData=colData_human, ~time)
dds_human <- DESeq(pre_dds_human)

##### all compared with 0 time point########
res_cglab_human_0_3 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","3h","0h"))
write.table(res_cglab_human_0_3, file="E:/projects/王秀丽/results/mouse/res_mouse_Ca_0_vs_3.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_3 <- read.delim("E:/projects/王秀丽/3-14/results/mouse/ca/res_mouse_Ca_0_vs_3.txt")
res_cglab_human_0_24 <- results(dds_human, cooksCutoff=FALSE, contrast = c("time","24h","0h"))
write.table(res_cglab_human_0_24, file="E:/projects/王秀丽/results/mouse/res_mouse_Ca_0_vs_24.txt", sep = "\t", quote = FALSE)
res_cglab_human_0_24 <- read.delim("E:/projects/王秀丽/3-14/results/mouse/ca/res_mouse_Ca_0_vs_24.txt")
#########Plotting of interpolated LFC #####################
all_data<-cbind.data.frame(res_cglab_human_0_3$log2FoldChange,res_cglab_human_0_3$padj,
                           res_cglab_human_0_24$log2FoldChange,res_cglab_human_0_24$padj)


rownames(all_data)<-rownames(res_cglab_human_0_24)
colnames(all_data)<-c("3", "padj3", "24", "padj24") #, "24c", "padj24c")

all_data<-all_data[complete.cases(all_data),]


all_data$`3`[all_data$padj3 > 0.01] <- 0
all_data$`24`[all_data$padj24 > 0.01] <- 0



LFC<-cbind.data.frame("0" = 0, all_data$`3`, all_data$`24`)
colnames(LFC)<-c("0",'3',"24")
rownames(LFC)<-rownames(all_data)

################## FILTERING
LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]

############ Plotting of interpolated LFC #################3
time<-c(0,3,24)
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

res_cglab_human_0_24 = na.omit(res_cglab_human_0_24)
res_cglab_human_0_24_up = nrow(res_cglab_human_0_24[(res_cglab_human_0_24$log2FoldChange>1.5 & res_cglab_human_0_24$padj < 0.01),])
res_cglab_human_0_24_down = nrow(res_cglab_human_0_24[(res_cglab_human_0_24$log2FoldChange<-1.5 & res_cglab_human_0_24$padj < 0.01),])


png('E:/projects/王秀丽/results/mouse/mouse_Ca_dyn.png', units="in", width=5, heigh=2.5, res=400)
ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
  geom_line() + xlab("Time (h)")+ylab("Log2 fold change") + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red",limits = c(-30,30),space = "Lab", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
        legend.text = element_text(size=14),legend.title = element_text(size=14),
        axis.title.x=element_blank(),
        axis.ticks.y = element_blank())+labs(colour = "L2FC\ngradient")+
  annotate("text", x = 3, y = 6, label = res_cglab_human_0_3_up)+
  annotate("text", x = 3, y = -6, label = res_cglab_human_0_3_down)+
  annotate("text", x = 23.5, y = 6, label = res_cglab_human_0_24_up)+
  annotate("text", x = 23.5, y = -6, label = res_cglab_human_0_24_down)+
  scale_y_continuous(breaks=seq(-30,30,5))+expand_limits(y=c(-30,30))+scale_x_continuous(breaks=c(0, 3, 24)) ### this is for changing ticks
ggsave('E:/projects/王秀丽/results/mouse/mouse_ca_dyn.pdf',width=5, heigh=2.5)

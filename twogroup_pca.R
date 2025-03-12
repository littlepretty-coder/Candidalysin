library(optparse)
if (!suppressWarnings(suppressMessages(require("maptools", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("maptools",character.only=T)
}

options(stringsAsFactors = F)
library(vegan)
library(maptools)
library(ggplot2)
library(ggrepel)
library(dplyr)


# 3. 读取输入文件

# 读取OTU表
dat <- log2(fpkm_filt+1)
mm_group = FPKM_group[37:72,]
mm_dat = dat[,mm_group$samples]
dat = mm_dat
groups = mm_group
genus<-t(dat)
#genus=decostand(genus,method = "hellinger")
otu_pca<- prcomp(genus,scal=F)
pc12 <- as.data.frame(otu_pca$x[,1:2])
pc12$samples<-rownames(pc12)

#pc12 <- as.data.frame(otu_pca$x[,1:2])
pc12$samples<-rownames(pc12)



colnames(groups)[3:4] = c('group1','group2')
groups$group1<-factor(groups$group1,levels = groups$group1[!duplicated(groups$group1)])
groups$group2<-factor(groups$group2,levels = groups$group2[!duplicated(groups$group2)])


pc <-summary(otu_pca)$importance[2,]*100

pc12<-merge(pc12,groups,by='samples')
pc12$group<-paste(pc12$group1,pc12$group2,sep = '')

colnames(pc12)[2:3]<-c('PC1','PC2')

mycol<-c('#E41A1C','#377EB8','#F0E68C')
myshape<-c(21:25)

ggplot(data = pc12,aes(PC1,PC2)) +
  #geom_text_repel(data = st,aes(RDA1,RDA2,label=row.names(st)),size=4)+#Show a Square
  geom_point(aes(fill=group2,shape=group1),size=3)+
  scale_shape_manual(values = myshape)+
  scale_fill_manual(values=mycol)+
  scale_color_manual(values = mycol)+
  guides(color = 'none',
         fill=guide_legend(override.aes = list(size=4,fill=mycol,shape=21)),
         size = guide_legend(override.aes = list(shape = myshape)))+
  geom_hline(yintercept=0,linetype=2) + 
  geom_vline(xintercept=0,linetype=2)+
  labs(x=paste0("PC1(",round(pc[1],2),"%",")"),y=paste0("PC2(",round(pc[2],2),"%)"))+
  theme_bw()+theme(panel.grid=element_blank())+
  stat_ellipse(aes(x = PC1, y =PC2, color = group2), linetype = 1, level = 0.95) 


ggsave('step11.twogroup.pca.pdf',width = 8,height = 6)



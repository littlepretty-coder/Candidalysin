#library(Seurat)
library(tidyverse)
#library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)

# new_row1 <- data.frame(Description = 'phosphatase activity', p.adjust = 1,fungi='1477_3h',time='3h')
# new_row2 <- data.frame(Description = 'phosphatase activity', p.adjust = 1,fungi='1477_12h',time='12h')
# 
# 
# new_row3 <- data.frame(Description = 'phosphatase activity', p.adjust = 1,fungi='2059_3h',time='3h')
# new_row4 <- data.frame(Description = 'phosphatase activity', p.adjust = 1,fungi='2059_12h',time='12h')

data1 = read.delim('E:/projects/王秀丽/results/human/3h/1477_BP.txt',sep='\t')[,c(2,6,9)]
data1$fungi = '1477_3h'
data1$time = '3h'
data2 = read.delim('E:/projects/王秀丽/results/human/3h/2174_BP.txt',sep='\t')[,c(2,6,9)]
data2$fungi = '2174_3h'
data2$time = '3h'
data3 = read.delim('E:/projects/王秀丽/results/human/3h/2059_BP.txt',sep='\t')[,c(2,6,9)]
data3$fungi = '2059_3h'
data3$time = '3h'
data11 = read.delim('E:/projects/王秀丽/results/human/6h/1477_BP.txt',sep='\t')[,c(2,6,9)]
data11$fungi = '1477_6h'
data11$time = '6h'
data22 = read.delim('E:/projects/王秀丽/results/human/6h/2174_BP.txt',sep='\t')[,c(2,6,9)]
data22$fungi = '2174_6h'
data22$time = '6h'
data33 = read.delim('E:/projects/王秀丽/results/human/6h/2059_BP.txt',sep='\t')[,c(2,6,9)]
data33$fungi = '2059_6h'
data33$time = '6h'

data44= read.delim('E:/projects/王秀丽/results/human/ca/ca_BP.txt',sep='\t')[,c(2,6,9)]
data44$fungi = 'ca_24h'
data44$time = '24h'

data4= read.delim('E:/projects/王秀丽/results/human/3h/ca_BP.txt',sep='\t')[,c(2,6,9)]
data4$fungi = 'ca_3h'
data4$time = '3h'

# data1 = rbind(data1,new_row1)
# data3 = rbind(data3,new_row3)
# data11 = rbind(data11,new_row2)
# data33 = rbind(data33,new_row4)


human_host = c('circadian regulation of gene expression',
               'regulation of transcription by RNA polymerase III',
               'regulation of neurogenesis',
               'positive regulation of neurogenesis',
               'skeletal muscle cell differentiation',
               'response to transforming growth factor beta',
               'response to temperature stimulus',
               'response to peptide hormone',
               'response to oxygen levels',
               'response to oxidative stress',
               'response to hypoxia',
               'response to glucocorticoid',
               'response to fibroblast growth factor',
               'regulation ofG1/S transition of mitotic cell cycle',
               'regulation of cyclin-dependent protein serine/threonine kinase activity',
               'positive regulation of mitotic cell cycle',
               'positive regulation of cell cycle G1/S phase transition',
               'G1/S transition of mitotic cell cycle')


human_host1 = c('circadian regulation of gene expression',
               'regulation of transcription by RNA polymerase III',
               'regulation of neurogenesis',
               'positive regulation of neurogenesis',
               'skeletal muscle cell differentiation',
               'response to transforming growth factor beta',
               'response to temperature stimulus',
               'response to peptide hormone',
               'response to oxygen levels',
               'response to oxidative stress',
               'regulation of response to oxidative stress',
               'regulation of cellular response to transforming growth factor beta stimulus',
               'regulation of cellular response to hypoxia',
               'intrinsic apoptotic signaling pathway in response to oxidative stress',
               'anterior/posterior pattern specification',
               'response to hypoxia',
               'response to glucocorticoid',
               'response to fibroblast growth factor',
               'regulation ofG1/S transition of mitotic cell cycle',
               'regulation of transcription from RNA polymerase ll promoter in response to hypoxia',
               'regulation of cyclin-dependent protein serine/threonine kinase activity',
               'positive regulation of mitotic cell cycle',
               'positive regulation of cell cycle G1/S phase transition',
               'G1/S transition of mitotic cell cycle')













all = rbind(data1,data2,data3,data11,data22,data33,data4,data44)

result <- all[all$Description %in% human_host1, ]
#result <- all[grepl(paste(mouse_host, collapse = "|"), all$Description), ]
#result <- all[grepl(paste(human_host, collapse = "|"), all$Description), ]
all = all[order(all$p.adjust),]
#result = rbind(result,all[1:100,])


#all = all[all$Description]
# all = rbind(data1[1:10,],data2[1:10,],data3[1:10,],data4[1:10,])
# all = rbind(all,data11[1:10,],data22[1:10,],data33[1:10,],data44[1:10,])
#all$p.adjust[all$p.adjust>0.05] = 0
result$group = sub("^(.*?)_.+", "\\1", result$fungi)
result$time =factor(result$time, levels = c("3h", "6h","24h"))
ggplot(data= result,aes(x=group, y=Description, size=Count,color = p.adjust)) +
  geom_point(alpha=0.5) +
  scale_x_discrete(limits = c("1477", "2174", "2059",'ca'))+
  scale_color_gradient(low = "red", high = "darkblue")+
  # geom_vline(xintercept = 3.5)+
  # geom_vline(xintercept = 7.5)+
  facet_grid(.~time)+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
# + scale_fill_continuous(limits = c(0.01,0.05), breaks = c(0.01, 0.04, 0.41, 0.42, 0.05), guide = guide_colorbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE))

ggsave('E:/projects/王秀丽/results/fungi/human/dotpot1.pdf',width=10,height=10)

# 
# 
# data4 = read.delim('E:/projects/王秀丽/results/fungi/mouse/3h/1477_BP.txt',sep='\t')[,c(2,6,9)]
# data4$fungi = '1477_3h'
# data4$time = '3h_mm'
# data5 = read.delim('E:/projects/王秀丽/results/fungi/mouse/3h/2174_BP.txt',sep='\t')[,c(2,6,9)]
# data5$fungi = '2174_3h'
# data5$time = '3h_mm'
# data6 = read.delim('E:/projects/王秀丽/results/fungi/mouse/3h/2059_BP.txt',sep='\t')[,c(2,6,9)]
# data6$fungi = '2059_3h'
# data6$time = '3h_mm'
# data44 = read.delim('E:/projects/王秀丽/results/fungi/mouse/12h/1477_BP.txt',sep='\t')[,c(2,6,9)]
# data44$fungi = '1477_12h'
# data44$time = '12h_mm'
# data55 = read.delim('E:/projects/王秀丽/results/fungi/mouse/12h/2174_BP.txt',sep='\t')[,c(2,6,9)]
# data55$fungi = '2174_12h'
# data55$time = '12h_mm'
# data66 = read.delim('E:/projects/王秀丽/results/fungi/mouse/12h/2059_BP.txt',sep='\t')[,c(2,6,9)]
# data66$fungi = '2059_12h'
# data66$time = '12h_mm'
# 
# all = rbind(data4,data5,data6,data44,data55,data66)
# all$p.adjust[all$p.adjust>0.05] = 0
# all$group = sub("^(.*?)_.+", "\\1", all$fungi)
# all$time =factor(all$time, levels = c("3h_mm", "12h_mm"))
# dev.new()
# ggplot(data= all,aes(x=group, y=Description, size=Count, color=p.adjust)) +
#   geom_point(alpha=0.5) +
#   scale_x_discrete(limits = c("1477", "2174", "2059"))+
#   scale_color_gradient(low = "red", high = "darkblue")+
#   # geom_vline(xintercept = 3.5)+
#   # geom_vline(xintercept = 7.5)+
#   facet_grid(.~time)+
#   theme_bw()+theme(panel.grid = element_blank(),
#                    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
# dev.off()
# # + scale_fill_continuous(limits = c(0.01,0.05), breaks = c(0.01, 0.04, 0.41, 0.42, 0.05), guide = guide_colorbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# ggsave('E:/projects/王秀丽/results/fungi/mouse/dotpot1.pdf',width=5,height=5)


# data1 = read.delim('E:/projects/王秀丽/results/human/3h/1477_BP.txt',sep='\t')[,c(2,6)]
# data1$fungi = '1477_3h'
# data1$time = '3h'
# # data2 = read.delim('E:/projects/王秀丽/results/human/3h/2174_BP.txt',sep='\t')[,c(2,6)]
# # data2$fungi = '2174_3h'
# # data2$time = '3h'
# data3 = read.delim('E:/projects/王秀丽/results/human/3h/2059_BP.txt',sep='\t')[,c(2,6)]
# data3$fungi = '2059_3h'
# data3$time = '3h'
# data4 = read.delim('E:/projects/王秀丽/results/human/3h/ca_BP.txt',sep='\t')[,c(2,6)]
# data4$fungi = 'ca_3h'
# data4$time = '3h'
# data11 = read.delim('E:/projects/王秀丽/results/human/6h/1477_BP.txt',sep='\t')[,c(2,6)]
# data11$fungi = '1477_6h'
# data11$time = '6h'
# data22 = read.delim('E:/projects/王秀丽/results/human/6h/2174_BP.txt',sep='\t')[,c(2,6)]
# data22$fungi = '2174_6h'
# data22$time = '6h'
# data33 = read.delim('E:/projects/王秀丽/results/human/6h/2059_BP.txt',sep='\t')[,c(2,6)]
# data33$fungi = '2059_6h'
# data33$time = '6h'
# data44 = read.delim('E:/projects/王秀丽/results/human/ca/ca_BP.txt',sep='\t')[,c(2,6)]
# data44$fungi = 'ca_24h'
# data44$time = '24h'
# #all = rbind(data1,data2,data3,data11,data22,data33)
# all = rbind(data1[1:10,],data3[1:10,],data4[1:10,])
# all = rbind(all,data11[1:10,],data22[1:10,],data33[1:10,],data44[1:10,])
# all$p.adjust[all$p.adjust>0.05] = 0
# all$group = sub("^(.*?)_.+", "\\1", all$fungi)
# all$time =factor(all$time, levels = c("3h", "6h", "24h"))
# ggplot(data= all,aes(x=group, y=Description, size=0.5, color=p.adjust)) +
#   geom_point(alpha=0.5) +
#   scale_x_discrete(limits = c("1477", "2174", "2059","ca"))+
#   scale_color_gradient(low = "red", high = "darkblue")+
#   # geom_vline(xintercept = 4.5)+
#   # geom_vline(xintercept = 7.5)+
#   facet_grid(.~time)+
#   theme_bw()+theme(panel.grid = element_blank(),
#                    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
# # + scale_fill_continuous(limits = c(0.01,0.05), breaks = c(0.01, 0.04, 0.41, 0.42, 0.05), guide = guide_colorbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# ggsave('E:/projects/王秀丽/results/human/dotpot1.pdf',width=10,height=10)

# data1 = read.delim('E:/projects/王秀丽/results/mouse/3h/1477_BP.txt',sep='\t')[,c(2,6)]
# data1$fungi = '1477_3h'
# data1$time = '3h'
# data2 = read.delim('E:/projects/王秀丽/results/mouse/3h/2174_BP.txt',sep='\t')[,c(2,6)]
# data2$fungi = '2174_3h'
# data2$time = '3h'
# data3 = read.delim('E:/projects/王秀丽/results/mouse/3h/2059_BP.txt',sep='\t')[,c(2,6)]
# data3$fungi = '2059_3h'
# data3$time = '3h'
# data4 = read.delim('E:/projects/王秀丽/results/mouse/3h/ca_BP.txt',sep='\t')[,c(2,6)]
# data4$fungi = 'ca_3h'
# data4$time = '3h'
# data11 = read.delim('E:/projects/王秀丽/results/mouse/12h/1477_BP.txt',sep='\t')[,c(2,6)]
# data11$fungi = '1477_12h'
# data11$time = '12h'
# data22 = read.delim('E:/projects/王秀丽/results/mouse/12h/2174_BP.txt',sep='\t')[,c(2,6)]
# data22$fungi = '2174_12h'
# data22$time = '12h'
# data33 = read.delim('E:/projects/王秀丽/results/mouse/12h/2059_BP.txt',sep='\t')[,c(2,6)]
# data33$fungi = '2059_12h'
# data33$time = '12h'
# data44 = read.delim('E:/projects/王秀丽/results/mouse/ca/ca_BP.txt',sep='\t')[,c(2,6)]
# data44$fungi = 'ca_24h'
# data44$time = '24h'
# 
# #all = rbind(data4,data5,data6,data44,data55,data66)
# #all = rbind(data1,data2,data3,data11,data22,data33)
# all = rbind(data1[1:10,],data2[1:10,],data3[1:10,],data4[1:10,])
# all = rbind(all,data11[1:10,],data22[1:10,],data33[1:10,],data44[1:10,])
# all$p.adjust[all$p.adjust>0.05] = 0
# all$group = sub("^(.*?)_.+", "\\1", all$fungi)
# all$time =factor(all$time, levels = c("3h", "12h", "24h"))
# #c("1477_3h", "2174_3h", "2059_3h","ca_3h","1477_12h", "2174_12h", "2059_12h","ca_24h"))+
# ggplot(data= all,aes(x=group, y=Description, size=0.5, color=p.adjust)) +
#   geom_point(alpha=0.5) +
#   scale_x_discrete(limits = c("1477", "2174", "2059","ca"))+
#   scale_color_gradient(low = "red", high = "darkblue")+
#   # geom_vline(xintercept = 3.5)+
#   # geom_vline(xintercept = 7.5)+
#   facet_grid(.~time)+
#   theme_bw()+theme(panel.grid = element_blank(),
#                    axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
#  # + scale_fill_continuous(limits = c(0.01,0.05), breaks = c(0.01, 0.04, 0.41, 0.42, 0.05), guide = guide_colorbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE))
# 
# ggsave('E:/projects/王秀丽/results/mouse/dotpot1.pdf',width=10,height=10)




data1 = read.delim('E:/projects/王秀丽/3-14/results/human/3h/1477_BP.txt',sep='\t')[,c(2,6)]
data1$fungi = '1477_3h'
data1$time = '3h'
data2 = read.delim('E:/projects/王秀丽/3-14/results/human/3h/2174_BP.txt',sep='\t')[,c(2,6)]
data2$fungi = '2174_3h'
data2$time = '3h'
data3 = read.delim('E:/projects/王秀丽/3-14/results/human/3h/2059_BP.txt',sep='\t')[,c(2,6)]
data3$fungi = '2059_3h'
data3$time = '3h'
data11 = read.delim('E:/projects/王秀丽/3-14/results/human/6h/1477_BP.txt',sep='\t')[,c(2,6)]
data11$fungi = '1477_6h'
data11$time = '6h'
data22 = read.delim('E:/projects/王秀丽/3-14/results/human/6h/2174_BP.txt',sep='\t')[,c(2,6)]
data22$fungi = '2174_6h'
data22$time = '6h'
data33 = read.delim('E:/projects/王秀丽/3-14/results/human/6h/2059_BP.txt',sep='\t')[,c(2,6)]
data33$fungi = '2059_6h'
data33$time = '6h'

data44= read.delim('E:/projects/王秀丽/3-14/results/human/ca/ca_BP.txt',sep='\t')[,c(2,6)]
data44$fungi = 'ca_24h'
data44$time = '24h'

data4= read.delim('E:/projects/王秀丽/3-14/results/human/3h/ca_BP.txt',sep='\t')[,c(2,6)]
data4$fungi = 'ca_3h'
data4$time = '3h'



all = rbind(data1,data3,data11,data22,data33,data4,data44)

#result <- all[all$Description %in% mouse_host, ]
#result <- all[grepl(paste(mouse_host, collapse = "|"), all$Description), ]
#result <- all[grepl(paste(human_host, collapse = "|"), all$Description), ]
all = all[all$p.adjust<0.05,]
result = all[1:100,]


#all = all[all$Description]
# all = rbind(data1[1:10,],data2[1:10,],data3[1:10,],data4[1:10,])
#all = rbind(all,data11[1:10,],data22[1:10,],data33[1:10,],data44[1:10,])
#all$p.adjust[all$p.adjust>0.05] = 0
result$group = sub("^(.*?)_.+", "\\1", result$fungi)
result$time =factor(result$time, levels = c("3h", "6h","24h"))
ggplot(data= result,aes(x=group, y=Description, size=0.5,color = p.adjust)) +
  geom_point(alpha=0.5) +
  scale_x_discrete(limits = c("1477", "2174", "2059",'ca'))+
  scale_color_gradient(low = "red", high = "darkblue")+
  # geom_vline(xintercept = 3.5)+
  # geom_vline(xintercept = 7.5)+
  facet_grid(.~time)+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))

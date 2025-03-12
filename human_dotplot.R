data1 = read.delim('E:/projects/王秀丽/results/human/kegg/kegg_3h_1477.txt',sep='\t')[,c(4,11,13,14)]
data1$fungi = '1477_3h'
data1$time = '3h'
data2 = read.delim('E:/projects/王秀丽/results/human/kegg/kegg_3h_2174.txt',sep='\t')[,c(4,11,13,14)]
data2$fungi = '2174_3h'
data2$time = '3h'
data3 = read.delim('E:/projects/王秀丽/results/human/kegg/kegg_3h_2059.txt',sep='\t')[,c(4,11,13,14)]
data3$fungi = '2059_3h'
data3$time = '3h'
data11 = read.delim('E:/projects/王秀丽/results/human/kegg/kegg_6h_1477.txt',sep='\t')[,c(4,11,13,14)]
data11$fungi = '1477_6h'
data11$time = '6h'
data22 = read.delim('E:/projects/王秀丽/results/human/kegg/kegg_6h_2174.txt',sep='\t')[,c(4,11,13,14)]
data22$fungi = '2174_6h'
data22$time = '6h'
data33 = read.delim('E:/projects/王秀丽/results/human/kegg/kegg_6h_2059.txt',sep='\t')[,c(4,11,13,14)]
data33$fungi = '2059_6h'
data33$time = '6h'
data44= read.delim('E:/projects/王秀丽/results/human/kegg/kegg_24h_ca.txt',sep='\t')[,c(4,11,13,14)]
data44$fungi = 'ca_24h'
data44$time = '24h'

data4= read.delim('E:/projects/王秀丽/results/human/kegg/kegg_3h_ca.txt',sep='\t')[,c(4,11,13,14)]
data4$fungi = 'ca_3h'
data4$time = '3h'


all = rbind(data1,data2,data3,data11,data22,data33,data4,data44)

all = all[order(all$p.adjust),]
threshold <- 0.05
# # 按 pathway 分组，筛选在 s, d, f 三组中 p 值都显著的通路
# significant_pathways <- aggregate(p.adjust ~ Description, data = all, 
#                                   function(x) all(x < threshold))
# # 筛选出所有 p 值都显著的通路
# significant_pathways <- significant_pathways[significant_pathways$p.adjust, "Description"]
# mouse0 = c(mouse,significant_pathways)
# result <- all[all$Description %in% mouse0, ]
#result <- all[grepl(paste(mouse_host, collapse = "|"), all$Description), ]
#result <- all[grepl(paste(mouse_host, collapse = "|"), all$Description), ]
library(ggplot2)
all = all[all$p.adjust<0.05,]
result = all
#result = all[grepl('12443',all$geneID),]
#all = all[all$Description]
# all = rbind(data1[1:10,],data2[1:10,],data3[1:10,],data4[1:10,])
# all = rbind(all,data11[1:10,],data22[1:10,],data33[1:10,],data44[1:10,])
# result$p.adjust[result$p.adjust>0.05] = 0
# result$Count[result$p.adjust == 0] <- 0
result$group = sub("^(.*?)_.+", "\\1", result$fungi)
result$time =factor(result$time, levels = c("3h", "6h","24h"))
ggplot(data= result,aes(x=group, y=Description, size=Count,,color = p.adjust)) +
  geom_point(alpha=0.5) +
  scale_x_discrete(limits = c("1477", "2174", "2059",'ca'))+
  scale_color_gradient(low = "red", high = "darkblue")+
  # geom_vline(xintercept = 3.5)+
  # geom_vline(xintercept = 7.5)+
  facet_grid(.~time)+
  theme_bw()+theme(panel.grid = element_blank(),
                   axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
write.table(all,'./kegg/kegg_human_all.txt',sep='\t',row.names = F)

# setwd("/disk191_3/AnalysisPipline/Pig/1_DNAseq/Chip_2/5_NJ")
# BiocManager::install("ggtree")
arg <- commandArgs(trailingOnly = TRUE)
tissue = arg[1]
print(tissue)
library(ggtree)
library(ggplot2) 

print(paste0("/disk201/zzy/1_GenomicSelection/8-TWASresearch/2.predictedExpression/3.Tkinship/Result/NJ/Pig_ELA_",tissue,"_neighbor.outtree"))
tree <- read.tree(paste0("/disk201/zzy/1_GenomicSelection/8-TWASresearch/2.predictedExpression/3.Tkinship/Result/NJ/Pig_ELA_",tissue,"_neighbor.outtree"))

# 读取分组信息
# group_file <- read.table("/disk201/zzy/1_GenomicSelection/8-TWASresearch/2.predictedExpression/3.Tkinship/Result/NJ/groupinfo.txt",header = F,row.names = 1)
# # 按类分组
# groupInfo <- split(row.names(group_file), group_file$V2)


group_file <- read.table("/disk201/zzy/1_GenomicSelection/8-TWASresearch/2.predictedExpression/3.Tkinship/Result/NJ/groupinfo.txt",header = F,row.names = 1)
# 按类分组
groupInfo <- split(row.names(group_file), group_file$V2)
# 将分组信息添加到树中
tree <- groupOTU(tree, groupInfo)


# 绘制进化树
# pdf
ggtree(tree, layout="circular", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right") +scale_color_manual(values =c("#EEAD0E","#6cd3f5","#EEB4B4","#4C6E84"))

ggsave(file=paste0("/disk201/zzy/1_GenomicSelection/8-TWASresearch/2.predictedExpression/3.Tkinship/Result/NJ/Pig_ELA_",tissue,"_neighbor.outtree.pdf"), height=8, width=8)
# rectangular
# 绘制进化树
ggtree(tree, layout="rectangular", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")+scale_color_manual(values =c("#EEAD0E","#6cd3f5","#EEB4B4","#4C6E84"))
# +scale_color_manual(values =c("#33a02c", '#AB274F',"#b2df8a","#009ACD","#EEAD0E","#6cd3f5","#EEB4B4","#4C6E84","#CD5E14","191970","#8B5A2B","#C276A0","#249A70","#698B22"))

ggsave(file=paste0("/disk201/zzy/1_GenomicSelection/8-TWASresearch/2.predictedExpression/3.Tkinship/Result/NJ/Pig_ELA_",tissue,"_neighbor.tree_rectangular.pdf"), height=10, width=10)

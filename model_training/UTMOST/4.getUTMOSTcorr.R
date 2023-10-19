library(data.table)
library(parallel)
geneinfo_path<-"path/gene_info.txt" # one column, gene list
est_file_path <- "path/UTMOST/Result/" # output byCTIMP
setwd("")
# ----------1.获取基因列表----------

geneinfo<-read.csv(geneinfo_path)
gene<- geneinfo[1,1]

##计算rsq的函数
rsq_create<- function(gene){
  library(data.table)
  # cat(gene)
  geneinfo_path<-"path/gene_info.txt" # one column, gene list
  est_file_path <- "path/UTMOST/Result/" # output byCTIMP
  # ----------2.read est ----------
  est_file <- paste0(est_file_path,gene,".est")
  est_df = as.data.frame(fread(est_file))
  tissue_vec = colnames(est_df)[4:ncol(est_df)]
  tissue_vec = as.character(sapply(tissue_vec, function(x) unlist(strsplit(x, "\\."))[1]))
  
  # ----------3.read cv.evaluation.RData and calculate correlation----------
  cv.evaluation.RData <- paste0(est_file_path,gene,".cv.evaluation.RData")
  load(cv.evaluation.RData)

  N.tissue = length(multi_res_test[[1]])
    cv_r<-cv_p<-c()
    # 合并五倍交叉验证结果
    for (j in 1:N.tissue){
        tmp <- data.frame()
        for (i in 1:5){
                tmp = rbind(tmp,multi_res_test[[i]][[j]])
            }
        fit<-cor.test(tmp[,1],tmp[,2])
        cv_r[j]<-as.numeric(fit$estimate)^2
        cv_p[j]<-as.numeric(fit$p.value)
    }

  res_df <- data.frame("tissues"=tissue_vec,"rsq"=cv_r,"pv" = cv_p)
  ##output
  
  write.csv(res_df,paste0("path/UTMOST/Performance/GenePerformance/",gene,".rsq"),row.names = F,quote = F)
}
rsq_create(gene)

# 计算可用线程数，并设置并行使用线程数
no_cores <-30
# 初始化
cl <- makeCluster(no_cores)

parLapply(cl, geneinfo$gene_ensg,  rsq_create)



# ----------4.obtain the rsq for each tissue----------
splitTissue_ <- function(df){
  tmp <- as.vector(df)
  write.table(data.frame(tmp[4],tmp[2],tmp[3]),paste0("path/UTMOST/Performance/TissuePerformance/",tmp[1],".rsq"),row.names = F,col.names = F,quote = F,sep = ",",append = T)
}

splitTissue<- function(gene){
  rsq_df = as.data.frame(fread(paste0("path/UTMOST/Performance/GenePerformance/",gene,".rsq")))
  rsq_df$gene <- gene
  apply(rsq_df,1, splitTissue_)
}

lapply( geneinfo$gene_ensg,  splitTissue)







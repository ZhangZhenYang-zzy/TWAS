##整和UTMOST的结果
##需要猪的基因信息
##需要猪的位点SNP信息
##样品个数

library(data.table)
library(dplyr)

sampleSIZE <- as.data.frame(fread("/disk191/zzy/0_work/11_animalTWAS/paper/statistics/2_tissue/samplesizePIG.csv"))

getwd()
setwd("/disk191/zzy/.0_TWAS/.0_data/eqtlTRAIN/Pig/UTMOST/Table")
est_file_dir = "/disk191/zzy/.0_TWAS/.0_data/eqtlTRAIN/Pig/UTMOST/Result"
# est_file_list = list.files(est_file_dir)


weight_list = list()
weight_snps = list()
construction_list = list()
sample_list = list()

gene_vec = c()


est_file_list = list.files(est_file_dir, pattern = "*.est")

est_files = paste0(est_file_dir, "/", est_file_list)
idx=1

for (idx in 1:length(est_files)) {
  gene = unlist(strsplit(est_file_list[idx],"\\."))[1]
  est_file = est_files[idx]
  est_df = as.data.frame(fread(est_file))
  tissue_vec = colnames(est_df)[4:ncol(est_df)]
  tissue_vec = as.character(sapply(tissue_vec, function(x) unlist(strsplit(x, "\\."))[1]))
  
  for (tissue_idx in 1:length(tissue_vec)) {
    tissue = tissue_vec[tissue_idx]
    #IRdisplay::display_html(tissue)
    # weight table
    if (!tissue %in% names(weight_list)) {
      weight_list[[tissue]] = data.frame(rsid = est_df$SNP,
                                         gene = gene,
                                         weight = est_df[, tissue_idx + 3],
                                         ref_allele = est_df$REF.0.,
                                         eff_allele = est_df$ALT.1.
      ) %>% 
        filter(weight != 0)
      
      
    } else {
      tmp_df  = data.frame(rsid = est_df$SNP,
                           gene = gene,
                           weight = est_df[, tissue_idx +3],
                           ref_allele = est_df$REF.0.,
                           eff_allele = est_df$ALT.1.
      ) %>% 
        filter(weight != 0)
      
      weight_list[[tissue]] = rbind(weight_list[[tissue]], tmp_df)
    }
    
    # extra table
    if (!tissue %in% names(weight_snps)) {
      weight_snps[[tissue]] = data.frame(gene = gene,
                                         genename = gene,
                                         pred.perf.R2 = 0.5,
                                         n.snps.in.model = nrow(est_df),
                                         pred.perf.pval = 1e-10,
                                         pred.perf.qval = 1e-10
      )
      
    } else {
      tmp_df  = data.frame(gene = gene,
                           genename = gene,
                           pred.perf.R2 = 0.5,
                           n.snps.in.model = nrow(est_df),
                           pred.perf.pval = 1e-10,
                           pred.perf.qval = 1e-10
      )
      
      weight_snps[[tissue]] = rbind(weight_snps[[tissue]], tmp_df)
      
      
      
    }
    
    # construction table
    construction_list[[tissue]] = data.frame(chr = 1:18, cv.seed = rep(1, 18))
    
    # sample table
    sample_list[[tissue]] = data.frame(n.samples = sampleSIZE[sampleSIZE$Tissue==tissue,"SampleSize"])
    
  }
}

# output 
output_dir <- "/disk191/zzy/.0_TWAS/.0_data/eqtlTRAIN/Pig/UTMOST/Table"
for (tissue in names(weight_list)) {
    weight_file = paste0(output_dir, "/", tissue, ".weight.txt")
    extra_file = paste0(output_dir, "/", tissue, ".extra.txt")
    construction_file = paste0(output_dir, "/", tissue, ".construction.txt")
    sample_file = paste0(output_dir, "/", tissue, ".sample.txt")
    write.table(weight_list[[tissue]], file = weight_file, quote = F, row.names = F)
    write.table(weight_snps[[tissue]], file = extra_file, quote = F, row.names = F)
    write.table(construction_list[[tissue]], file = construction_file, quote = F, row.names = F)
    write.table(sample_list[[tissue]], file = sample_file, quote = F, row.names = F)
}

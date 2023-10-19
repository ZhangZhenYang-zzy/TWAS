args=commandArgs(T)
tissue = args[1]

tissues = c(tissue)
print(tissue)


library(data.table)
out_tissue <- c()
geneCount <- c()

for (tissue in tissues) {

	genelist <- list.files(path=paste0("/disk191_1/zzy/TWAS/.0_data/pigEQTL/FUSION/",tissue,"/"),pattern="*.RDat")
	out_tissue <- c(out_tissue,tissue)
	geneCount <- c(geneCount,length(genelist))
    genenames <- genelist[1]
    # merge all the gene result

    top1 <- data.frame()
    blup <- data.frame()
    bslmm <- data.frame()

	for (genenames in genelist) {
		genename <- gsub(".wgt.RDat","",genenames)
		# print(genename)
		load(paste0("/disk191_1/zzy/TWAS/.0_data/pigEQTL/FUSION/",tissue,"/",genename,".wgt.RDat"))
        # 获取权重  rsid gene weight ref_allele eff_allele
        # 三个model
        colnames(snps) = c("chr","chr_pos","V3","position","ref_allele","eff_allele")
        df = merge(wgt.matrix,snps,by.x=0,by.y="chr_pos")
        df["gene"] = genename

# .weight.txt
        top1_t = df[,c("Row.names","gene","top1","ref_allele","eff_allele")]

        blup_t = df[,c("Row.names","gene","blup","ref_allele","eff_allele")]

        bslmm_t = df[,c("Row.names","gene","bslmm","ref_allele","eff_allele")]

        top1 <- rbind(top1,top1_t)
        blup <- rbind(blup,blup_t)
        bslmm <- rbind(bslmm,bslmm_t)
	}
    colnames(top1) <- c("rsid","gene","weight","ref_allele","eff_allele")
    colnames(blup) <- c("rsid","gene","weight","ref_allele","eff_allele")
    colnames(bslmm) <- c("rsid","gene","weight","ref_allele","eff_allele")

    extrafile = paste0("/disk201/zzy/1_GenomicSelection/8-TWASresearch/1.eQTLsummary/2.TissueSummary/FUSION/TissuePerformance/",tissue,".res")
    df_extra = read.csv(extrafile,sep="\t")
    # head(df_extra)

    top1_extra = data.frame(
        gene=df_extra$id,
        genename=df_extra$id,
        pred.perf.R2=df_extra$top1.adj.r2,
        n.snps.in.model=df_extra$nsnps,
        pred.perf.pval=df_extra$top1.pv,
        pred.perf.qval=df_extra$top1.pv
    )
    blup_extra = data.frame(
        gene=df_extra$id,
        genename=df_extra$id,
        pred.perf.R2=df_extra$blup.adj.r2,
        n.snps.in.model=df_extra$nsnps,
        pred.perf.pval=df_extra$blup.pv,
        pred.perf.qval=df_extra$blup.pv
    )
    bslmm_extra = data.frame(
        gene=df_extra$id,
        genename=df_extra$id,
        pred.perf.R2=df_extra$bslmm.adj.r2,
        n.snps.in.model=df_extra$nsnps,
        pred.perf.pval=df_extra$bslmm.pv,
        pred.perf.qval=df_extra$bslmm.pv
    )
    construction = data.frame(
        chr=1:18,
        cv.seed=rep(1, 18)
    )
    sampleSIZE <- as.data.frame(fread("/disk191/zzy/0_work/11_animalTWAS/paper/statistics/2_tissue/samplesizePIG.csv"))

    sample_tissue = data.frame(n.samples = sampleSIZE[sampleSIZE$Tissue==tissue,"SampleSize"])
    output_dir <- "Table/"
    weight_file1 = paste0(output_dir, "/", tissue, ".top1.weight.txt")
    extra_file1 = paste0(output_dir, "/", tissue, ".top1.extra.txt")

    weight_file2 = paste0(output_dir, "/", tissue, ".blup.weight.txt")
    extra_file2 = paste0(output_dir, "/", tissue, ".blup.extra.txt")

    weight_file3 = paste0(output_dir, "/", tissue, ".bslmm.weight.txt")
    extra_file3 = paste0(output_dir, "/", tissue, ".bslmm.extra.txt")


    construction_file = paste0(output_dir, "/", tissue, ".construction.txt")
    sample_file = paste0(output_dir, "/", tissue, ".sample.txt")

    write.table(na.omit(top1), file = weight_file1, quote = F, row.names = F)
    write.table(top1_extra, file = extra_file1, quote = F, row.names = F)

    write.table(na.omit(blup), file = weight_file2, quote = F, row.names = F)
    write.table(blup_extra, file = extra_file2, quote = F, row.names = F)

    write.table(na.omit(bslmm), file = weight_file3, quote = F, row.names = F)
    write.table(bslmm_extra, file = extra_file3, quote = F, row.names = F)

    write.table(construction, file = construction_file, quote = F, row.names = F)
    write.table(sample_tissue, file = sample_file, quote = F, row.names = F)

    # /disk201/zzy/1_GenomicSelection/8-TWASresearch/1.eQTLsummary/2.TissueSummary/FUSION/TissuePerformance/Adipose.res

}


# out <- data.frame("GeneList"=out_tissue,"geneCount"=geneCount)

# write.csv(out,paste0("Fusion/",tissue,".count.csv"),quote=F)

#!/usr/bin/Rscript
# GO KEGG GSEA.
start.time <- Sys.time()
suppressMessages(library("optparse"))

option_list = list(
make_option("--twas", action="store", default=NA, type='character',
		help="Path to genome-wide TWAS results [required]"),
make_option("--sig_p", action="store", default=NA, type='numeric',
		help="Threshold for filter result"),
make_option("--species", action="store", default=NA, type='character',
		help="the specie input"),
make_option("--output", action="store", default=NA, type='character',
		help="Path to save output [required]"),
make_option("--width", action="store", default=8, type='numeric',
		help="width of plot [optional]"),
make_option("--height", action="store", default=5, type='numeric',
		help="height of plot [optional]"),
make_option("--res", action="store", default=300, type='numeric',
		help="height of plot [optional]")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load libraries
suppressMessages(library(data.table))


suppressMessages(library(clusterProfiler))
suppressMessages(library(ggplot2))
suppressMessages(library(enrichplot))
suppressMessages(library(KEGG.db))
suppressMessages(library(KEGGREST))

# organism=c("org.Hs.eg.db","org.Ss.eg.db","org.Bt.eg.db")
# organism=c("org.Ss.eg.db","org.Bt.eg.db")
# BiocManager::install("org.Bt.eg.db")

twas<-data.frame(fread(opt$twas))
twas <- twas[!is.na(twas$GeneID),]
twas <- twas[twas$GeneID!="",]


suppressMessages(require(DOSE))

##Species 
if (opt$species == "pig"){
    organism = "org.Ss.eg.db"
    organismKEGG = "ssc"
    suppressMessages(library("org.Ss.eg.db"))
}
if (opt$species == "human"){
    organism = "org.Hs.eg.db"
    organismKEGG = "hsa"
    suppressMessages(library("org.Hs.eg.db"))
}
if (opt$species == "cow"){
    organism = "org.Bt.eg.db"
    organismKEGG = "bta"
    suppressMessages(library("org.Bt.eg.db"))
}

# library( "org.Hs.eg.db")
if(!is.na(opt$sig_p)){
    ##提取twas的genelist
	thresh <- opt$sig_p/nrow(twas)
	# thresh <- 0.05/nrow(twas)
    # twas$GeneName
    # 如果品种是人 需要转换一下id 去除版本号
    if (opt$species == "human") {
       twas$GeneID <- gsub("\\..*", "",  twas$GeneID)
    }
    # 转换id
    if(opt$species == "pig") {
    gene.ENTREZID <- bitr(geneID = twas$GeneName, 
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = organism)
    original_gene_list <- twas[match(gene.ENTREZID$SYMBOL,twas$GeneName),"TWAS.P"] 
    names(original_gene_list)  <- gene.ENTREZID$ENTREZID

    }else{
            gene.ENTREZID <- bitr(geneID = twas$GeneID, 
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = organism)
    original_gene_list <- twas[match(gene.ENTREZID$ENSEMBL,twas$GeneID),"TWAS.P"] 
    names(original_gene_list)  <- gene.ENTREZID$ENTREZID
    }


    # print(dim(df2))
    # print(dim(dedup_ids))

    gene_list<-na.omit(original_gene_list)
    gene_list<- -log10(gene_list)

    gene_list = sort(gene_list, decreasing = TRUE)
    # print(gene_list)

    # go分析
    gse <- gseGO(geneList=gene_list, 
                ont ="ALL", 
                keyType = "ENTREZID", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 1, 
                verbose = TRUE, 
                OrgDb = organism, 
                scoreType = "pos",
                pAdjustMethod = "BH")
    # gse <- gseGO(geneList=gene_list, 
    #             ont ="ALL", 
    #             keyType = "ENTREZID", 
    #             minGSSize = 3, 
    #             maxGSSize = 800, 
    #             pvalueCutoff = 0.99, 
    #             verbose = TRUE, 
    #             OrgDb = organism, 
    #             scoreType = "pos",
    #             pAdjustMethod = "none")
    result<- data.frame(
        ONTOLOGY = gse$ONTOLOGY,
        ID = gse$ID,
        Description = gse$Description,
        pvalue = gse$pvalue,
        p.adjust = gse$p.adjust,
        qvalues =gse$qvalues,
        setSize = gse$setSize,
        enrichmentScore =  gse$enrichmentScore)

    if(nrow(result) > 0 ){
        write.csv(result,paste0(opt$output,'.GO.csv'), row.names=FALSE,quote=FALSE)
        # 画图第一版
        # pdf(paste0(opt$output,'.GO.pdf'),  width = opt$width, height = opt$height)
        # print(dotplot(gse, showCategory=10))
        # dev.off()
        # 画图第二版
        result<- result[order(result$pvalue),]
        result<- result[!is.na(result$pvalue),]
        result_plot1 <- result[result$ONTOLOGY == "BP",]
        result_plot1 <- result_plot1[1:10,]
        result_plot2 <- result[result$ONTOLOGY == "MF",]
        result_plot2 <- result_plot2[1:10,]
        result_plot3 <- result[result$ONTOLOGY == "CC",]
        result_plot3 <- result_plot3[1:10,]
        result_plot <- rbind(result_plot1,result_plot2,result_plot3)
        pdf(paste0(opt$output,'.GO.pdf'),  width = opt$width, height = opt$height)
        # print(dotplot(gse, showCategory=10)+facet_grid(.~ONTOLOGY))
        p<- ggplot(result_plot,aes(y=enrichmentScore,x=Description,fill=pvalue)) + 
          geom_bar(stat="identity",position = "dodge") +
          facet_grid(ONTOLOGY~.,scales = "free",space = "free") + 
          coord_flip() + 
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5),
                strip.text.y = element_text(size = 14),
                legend.position="right",
                legend.title = element_text(size=18),
                legend.text = element_text(size=14),
                # axis.text.x = element_text(size=14),
                # axis.text.y = element_text(size=18),
                axis.title.x = element_text(size=14),
                axis.title.y = element_text(size=14))

        print(p)
        dev.off()
        # +facet_grid(.~ONTOLOGY)
    }else{
        print("GO NOT OK")
    }


# 转换一下id
    # ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

    # dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
    # df2 = twas[twas$GeneName %in% dedup_ids$SYMBOL,]
    # # print(dim(df2))
    # # print(dim(dedup_ids))
    # df2[match(dedup_ids$SYMBOL,df2$GeneName),"ENTREZID"] <- dedup_ids$ENTREZID
    # df2 <- df2[!duplicated(df2["ENTREZID"]),]
    # kegg_gene_list <- df2$TWAS.P
    # names(kegg_gene_list) <- df2$ENTREZID
    # kegg_gene_list<-na.omit(kegg_gene_list)
    # kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# print(gene_list)
    kk <- gseKEGG(geneList    = gene_list,
                minGSSize    = 3,
                maxGSSize    = 800,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                organism = organismKEGG,
                keyType = "ncbi-geneid",
                scoreType = "pos",use_internal_data =T)


# kk$Pathway<- apply(kk,1,function(x) keggGet(x["ID"])[[1]]$NAME)
# print(kk)
   tmp <- lapply(kk$ID,function(x) keggGet(x[1])[[1]]$NAME)
   tmp2 <- lapply(kk$core_enrichment,function(x) length(unlist(strsplit(x,split="/"))))
    # print(kk)
    result<- data.frame(
        ID = kk$ID,
        Description = unlist(tmp),
        enrichmentScore=kk$enrichmentScore,
        pvalue = kk$pvalue,
        p.adjust = kk$p.adjust,
        qvalues = kk$qvalues,
        core_enrichment=kk$core_enrichment,
        count = unlist(tmp2)
        )

# bta hsa
# 气泡图
    if(nrow(result) > 0 ){
        write.csv(result,paste0(opt$output,'.KEGG.csv'), row.names=FALSE,quote=FALSE)
        # # 画图
        # pdf(paste0(opt$output,'.KEGG.pdf'),  width = 5, height = opt$height)
        # print(dotplot(kk, showCategory = 10))

        # dev.off()
        ###画图
        # 按照显著性水平排序
        result <- result[order(result$pvalue),]

        # 取前20个通路绘制气泡图
        top20 <- head(result, 20)
        top20$Description <- apply(top20,1,function(x) unlist(strsplit(x["Description"],split=" - "))[1])
        # pdf(paste0(opt$output,'.KEGG.pdf'),  width = 5, height = opt$height)
        p <- ggplot(data = top20,mapping = aes(x = enrichmentScore,y = Description))+
        geom_point(aes(color= -log10(pvalue),size = count)) +
        scale_colour_gradient(high = '#d22d46',low = '#54ab77') +
        theme_bw()+
        labs(title = "KEGG Pathway Enrichment(top20)",
            x = 'enrichmentScore',
            y = 'Description')
        ggsave(paste0(opt$output,'.KEGG.pdf'),p,width = opt$width,height = opt$height)
        # dev.off()
        # +facet_grid(.~ONTOLOGY)
    }else{
         print("KEGG NOT OK")
    }

}

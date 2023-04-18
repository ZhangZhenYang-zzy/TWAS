import pandas as pd 
import threading
import os
import polars as pl 
import sys


startgene = sys.argv[1]
endgene = sys.argv[2]
# 基因型临时目录
# 
tmp_dir = "tmp_doseage"
if not os.path.exists(tmp_dir):
	os.mkdir(tmp_dir)

def ctimp(geneid,CHR,pos):
	genename = geneid
	print(genename)
	CHR = CHR
	pos = pos
	p1 = pos - 1000000
	if p1 <0:
		p1 = 0
	p2 = pos + 1000000
	# recode A
	cmd = "plink --cow --bfile cattle --make-bed --recode A --out {0}/{1}.tmp  --chr {2} --from-bp {3} --to-bp {4}".format(tmp_dir,genename,CHR,p1,p2)
	os.system(cmd)
	bim = pl.read_csv("{0}/{1}.tmp.bim".format(tmp_dir,genename),sep="\t",has_header=False)
	bim = bim[["column_2","column_5","column_6"]]
	bim = bim.rename({"column_2":"SNP","column_5":"REF.0.","column_6":"ALT.1."})
	INFO = "{0}/{1}.info".format(tmp_dir,genename)
	# bim.write_csv(INFO,sep="\t")
	bim.to_pandas().to_csv(INFO,sep="\t",index=None)
	doseageData = pl.read_csv("{0}/{1}.tmp.raw".format(tmp_dir,genename),sep=" ")
	# 去除多余的列
	doseageData = doseageData.drop(['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
	X = "{0}/{1}.X".format(tmp_dir,genename)
	# doseageData.write_csv(X,sep=" ",has_header=False)
	doseageData.to_pandas().to_csv(X,sep=" ",header=None,index=None)
	Y = "%s"%genename
	output_path = "Result"
	if not os.path.exists(output_path):
		os.mkdir(output_path)
	output_prefix = genename
	ntune=5 
	cmd2 = "Rscript CTIMP/main.R {0} {1} {2} 5 {3} {4}".format(X,INFO,Y,output_prefix,output_path)
	os.system(cmd2)
	os.system("rm {0}/{1}.tmp* & ".format(tmp_dir,genename))

gene_df = pd.read_csv("Gene_list_cattle.txt",sep=" ")
batch_start = int(startgene)
batch_end = int(endgene)
ts = []
j = 0
for i in range(batch_start,batch_end):
	if j<5:
		if os.path.exists("Result/{0}.est".format(gene_df.loc[i,"gene_id"])):
			print("Result/{0}.est".format(gene_df.loc[i,"gene_id"]),"exists")
		# if os.path.exists("tmp_doseage/{0}.X".format(gene_df.loc[i,"gene_id"])):
		# 	print("tmp_doseage/{0}.X".format(gene_df.loc[i,"gene_id"]),"exists")
			continue
		else:
			th = threading.Thread(target=ctimp, args=[gene_df.loc[i,"gene_id"],gene_df.loc[i,"#Chr"],gene_df.loc[i,"end"]])
			ts.append(th)
			j = j + 1
	else:
		for t in ts:
			t.start()
		for t in ts:
			t.join()
		ts = []
		j = 0

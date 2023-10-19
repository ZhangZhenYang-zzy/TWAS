import pandas as pd 
import threading
import os
import polars as pl 
import sys


batch_start = sys.argv[1]
batch_end = sys.argv[2]


# a dir store the genotype
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
	# extract the region(TSS-1MB & TSS + 1MB).
	cmd = "plink --bfile pig --make-bed --recode A --out {0}/{1}.tmp  --chr {2} --from-bp {3} --to-bp {4}".format(tmp_dir,genename,CHR,p1,p2)
	os.system(cmd)
	bim = pl.read_csv("{0}/{1}.tmp.bim".format(tmp_dir,genename),sep="\t",has_header=False)
	bim = bim[["column_2","column_5","column_6"]]
	bim = bim.rename({"column_2":"SNP","column_5":"REF.0.","column_6":"ALT.1."})
	INFO = "{0}/{1}.info".format(tmp_dir,genename)
	bim.write_csv(INFO,sep="\t")
	doseageData = pl.read_csv("{0}/{1}.tmp.raw".format(tmp_dir,genename),sep=" ")
	doseageData = doseageData.drop(['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
	X = "{0}/{1}.X".format(tmp_dir,genename)
	doseageData.write_csv(X,sep=" ",has_header=False)
	Y = "path/2_expression_eachtissue/%s"%genename  #TODO result by 2_getExpressionEachtissue.py
	output_path = "Result"
	if not os.path.exists(output_path):
		os.mkdir(output_path)
	output_prefix = genename
	ntune=5 
	cmd2 = "Rscript CTIMP/main.R {0} {1} {2} 5 {3} {4}".format(X,INFO,Y,output_prefix,output_path)
	os.system(cmd2)
	os.system("rm {0}/{1}.tmp*".format(tmp_dir,genename))

gene_df = pd.read_csv("path/Gene_list.txt",sep=" ") #TODO result by 1_getGenelist.py


ts = []
j = 0
for i in range(batch_start,batch_end):
	if j<10 :
		if os.path.exists("{0}/{1}.X".format(tmp_dir,gene_df.loc[i,"gene_id"])):
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
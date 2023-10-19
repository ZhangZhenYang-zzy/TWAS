import pandas as pd 
import os 

#pig
tissues = ['Heart','Morula','Uterus','Artery','Frontal_cortex','Small_intestine','Milk','Testis','Jejunum','Oocyte','Kidney','Synovial_membrane','Placenta','Lymph_node','Muscle','Blastocyst','Ovary','Cartilage','Large_intestine','Adipose','Brain','Embryo','Hypothalamus','Fetal_thymus','Duodenum','Blastomere','Blood','Pituitary','Spleen','Colon','Lung','Macrophage','Ileum','Liver']

#expression 
file_list = ["path/expression_tmm_inv/"+x+".expr_tmm_inv.bed.gz" for x in tissues]

genelist = pd.read_csv("path/Gene_list.txt",sep=" ")

for i in genelist.index:
	geneNAME = genelist.loc[i,"gene_id"]
	CHR = genelist.loc[i,"#Chr"]
	Start = genelist.loc[i,"start"]
	End = genelist.loc[i,"end"]
	print(geneNAME)
	if not os.path.exists(geneNAME):
		os.mkdir(geneNAME)
	#或者每个组织里该组织的表达量
	for tiss in tissues:
		tmp_gene = ['']
		x = "path/expression_tmm_inv/"+tiss+".expr_tmm_inv.bed.gz"
		header = os.popen("zcat %s|head -n 1"%x).read().strip().split("\t")
		tmp_gene = os.popen("tabix %s %s:%s-%s"%(x,str(CHR),str(Start),str(End))).read().strip().split("\t")
		if tmp_gene == ['']:    ##表示该组织没有这个基因
			pass
		else:
			try:
				res = pd.DataFrame({"0":header,"1":tmp_gene})
				res.iloc[4:,].to_csv("%s/%s.txt"%(geneNAME,tiss),index=0,header=0,sep=" ")
			except:
				print(tiss)
				print(header)
				print(tmp_gene)

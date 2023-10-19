import os
import pandas as pd
##获取基因列表

#pig
tissues = ['Heart','Morula','Uterus','Artery','Frontal_cortex','Small_intestine','Milk','Testis','Jejunum','Oocyte','Kidney','Synovial_membrane','Placenta','Lymph_node','Muscle','Blastocyst','Ovary','Cartilage','Large_intestine','Adipose','Brain','Embryo','Hypothalamus','Fetal_thymus','Duodenum','Blastomere','Blood','Pituitary','Spleen','Colon','Lung','Macrophage','Ileum','Liver']

file_list = ["path/expression_tmm_inv/"+x+".expr_tmm_inv.bed.gz" for x in tissues]

for file in file_list:
	os.system("zcat %s |awk '{print $1,$2,$3,$4}'>>ALLgene.txt"%file)
df = pd.read_csv("ALLgene.txt",sep=" ")
df_rmdup = df.drop_duplicates(keep="first")
df_rmdup = df_rmdup[df_rmdup["#Chr"]!="#Chr"]
df_rmdup.to_csv("Gene_list.txt",sep=" ",index=None)



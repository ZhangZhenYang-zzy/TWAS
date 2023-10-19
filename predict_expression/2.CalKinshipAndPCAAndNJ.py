import pandas as pd 
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os 
import matplotlib
matplotlib.use("nbagg")
#读取数据

# tissue = "Adipose"
tissue = sys.argv[1]

filename = "result/%s_predict.txt"%tissue ### TODO 记得修改输入文件的路径

df =pd.read_csv(filename,sep="\t")


##Keep the individuals
# idkeep =  pd.read_csv("GENOID_DMU",sep=" ",header=None)
# df = df[df["IID"].isin(idkeep[0])]  # 过滤个体


ids = df["IID"]
value_matrix = df.iloc[:,2:df.shape[1]].values
## Normalization
mean = np.mean(value_matrix)
std = np.std(value_matrix)
value_matrix = (value_matrix-mean)/std
kinship = np.dot(value_matrix,value_matrix.T)/value_matrix.shape[0]

ids_array = ids.values
indexname = np.triu_indices(ids_array.shape[0])
id1 = ids_array[indexname[0]]
id2 = ids_array[indexname[1]]
res = pd.DataFrame({
    'id1':id1,
    'id2':id2,
    'value':list(kinship[np.triu_indices(ids_array.shape[0])])
})


res.to_csv("Pig_ELA_TKinship_%s.txt"%tissue,index=None,sep="\t")  ### TODO Change the output name

###PCA
from sklearn.decomposition import PCA
meta = pd.read_csv("1.UsingID.meta",sep="\t")
meta["Breed"] = meta["Breed"].map({"D":"S1","TB":"S2"})

pca=PCA(n_components=2)
pca.fit(kinship)
pca_res = pca.transform(kinship)
pc_res = pd.DataFrame({
    "ID":ids_array,
    "PC1":pca_res[:,0],
    "PC2":pca_res[:,1]
})
pc_res2 = pd.merge(pc_res,meta[["ID1","Herd","Sex","Breed"]],left_on="ID",right_on="ID1")
pc_res2 = pc_res2[["ID","Herd","Sex","PC1","PC2","Breed"]]
pc_res2.to_csv("PCA/Pig_ELA_PCA_%s.txt"%tissue,index=None,sep="\t")  ### TODO 记得修改输出文件名

ratio=pca.explained_variance_ratio_

sns.scatterplot(x="PC1",y="PC2",data=pc_res2,hue="Herd",style="Breed")
plt.xlabel("PC1(%s%%)"%str(round(ratio[0]*100,2)))
plt.ylabel("PC2(%s%%)"%str(round(ratio[1]*100,2)))
plt.savefig("Result/PCA/Pig_ELA_PCA_%s.pdf"%tissue,dpi=300)

## NJtree
### 距离矩阵

mdist = 1 - kinship
# mdist
pc_res2["SeqID"] = ["%04d"%x for x in range(1,pc_res2.shape[0]+1)]
pc_res2["NewID"] = pc_res2["Herd"].str[2:4]+"-" + pc_res2["Breed"] + "-" + pc_res2["SeqID"]
# pc_res2["NewID"] = ["{:>15}".format(x) for x in pc_res2["NewID"].values]
mdist_df = pd.concat((pc_res2["NewID"],pd.DataFrame(mdist)),axis=1)
mdist_df
output = "Result/NJ/Pig_ELA_%s.mdist"%tissue
with open(output,"w") as f:
    f.write("  %d\n"%mdist.shape[0])
    for i in range(mdist.shape[0]):
        line = "{}\n".format("\t".join([str(x) for x in list(mdist_df.iloc[i,].values)]))
        f.write(line)
        
para = "Result/NJ/Pig_ELA_%s_neighbor.par"%tissue
out1 = "Result/NJ/Pig_ELA_%s_neighbor.outfile"%tissue
out2 = "Result/NJ/Pig_ELA_%s_neighbor.outtree"%tissue
with open(para,"w") as f:
    f.write("%s\n"%output)
    f.write("%s\n"%"L")
    f.write("%s\n"%"Y")

if not os.path.exists("tmp/%s"%tissue):
    os.mkdir("tmp/%s"%tissue)
os.chdir("tmp/%s"%tissue)

os.system("software/phylip-3.697/exe/neighbor < %s"%para)
os.system("mv outfile %s"%out1)
os.system("mv outtree %s"%out2)

### 画图
pc_res2["SubID"] = pc_res2["Herd"] + "-"+pc_res2["Breed"]

pc_res2[["NewID","SubID"]].to_csv("Result/NJ/groupinfo.txt",index=None,header=None,sep=" ")

os.system("Rscript ggtree.R %s"%tissue)


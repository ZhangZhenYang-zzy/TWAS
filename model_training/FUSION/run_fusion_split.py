import os
import sys

startgene = int(sys.argv[1])
endgene = int(sys.argv[2])
print(startgene,endgene)


#pig
tissues = ['Heart','Morula','Uterus','Artery','Frontal_cortex','Small_intestine','Milk','Testis','Jejunum','Oocyte','Kidney','Synovial_membrane','Placenta','Lymph_node','Muscle','Blastocyst','Ovary','Cartilage','Large_intestine','Adipose','Brain','Embryo','Hypothalamus','Fetal_thymus','Duodenum','Blastomere','Blood','Pituitary','Spleen','Colon','Lung','Macrophage','Ileum','Liver']

for tissue_each in tissues:
	end = os.popen("wc -l %s.expr_tmm_inv.bed"%tissue_each).read()
	final_end = end.split(" ")[0]
	os.system("bash run_fusion.sh %s %d %d & "%(tissue_each,startgene,endgene))

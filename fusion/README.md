# Functional Summary-based Imputation

FUSION is a suite of tools for performing transcriptome-wide and regulome-wide association studies (TWAS and RWAS). FUSION builds predictive models of the genetic component of a functional/molecular phenotype and predicts and tests that component for association with disease using GWAS summary statistics.

Please see [http://gusevlab.org/projects/fusion/](http://gusevlab.org/projects/fusion/) for documentation.

And we use parallel R package to make the FUSION.assoc_test.R can run in parallel
FUSION.profile_wgt.modified.R: 在原来的FUSION中，计算的RSQ是使用的lm function中的adj.r.sq，这会导致和其他方法的结果无法比较，在这个脚本中，我们增加了校正，计算出未校正的RSQ

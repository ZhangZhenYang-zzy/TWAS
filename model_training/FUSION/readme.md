# Training model for FUSION

The code come from https://github.com/gusevlab/fusion_twas

run_fusion_split.py:  using this code, user could run the FUSION.compute_weights.R in parallel.

FUSION.profile_wgt.modified.R: In the original code, the FUSION.profile_wgt.R provides the adjust rsq for correltion  coefficient. We modified the code, and provides the pvalue by 1-(1-adjust.rsq)*(N.tot-2)/(N.tot-1). 

changeTable.R: change the FUSION results for convert it th db

make_sqlite_db.py: change the result to db. python3 make_sqlite_db.py --output ${tissue}.${model}.db --betas ${tissue}.${model}.weight.txt --results ${tissue}.${model}.extra.txt --construction ${tissue}.construction.txt --meta ${tissue}.sample.txt . All the input files come from the changeTable.R.


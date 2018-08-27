# EPEM-GMM

## To run the entropy penalized EM algorithm with the code currently provided, you must have your input data in a specific format
## See SimData_R4_N400_K4_p2_q2.csv for an example
## Column names must be:
#  RepNum    GeneID     V1.timepoint1     V1.timepoint2      ...      V1.timepointT


### The user must specificy the fixed-effects polynomial order (p) and the random-effects polynomial order (r)
### You must also choose whether you want to average over your replicates (if any) by using R=1.
### If you have 4 replicates and would like to run the model without averaging over replicates, then specify R=4, for example.
### However, if you would like to run a simpler model without inclusion of a replicate-speficif random effect, then specify R=1.

### Run program is titled Run_EPEM.R
### Program with functions is titled EM_GMM_PenLik.R

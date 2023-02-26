library(dplyr)
library(data.table)

test_summ = fread("t2d.test.summaries", sep=" ") #Direction * Stat
full_summ = fread("t2d.summaries", sep=" ")

beta_test = test_summ$Direction * test_summ$Stat
beta_full = full_summ$Direction * full_summ$Stat
n = full_summ$n

#ref = fread("reference___") 
#X_j * X_l = n * corr(snp_j, snp_l)
XX = 

R = (t(matrix(beta_test)) %*% XX %*% matrix(beta_full)) / sqrt(n * t(matrix(beta_test)) %*% XX %*% matrix(beta_test))

print(paste0("R squared = ", R^2))
#########################################
##
##    Script for calculating
##    K-L Divergence scores
##
#########################################

rm(list=ls())  # clear global environment

# Load Libraries
library(yada)
library(doParallel)
library(dplyr)

registerDoParallel(detectCores())  # use multiple cores

###### Data Munge ######

# Reformat test data into collapsed(_c) and uncollapsed(_n) versions
var_info_c <- load_var_info("data/US_all_c_var_info.csv")
var_info_n <- load_var_info("data/US_all_var_info.csv")
test_samp_c <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_c,save_file=F)
test_samp_n <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_n,save_file=F)

# Extract US_88
kl_case_c <- test_samp_c %>% filter(SVAD_identifier=="US_0088")
kl_case_n <- test_samp_n %>% filter(SVAD_identifier=="US_0088")
known_age <- kl_case_c$agey

# Set up prior distribution and xcalc
xcalc <- seq(0,23,by=0.01)
th_x <- readRDS("results_univariate/solutionx_US_all_c.rds")
fprior <- calc_x_density(xcalc, th_x)  # prior distribution

# Store US_0088 and fprior for later use
saveRDS(kl_case_c, "data/US_0088.rds")
saveRDS(fprior, "data/priorx.rds")

###### K-L Divergence Calculations ######

# 18-var-c-cdep
y <- kl_case_c
mat <- y %>% select(max_M1,max_M2,max_PM2,man_M1,man_M2,
                    man_PM1,man_C,FGT_EF,HME_EF,RPE_EF,UDE_EF,CC_Oss,
                    ISPR_EF,ILIS_EF,FDL,TPB,HDL,HPB) %>%
  t() %>% as.matrix()
multi_cdep_c <- readRDS("results_eighteen_var_c/cdep_model_US_eighteen_var_c.rds")
post_multi_cdep_c <-
  calc_x_posterior(mat, th_x, multi_cdep_c, xcalc)
kl_div_multi_cdep_c <- calc_kl_div(post_multi_cdep_c, th_x)

# 18-var-c-cindep
multi_cindep_c <- readRDS("results_eighteen_var_c/cindep_model_US_eighteen_var_c.rds")
post_multi_cindep_c <-
  calc_x_posterior(mat, th_x, multi_cindep_c, xcalc)
kl_div_multi_cindep_c <- calc_kl_div(post_multi_cindep_c, th_x)

# 18-var-cdep
y <- kl_case_n
mat <- y %>% select(max_M1,max_M2,max_PM2,man_M1,man_M2,
                    man_PM1,man_C,FGT_EF,HME_EF,RPE_EF,UDE_EF,CC_Oss,
                    ISPR_EF,ILIS_EF,FDL,TPB,HDL,HPB) %>%
  t() %>% as.matrix()
multi_cdep <- readRDS("results_eighteen_var/cdep_model_US_eighteen_var.rds")
post_multi_cdep <-
  calc_x_posterior(mat, th_x, multi_cdep, xcalc)
kl_div_multi_cdep <- calc_kl_div(post_multi_cdep, th_x)

# 18-var-cindep
multi_cindep <- readRDS("results_eighteen_var/cindep_model_US_eighteen_var.rds")
post_multi_cindep <-
  calc_x_posterior(mat, th_x, multi_cindep, xcalc)
kl_div_multi_cindep <- calc_kl_div(post_multi_cindep, th_x)

# lb-cdep
y <- kl_case_c
mat <- y %>% select(FDL,FMSB,FDB,TDL,TPB,TMSB,TDB,FBDL,HDL,HPB,
                    HMSB,HDB,RDL,RPB,RMSB,RDB,UDL,UMSB) %>%
  t() %>% as.matrix()
lb_cdep <- readRDS("results_lb/cdep_model_US_lb.rds")
post_lb_cdep <-
  calc_x_posterior(mat, th_x, lb_cdep, xcalc)
kl_div_lb_cdep <- calc_kl_div(post_lb_cdep, th_x)

# lb-cindep
lb_cindep <- readRDS("results_lb/cindep_model_US_lb.rds")
post_lb_cindep <-
  calc_x_posterior(mat, th_x, lb_cindep, xcalc)
kl_div_lb_cindep <- calc_kl_div(post_lb_cindep, th_x)

# ef_oss_c-cdep
y <- kl_case_c
mat <- y %>% select(FH_EF,FGT_EF,FLT_EF,FDE_EF,TPE_EF,TDE_EF,FBPE_EF,
                    FBDE_EF,HH_Oss,HGT_Oss,HLT_Oss,HPE_EF,HC_Oss,HT_Oss,
                    HLE_Oss,HDE_EF,HME_EF,RPE_EF,RDE_EF,UPE_EF,UDE_EF,
                    CT_EF,CC_Oss,TC_Oss,ISPR_EF,ILIS_EF,PC_Oss,IC_EF) %>%
  t() %>% as.matrix()
ef_oss_cdep <- readRDS("results_ef_oss_c/cdep_model_US_ef_oss_c.rds")
post_ef_oss_cdep <-
  calc_x_posterior(mat, th_x, ef_oss_cdep, xcalc)
kl_div_ef_oss_cdep <- calc_kl_div(post_ef_oss_cdep, th_x)

# ef_oss_c-cindep
ef_oss_cindep <- readRDS("results_ef_oss_c/cindep_model_US_ef_oss_c.rds")
post_ef_oss_cindep <-
  calc_x_posterior(mat, th_x, ef_oss_cindep, xcalc)
kl_div_ef_oss_cindep <- calc_kl_div(post_ef_oss_cindep, th_x)

# dent-cdep
y <- kl_case_c
mat <- y %>% select(max_M1,max_M2,max_M3,max_PM1,max_PM2,max_C,
                    max_I1 ,max_I2,man_M1,man_M2,man_M3,man_PM1,
                    man_PM2,man_C,man_I1,man_I2) %>%
  t() %>% as.matrix()
dent_cdep <- readRDS("results_dent/cdep_model_US_dent.rds")
post_dent_cdep <-
  calc_x_posterior(mat, th_x, dent_cdep, xcalc)
kl_div_dent_cdep <- calc_kl_div(post_dent_cdep, th_x)

# dent-cindep
dent_cindep <- readRDS("results_dent/cindep_model_US_dent.rds")
post_dent_cindep <-
  calc_x_posterior(mat, th_x, dent_cindep, xcalc)
kl_div_dent_cindep <- calc_kl_div(post_dent_cindep, th_x)

# prox_dist_c-cdep
y <- kl_case_c
mat <- y %>% select(FH_EF,FDE_EF,TPE_EF,TDE_EF,FBPE_EF,FBDE_EF,
                    HH_Oss,HPE_EF,HDE_EF,RPE_EF,RDE_EF,UPE_EF,UDE_EF) %>%
  t() %>% as.matrix()
pd_c_cdep <- readRDS("results_prox_dist_c/cdep_model_US_prox_dist_c.rds")
post_pd_c_cdep <-
  calc_x_posterior(mat, th_x, pd_c_cdep, xcalc)
kl_div_pd_c_cdep <- calc_kl_div(post_pd_c_cdep, th_x)

# prox_dist_c-cindep
pd_c_cindep <- readRDS("results_prox_dist_c/cindep_model_US_prox_dist_c.rds")
post_pd_c_cindep <-
  calc_x_posterior(mat, th_x, pd_c_cindep, xcalc)
kl_div_pd_c_cindep <- calc_kl_div(post_pd_c_cindep, th_x)

# prox_dist-cdep
y <- kl_case_c
mat <- y %>% select(FH_EF,FDE_EF,TPE_EF,TDE_EF,FBPE_EF,FBDE_EF,
                    HH_Oss,HPE_EF,HDE_EF,RPE_EF,RDE_EF,UPE_EF,UDE_EF) %>%
  t() %>% as.matrix()
pd_cdep <- readRDS("results_prox_dist/cdep_model_US_prox_dist.rds")
post_pd_cdep <-
  calc_x_posterior(mat, th_x, pd_cdep, xcalc)
kl_div_pd_cdep <- calc_kl_div(post_pd_cdep, th_x)

# prox_dist-cindep
pd_cindep <- readRDS("results_prox_dist/cindep_model_US_prox_dist.rds")
post_pd_cindep <-
  calc_x_posterior(mat, th_x, pd_cindep, xcalc)
kl_div_pd_cindep <- calc_kl_div(post_pd_cindep, th_x)

# rdl-hetero
y <- kl_case_n
mat <- y %>% select(RDL) %>% t() %>% as.matrix()
rdl_hetero <- readRDS("results_univariate/solutiony_US_all_c_cont_k_13_RDL_pow_law_lin_pos_int.rds")
post_rdl_hetero <-
  calc_x_posterior(mat, th_x, rdl_hetero, xcalc)
kl_div_rdl_hetero <- calc_kl_div(post_rdl_hetero, th_x)

# rdl-homo
rdl_homo <- readRDS("results_univariate/solutiony_US_all_c_cont_k_13_RDL_pow_law_const.rds")
post_rdl_homo <-
  calc_x_posterior(mat, th_x, rdl_homo, xcalc)
kl_div_rdl_homo <- calc_kl_div(post_rdl_homo, th_x)

# fdl-hetero
y <- kl_case_n
mat <- y %>% select(FDL) %>% t() %>% as.matrix()
fdl_hetero <- readRDS("results_univariate/solutiony_US_all_c_cont_k_1_FDL_pow_law_lin_pos_int.rds")
post_fdl_hetero <-
  calc_x_posterior(mat, th_x, fdl_hetero, xcalc)
kl_div_fdl_hetero <- calc_kl_div(post_fdl_hetero, th_x)

# fdl-homo
fdl_homo <- readRDS("results_univariate/solutiony_US_all_c_cont_k_1_FDL_pow_law_const.rds")
post_fdl_homo <-
  calc_x_posterior(mat, th_x, fdl_homo, xcalc)
kl_div_fdl_homo <- calc_kl_div(post_fdl_homo, th_x)

# man_PM2-hetero
man_PM2_hetero <- readRDS("results_univariate/solutiony_US_all_c_ord_j_13_man_PM2_lin_ord_lin_pos_int.rds")
mat <- y %>% select(man_PM2) %>% t() %>% as.matrix()
post_man_PM2_hetero <-
  calc_x_posterior(mat, th_x, man_PM2_hetero, xcalc)
kl_div_man_PM2_hetero <- calc_kl_div(post_man_PM2_hetero, th_x)

# man_PM2-homo
man_PM2_homo <- readRDS("results_univariate/solutiony_US_all_c_ord_j_13_man_PM2_lin_ord_const.rds")
post_man_PM2_homo <-
  calc_x_posterior(mat, th_x, man_PM2_homo, xcalc)
kl_div_man_PM2_homo <- calc_kl_div(post_man_PM2_homo, th_x)

# pc_oss-hetero
pc_oss_hetero <- readRDS("results_univariate/solutiony_US_all_ord_j_27_PC_Oss_log_ord_lin_pos_int.rds")
mat <- y %>% select(PC_Oss) %>% t() %>% as.matrix()
post_pc_oss_hetero <-
  calc_x_posterior(mat, th_x, pc_oss_hetero, xcalc)
kl_div_pc_oss_hetero <- calc_kl_div(post_pc_oss_hetero, th_x)

# pc_oss-homo
pc_oss_homo <- readRDS("results_univariate/solutiony_US_all_ord_j_27_PC_Oss_log_ord_const.rds")
post_pc_oss_homo <-
  calc_x_posterior(mat, th_x, pc_oss_homo, xcalc)
kl_div_pc_oss_homo <- calc_kl_div(post_pc_oss_homo, th_x)

# Store all posterior densities and kl_div values in a single list
post_list <- list(post=list(post_multi_cdep_c,post_multi_cindep_c,
                            post_multi_cdep,post_multi_cindep,
                            post_lb_cdep, post_lb_cindep,
                            post_ef_oss_cdep, post_ef_oss_cindep,
                            post_dent_cdep, post_dent_cindep,
                            post_rdl_homo,post_rdl_hetero,
                            post_fdl_homo,post_fdl_hetero,
                            post_man_PM2_homo,post_man_PM2_hetero,
                            post_pc_oss_homo,post_pc_oss_hetero),
                  kl_div=list(kl_div_multi_cdep_c,kl_div_multi_cindep_c,
                              kl_div_multi_cdep,kl_div_multi_cindep,
                              kl_div_lb_cdep, kl_div_lb_cindep,
                              kl_div_ef_oss_cdep, kl_div_ef_oss_cindep,
                              kl_div_dent_cdep, kl_div_dent_cindep,
                              kl_div_rdl_homo,kl_div_rdl_hetero,
                              kl_div_fdl_homo,kl_div_fdl_hetero,
                              kl_div_man_PM2_homo,kl_div_man_PM2_hetero,
                              kl_div_pc_oss_homo,kl_div_pc_oss_hetero))
saveRDS(post_list, "data/US_88_list.rds")  # store in data folder

# Format the kl-div values into a dataframe
kl_df <- data.frame(matrix(nrow=length(post_list$kl_div), ncol=3))
kl_df[,1] <- c(rep("multi",4), rep("lb",2),
               rep("ef_oss",2), rep("dent",2),
               rep("rdl",2), rep("fdl",2),
               rep("man_PM2",2), rep("pc_oss",2))
kl_df[,2] <- c("cdep_c","cindep_c",
               rep(c("cdep","cindep"),4),
               rep(c("homo","hetero"),4))

for(i in 1:length(post_list$kl_div)) {
  kl_df[i, 3] <- post_list$kl_div[[i]]
  
}

names(kl_df) <- c("var","model","kl_div")

write.csv(kl_df, "data/kl_div_table.csv", row.names=F)









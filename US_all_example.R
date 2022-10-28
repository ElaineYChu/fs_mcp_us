#########################################
##
##    Example Pipeline for Univariate 
##  and Multivariate Model Optimization, 
##     Selection, and Prediction
##
#########################################


# Load libraries
library(yada)
library(doParallel)

# Clear the workspace
rm(list=ls())

base_seed <- 719921  # generated using random.org
set.seed(base_seed)  # set seed for reproducibility
registerDoParallel(cl=detectCores())  # prepare system for parallel processing

# Setting analysis to 'US_all' uses the /data/US_all_var_info.csv variable
# information file.
analysis <- 'US_all'
# Output fits and other files to this directory:
# TODO: considered renaming the following variable to data_dir for consistency
#       with previous scripts, though there is some potential for confusion
#       since there is, in fact a /data folder in the repository. For similar
#       reasons, condsider naming analysis to analysis_name.
target_dir <- 'results_US_all_example'

if (!file.exists(target_dir)) {
  dir.create(target_dir)
}

test_data_file_path <- 'data/SVAD_US_test.csv'
train_data_file_path <- 'data/SVAD_US_train.csv'
var_info_path <- paste0('data/', analysis, '_var_info.csv')
var_info <-  yada::load_var_info(var_info_path)
test_data  <- reformat_mcp_data(data_file_path=test_data_file_path,
                                var_info=var_info, 
                                save_file=T)
 
train_data <- reformat_mcp_data(data_file_path=train_data_file_path,
                                var_info=var_info, 
                                save_file=T)
main_problem <- generate_main_problem(data_file=train_data, 
                                      var_info=var_info)
save_problem(data_dir=target_dir,
             analysis_name=analysis,
             problem=main_problem,
             is_folds=F)

# Do the univariate optimizations

mean_specs_ord <- c(rep("pow_law_ord",2),rep("log_ord",2),rep("lin_ord",2))
noise_specs_ord <- rep(c("const","lin_pos_int"),3)
mean_specs_cont <- c(rep("pow_law",2))
noise_specs_cont <- c("const","lin_pos_int")

# Set up for reproducibility
ord_prob_list <- build_univariate_ord_problems(data_dir=target_dir,
                                               analysis_name=analysis,
                                               mean_specs=mean_specs_ord,
                                               noise_specs=noise_specs_ord,
                                               add_folds=F)

seed_vect <- sample.int(1000000, 
                        length(ord_prob_list), 
                        replace=F)  # vector of sampled numbers
ord_success <-
  foreach::foreach(i=1:length(ord_prob_list), .combine=cbind) %dopar% {
    yada::solve_ord_problem(data_dir=target_dir,
                            analysis_name=analysis,
                            ord_prob=ord_prob_list[[i]],
                            anneal_seed=seed_vect[i])
  }
 
cont_prob_list <- build_univariate_cont_problems(data_dir=target_dir,
                                                 analysis_name=analysis,
                                                mean_specs=mean_specs_cont,
                                                 noise_specs=noise_specs_cont,
                                                 add_folds=F)
cont_success <-
  foreach::foreach(i=1:length(cont_prob_list), .combine=cbind) %dopar% {
    yada::solve_cont_problem(data_dir=target_dir,
                             analysis_name=analysis,
                             cont_prob=cont_prob_list[[i]])
  }


ord_models <- build_model_vec(mean_models=unique(mean_specs_ord),
                              noise_models=unique(noise_specs_ord))

cont_models <- build_model_vec(mean_models=unique(mean_specs_cont),
                               noise_models=unique(noise_specs_ord))

# Univariate Model Selection
eval_data <- evaluate_univariate_models(data_dir=target_dir,
                                        analysis_name=analysis,
                                        eval_type="aic",
                                        ord_models=ord_models,
                                        cont_models=cont_models,
                                        cand_tol=0.05,
                                        scale_exp_min=0.01,
                                        beta2_max=5,
                                        save_file=T)

# Calculate prior parameterization on x
th_x <- calc_x_prior(data_dir=target_dir,
                     analysis_name=analysis,
                     prior_type="weibull",
                     offset=0.002,
                     seed=719921,
                     fold=NA,
                     save_file=T)

# Univariate Batch Calculations
# No seed is needed for the batch calc because xcalc is input. However, it must
# be explicitly set to NA in the input call.
univariate_batch_calc(data_dir=target_dir,
                      analysis_name=analysis,
                      test_samp=test_data, 
                      demo_cols=1:3,
                      input_cols=4:ncol(test_data),
                      ci_type="hdi", th_x=th_x, 
                      xcalc=seq(0,23,by=0.01),
                      seed=NA,save_file=T)

best_params <- get_best_univariate_params(data_dir=target_dir,
                                          analysis_name=analysis,
                                          save_file=T)

# Generate Ordinal Credible Interval Table example
max_M1_df <- generate_ord_ci(data_dir=data_dir,
                             analysis_name=analysis_name,
                             var_name="max_M1",  # CHANGE ME for different ordinal variable
                             th_x=th_x,
                             point_est="xmean",  # CHANGE ME if different central measure is desired
                             ci_type="hdi",
                             xcalc=seq(0,23,by=0.01), 
                             save_file=T)

# Generate the Main Conditionally Independent (cindep) multivariate model
cindep_model <- build_cindep_model(data_dir=target_dir,
                                   analysis_name=analysis,
                                   fold=NA,
                                   calc_se=T,
                                   allow_corner=T,
                                   remove_var=T,
                                   save_file=T)

# Use cindep model specifications to initialize Conditionally Dependent (cdep)
# Model Specifications (mod_spec)
cdep_mod_spec <- cindep_model$mod_spec

## Define assumptions for among-variable dependence - CHANGE ME
# J <- cdep_mod_spec$J  # number of ordinal variables
# K <- cdep_mod_spec$K  # number of continuous variables
# cdep_mod_spec$cdep_groups <- 1:(J+K)  # all variables are their own "group"

# THE FOLLOWING CODE IS DEMONSTRATION SPECIFIC, WHERE WE ARE DEFINING FOUR (4)
# GROUPS OF VARIABLES THAT WE ASSUME HAVE GREATER AMONG-GROUP CORRELATION. 
J <- cdep_mod_spec$J  # number of ordinal variables
K <- cdep_mod_spec$K  # number of continuous variables

ind_EF <- which(endsWith(main_problem$var_names,"EF"))  # EF variables
ind_Oss <- which(endsWith(main_problem$var_names,"Oss"))  # Oss variables
ind_Dent <- which(startsWith(main_problem$var_names,"ma"))  # dental variables
ind_LB <- J + (1:K)  # long bones

cdep_groups <- rep(NA,J+K)  # initialize empty vector for all variables
cdep_groups[ind_EF] <- 1  # all EF variables are in group 1
cdep_groups[ind_Oss] <- 2  # all Oss variables are in group 2
cdep_groups[ind_Dent] <- 3  # all dental variables are in group 3
cdep_groups[ind_LB] <- 4  # all long bones are in group 4

cdep_mod_spec$cdep_groups <- cdep_groups  # add cdep_groups to mod_spec

## Define the model as conditionally dependent
cdep_mod_spec$cdep_spec <- "dep"

# Initialize optimization progress and cdep multivariate model final file paths
hjk_progress_file <- build_file_path(data_dir=target_dir,
                                     analysis_name=analysis,
                                     file_type="hjk_progress",
                                     fold=NA)  # fold is NA because using AIC
cdep_save_file <- build_file_path(data_dir=data_dir,
                                  analysis_name=analysis,
                                  file_type="cdep_model",
                                  fold=NA)  # fold is NA because using AIC

# Generate the Main Conditionally Dependent (cdep) multivariate model
# The optimization initiated by the following line takes a vary long time
# (potentially weeks or even months). Consequently, we provide the final
# solution files for each case in this github repository.
cdep_model <- fit_multivariate(x=main_problem$x,
                               Y=main_problem$Y,
                               mod_spec=cdep_mod_spec,
                               cindep_model=cindep_model,
                               prog_file=hjk_progress_file,
                               save_file=cdep_save_file,
                               hjk_control=list(info=T))


# Multivariate Model Selection (coditionally independent versus dependent)
evaluate_multivariate_models(data_dir=target_dir,
                             analysis_name=analysis,
                             eval_type="aic")  # CHANGE ME if using CV

# Multivariate Batch Calculations
multivariate_batch_calc(data_dir=data_dir,
                        analysis_name=analysis_name,
                        test_samp=test_data,
                        demo_cols=1:3,
                        model_type="cindep", ci_type="hdi",  # CHANGE ME if for "cdep" model
                        th_x=th_x, xcalc=seq(0,23,by=0.01),
                        seed=NA, save_file=T)

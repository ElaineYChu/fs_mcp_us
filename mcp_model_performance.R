#########################################
##
##    Script for calculating 
##  MCP Model Performance Metrics
##
#########################################

## Don't forget to set your working directory!

rm(list=ls())  # clear global environment

# Load libraries
library(stringr)
library(yada)
library(doParallel)
library(dplyr)

registerDoParallel(detectCores())  # use multiple cores

######## Calculating Testing Accuracy and RMSE ########

pred_files <- list.files("test-predictions")
# reorder for univariate first, remove sex models
pred_files <- pred_files[c(1,16:44,49:92,2:15)]  
model_perf <- data.frame(matrix(nrow=length(pred_files), ncol=3))  # init df
names(model_perf) <- c("model","% Accuracy","RMSE")  # column names

# Loop through each prediction file and calculate testing accuracy and RMSE
for(i in 1:length(pred_files)) {
  df <- read.csv(paste0("test-predictions/",pred_files[i])) %>% 
    na.omit()
  
  model_perf[i,"model"] <- str_remove(pred_files[i],"_test_predictions.csv")
  model_perf[i,"% Accuracy"] <- round(length(which(df$agey <= df$upper95 & 
                                                     df$agey >= df$lower95))/nrow(df),2)
  model_perf[i,"RMSE"] <- round(Metrics::rmse(df$agey, df$xmean),3)
  
}

# Store object in /data folder
write.csv(model_perf, "data/mcp_model_performance.csv", row.names=F)

######## Test Mean Negative Log Posterior ########

# Import script to identify rows with all defined columns == NA
source("idx_all_na.R")

# For each univariate and multivariate model:
# 1. Calculate the posterior density for each test individual
# 2. Identify the density value for when xcalc==agey
# 3. Take the log of density value
# 4. Add to sum 
# 5. Take the -mean(sum) / N
# 6. Store that value and N to corresponding model
# 7. Use that value to order models for "performance" 

var_info <- load_var_info("data/US_all_var_info.csv")
var_info_c <- load_var_info("data/US_all_c_var_info.csv")
var_info_dent <- load_var_info("data/US_dent_var_info.csv")
var_info_18 <- load_var_info("data/US_eighteen_var_var_info.csv")
var_info_18_c <- load_var_info("data/US_eighteen_var_c_var_info.csv")
var_info_ef_oss <- load_var_info("data/US_ef_oss_var_info.csv")
var_info_lb <- load_var_info("data/US_lb_var_info.csv")
var_info_pd <- load_var_info("data/US_prox_dist_var_info.csv")
var_info_pd_c <- load_var_info("data/US_prox_dist_c_var_info.csv")
xcalc <- seq(0,23,by=0.01)
th_x <- readRDS("results_univariate/solutionx_US_all_c.rds")

models <- model_perf$model  # same model order as model_perf
# define data directories for each model
res_folder <- c(rep("results_univariate/",74),
                rep(c("results_dent/","results_ef_oss/",
                      "results_eighteen_var_c/","results_eighteen_var/",
                      "results_lb/","results_prox_dist_c/",
                      "results_prox_dist/"),2))

# Initialize empty data frame for storing
final <- data.frame(matrix(nrow=length(models), ncol=3))
names(final) <- c("model","N","tmnlp")
final$model <- models

# Loop through each model
for(i in 1:length(models)) {
  if(grepl("cdep|cindep",models[i])) {  # if multivariate
    model <- readRDS(paste0(res_folder[i],models[i],".rds"))
  } else {
    if(grepl("US_all_c",models[i])) {
      var_name <- str_remove(models[i],"_US_all_c")
      
      model <- load_best_univariate_model(data_dir=res_folder[i],
                                          analysis_name="US_all_c",
                                          var_name=var_name)
    } else {
      var_name <- str_remove(models[i],"_US_all")
      
      model <- load_best_univariate_model(data_dir=res_folder[i],
                                          analysis_name="US_all",
                                          var_name=var_name)
    }
  }
  
  # Reformat data per var info file
  var <- models[i]
  if(grepl("dent", var)) {
    test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_dent)
    var <- paste(var_info_dent$Variable[-(1:3)])
  } else {
    if(grepl("eighteen_var", var)) {
      if(grepl("_c",var)) {
        test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_18_c)
        var <- paste(var_info_18_c$Variable[-(1:3)])
      } else {
        test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_18)
        var <- paste(var_info_18$Variable[-(1:3)])
      }
    } else {
      if(grepl("ef_oss", var)) {
        test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_ef_oss)
        var <- paste(var_info_ef_oss$Variable[-(1:3)])
      } else {
        if(grepl("lb", var)) {
          test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_lb)
          var <- paste(var_info_lb$Variable[-(1:3)])
        } else {
          if(grepl("prox_dist", var)) {
            if(grepl("_c", var)) {
              test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_pd_c)
              var <- paste(var_info_pd_c$Variable[-(1:3)])
            } else {
              test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_pd)
              var <- paste(var_info_pd$Variable[-(1:3)])
            }
          } else {
            if(grepl("_c", var)) {
              test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info_c)
              var <- str_remove(var, "_US_all_c")
            } else {
              test <- reformat_mcp_data("data/SVAD_US_test.csv",var_info)
              var <- str_remove(var, "_US_all")
            }
          }
        }
      }
    }
  }
  
  # Initialize for actual TMNLP calculations
  df <- test %>% select(agey, all_of(var))
  
  if(ncol(df)==2) {
    df <- na.omit(df)
  } else {
    df <- df[-idx_all_na(df,2:ncol(df)),]
  }
  val_vec <- c()
  
  print(paste0("Starting analysis for: ", models[i]))
  
  for(j in 1:nrow(df)) {
    age <- round(df[j,"agey"],2)
    post <- calc_x_posterior(t(df[j,-1]), th_x, model, xcalc)
    val <- post$density[round(post$x,2)==age]
    
    val_vec <- c(val_vec, log(val))
  }
  
  if(length(val_vec) != nrow(df)) {
    stop("val_vec is not equal to number of individuals")
  } else {
    N <- final[i,"N"] <- length(val_vec)
  }
  
  final[i,"tmnlp"] <- round(-sum(val_vec, na.rm=T) / N, 4)
  
  
  write.csv(final, "data/tmnlp_progress.csv", row.names=F)
  
}

write.csv(final, "data/tmnlp_final.csv", row.names=F)

# combine model performance metrics after TMNLP is finished
full_model_perf <- full_join(final, model_perf, by="model")

write.csv(full_model_perf, 
          "data/full_mcp_model_performance.csv", row.names=F)






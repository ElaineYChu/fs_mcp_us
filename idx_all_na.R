###########################################
##
##     Function for identifying row
##       index with all NA using
##          specified columns
##
##        By: Elaine Y. Chu
##     Completed on: August 29, 2021
##
###########################################

#' @title All NA Row Indices
#' 
#' @description Generic function that identifies dataframe rows with all NA
#' 
#' @param df Dataframe of interest
#' @param columns Vector of column names or numbers of interest
#' 
#' @returns A vector of row numbers where all selected columns are NA

idx_all_na <- function(df, columns) {
  tmp_df <- df[c(columns)]
  
  idx <- which(rowSums(tmp_df, na.rm=TRUE)==0)
  
  return(idx)
}
















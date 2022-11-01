# Overview
This github repository provides details on the data generation (including visualizations and other validation parameters) process for the
following article (journal assigned ID: forensicsci-1982685):

>Stull et al. 2023 -- Subadult Age Estimation using the Mixed Cumulative Probit
and a Contemporary United States Population. Forensic Sciences. Special Issue:
Estimating Age in Forensic Anthropology. Special Issue editors: Dr. Kanya Godde
and Dr. Rebecca Taylor.

To cut straight to the chase (i.e., running some code), see the section
[Quicksteps](#quicksteps) below, which does assume some familiarity with R and the command
line / terminal.

The optimizations and other calculations that form the basis for the models and results presented in this
publication take a long time (many months). Consequently, our pipeline consists
of a range of scripts and .Rmd files run separately from each other on multiple
computers. This entailed moving and renaming some files. In the following
paragraphs we describe all the steps in detail. Relevant scripts and other
files are included in this github repository.

# Generation of training and test datasets
<!--Do we still need this if we provide it in the sample script? **I think so, yes, if we want to provide a full and accurate account of the pipeline**-->
The very first step in our pipeline is the creation of distinct training and
test datasets. This is accomplished with the script [train-test.R](train-test.R). 
The input to the script is the file "data/SVAD_US.csv" and the two output files are:
data/SVAD_US_train.csv  
data/SVAD_US_test.csv.

# Maximum likelihood optimizations and model selection
<!--TODO: change the script variable set name - Done-->

The next step is to do maximum likelihood fits and model selection for each model based on a combination of mean and noise response. A separate
fit is done for each of the six distinct possible models for each ordinal variable and
each of two distinct possible models for each continuous variable. Furthermore, a
separate set of fits is needed for each of nine distinct variables sets or combinations:
1. dentition ("US_dent"), 
2. epiphyseal fusion and ossification ("US_ef_oss"), 
3. eighteen-variable uncollapsed ("US_eighteen_var"),
4. eighteen-variable collapsed ("US_eighteen_var_c"), 
5. long bones ("US_lb"), 
6. proximal dist uncollapsed ("US_prox_dist"), 
7. proximal dist collapsed ("US_prox_dist_c"), 
8. all variables uncollapsed ("US_all"), 
9. all variables collapsed ("US_all_c". 
Each of these variable sets has a corresponding variable 
information ("var_info") file in the /data folder.


For additional context and details on each of the steps done to fit the models, please see
the vignette at https://rpubs.com/elainechu/mcp_vignette. In addition, we
provide an example script for one set of variables, the all variables uncollapsed. The script
is named "US_all_example.R". To avoid overwriting the files used
for subsequent steps, we output the files to the folder
"results_US_all_example". If the folder does not already exist, it is
first created.

The first parts of the script involve preprocessing to prepare for the next
step of univariate optimizations. The following files are generated:

/data/SVAD_US_train_reformatted.csv
/data/SVAD_US_test_reformatted.csv
/results_US_all_example/problem_US_all.rds

The two reformatted files in the /data folder were written over when we
processed each var_info file, unless `save_file=F` was specified in the header of the corresponding chunk of R code. 
The missing data results, **Figure 1** in the manuscript, use the results of the "all variables uncollapsed" model.


After preprocessing, we performed univariate model optimizations. This step
generates a number of .rds files with names beginning with *solutiony*. These contain the
univariate fits for each variable and model choice (mean and noise combinations).

Next, the AIC criterion is used to choose among univariate model
parameterizations. This yields the following model evaluation file:

/results_US_all_example/eval_data_US_all.rds

The results of AIC criterion model selection also yield the following 
parameterization file: 

/results_US_all_example/US_all_univariate_model_parameters.rds

This file was used to generate **Figures 5 and 6** in the manuscript.

Next, the univariate models chosen by the AIC criterion in the previous step
are combined to yield the conditionally independent model, which is saved in
the following file:

/results_US_all_example/cindep_model_US_all.rds

Next, the optimization is done for the multivariate conditionally dependent
model, for which the conditionally independent model provides the starting
parameter vector for the maximum likelihood maximization. The resulting file(s) 
hold the maximum likelihood solutions for the multivariate optimizations. 
Conditional dependence optimizations take a great deal of time to run. 
Since it is unrealistic that all users will want to run the optimizations 
themselves, we provide the resulting files in the directories that start 
with *results* in this repository. 
The multivariate optimization in the sample code yields the following file:

/results_US_all_example/cdep_model_US_all.rds

The specific optimization results may vary slightly if code is rerun with
different random number seeds since the optimizations utilize randomness at
key points. This applies also to the optimization that yields the prior 
parameterization on x (th_x, see below). In
addition, we have noticed occasional differences in the behavior of the random
number seeds and concluded that the following three factors may influence the
behavior: (1) operating system, (2) R environment (e.g., Rstudio versus base
R), and (3) whether one is running a script or a R markdown file. It's not
obvious that the random number seed functionality is "to blame" per se as there
may be nuanced interactions with such things as floating point number
representations (or even non-nuanced non-interactions that we just didn't
successfully identify in our sleuthing). All-in-all, these discrepancies have
been rare, but we have not always been able to pin them down to specific causes.

The next step is to create model selection results for the multivariate models,
again using AIC as the criterion. This does not inherently yield any file(s) 
at present, but reports the AIC values for both conditionally independent 
(cindep) and conditionally dependent (cdep) multivariate models.

# Intermediate calculations
A Weibull offset mixture is applied to the training dataset x-variable as a 
parameterization for the prior on x. A single parameterization was conducted on 
the full training set. This prior parameterization yields the 
following file:

/results_US_all_example/solutionx_US_all.rds

This file is used for new predictions using the optimized and selected models. 
While we have chosen this approach for the present analysis, we recognize that 
one of the complexities of growth is that not all variables are available across
the same age range, and therefore may only be appropriate in certain 
contexts. For example, diaphyseal long bone lengths are not measured after 
epiphyseal fusion is complete. Therefore, most (if not all) individuals in our 
sample after the age of 15 do not have diaphyseal lengths to use in age 
estimation. A prior on x that ranges from 0 to 21 may therefore not be the most 
appropriate prior distribution for long bone models. Put somewhat differently, 
that an observation is missing/non-missing may itself be informative, and one
way to accommodate this is through the use of different priors. Since a too
informative can over-influence the posterior, however, we chose to use the same
prior for all models. This choice does influence the performance metrics.

As previous mentioned, the same order of operations found in 
[US_all_example.R](US_all_example.R) was used to optimize and select different 
combinations of variables into nine multivariate models, resulting in the following 
folders containing univariate and multivariate model optimizations and selection
results: 

/results_dent  
/results_ef_oss 
/results_eighteen_var 
/results_eighteen_var_c
/results_lb
/results_prox_dist
/results_prox_dist_c
/results_univariate
(/results_US_all_example *example script only, not actually in repo)

# Results generation
The testing sample is used to evaluate performance of each univariate and 
multivariate model. Using the model selection results for univariate models 
(/results_univariate/eval_data_US_all.rds), the *solutiony* file containing 
the selected model is imported. Next, each individual in the testing 
sample is run through the selected *solutiony* model to calculate a posterior 
probability distribution on x. The posterior probability distribution is then 
analyzed for measures of centrality (mean, mode, median), 
and 95% and 99% Credible Intervals (CrI) using highest 
posterior interval (hdi). This process yields one *_test_predictions.csv* file 
per univariate model. For example: 

/results_US_all_example/FDL_US_all_test_predictions.csv  
/results_US_all_example/max_M1_US_all_test_predictions.csv

The same process was followed for both conditionally independent (cindep) 
and conditionally dependent (cdep) multivariate models using the testing 
dataset. This process yields one *_test_predictions.csv* file per 
multivariate model. For example:

/results_US_all_example/cindep_US_all_test_predictions.csv
/results_US_all_example/cdep_US_all_test_predictions.csv

**NOTE: ** For file organization purposes, all *_test_predictions.csv* files 
across multiple /results_* folders were copied into the following separate folder:

[/test-predictions](test-predictions)

These files were used to generate **Figures 7, 8, 13, 14, 15 and 17** in the manuscript.

From each file in /test-predictions, model performance metrics such as 
testing accuracy and root-mean-squared-error (RMSE) are calculated and compiled 
into a single file:

/data/mcp_model_performance.csv

A third model performance metric, [Negative] Test Mean Log Posterior (TMLP), is 
not yet integrated into the same process as above, because it takes a 
substantially longer amount of time to calculate per model. 
TMLP was calculated for each selected univariate model and both cindep and 
cdep multivariate models and stored in the following file:

/data/tmlp_final.csv

Results from both the model performance and TMLP steps were combined into a 
single dataframe and output to the following file:

/data/full_mcp_model_performance.csv

A [standalone script](mcp_model_performance.R) is provided to demonstrate this
process. 

This file was used for generating **Table 7** in the manuscript.

Ordinal stage credible interval (CrI) tables were also generated for each ordinal univariate model. 
They consist of the defined point estimate (mean), 95% and 99% CrIs for 
each ordinal stage. This process yields one *ordinal_ci_* file 
per ordinal univariate model. For example:

/results_US_all_example/ordinal_ci_US_all_max_M1.rds  
/results_US_all_example/ordinal_ci_US_all_HME_EF.rds

Each .rds file was also saved as a .csv file and stored in a new folder for 
file organization purposes. For example:

/ci-tables/ordinal_ci_US_all_max_M1.csv
/ci-tables/ordinal_ci_US_all_HME_EF.csv

These files were used to generate **Figures 9, 10, and 11** in the manuscript.

To compare the point estimates and 95% CrIs between uncollapsed and collapsed epiphyseal 
fusion variables, *ordinal_ci_* files were generated for proximal and distal 
epiphyseal fusion variables. The intermediate fusing stages of the uncollapsed seven-stage 
scoring method (0, 1, 1/2, 2, 2/3, 3, 4) were combined into a single point estimate using 
the mean, and the 95% CrI was expanded to include the minimum and maximum ages 
among the three collapsed stages (1/2, 2, and 2/3). The resulting files are stored as: 

/data/collapsed_ci.rds
/data/uncollapsed_ci.rds

These files were used to generate **Figure 12** in the manuscript.

Kullback-Leibler Divergence is used to describe the amount of "information" that is 
gained by a model's posterior distribution compared to the prior distribution. 
This information is referred to as "bits" of information. To demonstrate the 
differences between model information, a single individual ("US_0088") was 
selected because they had no missing information (i.e. missing data) for the 
univariate/multivariate models of choice. The process yielded the following 
files: 

/data/US_0088.rds
/data/fprior.rds
/data/US_88_list.rds
/data/kl_div_table.csv

A [standalone script](kl_div.R) is provided that demonstrates how the 
caluclations were performed.

These files were used to produce **Figure 16** in the manuscript.

Additional models were generated to generate sex-specific age estimates based on epiphyseal 
fusion univariate models (HPE_EF, TPE_EF) and compare them to the default pooled-sex models. The sex-specific models were 
left out of the model performance step, and are purely for demonstrative and 
discussion purposes. This process yielded the following files: 

/ci-tables/ordinal_ci_US_all_c_HPE_EF_sex-F.csv
/ci-tables/ordinal_ci_US_all_c_HPE_EF_sex-M.csv
/ci-tables/ordinal_ci_US_all_c_TPE_EF_sex-F.csv
/ci-tables/ordinal_ci_US_all_c_TPE_EF_sex-M.csv

These files were used to generate **Figure 18** in the manuscript.


# Quicksteps
**NOTE:** command line / terminal is assumed  

Clone the repository, change directory, and start R:

```
git clone https://github.com/ElaineYChu/fs_mcp_us
cd fs_mcp_us
R
```

Install yada (plus other dependencies as needed, though many users will
already have them).

As necessary:

```
install.packages('devtools')
install.packages('tidyverse')
install.packages('doParallel')
install.packages('magrittr')
install.packages('ggpubr')
install.packages('rmarkdown')
install.packages('prettydoc')
install.packages('ggplot2')
install.packages('Metrics')
```

The version of yada prior to article re-submission (on the dev branch; if
necessary, add force=TRUE):

```
library(devtools)
install_github('MichaelHoltonPrice/yada',
               ref='0cf71eb8b394ba607701b7b076e48530ab837fb4')
```

If desired, comment out the multivariate optimization near the end of the
script, which could take over a month. Then run the example script:

```
source('US_all_example.R')
```

The above script will complete all processes of model 1) optimization, 
2) selection, and 3) prediction. To calculate model performance metrics, 
run the corresponding script:

```
source('mcp_model_performance.R')
```

Next, to calculate K-L Divergence metrics (another form of model performance), 
run the following script:

```
source('kl_div.R')
```

The preceding steps should be enough to generate all files needed to create 
the tables and figures included in the manuscript.

Create the publication results by rendering the following R markdown file:

```
rmarkdown::render('visualizations.Rmd', output_file = 'visualizations.html')
```

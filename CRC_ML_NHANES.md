CRC_ML_NHANES
================
Minyang Xu
2025-03-24

# Import packages

``` r
library(haven)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(purrr)
library(ggplot2)
library(reshape2)
library(corrplot)
```

    ## corrplot 0.95 loaded

``` r
library(broom)
library(caret)
```

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
library(pROC)
```

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
library(tibble)
```

# Set Working Directory & Read Data

``` r
setwd("G:/My Drive/2024 Fall/GU2686/Repository/Dataset")
# "G:/My Drive/2024 Fall/GU2686/Repository/Dataset"
# "/Users/minyangmaxims_xu/Library/CloudStorage/GoogleDrive-mx2269@nyu.edu/My Drive/2024 Fall/GU2686/Repository/Dataset"
demorgraphic <- read_xpt("Demorgraphic.XPT")
totalcholoesterol <- read_xpt("Cholesterol - Total (P_TCHOL).XPT")
chromiunurine <- read_xpt("Chromium Urine.XPT")
completebloodcount <- read_xpt("Complete Blood Count with 5-Part Differential in Whole Blood (P_CBC).XPT")
medicalcondition <- read_xpt("Medical Conditions (P_MCQ).XPT")
smoking <- read_xpt("Smoking - Cigarette Use (P_SMQ).XPT")
bodymeasure <- read_xpt("Body Measure (P_BMX).XPT")
```

# Merge Data by SEQN

``` r
datasets <- list(
  demorgraphic = demorgraphic,
  totalcholoesterol = totalcholoesterol,
  chromiunurine = chromiunurine,
  completebloodcount = completebloodcount,
  medicalcondition = medicalcondition,
  smoking = smoking,
  bodymeasure = bodymeasure
)
# Merge data set and inner join by SEQN
merged_data <- reduce(datasets, ~ inner_join(.x, .y, by = "SEQN"))
head(merged_data)
```

    ## # A tibble: 6 × 153
    ##     SEQN SDDSRVYR RIDSTATR RIAGENDR RIDAGEYR RIDAGEMN RIDRETH1 RIDRETH3 RIDEXMON
    ##    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1 109266       66        2        2       29       NA        5        6        2
    ## 2 109273       66        2        1       36       NA        3        3        2
    ## 3 109274       66        2        1       68       NA        5        7        1
    ## 4 109290       66        2        2       68       NA        4        4        2
    ## 5 109295       66        2        2       54       NA        1        1        1
    ## 6 109300       66        2        2       54       NA        5        6        2
    ## # ℹ 144 more variables: DMDBORN4 <dbl>, DMDYRUSZ <dbl>, DMDEDUC2 <dbl>,
    ## #   DMDMARTZ <dbl>, RIDEXPRG <dbl>, SIALANG <dbl>, SIAPROXY <dbl>,
    ## #   SIAINTRP <dbl>, FIALANG <dbl>, FIAPROXY <dbl>, FIAINTRP <dbl>,
    ## #   MIALANG <dbl>, MIAPROXY <dbl>, MIAINTRP <dbl>, AIALANGA <dbl>,
    ## #   WTINTPRP <dbl>, WTMECPRP <dbl>, SDMVPSU <dbl>, SDMVSTRA <dbl>,
    ## #   INDFMPIR <dbl>, LBXTC <dbl>, LBDTCSI <dbl>, WTSAPRP <dbl>, URXUCM <dbl>,
    ## #   URDUCMLC <dbl>, LBXWBCSI <dbl>, LBXLYPCT <dbl>, LBXMOPCT <dbl>, …

``` r
rm(list = setdiff(ls(), "merged_data"))
```

# Select & Process Features

``` r
# Only keep selected predictors
selected_data <- merged_data %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, LBXTC, LBXWBCSI, LBXLYPCT, 
         LBXEOPCT, LBXRBCSI, LBXHGB, LBXPLTSI, MCQ080, MCQ220, SMQ020, 
         BMXBMI, MCQ230A, MCQ230B, MCQ230C)

# Recode outcome variable by combining 3 variables
selected_data <- selected_data %>%
  mutate(outcome = if_else(MCQ230A == 16 | MCQ230B == 16 | MCQ230C == 16, 1, 0))

selected_data$outcome[is.na(selected_data$outcome)] <- 0
selected_data <- selected_data %>% select(-MCQ230A, -MCQ230B, -MCQ230C)
```

# Data Balancing (Downsampling)

``` r
outcome_0 <- selected_data %>% filter(outcome == 0)
outcome_1 <- selected_data %>% filter(outcome == 1)
set.seed(42)
outcome_0_downsampled <- outcome_0 %>% slice_sample(n = 1200)
balanced_data <- bind_rows(outcome_0_downsampled, outcome_1) %>% sample_frac(1)
table(balanced_data$outcome)
```

    ## 
    ##    0    1 
    ## 1200   21

# Data Preprocessing

``` r
balanced_data <- balanced_data %>% select(-SEQN)
categorical_vars <- c("MCQ080", "MCQ220", "SMQ020")
balanced_data <- balanced_data %>%
  mutate(across(all_of(categorical_vars), as.factor))

balanced_data <- balanced_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

impute_mode <- function(x) {
  mode_val <- names(sort(table(x), decreasing = TRUE))[1]
  x[is.na(x)] <- mode_val
  return(x)
}
balanced_data <- balanced_data %>%
  mutate(across(where(is.factor), impute_mode))
```

# Univariate Analysis

``` r
library(dplyr)
library(openxlsx)
```

    ## Warning: package 'openxlsx' was built under R version 4.4.3

``` r
target_vars <- c(
  "SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", 
  "LBXTC", "LBXWBCSI", "LBXLYPCT", "LBXEOPCT", "LBXRBCSI", "LBXHGB", "LBXPLTSI", 
  "MCQ080", "MCQ220", "SMQ020", "BMXBMI", 
  "MCQ230A", "MCQ230B", "MCQ230C"
)

categorical_vars <- c("RIAGENDR", "RIDRETH1", "MCQ080", "MCQ220", "SMQ020", "MCQ230A", "MCQ230B", "MCQ230C")

summarize_variable <- function(data, var, is_cat = FALSE) {
  x <- data[[var]]
  
  if (is_cat) {
    x <- as.factor(x)
    freq_table <- table(x, useNA = "ifany")
    total <- sum(freq_table)
    categories <- names(freq_table)
    counts <- as.integer(freq_table)
    percentages <- round(100 * counts / total, 1)
    
    tibble(
      Variable = rep(var, length(categories)),
      Type = "Categorical",
      Class = categories,
      Count = counts,
      Percent = paste0(percentages, "%"),
      Mean = NA,
      Median = NA
    )
    
  } else {
    x <- as.numeric(x)
    tibble(
      Variable = var,
      Type = "Continuous",
      Class = NA,
      Count = sum(!is.na(x)),
      Percent = NA,
      Mean = round(mean(x, na.rm = TRUE), 2),
      Median = round(median(x, na.rm = TRUE), 2)
    )
  }
}

results <- bind_rows(lapply(target_vars, function(var) {
  is_cat <- var %in% categorical_vars
  summarize_variable(balanced_data, var, is_cat)
}))
```

# Train-Test Split (80/20)

``` r
# Split by 8 : 2
set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]
```

# Machine Learning Model

## Logistic Regression Model

``` r
logit_model <- glm(outcome ~ ., data = train_data, family = binomial)
summary(logit_model)
```

    ## 
    ## Call:
    ## glm(formula = outcome ~ ., family = binomial, data = train_data)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)  1.818e+00  5.113e+00   0.356   0.7221  
    ## RIDAGEYR     9.558e-03  2.573e-02   0.371   0.7103  
    ## RIAGENDR    -1.547e+00  8.074e-01  -1.916   0.0554 .
    ## RIDRETH1     2.085e-01  2.984e-01   0.699   0.4847  
    ## LBXTC        3.402e-04  7.844e-03   0.043   0.9654  
    ## LBXWBCSI    -1.796e-01  1.946e-01  -0.923   0.3561  
    ## LBXLYPCT    -1.035e-02  3.996e-02  -0.259   0.7957  
    ## LBXEOPCT    -2.581e-02  1.538e-01  -0.168   0.8667  
    ## LBXRBCSI    -1.558e+00  1.167e+00  -1.336   0.1816  
    ## LBXHGB       1.342e-01  3.753e-01   0.357   0.7207  
    ## LBXPLTSI     8.273e-03  5.918e-03   1.398   0.1621  
    ## MCQ0802      8.581e-01  7.874e-01   1.090   0.2758  
    ## MCQ2202     -2.143e+01  1.527e+03  -0.014   0.9888  
    ## MCQ2209     -1.957e+01  4.820e+04   0.000   0.9997  
    ## SMQ0202      6.355e-01  7.170e-01   0.886   0.3754  
    ## BMXBMI       6.633e-02  4.730e-02   1.402   0.1608  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 179.45  on 976  degrees of freedom
    ## Residual deviance:  78.53  on 961  degrees of freedom
    ## AIC: 110.53
    ## 
    ## Number of Fisher Scoring iterations: 21

``` r
logit_prob <- predict(logit_model, test_data, type = "response")
logit_pred <- ifelse(logit_prob > 0.5, 1, 0)
logit_pred <- factor(logit_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix <- confusionMatrix(logit_pred, test_data$outcome, positive = "1")
print(conf_matrix)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 237   1
    ##          1   4   2
    ##                                           
    ##                Accuracy : 0.9795          
    ##                  95% CI : (0.9528, 0.9933)
    ##     No Information Rate : 0.9877          
    ##     P-Value [Acc > NIR] : 0.9173          
    ##                                           
    ##                   Kappa : 0.4352          
    ##                                           
    ##  Mcnemar's Test P-Value : 0.3711          
    ##                                           
    ##             Sensitivity : 0.666667        
    ##             Specificity : 0.983402        
    ##          Pos Pred Value : 0.333333        
    ##          Neg Pred Value : 0.995798        
    ##              Prevalence : 0.012295        
    ##          Detection Rate : 0.008197        
    ##    Detection Prevalence : 0.024590        
    ##       Balanced Accuracy : 0.825035        
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_logit <- conf_matrix$overall["Accuracy"]
precision_logit <- conf_matrix$byClass["Precision"]
recall_logit <- conf_matrix$byClass["Sensitivity"]
specificity_logit <- conf_matrix$byClass["Specificity"]
f1_score_logit <- conf_matrix$byClass["F1"]

cat("\n📌 **Logistic Regression Performance Metrics**:\n")
```

    ## 
    ## 📌 **Logistic Regression Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_logit, 4), "\n")
```

    ## Accuracy: 0.9795

``` r
cat("Precision:", round(precision_logit, 4), "\n")
```

    ## Precision: 0.3333

``` r
cat("Recall:", round(recall_logit, 4), "\n")
```

    ## Recall: 0.6667

``` r
cat("Specificity:", round(specificity_logit,4), "\n")
```

    ## Specificity: 0.9834

``` r
cat("F1 Score:", round(f1_score_logit, 4), "\n")
```

    ## F1 Score: 0.4444

``` r
roc_curve <- roc(test_data$outcome, logit_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_value <- auc(roc_curve)
cat("AUC:", round(auc_value, 4), "\n")
```

    ## AUC: 0.9737

``` r
precision_logit <- conf_matrix$byClass["Precision"]
cat("Precision:", round(precision_logit, 4), "\n")
```

    ## Precision: 0.3333

``` r
ggplot() +
  geom_line(aes(x = roc_curve$specificities, y = roc_curve$sensitivities), color = "blue") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - Logistic Regression", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
feature_importance <- broom::tidy(logit_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(abs_coef = abs(estimate)) %>%
  arrange(desc(abs_coef))

# Display top features
print(feature_importance)
```

    ## # A tibble: 15 × 6
    ##    term       estimate   std.error statistic p.value  abs_coef
    ##    <chr>         <dbl>       <dbl>     <dbl>   <dbl>     <dbl>
    ##  1 MCQ2202  -21.4       1527.      -0.0140    0.989  21.4     
    ##  2 MCQ2209  -19.6      48196.      -0.000406  1.00   19.6     
    ##  3 LBXRBCSI  -1.56         1.17    -1.34      0.182   1.56    
    ##  4 RIAGENDR  -1.55         0.807   -1.92      0.0554  1.55    
    ##  5 MCQ0802    0.858        0.787    1.09      0.276   0.858   
    ##  6 SMQ0202    0.636        0.717    0.886     0.375   0.636   
    ##  7 RIDRETH1   0.209        0.298    0.699     0.485   0.209   
    ##  8 LBXWBCSI  -0.180        0.195   -0.923     0.356   0.180   
    ##  9 LBXHGB     0.134        0.375    0.357     0.721   0.134   
    ## 10 BMXBMI     0.0663       0.0473   1.40      0.161   0.0663  
    ## 11 LBXEOPCT  -0.0258       0.154   -0.168     0.867   0.0258  
    ## 12 LBXLYPCT  -0.0103       0.0400  -0.259     0.796   0.0103  
    ## 13 RIDAGEYR   0.00956      0.0257   0.371     0.710   0.00956 
    ## 14 LBXPLTSI   0.00827      0.00592  1.40      0.162   0.00827 
    ## 15 LBXTC      0.000340     0.00784  0.0434    0.965   0.000340

``` r
# Plot Feature Importance
feature_importance_logit <- ggplot(feature_importance, aes(x = reorder(term, abs_coef), y = abs_coef)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance - Logistic Regression",
       x = "Feature",
       y = "Absolute Coefficient Value") +
  theme_minimal()
```

## Random Forest

``` r
library(caret)
library(randomForest) 
```

    ## Warning: package 'randomForest' was built under R version 4.4.3

    ## randomForest 4.7-1.2

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(pROC)   
library(ggplot2)
balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]

set.seed(123)
# Building random forest model
rf_model <- randomForest(outcome ~ ., data = train_data, ntree = 500, mtry = sqrt(ncol(train_data) - 1), importance = TRUE)

rf_pred <- predict(rf_model, test_data)

rf_pred <- factor(rf_pred, levels = levels(test_data$outcome))

conf_matrix_rf <- confusionMatrix(rf_pred, test_data$outcome, positive = "1")
print(conf_matrix_rf)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
rf_prob <- predict(rf_model, test_data, type = "prob")[,2]
roc_curve_rf <- roc(test_data$outcome, rf_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_rf <- auc(roc_curve_rf)
cat("Random Forest AUC:", round(auc_rf, 4), "\n")
```

    ## Random Forest AUC: 0.9464

``` r
accuracy_rf <- conf_matrix_rf$overall["Accuracy"]
precision_rf <- conf_matrix_rf$byClass["Precision"]
recall_rf <- conf_matrix_rf$byClass["Sensitivity"]
specificity_rf <- conf_matrix_rf$byClass["Specificity"]
f1_score_rf <- conf_matrix_rf$byClass["F1"]
rf_prob <- predict(rf_model, test_data, type = "prob")[,2]
roc_curve_rf <- roc(test_data$outcome, rf_prob)
```

    ## Setting levels: control = 0, case = 1
    ## Setting direction: controls < cases

``` r
auc_rf <- auc(roc_curve_rf)

cat("\n **Random Forest Performance Metrics**:\n")
```

    ## 
    ##  **Random Forest Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_rf, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_rf, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall:", round(recall_rf, 4), "\n")
```

    ## Recall: 0

``` r
cat("Specificity:", round(specificity_rf, 4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_rf, 4), "\n")
```

    ## F1 Score: NA

``` r
cat("AUC:", round(auc_rf, 4), "\n")
```

    ## AUC: 0.9464

``` r
# Plotting ROC curve
ggplot() +
  geom_line(aes(x = roc_curve_rf$specificities, y = roc_curve_rf$sensitivities), color = "red") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - Random Forest", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Display feature importance
importance_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column(var = "Feature") %>%
  arrange(desc(MeanDecreaseGini))

importance_df_rf <- ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance - Random Forest", x = "Feature", y = "Importance (MeanDecreaseGini)") +
  theme_minimal()
```

## XGBoost

``` r
library(caret)
library(xgboost)
```

    ## Warning: package 'xgboost' was built under R version 4.4.3

    ## 
    ## Attaching package: 'xgboost'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     slice

``` r
library(pROC) 
library(ggplot2) 

balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]

# Transform to matrix
# 去掉 outcome 列
train_x <- model.matrix(outcome ~ . -1, data = train_data)  # -1 去掉 intercept
train_y <- as.numeric(as.character(train_data$outcome))

test_x <- model.matrix(outcome ~ . -1, data = test_data)
test_y <- as.numeric(as.character(test_data$outcome))


dtrain <- xgb.DMatrix(data = train_x, label = train_y)
dtest <- xgb.DMatrix(data = test_x, label = test_y)

xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc", 
  max_depth = 2,                 
  eta = 0.5,                      
  nthread = 2,                
  subsample = 0.5,               
  colsample_bytree = 0.5
)


# Train xgb model
set.seed(123)
xgb_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 60,   
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 10, 
  verbose = 1
)
```

    ## [1]  train-auc:0.961887  test-auc:0.956250 
    ## Multiple eval metrics are present. Will use test_auc for early stopping.
    ## Will train until test_auc hasn't improved in 10 rounds.
    ## 
    ## [2]  train-auc:0.968995  test-auc:0.948438 
    ## [3]  train-auc:0.970374  test-auc:0.948958 
    ## [4]  train-auc:0.977880  test-auc:0.948438 
    ## [5]  train-auc:0.977880  test-auc:0.948438 
    ## [6]  train-auc:0.985600  test-auc:0.946354 
    ## [7]  train-auc:0.987837  test-auc:0.943229 
    ## [8]  train-auc:0.991728  test-auc:0.950000 
    ## [9]  train-auc:0.992494  test-auc:0.937500 
    ## [10] train-auc:0.991544  test-auc:0.933333 
    ## [11] train-auc:0.994363  test-auc:0.926042 
    ## Stopping. Best iteration:
    ## [1]  train-auc:0.961887  test-auc:0.956250

``` r
xgb_prob <- predict(xgb_model, dtest)

xgb_pred <- ifelse(xgb_prob > 0.5, 1, 0)
xgb_pred <- factor(xgb_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix_xgb <- confusionMatrix(xgb_pred, test_data$outcome, positive = "1")
print(conf_matrix_xgb)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
# Showing performance matrics
accuracy_xgb <- conf_matrix_xgb$overall["Accuracy"]
precision_xgb <- conf_matrix_xgb$byClass["Precision"]
recall_xgb <- conf_matrix_xgb$byClass["Sensitivity"]
specificity_xgb <- conf_matrix_xgb$byClass["Specificity"]
f1_score_xgb <- conf_matrix_xgb$byClass["F1"]

cat("\n **XGBoost Performance Metrics**:\n")
```

    ## 
    ##  **XGBoost Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_xgb, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_xgb, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall (Sensitivity):", round(recall_xgb, 4), "\n")
```

    ## Recall (Sensitivity): 0

``` r
cat("Specificity:", round(specificity_xgb, 4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_xgb, 4), "\n")
```

    ## F1 Score: NA

``` r
roc_curve_xgb <- roc(test_data$outcome, xgb_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_xgb <- auc(roc_curve_xgb)
cat("AUC:", round(auc_xgb, 4), "\n")
```

    ## AUC: 0.9562

``` r
# Plotting ROC curve
ggplot() +
  geom_line(aes(x = roc_curve_xgb$specificities, y = roc_curve_xgb$sensitivities), color = "green") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - XGBoost", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
importance_matrix <- xgb.importance(feature_names = colnames(train_x), model = xgb_model)

print(importance_matrix)
```

    ##      Feature        Gain       Cover  Frequency
    ##       <char>       <num>       <num>      <num>
    ##  1:  MCQ2202 0.392149201 0.634036077 0.30434783
    ##  2: LBXEOPCT 0.169941227 0.069515381 0.17391304
    ##  3: LBXLYPCT 0.121531860 0.042992963 0.13043478
    ##  4: RIDAGEYR 0.089616985 0.113741931 0.08695652
    ##  5: LBXWBCSI 0.062232830 0.010520038 0.04347826
    ##  6: LBXRBCSI 0.056759745 0.009060854 0.04347826
    ##  7:   BMXBMI 0.043010260 0.010194150 0.04347826
    ##  8:   LBXHGB 0.032438236 0.016732400 0.08695652
    ##  9: LBXPLTSI 0.023558610 0.015142425 0.04347826
    ## 10:    LBXTC 0.008761047 0.078063781 0.04347826

``` r
importance_matrix_xgb <- ggplot(importance_matrix, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Feature Importance - XGBoost",
       x = "Feature",
       y = "Importance (Gain)") +
  theme_minimal()
```

## LASSO

``` r
library(glmnet) 
```

    ## Warning: package 'glmnet' was built under R version 4.4.3

    ## Loading required package: Matrix

    ## Loaded glmnet 4.1-8

``` r
library(caret)
library(pROC)   
library(ggplot2)

balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]

train_x <- as.matrix(train_data[, -which(names(train_data) == "outcome")])
train_y <- as.numeric(as.character(train_data$outcome))

test_x <- as.matrix(test_data[, -which(names(test_data) == "outcome")])
test_y <- as.numeric(as.character(test_data$outcome))

# Implementing cross validation bv n = 10
set.seed(123)
cv_lasso <- cv.glmnet(train_x, train_y, alpha = 1, family = "binomial", nfolds = 10)
best_lambda <- cv_lasso$lambda.min
cat("Best Lambda (Regularization Strength):", round(best_lambda, 6), "\n")
```

    ## Best Lambda (Regularization Strength): 0.005077

``` r
# Training LASSO model
lasso_model <- glmnet(train_x, train_y, alpha = 1, family = "binomial", lambda = best_lambda)

lasso_prob <- predict(lasso_model, newx = test_x, type = "response")

lasso_pred <- ifelse(lasso_prob > 0.5, 1, 0)
lasso_pred <- factor(lasso_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix_lasso <- confusionMatrix(lasso_pred, test_data$outcome, positive = "1")
print(conf_matrix_lasso)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_lasso <- conf_matrix_lasso$overall["Accuracy"]
precision_lasso <- conf_matrix_lasso$byClass["Precision"]
recall_lasso <- conf_matrix_lasso$byClass["Sensitivity"]
specificity_lasso <- conf_matrix_lasso$byClass["Specificity"]
f1_score_lasso <- conf_matrix_lasso$byClass["F1"]

cat("\n **LASSO Logistic Regression Performance Metrics**:\n")
```

    ## 
    ##  **LASSO Logistic Regression Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_lasso, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_lasso, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall (Sensitivity):", round(recall_lasso, 4), "\n")
```

    ## Recall (Sensitivity): 0

``` r
cat("Specificity:", round(specificity_lasso, 4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_lasso, 4), "\n")
```

    ## F1 Score: NA

``` r
roc_curve_lasso <- roc(test_data$outcome, lasso_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Warning in roc.default(test_data$outcome, lasso_prob): Deprecated use a matrix
    ## as predictor. Unexpected results may be produced, please pass a numeric vector.

    ## Setting direction: controls < cases

``` r
auc_lasso <- auc(roc_curve_lasso)
cat("AUC:", round(auc_lasso, 4), "\n")
```

    ## AUC: 0.951

``` r
# Plotting auc
ggplot() +
  geom_line(aes(x = roc_curve_lasso$specificities, y = roc_curve_lasso$sensitivities), color = "purple") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - LASSO Logistic Regression", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
library(glmnet)
library(ggplot2)
library(dplyr)

# Getting LASSO oefficients
lasso_coeff <- coef(lasso_model, s = best_lambda) 
lasso_coeff_df <- as.data.frame(as.matrix(lasso_coeff)) 
lasso_coeff_df$Feature <- rownames(lasso_coeff_df)
colnames(lasso_coeff_df) <- c("Coefficient", "Feature")

lasso_coeff_df <- lasso_coeff_df %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Importance = abs(Coefficient)) %>%
  arrange(desc(Importance))

# Feature importance
lasso_coeff_df <- ggplot(lasso_coeff_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "Feature Importance - LASSO Regression",
       x = "Feature",
       y = "Absolute Coefficient") +
  theme_minimal()
```

## Naive Bayes

``` r
library(e1071)
library(caret)
library(pROC) 
library(ggplot2)

balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]

# Training Naive Bayes
nb_model <- naiveBayes(outcome ~ ., data = train_data)

nb_prob <- predict(nb_model, test_data, type = "raw")[, 2]  # 获取属于类别 1 的概率
nb_pred <- ifelse(nb_prob > 0.5, 1, 0)
nb_pred <- factor(nb_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix_nb <- confusionMatrix(nb_pred, test_data$outcome, positive = "1")
print(conf_matrix_nb)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 234   4
    ##          1   6   0
    ##                                           
    ##                Accuracy : 0.959           
    ##                  95% CI : (0.9259, 0.9802)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.9974          
    ##                                           
    ##                   Kappa : -0.0201         
    ##                                           
    ##  Mcnemar's Test P-Value : 0.7518          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 0.97500         
    ##          Pos Pred Value : 0.00000         
    ##          Neg Pred Value : 0.98319         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.02459         
    ##       Balanced Accuracy : 0.48750         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_nb <- conf_matrix_nb$overall["Accuracy"]
precision_nb <- conf_matrix_nb$byClass["Precision"]
recall_nb <- conf_matrix_nb$byClass["Sensitivity"]
specificity_nb <- conf_matrix_nb$byClass["Specificity"]
f1_score_nb <- conf_matrix_nb$byClass["F1"]

cat("\n **Naïve Bayes Performance Metrics**:\n")
```

    ## 
    ##  **Naïve Bayes Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_nb, 4), "\n")
```

    ## Accuracy: 0.959

``` r
cat("Precision:", round(precision_nb, 4), "\n")
```

    ## Precision: 0

``` r
cat("Recall (Sensitivity):", round(recall_nb, 4), "\n")
```

    ## Recall (Sensitivity): 0

``` r
cat("Specificity:", round(specificity_nb, 4), "\n")
```

    ## Specificity: 0.975

``` r
cat("F1 Score:", round(f1_score_nb, 4), "\n")
```

    ## F1 Score: NaN

``` r
roc_curve_nb <- roc(test_data$outcome, nb_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_nb <- auc(roc_curve_nb)
cat("AUC:", round(auc_nb, 4), "\n")
```

    ## AUC: 0.9427

``` r
ggplot() +
  geom_line(aes(x = roc_curve_nb$specificities, y = roc_curve_nb$sensitivities), color = "orange") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - Naïve Bayes", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
feature_importance_nb <- data.frame(Feature = character(), Importance = numeric())

for (feature in names(nb_model$tables)) {
  prob_table <- nb_model$tables[[feature]]
  
  if (is.matrix(prob_table)) { 
    importance_score <- sum(abs(log(prob_table[,2] / prob_table[,1])), na.rm = TRUE)
    feature_importance_nb <- rbind(feature_importance_nb, data.frame(Feature = feature, Importance = importance_score))
  }
}

# Ranking feature importance
feature_importance_nb <- feature_importance_nb %>% arrange(desc(Importance))
print(feature_importance_nb)
```

    ##     Feature Importance
    ## 1    MCQ220        Inf
    ## 2  LBXRBCSI  4.3833537
    ## 3    LBXHGB  4.3438714
    ## 4     LBXTC  3.2429773
    ## 5    BMXBMI  2.7880734
    ## 6  LBXPLTSI  2.6757361
    ## 7  LBXWBCSI  2.6444149
    ## 8  LBXLYPCT  2.5785742
    ## 9  RIDAGEYR  2.3965374
    ## 10 RIAGENDR  2.1645353
    ## 11 RIDRETH1  2.0079563
    ## 12 LBXEOPCT  0.7545293
    ## 13   MCQ080  0.7414042
    ## 14   SMQ020  0.6553661

``` r
naivebayes_feature_importance <- ggplot(feature_importance_nb, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  labs(title = "Feature Importance - Naïve Bayes",
       x = "Feature",
       y = "Log Probability Contribution") +
  theme_minimal()
```

## Perfromance Metrics

``` r
performance_metrics <- data.frame(
  Model = c("Logistic Regression", "Random Forest", "Naïve Bayes", "XGBoost"),
  Accuracy = c(accuracy_logit, accuracy_rf, accuracy_nb, accuracy_xgb),
  Precision = c(precision_logit, precision_rf, precision_nb, precision_xgb),
  Specificity = c(specificity_logit, specificity_rf, specificity_nb, specificity_xgb),
  AUC = c(auc_value, auc_rf, auc_nb, auc_xgb)
)
print(performance_metrics)
```

    ##                 Model  Accuracy Precision Specificity       AUC
    ## 1 Logistic Regression 0.9795082 0.3333333   0.9834025 0.9737206
    ## 2       Random Forest 0.9836066        NA   1.0000000 0.9463542
    ## 3         Naïve Bayes 0.9590164 0.0000000   0.9750000 0.9427083
    ## 4             XGBoost 0.9836066        NA   1.0000000 0.9562500

# Visualize ROC Curve

``` r
library(ggplot2)
library(pROC)

roc_df_before <- rbind(
  data.frame(FPR = 1 - roc_curve$specificities, TPR = roc_curve$sensitivities, Model = "Logistic Regression"),
  data.frame(FPR = 1 - roc_curve_rf$specificities, TPR = roc_curve_rf$sensitivities, Model = "Random Forest"),
  data.frame(FPR = 1 - roc_curve_nb$specificities, TPR = roc_curve_nb$sensitivities, Model = "Naïve Bayes"),
  data.frame(FPR = 1 - roc_curve_xgb$specificities, TPR = roc_curve_xgb$sensitivities, Model = "XGBoost")
)

roc_plot_before <- ggplot(roc_df_before, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves - Before Feature Selection",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_minimal() +
  theme(legend.title = element_blank())
print(roc_plot_before)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

# Feature Selection

``` r
print(feature_importance_logit)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
print(importance_df_rf)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
print(importance_matrix_xgb)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

``` r
print(lasso_coeff_df)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

``` r
print(naivebayes_feature_importance)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-16-5.png)<!-- -->

``` r
library(gridExtra)
```

    ## Warning: package 'gridExtra' was built under R version 4.4.3

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:randomForest':
    ## 
    ##     combine

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(grid)
library(ggplot2)
combined_plot <- grid.arrange(
  feature_importance_logit, 
  importance_df_rf, 
  importance_matrix_xgb, 
  lasso_coeff_df, 
  ncol = 2
)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-16-6.png)<!-- -->

``` r
print(combined_plot)
```

    ## TableGrob (2 x 2) "arrange": 4 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (1-1,2-2) arrange gtable[layout]
    ## 3 3 (2-2,1-1) arrange gtable[layout]
    ## 4 4 (2-2,2-2) arrange gtable[layout]

# Remodelling after feature selection

## preprocessing

``` r
balanced_data <- balanced_data %>% select(-c(MCQ220, LBXHGB, LBXEOPCT, LBXLYPCT, SMQ020))
set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]
```

## Logistic Regression Model

``` r
logit_model <- glm(outcome ~ ., data = train_data, family = binomial)
summary(logit_model)
```

    ## 
    ## Call:
    ## glm(formula = outcome ~ ., family = binomial, data = train_data)
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -3.845700   4.188552  -0.918 0.358543    
    ## RIDAGEYR     0.077019   0.023229   3.316 0.000914 ***
    ## RIAGENDR    -1.052628   0.588404  -1.789 0.073622 .  
    ## RIDRETH1     0.146256   0.262072   0.558 0.576792    
    ## LBXTC        0.008688   0.007023   1.237 0.216010    
    ## LBXWBCSI    -0.092079   0.146160  -0.630 0.528703    
    ## LBXRBCSI    -1.384519   0.579467  -2.389 0.016881 *  
    ## LBXPLTSI     0.002992   0.003602   0.831 0.406098    
    ## MCQ0802     -0.301455   0.630648  -0.478 0.632644    
    ## BMXBMI       0.034475   0.039829   0.866 0.386717    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 171.45  on 976  degrees of freedom
    ## Residual deviance: 134.80  on 967  degrees of freedom
    ## AIC: 154.8
    ## 
    ## Number of Fisher Scoring iterations: 9

``` r
logit_prob <- predict(logit_model, test_data, type = "response")
logit_pred <- ifelse(logit_prob > 0.5, 1, 0)
logit_pred <- factor(logit_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix <- confusionMatrix(logit_pred, test_data$outcome, positive = "1")
print(conf_matrix)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_logit <- conf_matrix$overall["Accuracy"]
precision_logit <- conf_matrix$byClass["Precision"]
recall_logit <- conf_matrix$byClass["Sensitivity"]
specificity_logit <- conf_matrix$byClass["Specificity"]
f1_score_logit <- conf_matrix$byClass["F1"]

cat("\n **Logistic Regression Performance Metrics**:\n")
```

    ## 
    ##  **Logistic Regression Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_logit, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_logit, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall:", round(recall_logit, 4), "\n")
```

    ## Recall: 0

``` r
cat("Specificity:", round(specificity_logit,4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_logit, 4), "\n")
```

    ## F1 Score: NA

``` r
roc_curve <- roc(test_data$outcome, logit_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_value <- auc(roc_curve)
cat("AUC:", round(auc_value, 4), "\n")
```

    ## AUC: 0.7521

``` r
precision_logit <- conf_matrix$byClass["Precision"]
cat("Precision:", round(precision_logit, 4), "\n")
```

    ## Precision: NA

``` r
ggplot() +
  geom_line(aes(x = roc_curve$specificities, y = roc_curve$sensitivities), color = "blue") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - Logistic Regression", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
feature_importance <- broom::tidy(logit_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(abs_coef = abs(estimate)) %>%
  arrange(desc(abs_coef))

print(feature_importance)
```

    ## # A tibble: 9 × 6
    ##   term     estimate std.error statistic  p.value abs_coef
    ##   <chr>       <dbl>     <dbl>     <dbl>    <dbl>    <dbl>
    ## 1 LBXRBCSI -1.38      0.579      -2.39  0.0169    1.38   
    ## 2 RIAGENDR -1.05      0.588      -1.79  0.0736    1.05   
    ## 3 MCQ0802  -0.301     0.631      -0.478 0.633     0.301  
    ## 4 RIDRETH1  0.146     0.262       0.558 0.577     0.146  
    ## 5 LBXWBCSI -0.0921    0.146      -0.630 0.529     0.0921 
    ## 6 RIDAGEYR  0.0770    0.0232      3.32  0.000914  0.0770 
    ## 7 BMXBMI    0.0345    0.0398      0.866 0.387     0.0345 
    ## 8 LBXTC     0.00869   0.00702     1.24  0.216     0.00869
    ## 9 LBXPLTSI  0.00299   0.00360     0.831 0.406     0.00299

``` r
feature_importance_logit <- ggplot(feature_importance, aes(x = reorder(term, abs_coef), y = abs_coef)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance - Logistic Regression",
       x = "Feature",
       y = "Absolute Coefficient Value") +
  theme_minimal()
```

## Random Forest

``` r
library(caret)    
library(randomForest)
library(pROC)     
library(ggplot2)

balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]
set.seed(123)
rf_model <- randomForest(outcome ~ ., data = train_data, ntree = 500, mtry = sqrt(ncol(train_data) - 1), importance = TRUE)

rf_pred <- predict(rf_model, test_data)

rf_pred <- factor(rf_pred, levels = levels(test_data$outcome))

conf_matrix_rf <- confusionMatrix(rf_pred, test_data$outcome, positive = "1")
print(conf_matrix_rf)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
rf_prob <- predict(rf_model, test_data, type = "prob")[,2]
roc_curve_rf <- roc(test_data$outcome, rf_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_rf <- auc(roc_curve_rf)
cat("Random Forest AUC:", round(auc_rf, 4), "\n")
```

    ## Random Forest AUC: 0.5495

``` r
accuracy_rf <- conf_matrix_rf$overall["Accuracy"]
precision_rf <- conf_matrix_rf$byClass["Precision"]
recall_rf <- conf_matrix_rf$byClass["Sensitivity"]
specificity_rf <- conf_matrix_rf$byClass["Specificity"]
f1_score_rf <- conf_matrix_rf$byClass["F1"]
rf_prob <- predict(rf_model, test_data, type = "prob")[,2]
roc_curve_rf <- roc(test_data$outcome, rf_prob)
```

    ## Setting levels: control = 0, case = 1
    ## Setting direction: controls < cases

``` r
auc_rf <- auc(roc_curve_rf)

cat("\n **Random Forest Performance Metrics**:\n")
```

    ## 
    ##  **Random Forest Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_rf, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_rf, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall:", round(recall_rf, 4), "\n")
```

    ## Recall: 0

``` r
cat("Specificity:", round(specificity_rf, 4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_rf, 4), "\n")
```

    ## F1 Score: NA

``` r
cat("AUC:", round(auc_rf, 4), "\n")
```

    ## AUC: 0.5495

``` r
ggplot() +
  geom_line(aes(x = roc_curve_rf$specificities, y = roc_curve_rf$sensitivities), color = "red") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - Random Forest", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
importance_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column(var = "Feature") %>%
  arrange(desc(MeanDecreaseGini))

importance_df_rf <- ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance - Random Forest", x = "Feature", y = "Importance (MeanDecreaseGini)") +
  theme_minimal()
```

## XGBoost

``` r
library(caret)  
library(xgboost)
library(pROC)
library(ggplot2)

balanced_data$outcome <- as.factor(balanced_data$outcome)


set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]


train_x <- model.matrix(outcome ~ . -1, data = train_data)
train_y <- as.numeric(as.character(train_data$outcome))

test_x <- model.matrix(outcome ~ . -1, data = test_data)
test_y <- as.numeric(as.character(test_data$outcome))


dtrain <- xgb.DMatrix(data = train_x, label = train_y)
dtest <- xgb.DMatrix(data = test_x, label = test_y)

xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",         
  max_depth = 2,                
  eta = 0.5,                    
  nthread = 2,                      
  subsample = 0.5,              
  colsample_bytree = 0.5        
)


set.seed(123)
xgb_model <- xgb.train(
  params = xgb_params,
  data = dtrain,
  nrounds = 60,   
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 10,  
  verbose = 1
)
```

    ## [1]  train-auc:0.715564  test-auc:0.601042 
    ## Multiple eval metrics are present. Will use test_auc for early stopping.
    ## Will train until test_auc hasn't improved in 10 rounds.
    ## 
    ## [2]  train-auc:0.718903  test-auc:0.600000 
    ## [3]  train-auc:0.835386  test-auc:0.637500 
    ## [4]  train-auc:0.884589  test-auc:0.750521 
    ## [5]  train-auc:0.921232  test-auc:0.814583 
    ## [6]  train-auc:0.918597  test-auc:0.820833 
    ## [7]  train-auc:0.920404  test-auc:0.790104 
    ## [8]  train-auc:0.910355  test-auc:0.687500 
    ## [9]  train-auc:0.917004  test-auc:0.660937 
    ## [10] train-auc:0.956127  test-auc:0.531250 
    ## [11] train-auc:0.960355  test-auc:0.498958 
    ## [12] train-auc:0.955637  test-auc:0.616667 
    ## [13] train-auc:0.955055  test-auc:0.520833 
    ## [14] train-auc:0.970221  test-auc:0.566667 
    ## [15] train-auc:0.975245  test-auc:0.642708 
    ## [16] train-auc:0.981801  test-auc:0.604167 
    ## Stopping. Best iteration:
    ## [6]  train-auc:0.918597  test-auc:0.820833

``` r
xgb_prob <- predict(xgb_model, dtest)

xgb_pred <- ifelse(xgb_prob > 0.5, 1, 0)
xgb_pred <- factor(xgb_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix_xgb <- confusionMatrix(xgb_pred, test_data$outcome, positive = "1")
print(conf_matrix_xgb)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_xgb <- conf_matrix_xgb$overall["Accuracy"]
precision_xgb <- conf_matrix_xgb$byClass["Precision"]
recall_xgb <- conf_matrix_xgb$byClass["Sensitivity"]
specificity_xgb <- conf_matrix_xgb$byClass["Specificity"]
f1_score_xgb <- conf_matrix_xgb$byClass["F1"]

cat("\n **XGBoost Performance Metrics**:\n")
```

    ## 
    ##  **XGBoost Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_xgb, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_xgb, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall (Sensitivity):", round(recall_xgb, 4), "\n")
```

    ## Recall (Sensitivity): 0

``` r
cat("Specificity:", round(specificity_xgb, 4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_xgb, 4), "\n")
```

    ## F1 Score: NA

``` r
roc_curve_xgb <- roc(test_data$outcome, xgb_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_xgb <- auc(roc_curve_xgb)
cat("AUC:", round(auc_xgb, 4), "\n")
```

    ## AUC: 0.8208

``` r
ggplot() +
  geom_line(aes(x = roc_curve_xgb$specificities, y = roc_curve_xgb$sensitivities), color = "green") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - XGBoost", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
importance_matrix <- xgb.importance(feature_names = colnames(train_x), model = xgb_model)

print(importance_matrix)
```

    ##     Feature        Gain       Cover  Frequency
    ##      <char>       <num>       <num>      <num>
    ## 1: RIDAGEYR 0.263208029 0.433400179 0.21621622
    ## 2: LBXPLTSI 0.166086051 0.068254185 0.16216216
    ## 3: LBXWBCSI 0.161231723 0.053616484 0.13513514
    ## 4:    LBXTC 0.143393435 0.228512336 0.18918919
    ## 5: LBXRBCSI 0.131790015 0.156786008 0.13513514
    ## 6:   BMXBMI 0.100890581 0.034371098 0.08108108
    ## 7: RIDRETH1 0.016696173 0.015603458 0.02702703
    ## 8:  MCQ0802 0.012952800 0.005101574 0.02702703
    ## 9: RIAGENDR 0.003751192 0.004354678 0.02702703

``` r
importance_matrix_xgb <- ggplot(importance_matrix, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(title = "Feature Importance - XGBoost",
       x = "Feature",
       y = "Importance (Gain)") +
  theme_minimal()
```

## LASSO

``` r
library(glmnet) 
library(caret)  
library(pROC)   
library(ggplot2) 

balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]

train_x <- as.matrix(train_data[, -which(names(train_data) == "outcome")])
train_y <- as.numeric(as.character(train_data$outcome))

test_x <- as.matrix(test_data[, -which(names(test_data) == "outcome")])
test_y <- as.numeric(as.character(test_data$outcome))

set.seed(123)
cv_lasso <- cv.glmnet(train_x, train_y, alpha = 1, family = "binomial", nfolds = 10)
best_lambda <- cv_lasso$lambda.min
cat("Best Lambda (Regularization Strength):", round(best_lambda, 6), "\n")
```

    ## Best Lambda (Regularization Strength): 0.004826

``` r
lasso_model <- glmnet(train_x, train_y, alpha = 1, family = "binomial", lambda = best_lambda)

lasso_prob <- predict(lasso_model, newx = test_x, type = "response")

lasso_pred <- ifelse(lasso_prob > 0.5, 1, 0)
lasso_pred <- factor(lasso_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix_lasso <- confusionMatrix(lasso_pred, test_data$outcome, positive = "1")
print(conf_matrix_lasso)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 240   4
    ##          1   0   0
    ##                                           
    ##                Accuracy : 0.9836          
    ##                  95% CI : (0.9586, 0.9955)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.6288          
    ##                                           
    ##                   Kappa : 0               
    ##                                           
    ##  Mcnemar's Test P-Value : 0.1336          
    ##                                           
    ##             Sensitivity : 0.00000         
    ##             Specificity : 1.00000         
    ##          Pos Pred Value :     NaN         
    ##          Neg Pred Value : 0.98361         
    ##              Prevalence : 0.01639         
    ##          Detection Rate : 0.00000         
    ##    Detection Prevalence : 0.00000         
    ##       Balanced Accuracy : 0.50000         
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_lasso <- conf_matrix_lasso$overall["Accuracy"]
precision_lasso <- conf_matrix_lasso$byClass["Precision"]
recall_lasso <- conf_matrix_lasso$byClass["Sensitivity"]
specificity_lasso <- conf_matrix_lasso$byClass["Specificity"]
f1_score_lasso <- conf_matrix_lasso$byClass["F1"]

cat("\n **LASSO Logistic Regression Performance Metrics**:\n")
```

    ## 
    ##  **LASSO Logistic Regression Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_lasso, 4), "\n")
```

    ## Accuracy: 0.9836

``` r
cat("Precision:", round(precision_lasso, 4), "\n")
```

    ## Precision: NA

``` r
cat("Recall (Sensitivity):", round(recall_lasso, 4), "\n")
```

    ## Recall (Sensitivity): 0

``` r
cat("Specificity:", round(specificity_lasso, 4), "\n")
```

    ## Specificity: 1

``` r
cat("F1 Score:", round(f1_score_lasso, 4), "\n")
```

    ## F1 Score: NA

``` r
roc_curve_lasso <- roc(test_data$outcome, lasso_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Warning in roc.default(test_data$outcome, lasso_prob): Deprecated use a matrix
    ## as predictor. Unexpected results may be produced, please pass a numeric vector.

    ## Setting direction: controls < cases

``` r
auc_lasso <- auc(roc_curve_lasso)
cat("AUC:", round(auc_lasso, 4), "\n")
```

    ## AUC: 0.7875

``` r
ggplot() +
  geom_line(aes(x = roc_curve_lasso$specificities, y = roc_curve_lasso$sensitivities), color = "purple") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - LASSO Logistic Regression", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
library(glmnet)
library(ggplot2)
library(dplyr)

lasso_coeff <- coef(lasso_model, s = best_lambda) 
lasso_coeff_df <- as.data.frame(as.matrix(lasso_coeff))  
lasso_coeff_df$Feature <- rownames(lasso_coeff_df)
colnames(lasso_coeff_df) <- c("Coefficient", "Feature")

lasso_coeff_df <- lasso_coeff_df %>%
  filter(Feature != "(Intercept)") %>%
  mutate(Importance = abs(Coefficient)) %>%
  arrange(desc(Importance))

lasso_coeff_df <- ggplot(lasso_coeff_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(title = "Feature Importance - LASSO Regression",
       x = "Feature",
       y = "Absolute Coefficient") +
  theme_minimal()
```

## Naive Bayes

``` r
library(e1071)
library(caret)
library(pROC)  
library(ggplot2) 

balanced_data$outcome <- as.factor(balanced_data$outcome)

set.seed(123)
train_index <- createDataPartition(balanced_data$outcome, p = 0.8, list = FALSE)
train_data <- balanced_data[train_index, ]
test_data <- balanced_data[-train_index, ]

nb_model <- naiveBayes(outcome ~ ., data = train_data)

nb_prob <- predict(nb_model, test_data, type = "raw")[, 2]
nb_pred <- ifelse(nb_prob > 0.5, 1, 0)
nb_pred <- factor(nb_pred, levels = c(0, 1))
test_data$outcome <- factor(test_data$outcome, levels = c(0, 1))

conf_matrix_nb <- confusionMatrix(nb_pred, test_data$outcome, positive = "1")
print(conf_matrix_nb)
```

    ## Confusion Matrix and Statistics
    ## 
    ##           Reference
    ## Prediction   0   1
    ##          0 239   4
    ##          1   1   0
    ##                                           
    ##                Accuracy : 0.9795          
    ##                  95% CI : (0.9528, 0.9933)
    ##     No Information Rate : 0.9836          
    ##     P-Value [Acc > NIR] : 0.7864          
    ##                                           
    ##                   Kappa : -0.0066         
    ##                                           
    ##  Mcnemar's Test P-Value : 0.3711          
    ##                                           
    ##             Sensitivity : 0.000000        
    ##             Specificity : 0.995833        
    ##          Pos Pred Value : 0.000000        
    ##          Neg Pred Value : 0.983539        
    ##              Prevalence : 0.016393        
    ##          Detection Rate : 0.000000        
    ##    Detection Prevalence : 0.004098        
    ##       Balanced Accuracy : 0.497917        
    ##                                           
    ##        'Positive' Class : 1               
    ## 

``` r
accuracy_nb <- conf_matrix_nb$overall["Accuracy"]
precision_nb <- conf_matrix_nb$byClass["Precision"]
recall_nb <- conf_matrix_nb$byClass["Sensitivity"]
specificity_nb <- conf_matrix_nb$byClass["Specificity"]
f1_score_nb <- conf_matrix_nb$byClass["F1"]

cat("\n**Naïve Bayes Performance Metrics**:\n")
```

    ## 
    ## **Naïve Bayes Performance Metrics**:

``` r
cat("Accuracy:", round(accuracy_nb, 4), "\n")
```

    ## Accuracy: 0.9795

``` r
cat("Precision:", round(precision_nb, 4), "\n")
```

    ## Precision: 0

``` r
cat("Recall (Sensitivity):", round(recall_nb, 4), "\n")
```

    ## Recall (Sensitivity): 0

``` r
cat("Specificity:", round(specificity_nb, 4), "\n")
```

    ## Specificity: 0.9958

``` r
cat("F1 Score:", round(f1_score_nb, 4), "\n")
```

    ## F1 Score: NaN

``` r
roc_curve_nb <- roc(test_data$outcome, nb_prob)
```

    ## Setting levels: control = 0, case = 1

    ## Setting direction: controls < cases

``` r
auc_nb <- auc(roc_curve_nb)
cat("AUC:", round(auc_nb, 4), "\n")
```

    ## AUC: 0.6927

``` r
ggplot() +
  geom_line(aes(x = roc_curve_nb$specificities, y = roc_curve_nb$sensitivities), color = "orange") +
  geom_abline(linetype = "dashed") +
  labs(title = "ROC Curve - Naïve Bayes", x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal()
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
feature_importance_nb <- data.frame(Feature = character(), Importance = numeric())

for (feature in names(nb_model$tables)) {
  prob_table <- nb_model$tables[[feature]]
  
  if (is.matrix(prob_table)) {  
    importance_score <- sum(abs(log(prob_table[,2] / prob_table[,1])), na.rm = TRUE)
    
    feature_importance_nb <- rbind(feature_importance_nb, data.frame(Feature = feature, Importance = importance_score))
  }
}

feature_importance_nb <- feature_importance_nb %>% arrange(desc(Importance))

print(feature_importance_nb)
```

    ##    Feature Importance
    ## 1 LBXRBCSI  4.3833537
    ## 2    LBXTC  3.2429773
    ## 3   BMXBMI  2.7880734
    ## 4 LBXPLTSI  2.6757361
    ## 5 LBXWBCSI  2.6444149
    ## 6 RIDAGEYR  2.3965374
    ## 7 RIAGENDR  2.1645353
    ## 8 RIDRETH1  2.0079563
    ## 9   MCQ080  0.7414042

``` r
naivebayes_feature_importance <- ggplot(feature_importance_nb, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  labs(title = "Feature Importance - Naïve Bayes",
       x = "Feature",
       y = "Log Probability Contribution") +
  theme_minimal()
```

## Perfromance Metrics

``` r
performance_metrics <- data.frame(
  Model = c("Logistic Regression", "Random Forest", "Naïve Bayes", "XGBoost", "Lasso"),
  Accuracy = c(accuracy_logit, accuracy_rf, accuracy_nb, accuracy_xgb, accuracy_lasso),
  Precision = c(precision_logit, precision_rf, precision_nb, precision_xgb, precision_lasso),
  Specificity = c(specificity_logit, specificity_rf, specificity_nb, specificity_xgb, specificity_lasso),
  AUC = c(auc_value, auc_rf, auc_nb, auc_xgb, auc_lasso)
)

print(performance_metrics)
```

    ##                 Model  Accuracy Precision Specificity       AUC
    ## 1 Logistic Regression 0.9836066        NA   1.0000000 0.7520833
    ## 2       Random Forest 0.9836066        NA   1.0000000 0.5494792
    ## 3         Naïve Bayes 0.9795082         0   0.9958333 0.6927083
    ## 4             XGBoost 0.9836066        NA   1.0000000 0.8208333
    ## 5               Lasso 0.9836066        NA   1.0000000 0.7875000

# Visualize AUC

``` r
roc_df_after <- rbind(
  data.frame(FPR = 1 - roc_curve$specificities, TPR = roc_curve$sensitivities, Model = "Logistic Regression"),
  data.frame(FPR = 1 - roc_curve_rf$specificities, TPR = roc_curve_rf$sensitivities, Model = "Random Forest"),
  data.frame(FPR = 1 - roc_curve_nb$specificities, TPR = roc_curve_nb$sensitivities, Model = "Naïve Bayes"),
  data.frame(FPR = 1 - roc_curve_xgb$specificities, TPR = roc_curve_xgb$sensitivities, Model = "XGBoost"),
  data.frame(FPR = 1 - roc_curve_lasso$specificities, TPR = roc_curve_lasso$sensitivities, Model = "LASSO")
)
roc_plot_after <- ggplot(roc_df_after, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.1) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves - After Feature Selection",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  theme_minimal() +
  theme(legend.title = element_blank())
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
print(roc_plot_after)
```

![](CRC_ML_NHANES_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

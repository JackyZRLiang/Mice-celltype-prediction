### best model
library(tidyverse)
library(tidymodels)
library(randomForest)
library(skimr)
library(glmnet)
library(leaps)
library(xgboost)


tr <- read_csv("mouse_tr.csv") %>%
  mutate(celltype = factor(celltype))
ts <- read_csv("mouse_ts_mask.csv")

grid <- 10^seq(10,-2,length=100)
x <- model.matrix(celltype~., tr[,-1])[,-1]
y <- tr$celltype


#### lasso
lasso.mod <- glmnet(x, y, alpha=1, lambda=grid, family = "multinomial")

co_1 <- as.matrix(lasso.mod$beta$cardiomyocyte)
co_1 <- cbind(var=rownames(co_1), co_1)

co_2 <- as.matrix(lasso.mod$beta$endothelium)
co_2 <- cbind(var=rownames(co_2), co_2)

co_3 <- as.matrix(lasso.mod$beta$erythroid)
co_3 <- cbind(var=rownames(co_3), co_3)

co_4 <- as.matrix(lasso.mod$beta$forebrain)
co_4 <- cbind(var=rownames(co_4), co_4)

co_5 <- as.matrix(lasso.mod$beta$gut)
co_5 <- cbind(var=rownames(co_5), co_5)

co_6 <- as.matrix(lasso.mod$beta$mesenchyme)
co_6 <- cbind(var=rownames(co_6), co_6)

co_7 <- as.matrix(lasso.mod$beta$mid_hindbrain)
co_7 <- cbind(var=rownames(co_7), co_7)

co_8 <- as.matrix(lasso.mod$beta$neural_crest)
co_8 <- cbind(var=rownames(co_8), co_8)

co_9 <- as.matrix(lasso.mod$beta$somitic_mesoderm)
co_9 <- cbind(var=rownames(co_9), co_9)

co_10 <- as.matrix(lasso.mod$beta$spinalCord)
co_10 <- cbind(var=rownames(co_10), co_10)


trans_1 <- co_1 %>% as_tibble() %>%
  gather(nval, co_1, -var) %>%
  mutate(coeffs=as.numeric(co_1)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_2 <- co_2 %>% as_tibble() %>%
  gather(nval, co_2, -var) %>%
  mutate(coeffs=as.numeric(co_2)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_3 <- co_3 %>% as_tibble() %>%
  gather(nval, co_3, -var) %>%
  mutate(coeffs=as.numeric(co_3)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_4 <- co_4 %>% as_tibble() %>%
  gather(nval, co_4, -var) %>%
  mutate(coeffs=as.numeric(co_4)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_5 <- co_5 %>% as_tibble() %>%
  gather(nval, co_5, -var) %>%
  mutate(coeffs=as.numeric(co_5)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_6 <- co_6 %>% as_tibble() %>%
  gather(nval, co_6, -var) %>%
  mutate(coeffs=as.numeric(co_6)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_7 <- co_7 %>% as_tibble() %>%
  gather(nval, co_7, -var) %>%
  mutate(coeffs=as.numeric(co_7)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_8 <- co_8 %>% as_tibble() %>%
  gather(nval, co_8, -var) %>%
  mutate(coeffs=as.numeric(co_8)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_9 <- co_9 %>% as_tibble() %>%
  gather(nval, co_9, -var) %>%
  mutate(coeffs=as.numeric(co_9)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))

trans_10 <- co_10 %>% as_tibble() %>%
  gather(nval, co_10, -var) %>%
  mutate(coeffs=as.numeric(co_10)) %>%
  mutate(lambda=rep(lasso.mod$lambda, rep(1003, 100)))



lasso_v_1 <- trans_1 %>%
  filter(co_1 != 0) %>%
  select(var) %>%
  unique()

lasso_v_2 <- trans_2 %>%
  filter(co_2 != 0) %>%
  select(var) %>%
  unique()

lasso_v_3 <- trans_3 %>%
  filter(co_3 != 0) %>%
  select(var) %>%
  unique()

lasso_v_4 <- trans_4 %>%
  filter(co_4 != 0) %>%
  select(var) %>%
  unique()

lasso_v_5 <- trans_5 %>%
  filter(co_5 != 0) %>%
  select(var) %>%
  unique()

lasso_v_6 <- trans_6 %>%
  filter(co_6 != 0) %>%
  select(var) %>%
  unique()

lasso_v_7 <- trans_7 %>%
  filter(co_7 != 0) %>%
  select(var) %>%
  unique()

lasso_v_8 <- trans_8 %>%
  filter(co_8 != 0) %>%
  select(var) %>%
  unique()

lasso_v_9 <- trans_9 %>%
  filter(co_9 != 0) %>%
  select(var) %>%
  unique()

lasso_v_10 <- trans_10 %>%
  filter(co_10 != 0) %>%
  select(var) %>%
  unique()


lasso_var_all <- rbind(lasso_v_1,lasso_v_2,lasso_v_3,lasso_v_4,lasso_v_5,
                       lasso_v_6,lasso_v_7,lasso_v_8,lasso_v_9,lasso_v_10) %>%
  unique()



tr_lasso <- tr %>%
  select(location, lasso_var_all[[1]], celltype)
lasso_variables <- colnames(tr_lasso)
lasso_variables <- lasso_variables[-c(1,220,221,222,223,224,225,226,227,228,1005)]
library(MASS)

outcome <- "celltype"
fmla <- as.formula(paste( paste('as.factor(', outcome, ')', sep = ''),
                          paste(c(lasso_variables), collapse="+"), sep="~"))

lda_model <- lda(fmla, data = tr)
pred <- predict(lda_model, newdata=ts, type="class")$class

mouse_pred <- ts
mouse_pred$celltype <- pred

write_csv(mouse_pred[,c(1, 1005)], file="mouse_pre_lasso_lda2.csv")
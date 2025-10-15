#load('msd.RData')
#covarfinal1<- read.csv(file='covarfinal2.csv')
#endpoint_crc<- read.csv(file='endpoint_crc1.csv')

#msd <- persd
#persd <- merge(covarfinal1,msd,by.x = "f.eid",by.y = "id")
#persd <- merge(endpoint_crc,persd,by.x = "id",by.y = "f.eid")

# set.seed(123)
# crc_1 <- persd[persd$outcome == 1, ]
# crc_0 <- persd[persd$outcome == 0, ]
# 
# sample_crc_1 <- crc_1[sample(nrow(crc_1), 50), ]
# sample_crc_0 <- crc_0[sample(nrow(crc_0), 50), ]
#sample_df <- rbind(sample_crc_1, sample_crc_0)
################CRC LASSO#########
library(survival)
library(glmnet)
library(data.table)
library(dplyr)
sample_df = data.table::fread("/data/cox_example.csv")
persd = sample_df

persd$outcome <- factor(persd$outcome,levels = c(0,1))  
persd$f.31.0.0 <- factor(persd$f.31.0.0,levels = c(0,1))  
persd$censor <- factor(persd$censor, labels = c(1, 2))
persd$ethnic2 <- factor(persd$ethnic2,
                        levels = c(0, 1, 999),
                        labels = c("0", "1", "999"))
persd$fasting_status <- factor(persd$fasting_status,
                               levels = c(0, 1, 999),
                               labels = c("0", "1", "999"))
persd$smoking_staus <- factor(persd$smoking_staus,
                              levels = c(0, 1, 2, 3, 4, 999),
                              labels = c("0", "1", "2", "3", "4", "999"))
persd$drinking_freq <- factor(persd$drinking_freq,
                              levels = c(1, 2, 3, 4, 5, 999),
                              labels = c("1", "2", "3", "5", "6", "999"))
persd$family_history <- factor(persd$family_history,
                               levels = c(0, 1, 999),
                               labels = c("0", "1", "999"))
persd$bmi.class <- factor(persd$bmi.class,levels = c(1,2,3))  


set.seed(2022)
res = fread('/data/perlogcrcM1.csv') #the association between metabolites and CRC in the model 1
res$fdr <- as.numeric(res$fdr)
metid <- res$Metabolites.ID[res$fdr<0.05]

lasso <- persd %>% select(c("id",metid))
x <- as.matrix(lasso[,1:length(lasso)])

time <- as.double(persd$time)
outcome <- as.double(persd$outcome)

surv <- Surv(time,outcome)

alpha1_fit <- glmnet(x,surv,alpha=1,family="cox")  
plot(alpha1_fit,xvar="lambda",label=TRUE)

cvfit <- cv.glmnet(x,surv,type.measure = "deviance",alpha=1,family="cox",nfolds = 10)
plot(cvfit)
print(cvfit)

l.coef1 <- coef(cvfit$glmnet.fit,s=cvfit$lambda.min,exact = F)
l.coef1 <- coef(cvfit$glmnet.fit,s=cvfit$lambda.1se ,exact = F)

get_coe <- function(the_fit,the_lamb){
  Coefficients <- coef(the_fit,s = the_lamb)
  Active.Index <- which(Coefficients!=0)
  Active.Coefficients <- Coefficients[Active.Index]
  re <- data.frame(rownames(Coefficients)[Active.Index],Active.Coefficients)
  re <- data.table('var_names' = rownames(Coefficients)[Active.Index],
                   'coef' = Active.Coefficients)
  
  re$expcoef <- exp(re$coef)
  return(re[order(expcoef)])
}

z <- get_coe(cvfit,cvfit$lambda.min)
print(z)
#round(z$coef,6)


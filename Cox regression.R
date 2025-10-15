#Cox regression
#use CRC as example
################## CRC Model1 ###########
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
library(dplyr)
library(survival)

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

res1 <- data.frame()
for (i in 42:290) {
  tryCatch({
    result <- coxph(Surv(time, outcome) ~ persd[[i]] + f.21003.0.0 + f.31.0.0 + censor + ethnic2,
                    data = persd, id = id)
    summary_result <- summary(result)
    
    out <- data.frame(
      Metabolites.ID = names(persd)[i],
      HR = round(exp(coef(result))[1], 4),
      beta = round(coef(result)[1], 4),
      CI1 = round(exp(confint(result))[1, 1], 4),
      CI2 = round(exp(confint(result))[1, 2], 4),
      p = round(summary(result)$coefficients[1, 6], 6),  
      stringsAsFactors = FALSE
    )
    res1 <- rbind(res1, out)
  }, error = function(e) {
    cat("error in", names(persd)[i], ":", conditionMessage(e), "\n")
  })
}
p <- res1$p
res1$fdr <- round(p.adjust(p,method = 'fdr',n=length(p)),5)
res1$CI <- paste("(",res1$CI1,",",res1$CI2,")",sep = "")
head(res1)
write.csv(res1,"crc_model1.csv",row.names = F)
################## CRC Model2 ###########
res2 <- data.frame()
for (i in 42:290) {
  tryCatch({
    result <- coxph(Surv(time, outcome) ~ persd[[i]] + f.21003.0.0 + f.31.0.0 + censor + ethnic2 + Townsend + fasting_status + smoking_staus + drinking_freq + family_history,
                    data = persd, id = id)
    summary_result <- summary(result)
    out <- data.frame(
      Metabolites.ID = names(persd)[i],
      HR = round(exp(coef(result))[1], 4),
      beta = round(coef(result)[1], 4),
      CI1 = round(exp(confint(result))[1, 1], 4),
      CI2 = round(exp(confint(result))[1, 2], 4),
      p = round(summary(result)$coefficients[1, 6], 6),  
      stringsAsFactors = FALSE
    )
    res2 <- rbind(res2, out)
  }, error = function(e) {
    cat("Warning", names(persd)[i], ":", conditionMessage(e), "\n")
  })
}
p <- res2$p
res2$fdr <- round(p.adjust(p,method = 'fdr',n=length(p)),5)
res2$CI <- paste("(",res2$CI1,",",res2$CI2,")",sep = "")
head(res2)
################## CRC Model3 ###########
res3 <- data.frame()
for (i in 42:290) {
  tryCatch({
    result <- coxph(Surv(time, outcome) ~ persd[[i]] + f.21003.0.0 + f.31.0.0 + censor + ethnic2 + Townsend + fasting_status + smoking_staus + drinking_freq + family_history + bmi.class,
                    data = persd, id = id)
    
    summary_result <- summary(result)
    out <- data.frame(
      Metabolites.ID = names(persd)[i],
      HR = round(exp(coef(result))[1], 4),
      beta = round(coef(result)[1], 4),
      CI1 = round(exp(confint(result))[1, 1], 4),
      CI2 = round(exp(confint(result))[1, 2], 4),
      p = round(summary(result)$coefficients[1, 6], 6),  
      stringsAsFactors = FALSE
    )
    res3 <- rbind(res3, out)
  }, error = function(e) {
    cat("Warning", names(persd)[i], ":", conditionMessage(e), "\n")
  })
}
p <- res3$p
res3$fdr <- round(p.adjust(p,method = 'fdr',n=length(p)),5)
res3$CI <- paste("(",res3$CI1,",",res3$CI2,")",sep = "")
head(res3)
################## CRC Model4 ###########
res4 <- data.frame()
for (i in 42:290) {
  tryCatch({
    result <- coxph(Surv(time, outcome) ~ persd[[i]] + f.21003.0.0 + f.31.0.0 + censor + ethnic2 + Townsend + fasting_status + smoking_staus + drinking_freq + family_history + bmi.class + MDscore,
                    data = persd, id = id)
    summary_result <- summary(result)
    out <- data.frame(
      Metabolites.ID = names(persd)[i],
      HR = round(exp(coef(result))[1], 4),
      beta = round(coef(result)[1], 4),
      CI1 = round(exp(confint(result))[1, 1], 4),
      CI2 = round(exp(confint(result))[1, 2], 4),
      p = round(summary(result)$coefficients[1, 6], 6),  
      stringsAsFactors = FALSE
    )
    res4 <- rbind(res4, out)
  }, error = function(e) {
    cat("Warning", names(persd)[i], ":", conditionMessage(e), "\n")
  })
}

p <- res4$p
res4$fdr <- round(p.adjust(p,method = 'fdr',n=length(p)),5)
res4$CI <- paste("(",res4$CI1,",",res4$CI2,")",sep = "")
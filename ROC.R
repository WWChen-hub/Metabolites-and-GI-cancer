#install.packages("survminer")
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
sample_df = data.table::fread("/data/cox_example.csv")
persd = sample_df

library(broom)
library(survival)
library(survminer)
library(timeROC)
library(pROC)  # DeLong
library(pec)

z <- fread('/Volumes/ww/1 project/GI 代谢组/CEBP/2025review/lasso-CRC.csv',header = T)

score<- persd %>% dplyr::select(z$var_names)
sss <- score %>% as.data.frame()

for (i in 1:length(score)) {
  sss[[i]] <- as.numeric(z$coef[i])*as.numeric(score[[i]])
  sss$score <- rowSums(sss[1:i])
}
sss <- cbind(persd$id,sss)
names(sss)[match("persd$id", names(sss))] <- "id"

m <- merge(persd, sss, by = "id")
m$tertile <- cut(m$score, breaks = quantile(m$score, probs = seq(0, 1, by = 1/3)), 
                 labels = c("T1", "T2", "T3"), include.lowest = TRUE)
m$time <- m$time / 365

#food_roc <- fread('food_roc1.csv',header = T)
#food_roc$id = rownames(food_roc)
#food_roc = food_roc[food_roc$id%in% sample_df$id,]
#write.csv(food_roc,"food_roc.csv",row.names = F)

food_roc = data.table::fread("/data/food_roc.csv")
m <- merge(m, food_roc, by = "id")

cox_model <- coxph(Surv(time, outcome) ~ f.21003.0.0 + f.31.0.0 + smoking_staus + drinking_freq + family_history + f.21001.0.0 + MET2 + food_vegetable + food_fruit + food_meat, data = m)

#####
cox_results <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)

cox_table <- cox_results %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  dplyr::rename(
    Variable = term,
    HR = estimate,
    `95%CI_lower` = conf.low,
    `95%CI_upper` = conf.high,
    P_value = p.value
  ) %>%
  mutate(
    `95% CI` = paste0(round(`95%CI_lower`, 6), "-", round(`95%CI_upper`, 6)),
    HR = round(HR, 6),
    P_value = round(P_value, 6)
  ) %>%
  dplyr::select(Variable, HR, `95% CI`, P_value)


m$time <- m$time*365

m_complete <- na.omit(m[, c("time", "outcome", "score", "f.21003.0.0", "f.31.0.0", 
                            "f.21001.0.0")])

time_point <- 5 * 365

# 1. score
roc_score <- timeROC(
  T = m_complete$time,
  delta = m_complete$outcome,
  marker = m_complete$score,
  cause = 1,
  weighting = "marginal",
  times = time_point,
  ROC = TRUE,
  iid = FALSE
)

# 2. clinical
cox_clinical <- coxph(Surv(time, outcome) ~ f.21003.0.0 + f.31.0.0 + 
                        f.21001.0.0, 
                      data = m_complete)
clinical_score <- predict(cox_clinical, type = "lp")

roc_clinical <- timeROC(
  T = m_complete$time,
  delta = m_complete$outcome,
  marker = clinical_score,
  cause = 1,
  weighting = "marginal",
  times = time_point,
  ROC = TRUE,
  iid = FALSE
)

# 3. score + clinical
cox_combined <- coxph(Surv(time, outcome) ~ score + f.21003.0.0 + f.31.0.0 + 
                        f.21001.0.0, 
                      data = m_complete)
combined_score <- predict(cox_combined, type = "lp")

roc_combined <- timeROC(
  T = m_complete$time,
  delta = m_complete$outcome,
  marker = combined_score,
  cause = 1,
  weighting = "marginal",
  times = time_point,
  ROC = TRUE,
  iid = FALSE
)


auc_score <- round(roc_score$AUC[2], 3)
auc_clinical <- round(roc_clinical$AUC[2], 3)
auc_combined <- round(roc_combined$AUC[2], 3)


roc1_simple <- roc(m_complete$outcome, m_complete$score)
roc2_simple <- roc(m_complete$outcome, clinical_score)
roc3_simple <- roc(m_complete$outcome, combined_score)

delong_test1 <- roc.test(roc1_simple, roc2_simple, method = "delong")
delong_test2 <- roc.test(roc1_simple, roc3_simple, method = "delong")
delong_test3 <- roc.test(roc2_simple, roc3_simple, method = "delong")

cat("=== Time-dependent AUC Values ===\n")
cat("Metabolites only: AUC =", auc_score, "\n")
cat("Clinical variables only: AUC =", auc_clinical, "\n")
cat("Metabolites + Clinical variables: AUC =", auc_combined, "\n")

cat("\n=== DeLong Test Results ===\n")
cat("Metabolites vs Clinical: p =", round(delong_test1$p.value, 4), "\n")
cat("Metabolites vs Combined: p =", round(delong_test2$p.value, 4), "\n")
cat("Clinical vs Combined: p =", round(delong_test3$p.value, 4), "\n")

plot(roc_score, time = time_point, col = "#BC3C29FF", title = FALSE, lwd = 2)
plot(roc_clinical, time = time_point, add = TRUE, col = "#0000CD", lwd = 2)
plot(roc_combined, time = time_point, add = TRUE, col = "#2E8B57", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "gray")


legend("bottomright", 
       c(paste0("Metabolites only (C-statistics = ", auc_score, ")"),
         paste0("Clinical variables only (C-statistics = ", auc_clinical, ")"),
         paste0("Combined (C-statistics = ", auc_combined, ")")),
       col = c("#BC3C29FF", "#0000CD", "#2E8B57"), 
       lty = 1, lwd = 2, cex = 0.8, bty = "n")

text(x = 0.5, y = 0.20, cex = 0.85,
     labels = paste0("DeLong Test P-values:\n",
                     "Metabolites vs Clinical: ", 
                     ifelse(delong_test1$p.value < 0.001, "P < 0.001", 
                            paste0("P = ", round(delong_test1$p.value, 4))), "\n",
                     "Metabolites vs Combined: ", 
                     ifelse(delong_test2$p.value < 0.001, "P < 0.001", 
                            paste0("P = ", round(delong_test2$p.value, 4))), "\n",
                     "Clinical vs Combined: ", 
                     ifelse(delong_test3$p.value < 0.001, "P < 0.001", 
                            paste0("P = ", round(delong_test3$p.value, 4)))),
     adj = 0, col = "black")

title("5-Year Time-dependent ROC Curves for CRC Prediction")

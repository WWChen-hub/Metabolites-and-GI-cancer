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
########## metabolic scores and CRC#############
library(survival)
#library(glmnet)
library(data.table)
library(dplyr)
library(ggplot2)
library(survminer)

sample_df = data.table::fread("/data/cox_example.csv")
persd = sample_df
z <- fread('/data/lasso-CRC.csv',header = T)

score<- persd %>% select(z$var_names)
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
cox_model <- coxph(Surv(time, outcome) ~ tertile + f.21003.0.0 + f.31.0.0 + censor + ethnic2 + Townsend + fasting_status + smoking_staus + drinking_freq + family_history + bmi.class + MDscore + MET2, data = m)
cox_summary <- summary(cox_model)

hr_t2 <- round(cox_summary$coefficients[1, 2], 2)  
ci_t2_lower <- round(cox_summary$conf.int[1, 3], 2) 
ci_t2_upper <- round(cox_summary$conf.int[1, 4], 2)  

hr_t3 <- round(cox_summary$coefficients[2, 2], 2)  
ci_t3_lower <- round(cox_summary$conf.int[2, 3], 2) 
ci_t3_upper <- round(cox_summary$conf.int[2, 4], 2) 

legend_labels <- c(
  paste0("T1: HR (Ref)"),
  paste0("T2: HR ", hr_t2, " 95%CI: ", ci_t2_lower, "-", ci_t2_upper),
  paste0("T3: HR ", hr_t3, " 95%CI: ", ci_t3_lower, "-", ci_t3_upper)
)

survival_obj <- Surv(time = m$time, event = m$outcome)
km_fit <- survfit(survival_obj ~ tertile, data = m)


km_plot <- ggsurvplot(km_fit, 
                      data = m, 
                      fun = "event",
                      conf.int = TRUE,
                      risk.table = FALSE,
                      censor = FALSE,
                      xlab = "Time (Years)", 
                      ylab = "Cumulative Probability of EC",
                      palette = c("green4", "blue4", "red3"),
                      ylim = c(0, 0.06),
                      break.x.by = 2.5,
                      legend.labs = legend_labels,  # 修改图例标签
                      legend.title = "")  # 移除图例标题


summary_km_5_years <- summary(km_fit, times = 5)
summary_km_10_years <- summary(km_fit, times = 10)


survival_t1_5 <- summary_km_5_years$surv[1]
survival_t2_5 <- summary_km_5_years$surv[2]
survival_t3_5 <- summary_km_5_years$surv[3]

cumulative_incidence_t1_5 <- 1 - survival_t1_5
cumulative_incidence_t2_5 <- 1 - survival_t2_5
cumulative_incidence_t3_5 <- 1 - survival_t3_5


survival_t1_10 <- summary_km_10_years$surv[1]
survival_t2_10 <- summary_km_10_years$surv[2]
survival_t3_10 <- summary_km_10_years$surv[3]

cumulative_incidence_t1_10 <- 1 - survival_t1_10
cumulative_incidence_t2_10 <- 1 - survival_t2_10
cumulative_incidence_t3_10 <- 1 - survival_t3_10


km_plot$plot <- km_plot$plot +
  geom_vline(xintercept = c(5, 10), linetype = "dashed", alpha = 0.5, color = "gray50") +

  annotate("text", x = 5.2, y = 0.04, 
           label = paste0("T1: ", round(cumulative_incidence_t1_5, 4)),
           color = "green4", size = 3, hjust = 0, vjust = 0) +
  annotate("text", x = 5.2, y = 0.038, 
           label = paste0("T2: ", round(cumulative_incidence_t2_5, 4)),
           color = "blue4", size = 3, hjust = 0, vjust = 0) +
  annotate("text", x = 5.2, y = 0.036, 
           label = paste0("T3: ", round(cumulative_incidence_t3_5, 4)),
           color = "red3", size = 3, hjust = 0, vjust = 0) +
  annotate("text", x = 10.2, y = 0.04, 
           label = paste0("T1: ", round(cumulative_incidence_t1_10, 4)),
           color = "green4", size = 3, hjust = 0, vjust = 0) +
  annotate("text", x = 10.2, y = 0.038, 
           label = paste0("T2: ", round(cumulative_incidence_t2_10, 4)),
           color = "blue4", size = 3, hjust = 0, vjust = 0) +
  annotate("text", x = 10.2, y = 0.036, 
           label = paste0("T3: ", round(cumulative_incidence_t3_10, 4)),
           color = "red3", size = 3, hjust = 0, vjust = 0) +

  annotate("text", x = 5, y = 0.055, label = "5 Years", hjust = -0.1, size = 3, fontface = "bold") +
  annotate("text", x = 10, y = 0.055, label = "10 Years", hjust = -0.1, size = 3, fontface = "bold") +

  theme(legend.position = c(0.21, 0.8), 
        legend.background = element_rect(fill = "white", color = "white"),
        legend.key = element_rect(fill = "white"))

cat("T2 vs T1: HR =", hr_t2, "95%CI:", ci_t2_lower, "-", ci_t2_upper, "\n")
cat("T3 vs T1: HR =", hr_t3, "95%CI:", ci_t3_lower, "-", ci_t3_upper, "\n")

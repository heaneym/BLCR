data_full <- read.csv('./eleanor_data_full.csv')

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



# Scatter plot of vehicle EPO vs. LPS EPO

plot(
  data_full$EPO, 
  data_full$EPO_LPS, 
  col = cbPalette[2:8][factor(data_full$riskgroup)],
  pch = 16,
  xlab = 'EPO Vehicle',
  ylab = 'EPO LPS',
  main = 'Scatter Plot of LPS EPO vs. Vehicle EPO' 
)

legend(
  'topright',
  legend = c('Controls', 'Cerebral Palsy', 'High Risk'),
  pch = 16,
  col = cbPalette[2:8][factor(levels(factor(data_full$riskgroup)))]
)


# Creating a bar chart of percentage increase in the value of vehicle by LPS for each of the groups

data_full$EPO_percentage_increase <- 100*(data_full$EPO_LPS - data_full$EPO)/(data_full$EPO) 

CN_percentage_increase_EPO_mean <- mean(data_full$EPO_percentage_increase[data_full$riskgroup == 'CN'], na.rm = TRUE)
HR_percentage_increase_EPO_mean <- mean(data_full$EPO_percentage_increase[data_full$riskgroup == 'HR'], na.rm = TRUE)
CP_percentage_increase_EPO_mean <- mean(data_full$EPO_percentage_increase[data_full$riskgroup == 'CP'], na.rm = TRUE)
  
CN_percentage_increase_EPO_ci <- t.test(data_full$EPO_percentage_increase[data_full$riskgroup == 'CN'])$conf.int
HR_percentage_increase_EPO_ci <- t.test(data_full$EPO_percentage_increase[data_full$riskgroup == 'HR'])$conf.int
CP_percentage_increase_EPO_ci <- t.test(data_full$EPO_percentage_increase[data_full$riskgroup == 'CP'])$conf.int

means <- c(CN_percentage_increase_EPO_mean, HR_percentage_increase_EPO_mean, CP_percentage_increase_EPO_mean)
ci <- rbind(CN_percentage_increase_EPO_ci, HR_percentage_increase_EPO_ci, CP_percentage_increase_EPO_ci)

bp <- barplot(
  names = c('Controls','High Risk','Cerebral Palsy'),
  height = log(means),
  ylim = c(min(ci[,2]), max(ci[,2]))*1.2
)
arrows(x0 = bp, y0 = ci[,1],
       x1 = bp, y1 = ci[,2],
       angle = 90, code = 3, length = 0.05, lwd = 1.5)

boxplot(EPO_percentage_increase ~ riskgroup, data = data_full, ylim = c(-100,100))

# Assessing differences using Kruskal-Wallis
kw_test_EPO_percentage_increase <- kruskal.test(EPO_percentage_increase ~ riskgroup, data = data_full)
# Applying ANOVA test
anova_test_EPO_percentage_increase <- aov(EPO_percentage_increase ~ riskgroup, data = data_full)

# Kruskal-Wallis gives a non-significant p, but anova gives almost significant one (likely an artifact of the extreme outliers)

# Looking at the log ratio values instead of the percentage increase

boxplot(EPO_diff_log2 ~ riskgroup, data = data_full)

# Assessing differences using Kruskal-Wallis
kw_test_EPO_log_ratio <- kruskal.test(EPO_diff_log2 ~ riskgroup, data = data_full)
# Applying ANOVA test
anova_test_EPO_log_ratio <- aov(EPO_diff_log2 ~ riskgroup, data = data_full)

#Similar results found for log ratio







# Looking at the CSHQ total scores now and testing for differences between groups

boxplot(CSHQ_total ~ riskgroup, data = data_full)

# Formal test with Kruskal-Wallis

kw_test_CSHQ_scores <- kruskal.test(CSHQ_total ~ riskgroup, data = data_full)

# kw test gives a significant p value, looking at the post-hoc Dunn's test to see where these differences lie


dunn_test_CSHQ <- dunnTest(CSHQ_total ~ riskgroup, data = data_full)

#This gives a significant difference between the CP group and the HR group. 
tapply(data_full$CSHQ_total, data_full$riskgroup, mean, na.rm = TRUE)
tapply(data_full$CSHQ_total, data_full$riskgroup, median, na.rm = TRUE)
# CSHQ total is significantly greater in the CP group than in the HR group

















library(BayesLCA)
library(mcclust.ext)
library(pheatmap)
library(misty)
library(mclust)
library(viridis)
library(mclust)
library(FSA)
library(car)
library(mice)
source('./LCR_Gibbs.R')

ord <- c("BR1", "BR2", "BR3", "BR5",
         "SOD1",
         "SD1", "SD2", "SD3",
         "SA1", "SA2", "SA3", "SA4",
         "NW1", "NW2", "NW3",
         "P1", "P2", "P3", "P4", "P5", "P6", "P7",
         "SDB1", "SDB2", "SDB3",
         "DS1", "DS2", "DS3", "DS4", "DS5", "DS6", "DS7", "DS8")




CSHQ_LCR_varsel <- readRDS('./4 group BLCR model with variable selection (CSHQ)/CSHQ_LCR_varsel.RDS')

#Focusing on the Firefly Dataset first

firefly_control_raw <- read.csv("./Eleanor_CSHQ_data/Firefly_rawscores_CN.csv",
                                stringsAsFactors = FALSE)
firefly_control_raw_totals <- read.csv("./Eleanor_CSHQ_data/CSHQ_data_Eleanor_firefly_controls.csv")

data_firefly_control <- firefly_control_raw[2:21, c(1, 6:23, 25:39)]
colnames(data_firefly_control) <- c('ID', 'BR1', 'BR2', 'SOD1', 'SD2', 'SD3', 'DS1', 'DS7', 'DS8',
                                    'BR3', 'BR5', 'SA1', 'SA3', 'SD1', 'SA2', 'SA4', 'NW1', 
                                    'NW2', 'NW3', 'P2', 'P3', 'P4', 'P1', 'P5', 'P7', 'P6', 
                                    'SDB1', 'SDB2', 'SDB3', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6')

data_firefly_control[, -1] <- lapply(data_firefly_control[, -1], function(x) {
  as.numeric(trimws(as.character(x)))
})

data_firefly_control[data_firefly_control == ''] <- NA

#data_firefly_control[,2:7] <- 4 - data_firefly_control[,2:7]
data_firefly_control[,c(8,9)] <- data_firefly_control[,c(8,9)] + 1

data_firefly_control <- data_firefly_control[, c("ID", ord)]


#adding a column to indicate the original dataset

data_firefly_control$dataset <- 'firefly_control'


Y_firefly_control <- data_firefly_control[,-c(1,35)]


firefly_NE_raw <- read.csv('./Eleanor_CSHQ_data/Firefly_rawscores_Pt.csv')

firefly_NE_raw_totals <- read.csv('./Eleanor_CSHQ_data/CSHQ_data_Eleanor_firefly_NE.csv')

data_firefly_NE <- firefly_NE_raw[2:35, c(1, 6:23, 25:39)]

colnames(data_firefly_NE) <- c('ID', 'BR1', 'BR2', 'SOD1', 'SD2', 'SD3', 'DS1', 'DS7', 'DS8',
                               'BR3', 'BR5', 'SA1', 'SA3', 'SD1', 'SA2', 'SA4', 'NW1', 
                               'NW2', 'NW3', 'P2', 'P3', 'P4', 'P1', 'P5', 'P7', 'P6', 
                               'SDB1', 'SDB2', 'SDB3', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6')

data_firefly_NE[, -1] <- lapply(data_firefly_NE[, -1], function(x) {
  as.numeric(trimws(as.character(x)))
})

data_firefly_NE[data_firefly_NE == ''] <- NA

#data_firefly_NE[,2:7] <- 4 - data_firefly_NE[,2:7]
data_firefly_NE[,c(8,9)] <- data_firefly_NE[,c(8,9)] + 1




data_firefly_NE <- data_firefly_NE[, c("ID", ord)]



#adding a column to indicate the original dataset

data_firefly_NE$dataset <- 'firefly_NE'

#There is one observation (row 34) that has a value of 4 in the item DS6 - Going to assume this a typo and should be a 3

data_firefly_NE[which(data_firefly_NE$DS8 == 4), 'DS8'] <- 3

Y_firefly_NE <- data_firefly_NE[,-c(1,35)]

data_firefly_full <- rbind(data_firefly_control, data_firefly_NE)

Y_firefly_full <- rbind(Y_firefly_control, Y_firefly_NE)


#Importing the cytokine data

# firefly_cytokine_raw <- read.csv('./Eleanor_CSHQ_data/Cytokines Firefly neonates and Follow-up 27.02.csv')





#Working with the starfish data now

starfish_preterm_raw <- read.csv("./Eleanor_CSHQ_data/starfish_preterm_followup_raw_CSHQ.csv")



data_starfish_preterm <- starfish_preterm_raw[,c(1, 10:15, 18, 19, 22:46)]

colnames(data_starfish_preterm) <- c('ID', 'BR1', 'BR2', 'SOD1', 'SD2', 'SD3', 'DS1', 'DS7', 'DS8',
                                     'BR3', 'BR5', 'SA1', 'SA3', 'SD1', 'SA2', 'SA4', 'NW1', 
                                     'NW2', 'NW3', 'P2', 'P3', 'P4', 'P1', 'P5', 'P7', 'P6', 
                                     'SDB1', 'SDB2', 'SDB3', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6')

data_starfish_preterm[, -1] <- lapply(data_starfish_preterm[, -1], function(x) {
  as.numeric(trimws(as.character(x)))
})

data_starfish_preterm[data_starfish_preterm == ''] <- NA

#data_starfish_preterm[,2:7] <- 4 - data_starfish_preterm[,2:7]
data_starfish_preterm[,c(8,9)] <- data_starfish_preterm[,c(8,9)] + 1



data_starfish_preterm <- data_starfish_preterm[, c("ID", ord)]

#There's a row containing values of 4 here - not sure whats going on exactly

invalid_response_nrow <- sum(rowSums(data_starfish_preterm[,-1] == 4, na.rm = TRUE) > 0)

#For the moment we'll remove this row (row 17 - ID S17)

data_starfish_preterm <- data_starfish_preterm[-which(data_starfish_preterm$ID == 'S17'),]




#adding a column to indicate the original dataset

data_starfish_preterm$dataset <- 'starfish_preterm' 


Y_starfish_preterm <- data_starfish_preterm[,-c(1,35)]



#Serenity Data

serenity_raw <- read.csv('./Eleanor_CSHQ_data/SERENITY_raw_data.csv')

serenity_raw_totals <- read.csv('./Eleanor_CSHQ_data/CSHQ_data_Eleanor_serenity.csv')

data_serenity <- serenity_raw[3:21,c(1,6:23, 25:39)]

colnames(data_serenity) <- c('ID', 'BR1', 'BR2', 'SOD1', 'SD2', 'SD3', 'DS1', 'DS7', 'DS8',
                                     'BR3', 'BR5', 'SA1', 'SA3', 'SD1', 'SA2', 'SA4', 'NW1', 
                                     'NW2', 'NW3', 'P2', 'P3', 'P4', 'P1', 'P5', 'P7', 'P6', 
                                     'SDB1', 'SDB2', 'SDB3', 'DS2', 'DS3', 'DS4', 'DS5', 'DS6')


data_serenity[, -1] <- lapply(data_serenity[, -1], function(x) {
  as.numeric(trimws(as.character(x)))
})




data_serenity[data_serenity == ''] <- NA

#data_serenity[,2:7] <- 4 - data_serenity[,2:7]
data_serenity[,c(8,9)] <- data_serenity[,c(8,9)] + 1




data_serenity <- data_serenity[, c("ID", ord)]



data_serenity$dataset <- 'serenity' 


Y_serenity <- data_serenity[,-c(1,35)]





#Creating a combined dataframe of firefly, starfish and serenity as well as the corresponding 
#response matrix Y.


data_full <- rbind(data_firefly_full, data_starfish_preterm, data_serenity)

Y_full <- rbind(Y_firefly_full, Y_starfish_preterm, Y_serenity)


#Removing any rows of entirely NAs (Specifically for the Survey responses)

data_full <- data_full[rowSums(is.na(Y_full)) != ncol(Y_full), ]

Y_full <- Y_full[rowSums(is.na(Y_full)) != ncol(Y_full), ]

row.names(Y_full) <- 1:nrow(Y_full)

row.names(data_full) <- 1:nrow(data_full)


#We need to have the columns in a consistent order with the previous analysis

data_full <- data_full[, c("ID", ord, "dataset")]

Y_full <- Y_full[, ord]




#Looking at summaries first

#Histogram of age distributions between cohorts
x_lower <- min(c(as.numeric(firefly_control_raw_totals$X.1[3:22]),
                 as.numeric(firefly_NE_raw_totals$X.1[3:31]),
                 starfish_preterm_raw$AGE..years.*12,
                 as.numeric(serenity_raw_totals$X.1[3:30])), na.rm = TRUE)
x_upper <- max(c(as.numeric(firefly_control_raw_totals$X.1[3:22]),
                 as.numeric(firefly_NE_raw_totals$X.1[3:31]),
                 starfish_preterm_raw$AGE..years.*12,
                 as.numeric(serenity_raw_totals$X.1[3:30])), na.rm = TRUE)

pdf('./Eleanor_CSHQ_data/age_hist1.pdf', height = 8, width = 10)
png('./Eleanor_CSHQ_data/age_hist1.png', width = 1000, height = 800, units = "px")
par(mfrow=c(2, 2))

#Assuming that serenity contains age in months
hist(as.numeric(firefly_control_raw_totals$X.1[3:22]), main="Histogram of Age for Firefly Controls", xlab="Age (months)", ylab="Frequency", col="lightblue", xlim = c(x_lower,x_upper), probability = TRUE)
hist(as.numeric(firefly_NE_raw_totals$X.1[3:31]), main="Histogram of Age for Firefly NE", xlab="Age (months)", ylab="Frequency", col="lightgreen", xlim = c(x_lower,x_upper), probability = TRUE)
hist(starfish_preterm_raw$AGE..years.*12, main="Histogram of Age for Starfish", xlab="Age (months)", ylab="Frequency", col="orange", xlim = c(x_lower,x_upper), probability = TRUE)
hist(as.numeric(serenity_raw_totals$X.1[3:30]), main="Histogram of Age for Serenity", xlab="Age (months)", ylab="Frequency", col="purple", xlim = c(x_lower,x_upper), probability = TRUE)


dev.off()
dev.off()

pdf('./Eleanor_CSHQ_data/age_hist2.pdf', height = 8, width = 10)
png('./Eleanor_CSHQ_data/age_hist2.png', width = 1000, height = 800, units = "px")
par(mfrow=c(2, 2))

x_lower <- min(c(as.numeric(firefly_control_raw_totals$X.1[3:22]),
                 as.numeric(firefly_NE_raw_totals$X.1[3:31]),
                 starfish_preterm_raw$AGE..years.*12,
                 as.numeric(serenity_raw_totals$X.1[3:30])*0.230137), na.rm = TRUE)
x_upper <- max(c(as.numeric(firefly_control_raw_totals$X.1[3:22]),
                 as.numeric(firefly_NE_raw_totals$X.1[3:31]),
                 starfish_preterm_raw$AGE..years.*12,
                 as.numeric(serenity_raw_totals$X.1[3:30])*0.230137), na.rm = TRUE)

#Assuming that serenity contains age in weeks - TBC?
hist(as.numeric(firefly_control_raw_totals$X.1[3:22]), main="Histogram of Age for Firefly Controls", xlab="Age (months)", ylab="Frequency", col="lightblue", xlim = c(x_lower,x_upper), probability = TRUE)
hist(as.numeric(firefly_NE_raw_totals$X.1[3:31]), main="Histogram of Age for Firefly NE", xlab="Age (months)", ylab="Frequency", col="lightgreen", xlim = c(x_lower,x_upper), probability = TRUE)
hist(starfish_preterm_raw$AGE..years.*12, main="Histogram of Age for Starfish", xlab="Age (months)", ylab="Frequency", col="orange", xlim = c(x_lower,x_upper), probability = TRUE)
hist(as.numeric(serenity_raw_totals$X.1[3:30])*0.230137, main="Histogram of Age for Serenity", xlab="Age (months)", ylab="Frequency", col="purple", xlim = c(x_lower,x_upper), probability = TRUE)

dev.off()
dev.off()


# Missingness

#Creating a heatmap of missingness

Y_missing <- na.indicator(Y_full)

png('./Eleanor_CSHQ_data/missingness_heatmap.png', width = 1000, height = 800, units = "px")
pheatmap(1-Y_missing[,34:66],
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = c('red','yellow'), 
         gaps_row = c(20,54,75),
         gaps_col = c(4,5,8,12,15,22,25,33),
         labels_row = rep('', nrow(Y_missing)),
         labels_col = colnames(Y_full),
         main = 'Missing Indicator Heatmap')
dev.off()


#Looking at item-level summaries

firefly_control_item_mean <- colMeans(Y_full[data_full$dataset == 'firefly_control',], na.rm = TRUE)
firefly_NE_item_mean <- colMeans(Y_full[data_full$dataset == 'firefly_NE',], na.rm = TRUE)


#Plotting a heatmap of the CSHQ responses
png('./Eleanor_CSHQ_data/CSHQ_response_heatmap.png', width = 1000, height = 800, units = "px")
pheatmap(Y_full, cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = viridis(3), 
         gaps_row = c(20,54,75),
         gaps_col = c(4,5,8,12,15,22,25,33),
         labels_row = rep('', nrow(Y_full)),
         labels_col = colnames(Y_full),
         main = 'CSHQ Response Heatmap')

dev.off()

#Creating a matrix to plot a heatmap for just the cohort means and sds

firefly_control_item_mean <- colMeans(Y_full[data_full$dataset == 'firefly_control',], na.rm = TRUE)
firefly_NE_item_mean <- colMeans(Y_full[data_full$dataset == 'firefly_NE',], na.rm = TRUE)
starfish_item_mean <- colMeans(Y_full[data_full$dataset == 'starfish_preterm',], na.rm = TRUE)
serenity_item_mean <- colMeans(Y_full[data_full$dataset == 'serenity',], na.rm = TRUE)

cohorts_item_mean_matrix <- rbind(firefly_control_item_mean, firefly_NE_item_mean, starfish_item_mean, serenity_item_mean) 

#Adding in NA for the NW1 for the control cohort since almost all were missing this item
cohorts_item_mean_matrix[1,13] <- NA

firefly_control_item_sd <- apply(Y_full[data_full$dataset == 'firefly_control',], 2, sd, na.rm = TRUE)
firefly_NE_item_sd <- apply(Y_full[data_full$dataset == 'firefly_NE',], 2, sd, na.rm = TRUE)
starfish_item_sd <- apply(Y_full[data_full$dataset == 'starfish_preterm',], 2, sd, na.rm = TRUE)
serenity_item_sd <- apply(Y_full[data_full$dataset == 'serenity',], 2, sd, na.rm = TRUE)

cohorts_item_sd_matrix <- rbind(firefly_control_item_sd, firefly_NE_item_sd, starfish_item_sd, serenity_item_sd) 

cell_labels <- round(cohorts_item_mean_matrix, 2)


png('./Eleanor_CSHQ_data/CSHQ_cohort_mean_response_heatmap.png', width = 14, height = 3, units = "in", res = 600)

pheatmap(cohorts_item_mean_matrix,
         display_numbers = cell_labels,          
         cluster_rows = FALSE,                   
         cluster_cols = FALSE,                    
         color = colorRampPalette(c("green", "white", "red"))(100),
         #color = viridis(100),
         main = "Mean CSHQ item scores by cohort",
         fontsize_number = 9,
         fontcolor_number = "gray30",
         angle_col = 45,
         border_color = "grey60",
         cellwidth = 25,
         cellheight = 35,
         labels_row = c('firefly controls', 'firefly NE', 'starfish', 'serenity'))

dev.off()



# 
# #Creating a similar heatmap but summing to get subscale totals rather than items - We use rowMeans to account for differing numbers of items in subscales
# cohorts_subscale_mean_matrix <- cbind(
#   rowMeans(cohorts_item_mean_matrix[,1:4], na.rm = TRUE),
#   cohorts_item_mean_matrix[,5],
#   rowMeans(cohorts_item_mean_matrix[,6:8], na.rm = TRUE),
#   rowMeans(cohorts_item_mean_matrix[,9:12], na.rm = TRUE),
#   rowMeans(cohorts_item_mean_matrix[,13:15], na.rm = TRUE),
#   rowMeans(cohorts_item_mean_matrix[,16:22], na.rm = TRUE),
#   rowMeans(cohorts_item_mean_matrix[,23:25], na.rm = TRUE),
#   rowMeans(cohorts_item_mean_matrix[,26:33], na.rm = TRUE)
# )
# colnames(cohorts_subscale_mean_matrix) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'SDB', 'DS')
# 
# 
# #We can add in the Champion NE group here since we have the subscale totals for this cohort
# 
# champion_NE_subscales <- read.csv("./Eleanor_CSHQ_data/CSHQ_data_Eleanor_champion_NE.csv")
# #getting only the IDs and subscale totals
# champion_NE_subscales_only <- champion_NE_subscales[3:55,c(1,3:10)]
# #Removing empty rows
# champion_NE_subscales_only <- champion_NE_subscales_only[-which(is.na(champion_NE_subscales_only$Bed.time.res)),]
# colnames(champion_NE_subscales_only) <- c('ID', 'BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'SDB', 'DS')
# #relevelling DS column due to difference in scoring for 2 items
# champion_NE_subscales_only$DS <- champion_NE_subscales_only$DS + 2
# #There is a 0 in one of the entries - we can assume missingness and impute?
# champion_NE_subscales_only$SD[which(champion_NE_subscales_only$SD == 0)] <- NA
# champion_imp <- mice(champion_NE_subscales_only, method = 'pmm', seed = 123)
# champion_NE_subscales_only <- complete(champion_imp)
# 
# cohorts_subscale_mean_matrix <- rbind(cohorts_subscale_mean_matrix , colMeans(champion_NE_subscales_only[,2:9])/(c(4,1,3,4,3,7,3,8)))
# 
# cell_labels_subtotals <- round(cohorts_subscale_mean_matrix, 2)
# 
# 
# png('./Eleanor_CSHQ_data/CSHQ_cohort_mean_subtotal_heatmap.png', width = 5, height = 4, units = "in", res = 600)
# 
# pheatmap(cohorts_subscale_mean_matrix,
#          display_numbers = cell_labels_subtotals,          
#          cluster_rows = FALSE,                   
#          cluster_cols = FALSE,                    
#          color = colorRampPalette(c("green", "white", "red"))(100),
#          #color = viridis(100),
#          main = "Mean CSHQ item scores by cohort",
#          fontsize_number = 9,
#          fontcolor_number = "gray30",
#          angle_col = 45,
#          border_color = "grey60",
#          cellwidth = 25,
#          cellheight = 35,
#          labels_row = c('firefly controls', 'firefly NE', 'starfish', 'serenity', 'champion NE'))
# 
# dev.off()



# Computing proportions above threshold of 41 
# (range of possible proportions given missingness)

# We'll drop the rows with severe amounts of missingness

data_full_clean <- data_full[rowSums(is.na(Y_full)) <= 16,]
Y_full_clean <- Y_full[rowSums(is.na(Y_full)) <= 16,]


row_score_full_sample <- rowSums(Y_full_clean, na.rm = TRUE)
row_missing_full_sample <- rowSums(is.na(Y_full_clean))

row_min_possible_full_sample <- row_score_full_sample + row_missing_full_sample*1
row_max_possible_full_sample <- row_score_full_sample + row_missing_full_sample*3
sum(row_min_possible_full_sample>41)/nrow(Y_full_clean)
sum(row_max_possible_full_sample>41)/nrow(Y_full_clean)




row_score_firefly_control <- rowSums(Y_full_clean[which(data_full_clean$dataset == 'firefly_control'),], na.rm = TRUE)
row_missing_firefly_control <- rowSums(is.na(Y_full_clean[which(data_full_clean$dataset == 'firefly_control'),]))

row_min_possible_firefly_control <- row_score_firefly_control + row_missing_firefly_control*1
row_max_possible_firefly_control <- row_score_firefly_control + row_missing_firefly_control*3
sum(row_min_possible_firefly_control>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'firefly_control'),])
sum(row_max_possible_firefly_control>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'firefly_control'),])



row_score_firefly_NE <- rowSums(Y_full_clean[which(data_full_clean$dataset == 'firefly_NE'),], na.rm = TRUE)
row_missing_firefly_NE <- rowSums(is.na(Y_full_clean[which(data_full_clean$dataset == 'firefly_NE'),]))

row_min_possible_firefly_NE <- row_score_firefly_NE + row_missing_firefly_NE*1
row_max_possible_firefly_NE <- row_score_firefly_NE + row_missing_firefly_NE*3
sum(row_min_possible_firefly_NE>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'firefly_NE'),])
sum(row_max_possible_firefly_NE>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'firefly_NE'),])




row_score_starfish <- rowSums(Y_full_clean[which(data_full_clean$dataset == 'starfish_preterm'),], na.rm = TRUE)
row_missing_starfish <- rowSums(is.na(Y_full_clean[which(data_full_clean$dataset == 'starfish_preterm'),]))

row_min_possible_starfish <- row_score_starfish + row_missing_starfish*1
row_max_possible_starfish <- row_score_starfish + row_missing_starfish*3
sum(row_min_possible_starfish>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'starfish_preterm'),])
sum(row_max_possible_starfish>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'starfish_preterm'),])






row_score_serenity <- rowSums(Y_full_clean[which(data_full_clean$dataset == 'serenity'),], na.rm = TRUE)
row_missing_serenity <- rowSums(is.na(Y_full_clean[which(data_full_clean$dataset == 'serenity'),]))

row_min_possible_serenity <- row_score_serenity + row_missing_serenity*1
row_max_possible_serenity <- row_score_serenity + row_missing_serenity*3
sum(row_min_possible_serenity>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'serenity'),])
sum(row_max_possible_serenity>41)/nrow(Y_full_clean[which(data_full_clean$dataset == 'serenity'),])






#Conducting formal tests for score differences between cohorts

#Dropping the NW1 column for this analysis due to mostly missing values in control cohort 

Y_full_clean_minus_NW1 <- Y_full_clean[,-13]

#For the differences in total scores, we can include the champion cohorts as well, since we have subscale totals for some and CSHQ totals for all

champion_NE_subscales <- read.csv("./Eleanor_CSHQ_data/CSHQ_data_Eleanor_champion_NE.csv")
champion_CP_scores <- read.csv('./Eleanor_CSHQ_data/CSHQ_data_Eleanor_champion_CP.csv')
champion_CP_scores <- champion_CP_scores[3:21,]
champion_CP_scores_only <- champion_CP_scores[,c(1,4)]
colnames(champion_CP_scores_only) <- c('ID', 'total')
#getting only the IDs and subscale totals
champion_NE_subscales_only <- champion_NE_subscales[3:55,c(1,3:10)]
#Removing empty rows
champion_NE_subscales_only <- champion_NE_subscales_only[-which(is.na(champion_NE_subscales_only$Bed.time.res)),]
colnames(champion_NE_subscales_only) <- c('ID', 'BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'SDB', 'DS')
#relevelling DS column due to difference in scoring for 2 items
champion_NE_subscales_only$DS <- champion_NE_subscales_only$DS + 2
#There is a 0 in one of the entries - we can assume missingness and impute?
champion_NE_subscales_only$SD[which(champion_NE_subscales_only$SD == 0)] <- NA
champion_imp <- mice(champion_NE_subscales_only, method = 'pmm', seed = 123)
champion_NE_subscales_only <- complete(champion_imp)
#We'll need to deal somehow with the fact that the NW1 item is not included in the other cohort subtotals?
#For the moment, we'll split the NW subscale equally across items for each observation and remove the contribution
#This makes an assumption that each item in NW contributes equally to the subtotal for each observation in champion NE
#Not necessarily a great assumption, but it's only one item so shouldn't bias results super heavily 
#There are 3 items under NW. We relevel as
champion_NE_subscales_only$NW <- champion_NE_subscales_only$NW - champion_NE_subscales_only$NW/3


#There are two items that are repeated in the BR and SA subscale, we drop this from BR using a similar approach
champion_NE_subscales_only$BR <- champion_NE_subscales_only$BR - 2 * champion_NE_subscales_only$BR / 6


#Dealing with the champion CP data similarly
#Firstly relevelling to adjust for the question scoring differences
champion_CP_scores_only$total <- as.numeric(champion_CP_scores_only$total) + 2

#Secondly dealing with the removal of the NW1 column. We'll carry out a similar approcimation as the other one by just
# splitting scores across all items and subtracting this value. We have

champion_CP_scores_only$total <- champion_CP_scores_only$total - champion_CP_scores_only$total/33

  

#For the moment we'll use MICE to impute remaining missing values 

imputation_frame <- cbind(Y_full_clean_minus_NW1, data_full_clean$dataset)
imp <- mice(imputation_frame, m = 5, method = "pmm", seed = 123)

completed_Y_clean <- complete(imp)


#We'll compare subscale totals for the moment

BR_Y <- completed_Y_clean[,1:4]
SOD_Y <- completed_Y_clean[,5]
SD_Y <- completed_Y_clean[,6:8]
SA_Y <- completed_Y_clean[,9:12]
NW_Y <- completed_Y_clean[,13:14]
P_Y <- completed_Y_clean[,15:21]
SDB_Y <- completed_Y_clean[,22:24]
DS_Y <- completed_Y_clean[,25:32]

BR_total <- c(rowSums(BR_Y), champion_NE_subscales_only$BR)
SOD_total <- c(SOD_Y, champion_NE_subscales_only$SOD)
SD_total <- c(rowSums(SD_Y), champion_NE_subscales_only$SD)
SA_total <- c(rowSums(SA_Y), champion_NE_subscales_only$SA)
NW_total <- c(rowSums(NW_Y), champion_NE_subscales_only$NW)
P_total <- c(rowSums(P_Y), champion_NE_subscales_only$P)
SDB_total <- c(rowSums(SDB_Y), champion_NE_subscales_only$SDB)
DS_total <- c(rowSums(DS_Y), champion_NE_subscales_only$DS)
CSHQ_total <- rowSums(completed_Y_clean[,1:32])




# subscale_total_frame <- cbind(BR_total, SOD_total, 
#                               SD_total, SA_total,
#                               NW_total, P_total,
#                               SDB_total, DS_total,
#                               CSHQ_total,
#                               data_full_clean$dataset)

subscale_total_frame_with_champion_NE <- cbind(BR_total, SOD_total, 
                                               SD_total, SA_total, 
                                               NW_total, P_total,
                                               SDB_total, DS_total, 
                                               c(data_full_clean$dataset, rep('champion_NE', nrow(champion_NE_subscales_only))))

subscale_heatmap_mat <- apply(subscale_total_frame_with_champion_NE[,1:8], 2, as.numeric)

subscale_heatmap_mat_champion_NE <- subscale_heatmap_mat[which(subscale_total_frame_with_champion_NE[,9] == 'champion_NE'),]
subscale_heatmap_mat_firefly_control <- subscale_heatmap_mat[which(subscale_total_frame_with_champion_NE[,9] == 'firefly_control'),]
subscale_heatmap_mat_firefly_NE <- subscale_heatmap_mat[which(subscale_total_frame_with_champion_NE[,9] == 'firefly_NE'),]
subscale_heatmap_mat_starfish <- subscale_heatmap_mat[which(subscale_total_frame_with_champion_NE[,9] == 'starfish_preterm'),]
subscale_heatmap_mat_serenity <- subscale_heatmap_mat[which(subscale_total_frame_with_champion_NE[,9] == 'serenity'),]

subscale_heatmap_mat_scaled <- rbind(
  colMeans(subscale_heatmap_mat_firefly_control)/c(4,1,3,4,2,7,3,8),
  colMeans(subscale_heatmap_mat_firefly_NE)/c(4,1,3,4,2,7,3,8),
  colMeans(subscale_heatmap_mat_starfish)/c(4,1,3,4,2,7,3,8),
  colMeans(subscale_heatmap_mat_serenity)/c(4,1,3,4,2,7,3,8),
  colMeans(subscale_heatmap_mat_champion_NE)/c(4,1,3,4,2,7,3,8)
)

cell_labels_subtotals <- round(subscale_heatmap_mat_scaled, 2)


png('./Eleanor_CSHQ_data/CSHQ_cohort_mean_subtotal_heatmap.png', width = 5, height = 4, units = "in", res = 600)

pheatmap(subscale_heatmap_mat_scaled,
         display_numbers = cell_labels_subtotals,          
         cluster_rows = FALSE,                   
         cluster_cols = FALSE,                    
         color = colorRampPalette(c("green", "white", "red"))(100),
         #color = viridis(100),
         main = "Mean CSHQ item scores by cohort",
         fontsize_number = 9,
         fontcolor_number = "gray30",
         angle_col = 45,
         border_color = "grey60",
         cellwidth = 25,
         cellheight = 35,
         labels_row = c('firefly controls', 'firefly NE', 'starfish', 'serenity', 'champion NE'))

dev.off()

#colnames(subscale_total_frame)[10] <- 'cohort'

colnames(subscale_total_frame_with_champion_NE)[9] <- 'cohort'

subscale_total_frame_with_champion_NE <- as.data.frame(subscale_total_frame_with_champion_NE)

cols_to_convert <- setdiff(names(subscale_total_frame_with_champion_NE), "cohort")
subscale_total_frame_with_champion_NE[cols_to_convert] <- lapply(subscale_total_frame_with_champion_NE[cols_to_convert], as.numeric)

#Creating a separate from for the CSHQ total scores
CSHQ_total_frame <- data.frame(
  CSHQ_total = c(rowSums(completed_Y_clean[,1:32]), rowSums(champion_NE_subscales_only[,2:9]), champion_CP_scores_only$total),
  cohort = c(data_full_clean$dataset, rep('champion_NE', nrow(champion_NE_subscales_only)), rep('champion_CP', nrow(champion_CP_scores_only)))
)

# We want to assess whether assumptions hold to justify an anova test

# Equality of variances between groups (Levene test)

# BR_levine <- leveneTest(BR_total ~ cohort, data = as.data.frame(subscale_total_frame))
# SOD_levine <- leveneTest(SOD_total ~ cohort, data = as.data.frame(subscale_total_frame))
# SD_levine <- leveneTest(SD_total ~ cohort, data = as.data.frame(subscale_total_frame))
# SA_levine <- leveneTest(SA_total ~ cohort, data = as.data.frame(subscale_total_frame))
# NW_levine <- leveneTest(NW_total ~ cohort, data = as.data.frame(subscale_total_frame))
# P_levine <- leveneTest(P_total ~ cohort, data = as.data.frame(subscale_total_frame))
# SDB_levine <- leveneTest(SDB_total ~ cohort, data = as.data.frame(subscale_total_frame))
# DS_levine <- leveneTest(DS_total ~ cohort, data = as.data.frame(subscale_total_frame))
# CSHQ_levine <- leveneTest(CSHQ_total ~ cohort, data = as.data.frame(subscale_total_frame))

BR_levine <- leveneTest(BR_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
SOD_levine <- leveneTest(SOD_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
SD_levine <- leveneTest(SD_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
SA_levine <- leveneTest(SA_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
NW_levine <- leveneTest(NW_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
P_levine <- leveneTest(P_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
SDB_levine <- leveneTest(SDB_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
DS_levine <- leveneTest(DS_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
CSHQ_levine <- leveneTest(CSHQ_total ~ cohort, data = as.data.frame(CSHQ_total_frame))

#It appears that they all have non-significant p-values other than the 
# SD (sleep duration), SOD (sleep onset delay) ,SDB (sleep disordered breathing) subscales. Hence we can apply ANOVA for all but these subscale.
# We'll use a Kruskal-Wallis test for those ones


kruskal_SDB <- kruskal.test(SDB_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_SD <- kruskal.test(SD_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_SOD <- kruskal.test(SOD_total ~ cohort, data = subscale_total_frame_with_champion_NE)


anova_BR <- aov(BR_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
anova_SA <- aov(SA_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
anova_NW <- aov(NW_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
anova_P <- aov(P_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
anova_DS <- aov(DS_total ~ cohort, data = as.data.frame(subscale_total_frame_with_champion_NE))
anova_CSHQ <- aov(CSHQ_total ~ cohort, data = CSHQ_total_frame)

#Checking normality of residuals assumptions for anova tests

residuals_BR <- anova_BR$residuals
residuals_SA <- anova_SA$residuals
residuals_NW <- anova_NW$residuals
residuals_P <- anova_P$residuals
residuals_DS <- anova_DS$residuals
residuals_CSHQ <- anova_CSHQ$residuals


qqnorm(residuals_BR, main = 'SOD Residuals')
qqline(residuals_BR)
shapiro.test(residuals_BR)


qqnorm(residuals_SA, main = 'SA Residuals')
qqline(residuals_SA)
shapiro.test(residuals_SA)

qqnorm(residuals_NW, main = 'NW Residuals')
qqline(residuals_NW)
shapiro.test(residuals_NW)

qqnorm(residuals_P, main = 'SDB Residuals')
qqline(residuals_P)
shapiro.test(residuals_P)

qqnorm(residuals_DS, main = 'DS Residuals')
qqline(residuals_DS)
shapiro.test(residuals_DS)

qqnorm(residuals_CSHQ, main = 'CSHQ Residuals')
qqline(residuals_CSHQ)
shapiro.test(residuals_CSHQ)

#Normality of residuals is violated for all according to QQ-plot and residuals

#Hence we run the Kruskal-Wallis for all of the subscales

kruskal_BR <- kruskal.test(BR_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_SA <- kruskal.test(SA_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_NW <- kruskal.test(NW_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_P <- kruskal.test(P_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_DS <- kruskal.test(DS_total ~ cohort, data = subscale_total_frame_with_champion_NE)
kruskal_CSHQ <- kruskal.test(CSHQ_total ~ cohort, data = CSHQ_total_frame)

#Adjustment for multiple tests
p_vals <- c(
  kruskal_BR$p.value,
  kruskal_SOD$p.value,
  kruskal_SD$p.value,
  kruskal_SA$p.value,
  kruskal_NW$p.value,
  kruskal_P$p.value,
  kruskal_SDB$p.value,
  kruskal_DS$p.value,
  kruskal_CSHQ$p.value
) 

#Significant difference in CSHQ totals (p ~ 0.03)


#Adjusting for multiple testing in the subscale total tests (not including CSHQ total score in this correction)

adjust_p_fdr <- rbind(colnames(subscale_total_frame_with_champion_NE)[1:8] , 
                      round(p.adjust(p_vals[1:8], method = 'fdr'),4))

# Post-correction, we have p<0.05 for the subscales of bedtime resistance, sleep onset delay, sleep duration, sleep disordered breathing


# We apply the Dunn test to identify where the differences lie
dunn_SOD <- dunnTest(SOD_total ~ cohort, data = subscale_total_frame_with_champion_NE, method = "bh")
dunn_SD <- dunnTest(SD_total ~ cohort, data = subscale_total_frame_with_champion_NE, method = "bh")
dunn_NW <- dunnTest(NW_total ~ cohort, data = subscale_total_frame_with_champion_NE, method = "bh")
dunn_SDB <- dunnTest(SDB_total ~ cohort, data = subscale_total_frame_with_champion_NE, method = "bh")

dunn_SOD$res
dunn_SD$res
dunn_NW$res
dunn_SDB$res


# Bedtime Resistance results
# significant difference between champion NE and firefly control, champion NE and firefly NE, champion NE and serenity, champion NE and starfish preterm
# In particular, the BR subscale totals are significantly greater for champion NE than any other group

# Sleep Onset Delay results
# Significant difference between champion NE and serenity, firefly NE and serenity, serenity and starfish preterm
# In particular, serenity had a higher mean SOD subtotal when compared with these other groups


# Sleep Duration results
# Significant difference between champion NE and serenity, firefly NE and serenity, serenity and starfish preterm
# serenity has a significantly greater subscale total than these other subscales

#Sleep Duration results
#

# Sleep Disordered Breathing results
# Significant difference between champion NE and serenity, firefly control and serenity, firefly NE and serenity, serenity and starfish preterm
# serenity has a significantly elevated subscale total compared with these groups.

# CSHQ results
# post fdr correction, there is are no significant differences up to alpha = 0.05 significance - needs a futher look maybe?







#Analysis of cytokine data

#Reading in the data
firefly_NE_old_frame <- read.csv('./Eleanor_CSHQ_data/CSHQ_data_Eleanor_firefly_NE.csv')
firefly_control_old_frame <- read.csv('./Eleanor_CSHQ_data/CSHQ_data_Eleanor_firefly_controls.csv') 

firefly_NE_age_frame <- firefly_NE_old_frame[,c(1,3)]
firefly_control_age_frame <- firefly_control_old_frame[,c(1,3)]

firefly_cytokine_raw <- read.csv("./Eleanor_CSHQ_data/FIREFLY_cytokines.csv")

firefly_cytokine_raw_vehicle <- firefly_cytokine_raw[2:52, ]

firefly_sarnat <- read.csv('./Eleanor_CSHQ_data/FIREFLY_sarnat.csv')

firefly_sarnat <- firefly_sarnat[,2:4]

colnames(firefly_cytokine_raw_vehicle)[2] <- 'ID'

cytokine_X_firefly <- firefly_cytokine_raw_vehicle[,3:6]

#Merging cytokin frame with CSHQ frame, only retaining ID values that contain both
firefly_sleep_cytokine_merged <- merge(data_firefly_full, firefly_cytokine_raw_vehicle, by = "ID")

#we have 6 individuals lost due to merging

#We again remove the rows with large amounts of missingness, as well as the column of NW1
firefly_sleep_cytokine_merged <- firefly_sleep_cytokine_merged[
  -c(6, 13, 32),
  -c(which(names(firefly_sleep_cytokine_merged) == "NW1"), 36, 41)
]
#There is a single missing value in the cytokines, I'll remove that row for the moment
firefly_sleep_cytokine_merged <- firefly_sleep_cytokine_merged[-21,]

#Carrying out imputation for sleep survey responses as before
imp <- mice(firefly_sleep_cytokine_merged[,2:34], m = 5, method = "pmm", seed = 123)

firefly_sleep_cytokine_merged_completed <- firefly_sleep_cytokine_merged_completed <- cbind(
  ID = firefly_sleep_cytokine_merged[, 1],
  complete(imp),
  firefly_sleep_cytokine_merged[, 35:48]
)

colnames(firefly_sleep_cytokine_merged_completed) <- colnames(firefly_sleep_cytokine_merged_completed)

firefly_sleep_cytokine_merged_completed[,36:48] <- sapply(firefly_sleep_cytokine_merged_completed[,36:48],as.numeric)

ggpairs(firefly_sleep_cytokine_merged_completed[,34:38], columns = 2:5, aes(color = dataset))
ggsave('./Eleanor_CSHQ_data/firefly_cytokine_pairs_plot.png', width = 10, height = 10)

res <- rcorr(as.matrix(firefly_sleep_cytokine_merged_completed[,34:38][,2:5]), type = "pearson")

res$r   
res$P 

#It seems like there is a positive correlation between VEGF and IL8 (maybe to be expected?) Checking this using a rank correlation (spearman)
spearman_test_VEGF_IL8 <- cor.test(firefly_sleep_cytokine_merged_completed$VEGF, firefly_sleep_cytokine_merged_completed$IL.8, method = "spearman")
group_level_spearman_VEGF_IL8 <- by(firefly_sleep_cytokine_merged_completed, firefly_sleep_cytokine_merged_completed$dataset, 
                                    function(d) cor.test(d$VEGF, d$IL.8, method = "spearman"))
#So these markers seem to be positively correlated for the NE group, with insufficient evidence to say the same for the controls

standardised_X <- scale(firefly_sleep_cytokine_merged_completed[,35:38])

firefly_sleep_cytokine_merged_completed$EPO_z   <- standardised_X[,1]
firefly_sleep_cytokine_merged_completed$GMCSF_z <- standardised_X[,2]
firefly_sleep_cytokine_merged_completed$VEGF_z  <- standardised_X[,3]
firefly_sleep_cytokine_merged_completed$IL8_z   <- standardised_X[,4]


#Getting the columns in order so that we can compute the subscale totals/CSHQ totals
ord_no_NW1 <- setdiff(ord, "NW1")  

firefly_sleep_cytokine_merged_completed <- firefly_sleep_cytokine_merged_completed[
  , c("ID",
      ord_no_NW1,
      setdiff(names(firefly_sleep_cytokine_merged_completed),
              c("ID", ord_no_NW1)))   
]

firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_BR_total <- rowSums(firefly_sleep_cytokine_merged_completed[,2:5])
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_SOD_total <- firefly_sleep_cytokine_merged_completed[,6]
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_SD_total <- rowSums(firefly_sleep_cytokine_merged_completed[,7:8])
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_SA_total <- rowSums(firefly_sleep_cytokine_merged_completed[,9:12])
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_NW_total <- rowSums(firefly_sleep_cytokine_merged_completed[,13:15])
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_P_total <- rowSums(firefly_sleep_cytokine_merged_completed[,16:22])
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_SDB_total <- rowSums(firefly_sleep_cytokine_merged_completed[,23:25])
firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_DS_total <- rowSums(firefly_sleep_cytokine_merged_completed[,26:33])

firefly_sleep_cytokine_merged_completed$firefly_sleep_cytokine_CSHQ_total <-
  rowSums(
    firefly_sleep_cytokine_merged_completed[, grep("_total$", names(firefly_sleep_cytokine_merged_completed))]
  )



lm_BR_cytokine <- lm(firefly_sleep_cytokine_BR_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_SOD_cytokine <- lm(firefly_sleep_cytokine_SOD_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_SD_cytokine <- lm(firefly_sleep_cytokine_SD_total ~ dataset  + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_SA_cytokine <- lm(firefly_sleep_cytokine_SA_total ~ dataset  + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_NW_cytokine <- lm(firefly_sleep_cytokine_NW_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_P_cytokine <- lm(firefly_sleep_cytokine_P_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_SDB_cytokine <- lm(firefly_sleep_cytokine_SDB_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)
lm_DS_cytokine <- lm(firefly_sleep_cytokine_DS_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)

lm_CSHQ_cytokine <- lm(firefly_sleep_cytokine_CSHQ_total ~ dataset + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = firefly_sleep_cytokine_merged_completed)

#Applying FDR correction for each predictor across all of the subscales
p_vals_EPO <- c(
  summary(lm_BR_cytokine)$coefficients[3,4], 
  summary(lm_SOD_cytokine)$coefficients[3,4], 
  summary(lm_SD_cytokine)$coefficients[3,4],
  summary(lm_SA_cytokine)$coefficients[3,4],
  summary(lm_NW_cytokine)$coefficients[3,4],
  summary(lm_P_cytokine)$coefficients[3,4],
  summary(lm_SDB_cytokine)$coefficients[3,4],
  summary(lm_DS_cytokine)$coefficients[3,4]
)
p_vals_GMCSF <- c(
  summary(lm_BR_cytokine)$coefficients[4,4], 
  summary(lm_SOD_cytokine)$coefficients[4,4], 
  summary(lm_SD_cytokine)$coefficients[4,4],
  summary(lm_SA_cytokine)$coefficients[4,4],
  summary(lm_NW_cytokine)$coefficients[4,4],
  summary(lm_P_cytokine)$coefficients[4,4],
  summary(lm_SDB_cytokine)$coefficients[4,4],
  summary(lm_DS_cytokine)$coefficients[4,4]
)
p_vals_VEGF <- c(
  summary(lm_BR_cytokine)$coefficients[5,4], 
  summary(lm_SOD_cytokine)$coefficients[5,4], 
  summary(lm_SD_cytokine)$coefficients[5,4],
  summary(lm_SA_cytokine)$coefficients[5,4],
  summary(lm_NW_cytokine)$coefficients[5,4],
  summary(lm_P_cytokine)$coefficients[5,4],
  summary(lm_SDB_cytokine)$coefficients[5,4],
  summary(lm_DS_cytokine)$coefficients[5,4]
)
p_vals_IL8 <- c(
  summary(lm_BR_cytokine)$coefficients[6,4], 
  summary(lm_SOD_cytokine)$coefficients[6,4], 
  summary(lm_SD_cytokine)$coefficients[6,4],
  summary(lm_SA_cytokine)$coefficients[6,4],
  summary(lm_NW_cytokine)$coefficients[6,4],
  summary(lm_P_cytokine)$coefficients[6,4],
  summary(lm_SDB_cytokine)$coefficients[6,4],
  summary(lm_DS_cytokine)$coefficients[6,4]
)

adjust_p_fdr_EPO <- rbind(colnames(subscale_total_frame_with_champion_NE)[1:8] , 
                      round(p.adjust(p_vals_EPO, method = 'fdr'),4))
adjust_p_fdr_GMCSF <- rbind(colnames(subscale_total_frame_with_champion_NE)[1:8] , 
                          round(p.adjust(p_vals_GMCSF, method = 'fdr'),4))
adjust_p_fdr_VEGF <- rbind(colnames(subscale_total_frame_with_champion_NE)[1:8], 
                          round(p.adjust(p_vals_VEGF, method = 'fdr'),4))
adjust_p_fdr_IL8 <- rbind(colnames(subscale_total_frame_with_champion_NE)[1:8], 
                          round(p.adjust(p_vals_IL8, method = 'fdr'),4))


# Summary of statistically significant results 
# Following FDR correction for each predictor across subscale totals, there seems to be no significant predictive effect of any of the biomarkers for the subscale totals 
# at the alpha = 0.05 level. EPO seems to be predictive of parasomnias at the alpha = 0.1 level (p = 0.07, beta = 0.78).





#Creating plots for visualising the differences in subscale total and CSHQ total distributions


cohort_levels <- c("firefly_control", "firefly_NE", "champion_NE",
                   "starfish_preterm", "serenity", "champion_CP")
cohort_labels <- c(firefly_control  = "Firefly controls",
                   firefly_NE       = "Firefly NE",
                   champion_NE      = "Champion NE",
                   starfish_preterm = "Starfish (preterm)",
                   serenity         = "Serenity",
                   champion_CP      = "Champion CP")

subscale_labels <- c(BR_total  = "Bedtime resistance",
                     SOD_total = "Sleep onset delay",
                     SD_total  = "Sleep duration",
                     SA_total  = "Sleep anxiety",
                     NW_total  = "Night waking",
                     P_total   = "Parasomnias",
                     SDB_total = "Sleep-disordered breathing",
                     DS_total  = "Daytime sleepiness")


relabel_cohort <- function(x) {
  present <- cohort_levels[cohort_levels %in% unique(as.character(x))]
  factor(x, levels = present, labels = cohort_labels[present])
}

ang_x <- theme(axis.text.x = element_text(angle = 25, hjust = 1))


cshq_df <- CSHQ_total_frame %>%
  mutate(cohort = relabel_cohort(cohort))

kw_cshq_p <- signif(kruskal_CSHQ$p.value, 2)

p_cshq_total <- ggplot(cshq_df, aes(cohort, CSHQ_total, fill = cohort)) +
  geom_violin(alpha = 0.30, colour = NA, width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.5) +
  geom_hline(yintercept = 41, linetype = "dashed", colour = "red") +
  annotate("text", x = Inf, y = 41, label = "Clinical cut-off (41)",
           hjust = 1.05, vjust = -0.6, size = 3, colour = "red") +
  scale_fill_viridis_d(guide = "none") +
  labs(title = "CSHQ total score by cohort",
       subtitle = paste0("Kruskal\u2013Wallis p = ", kw_cshq_p),
       x = NULL, y = "CSHQ total score") +
  theme_minimal(base_size = 12) + ang_x

ggsave("./Eleanor_CSHQ_data/cshq_total_by_cohort.png", p_cshq_total, width = 8.5, height = 5.5, dpi = 600)







subscale_df <- as.data.frame(subscale_total_frame_with_champion_NE,
                             stringsAsFactors = FALSE)
names(subscale_df)[ncol(subscale_df)] <- "cohort"   

sig <- c("SOD_total", "SD_total", "NW_total", "SDB_total")
sub_labs <- c(SOD_total = "Sleep onset delay",
              SD_total  = "Sleep duration",
              NW_total  = "Night waking",
              SDB_total = "Sleep-disordered breathing")

long <- do.call(rbind, lapply(sig, function(s) {
  data.frame(cohort   = subscale_df$cohort,
             subscale = sub_labs[[s]],
             score    = as.numeric(as.character(subscale_df[[s]])),
             stringsAsFactors = FALSE)
}))
long$subscale <- factor(long$subscale, levels = unname(sub_labs))

nice <- c(firefly_control  = "Firefly controls",
          firefly_NE       = "Firefly NE",
          champion_NE      = "Champion NE",
          starfish_preterm = "Starfish (preterm)",
          serenity         = "Serenity")
long$cohort <- factor(nice[long$cohort], levels = nice[names(nice) %in% long$cohort])


p <- ggplot(long, aes(cohort, score, fill = cohort)) +
  geom_violin(alpha = 0.30, colour = NA, trim = FALSE) +              
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +       
  geom_jitter(width = 0.10, height = 0, size = 0.7, alpha = 0.35) +
  facet_wrap(~ subscale, scales = "free_y") +
  scale_fill_viridis_d(guide = "none") +
  labs(title = "CSHQ subscale totals by cohort",
       x = NULL, y = "Subscale total") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

print(p)


ggsave("./Eleanor_CSHQ_data/subscale_violin_plot.png",
       p, width = 7.5, height = 5, dpi = 600)












p_epo <- ggplot(firefly_sleep_cytokine_merged_completed,
                aes(EPO_z, firefly_sleep_cytokine_CSHQ_total)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +
  geom_point(aes(colour = dataset), size = 2, alpha = 0.85) +
  scale_colour_brewer(palette = 'Set1', name = "Cohort", labels = c(firefly_control = "Firefly controls", firefly_NE = "Firefly NE")) +
  labs(title = "EPO vs CSHQ total score",
       subtitle = "p = 0.04, beta = 2.58",
       x = "EPO (standardised)", y = "CSHQ total score") +
  theme_minimal(base_size = 12)

ggsave("./Eleanor_CSHQ_data/epo_vs_cshq_total.png",
       p_epo, width = 7.5, height = 5, dpi = 600)














































#Currently,the only issues with the data seem to be
# 1. value of 4 in one of the items for 1 observation in firefly NE
# 2. multiple items having a value of 4 in observation 17 in starfish
# 3. A few missing entries here and there


# We will proceed with a preliminary analysis, assuming that issue #1. is a typo (should be 3)
# and dropping the row with loads of 4s in the starfish dataset

Y_full <- Y_full[-71,]

#Looking first at the missingness 

M <-1*is.na(Y_full)

col_n_missing <-  colSums(is.na(Y_full))

row_n_missing <- rowSums(is.na(Y_full))



#Classifying observations using the parameters from the old LCA model

appended_theta_CSHQ_varsel <- CSHQ_LCR_varsel$itemprob

appended_theta_CSHQ_varsel[[19]] <- cbind(appended_theta_CSHQ_varsel[[19]], 0)
appended_theta_CSHQ_varsel[[21]] <- cbind(appended_theta_CSHQ_varsel[[21]], 0)
appended_theta_CSHQ_varsel[[24]] <- cbind(appended_theta_CSHQ_varsel[[24]], 0)
appended_theta_CSHQ_varsel[[25]] <- cbind(appended_theta_CSHQ_varsel[[25]], 0)


new_data_LCR_classification <- LCR_posterior_membership_prob(theta = appended_theta_CSHQ_varsel, Y = Y_full, pi = CSHQ_LCR_varsel$pi)

chisq.test(table(factor(max.col(new_data_LCR_classification), levels = 1:4)), p = CSHQ_LCR_varsel$pi)


#testing a 4 gruop LCA model for the new data for comparison
new_data_gibbs_LCA <- blca.gibbs(X = Y_full, G = 4, impute.missing = TRUE, verbose = TRUE)







#Attempting some kind of overfitted mixture to try and appropriately impute missing values etc. - might be able to run the collapsed sampler on the imputed dataset?

#Fitting with 10 groups, with mixture proportion prior being 1/10
overfitted_gibbs_LCA <- blca.gibbs(X = Y_full, G = 10, impute.missing = TRUE, verbose = TRUE, delta = 1/10)

overfitted_gibbs_LCA$classprob
overfitted_gibbs_LCA$Z
max.col(overfitted_gibbs_LCA$Z)
table(max.col(overfitted_gibbs_LCA$Z))

#From this it seems as though the last ~5 groups certainly empty out, with the other groups being potentially present.
#However there is a huge number of observations being assigned to group 1, with not many in the others, so not sure about this solution
#I will use multiple draws from the imputed samples and run several collapsed samplers to see about variable and model selection


run_blca_collapsed_on_imputations <- function(
    Y,
    missing_draws,
    m = 20,
    burn_in_impute = 0,
    thin_impute = 1,
    draw_ids = NULL,
    seed = NULL,
    replace = FALSE,
    verbose = TRUE,
    G = 3,
    formula = NULL,
    ncat = NULL,
    alpha = 1,
    beta = 1,
    delta = 1,
    start.vals = c("single"),
    counts.n = NULL,
    iter = 5000,
    burn.in = 500,
    thin = 1,
    G.sel = TRUE,
    var.sel = TRUE,
    post.hoc.run = TRUE,
    control.post.hoc = list(iter = 2000, burn.in = 500, thin = 1),
    var.prob.thresh = 0.75,
    n.gibbs = nrow(Y),
    only.gibbs = TRUE,
    G.max = 30,
    G.prior = dpois(1:G.max, lambda = 1) / sum(dpois(1:G.max, lambda = 1)),
    prob.inc = 0.5,
    hprior.model = FALSE,
    relabel = TRUE,
    blca_verbose = FALSE,
    verbose.update = 1000
) {
  
  if (!requireNamespace("BayesLCA", quietly = TRUE)) {
    stop("Package 'BayesLCA' is required.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  Y_is_df <- is.data.frame(Y)
  Y_mat <- if (Y_is_df) as.matrix(Y) else as.matrix(Y)
  
  if (length(dim(Y_mat)) != 2) {
    stop("Y must be coercible to a 2D matrix/data.frame.")
  }
  
  nr <- nrow(Y_mat)
  nc <- ncol(Y_mat)
  dn <- dimnames(Y_mat)
  
  mis_idx <- which(is.na(Y_mat), arr.ind = TRUE)
  n_mis <- nrow(mis_idx)
  if (ncol(missing_draws) == n_mis) {
    md <- as.matrix(missing_draws)
  } else if (nrow(missing_draws) == n_mis) {
    md <- t(as.matrix(missing_draws))
  } else {
    stop(sprintf(
      "missing_draws must have either ncol == %d or nrow == %d",
      n_mis, n_mis
    ))
  }
  
  n_samples <- nrow(md)
  
  if (is.null(draw_ids)) {
    if (burn_in_impute >= n_samples) {
      stop("burn_in_impute must be smaller than the number of imputation draws.")
    }
    
    keep <- seq(from = burn_in_impute + 1, to = n_samples, by = thin_impute)
    
    if (length(keep) < m && !replace) {
      stop("Not enough retained imputation draws to sample m draws without replacement.")
    }
    
    draw_ids <- sample(keep, size = m, replace = replace)
  } else {
    draw_ids <- as.integer(draw_ids)
    m <- length(draw_ids)
    
    if (any(draw_ids < 1 | draw_ids > n_samples)) {
      stop("All draw_ids must be between 1 and the number of imputation draws.")
    }
  }
  
  completed_data <- vector("list", m)
  fits <- vector("list", m)
  errors <- vector("list", m)
  
  for (k in seq_len(m)) {
    d <- draw_ids[k]

    Yk <- matrix(Y_mat, nrow = nr, ncol = nc, dimnames = dn)
    Yk[mis_idx] <- md[d, ]

    storage.mode(Yk) <- storage.mode(Y_mat)

    if (Y_is_df) {
      Yk_df <- as.data.frame(Yk, stringsAsFactors = FALSE)
      
      for (j in seq_along(Y)) {
        if (is.factor(Y[[j]])) {
          Yk_df[[j]] <- factor(Yk_df[[j]], levels = levels(Y[[j]]))
        } else if (is.integer(Y[[j]])) {
          Yk_df[[j]] <- as.integer(Yk_df[[j]])
        } else if (is.numeric(Y[[j]])) {
          Yk_df[[j]] <- as.numeric(Yk_df[[j]])
        }
      }
      
      X_in <- Yk_df
    } else {
      X_in <- Yk
    }
    
    completed_data[[k]] <- X_in
    
    if (verbose) {
      message(sprintf(
        "Running blca.collapsed on imputation %d/%d using imputation draw %d",
        k, m, d
      ))
    }
    
    fits[[k]] <- tryCatch(
      BayesLCA::blca.collapsed(
        X = X_in,
        G = G,
        formula = formula,
        ncat = ncat,
        alpha = alpha,
        beta = beta,
        delta = delta,
        start.vals = start.vals,
        counts.n = counts.n,
        iter = iter,
        burn.in = burn.in,
        thin = thin,
        G.sel = G.sel,
        var.sel = var.sel,
        post.hoc.run = post.hoc.run,
        control.post.hoc = control.post.hoc,
        var.prob.thresh = var.prob.thresh,
        n.gibbs = n.gibbs,
        only.gibbs = only.gibbs,
        G.max = G.max,
        G.prior = G.prior,
        prob.inc = prob.inc,
        hprior.model = hprior.model,
        relabel = relabel,
        verbose = blca_verbose,
        verbose.update = verbose.update
      ),
      error = function(e) {
        errors[[k]] <<- list(
          draw_id = d,
          message = conditionMessage(e),
          dim_X = dim(X_in),
          class_X = class(X_in)
        )
        NULL
      }
    )
  }
  
  out <- list(
    draw_ids = draw_ids,
    completed_data = completed_data,
    fits = fits,
    errors = errors,
    mis_idx = mis_idx
  )
  
  class(out) <- "repeated_blca_collapsed"
  out
}


get_G_posterior_matrix <- function(fits) {
  G_tabs <- lapply(fits, function(fit) {
    prop.table(table(fit$samples$G))
  })

  all_G <- sort(unique(as.integer(unlist(lapply(G_tabs, names)))))

  G_mat <- t(sapply(G_tabs, function(tab) {
    out <- setNames(rep(0, length(all_G)), all_G)
    out[names(tab)] <- as.numeric(tab)
    out
  }))
  
  rownames(G_mat) <- paste0("imp_", seq_along(fits))
  colnames(G_mat) <- paste0("G_", all_G)
  
  G_mat
}


get_varprob_matrix_from_mcmc <- function(fits, var_names = NULL, drop_failed = TRUE) {
  if (drop_failed) {
    ok <- !vapply(fits, is.null, logical(1))
    fits_use <- fits[ok]
  } else {
    fits_use <- fits
    ok <- rep(TRUE, length(fits))
  }
  
  if (length(fits_use) == 0) {
    stop("No non-NULL fits found.")
  }
  
  run_probs <- lapply(fits_use, function(fit) {
    vi <- fit$samples$var.ind
    
    if (is.null(vi)) {
      stop("A fit has no samples$var.ind component.")
    }

    vi <- vi[!vapply(vi, is.null, logical(1))]
    vi <- vi[vapply(vi, function(x) length(x) > 0, logical(1))]
    
    if (length(vi) == 0) {
      stop("A fit has no non-empty variable-indicator samples.")
    }
    

    vi_mats <- lapply(vi, function(x) {
      x <- as.matrix(x)
      if (is.null(colnames(x))) {
        stop("var.ind matrix has no column names.")
      }
      x
    })
    

    all_vars <- unique(unlist(lapply(vi_mats, colnames)))
    all_vars <- sort(all_vars)
    
    vi_aligned <- lapply(vi_mats, function(mat) {
      out <- matrix(0, nrow = nrow(mat), ncol = length(all_vars))
      colnames(out) <- all_vars
      out[, colnames(mat)] <- mat
      out
    })
    
    vi_all <- do.call(rbind, vi_aligned)
    colMeans(vi_all)
  })
  

  all_vars <- unique(unlist(lapply(run_probs, names)))
  all_vars <- sort(all_vars)
  
  prob_mat <- t(sapply(run_probs, function(v) {
    out <- setNames(rep(NA_real_, length(all_vars)), all_vars)
    out[names(v)] <- as.numeric(v)
    out
  }))
  
  rownames(prob_mat) <- paste0("imp_", which(ok))
  
  if (!is.null(var_names)) {

    if (all(var_names %in% colnames(prob_mat))) {
      prob_mat <- prob_mat[, var_names, drop = FALSE]
    }
  }
  
  prob_mat
}





missing_draws_num <- apply(overfitted_gibbs_LCA$samples$missing, 2, as.numeric)

blca_collapsed_imputation_run <- run_blca_collapsed_on_imputations(Y = Y_full, 
                                                                   missing_draws = missing_draws_num,
                                                                   m = 100)


G_post_across_runs <- get_G_posterior_matrix(blca_collapsed_imputation_run$fits)

G_post_mean_across_runs <- colMeans(G_post_across_runs)


varsel_mat <- get_varprob_matrix_from_mcmc(blca_collapsed_imputation_run$fits, 
                                           var_names = colnames(Y_full), drop_failed = TRUE)

varsel_mat_mean_across_runs <- colMeans(varsel_mat)


strong_vars <- sort(varsel_mat_mean_across_runs[varsel_mat_mean_across_runs >= 0.90], decreasing = TRUE)
moderate_vars <- sort(varsel_mat_mean_across_runs[varsel_mat_mean_across_runs >= 0.75 & varsel_mat_mean_across_runs < 0.90], decreasing = TRUE)
borderline_vars <- sort(varsel_mat_mean_across_runs[varsel_mat_mean_across_runs >= 0.50 & varsel_mat_mean_across_runs < 0.75], decreasing = TRUE)
weak_vars <- sort(varsel_mat_mean_across_runs[varsel_mat_mean_across_runs < 0.50], decreasing = TRUE)

strong_vars
moderate_vars
borderline_vars
weak_vars


varsel_sd <- apply(varsel_mat, 2, sd)

variable_summary <- data.frame(
  mean_inclusion = colMeans(varsel_mat),
  sd_inclusion = varsel_sd,
  min_inclusion = apply(varsel_mat, 2, min),
  max_inclusion = apply(varsel_mat, 2, max)
)

variable_summary <- variable_summary[order(-variable_summary$mean_inclusion), ]
variable_summary



#Running LCR with item selection on the 4 and 5 group model to see how it performs

#Creating a dummy X to consider the intercept-only model ignoring covariates firstly
X_null <- matrix(nrow = nrow(Y_full), ncol = 0)

LCA_4group_item_sel_new_set <- LCR_Gibbs(X = X_null, Y = Y_full, G = 4, 
                                         beta_prior_cov = diag(10^2,1),
                                         beta_prior_mean = rep(0,1),
                                         theta_hyperparam = 1, 
                                         clust_var_prior = 0.5, 
                                         item.sel = TRUE, 
                                         cov.sel = FALSE, 
                                         verbose = TRUE, 
                                         relabel = TRUE, 
                                         n_samples = 10000, 
                                         burnin = 1000, thinby = 10)

LCA_5group_item_sel_new_set <- LCR_Gibbs(X = X_null, Y = Y_full, G = 5, 
                                         beta_prior_cov = diag(10^2,1),
                                         beta_prior_mean = rep(0,1),
                                         theta_hyperparam = 1, 
                                         clust_var_prior = 0.5, 
                                         item.sel = TRUE, 
                                         cov.sel = FALSE, 
                                         verbose = TRUE, 
                                         relabel = TRUE, 
                                         n_samples = 10000, 
                                         burnin = 1000, thinby = 10)

# Checking the item selection results of the 4 group sampler run, and then compating with the 
# results of the repeated imputation approach

LCA_4group_item_sel_new_set$item_inclusion_prob


sum(1*(varsel_mat_mean_across_runs>0.5) == 1*(LCA_4group_item_sel_new_set$item_inclusion_prob>0.5))

# So we have agreement on selection of variables between the imputed datasets 
# and the collapsed output on the dataset with the missing entries for 31/33 of the items.
# The only disagreement is in the variables P5 and P7, which were borderline anyway.
# In the collapsed run, these had PIP of 0.41, 0.39 respectively
# In the imputation multiple sampler run they had 0.55 and 0.56 respectively.


# We want to possibly run a 5 group model
# as well as clustering point estimation using Wade and Gharamani





#Things to look at now as a preliminary look

# 1. size of the groups
# 2. Clustering point estimation via Wade and Gharahmani - try with collapsed across the imputations
# This might be straightforward enough given that the output is already there.
# 3. summary characteristics - for the moment just looking at subscale totals, full totals, 
# potentially making a profile plot or just mosaic plots





# 1. Proportion of the groups

LCA_4group_item_sel_new_set$pi
LCA_5group_item_sel_new_set$pi

#2. Computing clustering point estimates using Wade and Ghahramani - across all 100 imputation runs

CSHQ_collapsed_cluster_mat <- do.call(rbind,CSHQ_preliminary_collapsed_run$samples$labels)

CSHQ_collapsed_psm_mat <- comp.psm(CSHQ_collapsed_cluster_mat)

CSHQ_minVI_cluster <- minVI(psm = CSHQ_collapsed_psm_mat, method = 'greedy', start.cl = max.col(CSHQ_preliminary_collapsed_run$Z), suppress.comment = FALSE)


valid_runs <- list()

for (run in 1:100) {
  print(run)
  
  labels <- blca_collapsed_imputation_run$fits[[run]]$samples$labels
  
  if (is.null(labels) || !is.list(labels)) {
    message("Skipping run ", run, ": labels is NULL or not a list")
    next
  }
  
  run_mat <- tryCatch(
    do.call(rbind, labels),
    error = function(e) NULL
  )
  
  if (is.null(run_mat)) {
    message("Skipping run ", run, ": could not bind labels")
    next
  }
  
  valid_runs[[length(valid_runs) + 1]] <- run_mat
}

full_cluster_mat <- do.call(rbind, valid_runs)

full_psm_mat <- comp.psm(full_cluster_mat)


# For the moment, I'll just use as the initial clustering point estimate the posterior mode from the 4-group fit

full_minVI_cluster <- minVI(psm = full_psm_mat, method = 'greedy', start.cl = max.col(LCA_4group_item_sel_new_set$Z), suppress.comment = FALSE)

full_minVI_cluster$cl
table(full_minVI_cluster$cl)

# A 4 group solution seems pretty plausible at the moment - although the groups 3 and 4 have only a few observations

#4. Looking at summary properties of the 4 group solution

# Firstly cumulative totals

group_A_Y <- Y_full[which(full_minVI_cluster$cl == 1),]
group_B_Y <- Y_full[which(full_minVI_cluster$cl == 2),]
group_C_Y <- Y_full[which(full_minVI_cluster$cl == 3),]
group_D_Y <- Y_full[which(full_minVI_cluster$cl == 4),]




# GROUP A



# observed score, ignoring missing values
row_score_group_A <- rowSums(group_A_Y, na.rm = TRUE)

# number of missing items in each row
row_missing_group_A <- rowSums(is.na(group_A_Y))

# possible minimum and maximum total, since each item is between 1 and 3
row_min_possible_group_A <- row_score_group_A + row_missing_group_A * 1
row_max_possible_group_A <- row_score_group_A + row_missing_group_A * 3

row_summary_group_A <- data.frame(
  observed_score = row_score_group_A,
  missing_n = row_missing_group_A,
  min_possible = row_min_possible_group_A,
  max_possible = row_max_possible_group_A
)

row_summary_group_A

#Computing the mean of each of the item scores and summing these

mean_summary_score_group_A = sum(colMeans(group_A_Y, na.rm = TRUE))

#getting the mean of the minimum possible scores

mean_min_scores_group_A <- mean(row_summary_group_A$min_possible)

#getting the mean of the max possible scores

mean_max_scores_group_A <- mean(row_summary_group_A$max_possible)



# GROUP B


# observed score, ignoring missing values
row_score_group_B <- rowSums(group_B_Y, na.rm = TRUE)

# number of missing items in each row
row_missing_group_B <- rowSums(is.na(group_B_Y))

# possible minimum and maximum total, since each item is between 1 and 3
row_min_possible_group_B <- row_score_group_B + row_missing_group_B * 1
row_max_possible_group_B <- row_score_group_B + row_missing_group_B * 3

row_summary_group_B <- data.frame(
  observed_score = row_score_group_B,
  missing_n = row_missing_group_B,
  min_possible = row_min_possible_group_B,
  max_possible = row_max_possible_group_B
)

row_summary_group_B

#Computing the mean of each of the item scores and summing these

mean_summary_score_group_B = sum(colMeans(group_B_Y, na.rm = TRUE))

#getting the mean of the minimum possible scores

mean_min_scores_group_B <- mean(row_summary_group_B$min_possible)

#getting the mean of the max possible scores

mean_max_scores_group_B <- mean(row_summary_group_B$max_possible)



# GROUP C


# observed score, ignoring missing values
row_score_group_C <- rowSums(group_C_Y, na.rm = TRUE)

# number of missing items in each row
row_missing_group_C <- rowSums(is.na(group_C_Y))

# possible minimum and maximum total, since each item is between 1 and 3
row_min_possible_group_C <- row_score_group_C + row_missing_group_C * 1
row_max_possible_group_C <- row_score_group_C + row_missing_group_C * 3

row_summary_group_C <- data.frame(
  observed_score = row_score_group_C,
  missing_n = row_missing_group_C,
  min_possible = row_min_possible_group_C,
  max_possible = row_max_possible_group_C
)

row_summary_group_C

#Computing the mean of each of the item scores and summing these

mean_summary_score_group_C = sum(colMeans(group_C_Y, na.rm = TRUE))

#getting the mean of the minimum possible scores

mean_min_scores_group_C <- mean(row_summary_group_C$min_possible)

#getting the mean of the max possible scores

mean_max_scores_group_C <- mean(row_summary_group_C$max_possible)




#GROUP D


# observed score, ignoring missing values
row_score_group_D <- rowSums(group_D_Y, na.rm = TRUE)

# number of missing items in each row
row_missing_group_D <- rowSums(is.na(group_D_Y))

# possible minimum and maximum total, since each item is between 1 and 3
row_min_possible_group_D <- row_score_group_D + row_missing_group_D * 1
row_max_possible_group_D <- row_score_group_D + row_missing_group_D * 3

row_summary_group_D <- data.frame(
  observed_score = row_score_group_D,
  missing_n = row_missing_group_D,
  min_possible = row_min_possible_group_D,
  max_possible = row_max_possible_group_D
)

row_summary_group_D

#Computing the mean of each of the item scores and summing these

mean_summary_score_group_D = sum(colMeans(group_D_Y, na.rm = TRUE))

#getting the mean of the minimum possible scores

mean_min_scores_group_D <- mean(row_summary_group_D$min_possible)

#getting the mean of the max possible scores

mean_max_scores_group_D <- mean(row_summary_group_D$max_possible)


#Creating a table with the 4 groups summarised by cumulative totals

group_summary_table <- data.frame(
  group = c("A", "B", "C", "D"),
  n = c(nrow(group_A_Y), nrow(group_B_Y), nrow(group_C_Y), nrow(group_D_Y)),
  mean_summary_score = c(
    mean_summary_score_group_A,
    mean_summary_score_group_B,
    mean_summary_score_group_C,
    mean_summary_score_group_D
  ),
  mean_observed_score = c(
    mean(row_summary_group_A$observed_score),
    mean(row_summary_group_B$observed_score),
    mean(row_summary_group_C$observed_score),
    mean(row_summary_group_D$observed_score)
  ),
  mean_missing_n = c(
    mean(row_summary_group_A$missing_n),
    mean(row_summary_group_B$missing_n),
    mean(row_summary_group_C$missing_n),
    mean(row_summary_group_D$missing_n)
  ),
  mean_min_score = c(
    mean_min_scores_group_A,
    mean_min_scores_group_B,
    mean_min_scores_group_C,
    mean_min_scores_group_D
  ),
  mean_max_score = c(
    mean_max_scores_group_A,
    mean_max_scores_group_B,
    mean_max_scores_group_C,
    mean_max_scores_group_D
  )
)

group_summary_table




#Now subscale totals - currently working with the full datasets - no item selection


minBR <- 4
minSOD <- 1
minSD <- 3
minSA <- 4
minNW <- 3
minP <- 7
minSDB <- 3
minDS <- 8

Y_full_BR <- Y_full[,1:4]
Y_full_SOD <- Y_full[,5]
Y_full_SD <- Y_full[,6:8]
Y_full_SA <- Y_full[,9:12]
Y_full_NW <- Y_full[,13:15]
Y_full_P <- Y_full[,16:22]
Y_full_SDB <- Y_full[,23:25]
Y_full_DS <- Y_full[,26:33]


# Creating a function for the missing entries tha
# computes the mean of the observed entries on a given subscale, 
# and then imputes this value for the missing ones 
# (provided a certain number of the items are answered). 
# This makes the assumption that on a given subscale, the missing items are answered 
# similarly for an individual as the observed ones.

impute_item_score <- function(mat, min_prop = 0.75) {
  if (is.null(dim(mat))) return(mat)
  n_items <- ncol(mat)
  n_answered <- rowSums(!is.na(mat))
  total <- rowMeans(mat, na.rm = TRUE) * n_items
  ifelse(n_answered / n_items >= min_prop, total, NA)
}




  
# observed score, ignoring missing values
total_BR_observed <- impute_item_score(Y_full_BR)
total_SOD_observed <- Y_full_SOD   
total_SD_observed <- impute_item_score(Y_full_SD)
total_SA_observed <- impute_item_score(Y_full_SA)
total_NW_observed <- impute_item_score(Y_full_NW)
total_P_observed <- impute_item_score(Y_full_P)
total_SDB_observed <- impute_item_score(Y_full_SDB)
total_DS_observed <- impute_item_score(Y_full_DS)

# Relevelled observed score
total_BR_relevelled <- total_BR_observed  - minBR
total_SOD_relevelled <- total_SOD_observed - minSOD
total_SD_relevelled <- total_SD_observed  - minSD
total_SA_relevelled <- total_SA_observed  - minSA
total_NW_relevelled <- total_NW_observed  - minNW
total_P_relevelled <- total_P_observed   - minP
total_SDB_relevelled <- total_SDB_observed - minSDB
total_DS_relevelled <- total_DS_observed  - minDS


# Standardising subscale totals

total_BR_standardised  <- scale(total_BR_relevelled)
total_SOD_standardised <- scale(total_SOD_relevelled)
total_SD_standardised  <- scale(total_SD_relevelled)
total_SA_standardised  <- scale(total_SA_relevelled)
total_NW_standardised  <- scale(total_NW_relevelled)
total_P_standardised   <- scale(total_P_relevelled)
total_SDB_standardised <- scale(total_SDB_relevelled)  
total_DS_standardised  <- scale(total_DS_relevelled)

subscale_total_mat_CSHQ_full <- cbind(
  total_BR_relevelled, total_SOD_relevelled, total_SD_relevelled,
  total_SA_relevelled, total_NW_relevelled, total_P_relevelled,
  total_SDB_relevelled, total_DS_relevelled
)

subscale_total_mat_CSHQ_full_standardised <- cbind(
  total_BR_standardised, total_SOD_standardised, total_SD_standardised,
  total_SA_standardised, total_NW_standardised, total_P_standardised,
  total_SDB_standardised, total_DS_standardised
)


colnames(subscale_total_mat_CSHQ_full) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'SDB', 'DS')
colnames(subscale_total_mat_CSHQ_full_standardised) <- c('BR', 'SOD', 'SD', 'SA', 'NW', 'P', 'SDB', 'DS')



#Now getting the group-specific values

groups <- full_minVI_cluster$cl   

group_mean <- function(mat, g, k) {
  apply(mat[g == k, ], 2, mean, na.rm = TRUE)
}

CSHQ_LCA_profile_mat_full <- rbind(
  group_mean(subscale_total_mat_CSHQ_full, groups, 1),
  group_mean(subscale_total_mat_CSHQ_full, groups, 2),
  group_mean(subscale_total_mat_CSHQ_full, groups, 3),
  group_mean(subscale_total_mat_CSHQ_full, groups, 4)
)

CSHQ_LCA_profile_mat_full_standardised <- rbind(
  group_mean(subscale_total_mat_CSHQ_full_standardised, groups, 1),
  group_mean(subscale_total_mat_CSHQ_full_standardised, groups, 2),
  group_mean(subscale_total_mat_CSHQ_full_standardised, groups, 3),
  group_mean(subscale_total_mat_CSHQ_full_standardised, groups, 4)
)




subscale_names_full <- c("Bedtime Resistance", "Sleep Onset Delay",
                         "Sleep Duration", "Sleep Anxiety", "Night Waking",
                         "Parasomnias", "Sleep Disordered Breathing",
                         "Daytime Sleepiness")

# Desired display order for the plot - consistent with the old plot
subscale_order_plot <- c("Bedtime Resistance", "Daytime Sleepiness",
                         "Night Waking", "Parasomnias", "Sleep Anxiety",
                         "Sleep Duration", "Sleep Onset Delay",
                         "Sleep Disordered Breathing")

CSHQ_LCA_profile_df_full <- tibble(
  Subscale = factor(rep(subscale_names_full, times = 4),
                    levels = subscale_order_plot),
  Score    = c(t(CSHQ_LCA_profile_mat_full)),
  Group    = rep(c('G1', 'G2', 'G3', 'G4'), each = length(subscale_names_full))
)

CSHQ_LCA_profile_df_full_standardised <- tibble(
  Subscale = factor(rep(subscale_names_full, times = 4),
                    levels = subscale_order_plot),
  Score    = c(t(CSHQ_LCA_profile_mat_full_standardised)),
  Group    = rep(c('G1', 'G2', 'G3', 'G4'), each = length(subscale_names_full))
)

CSHQ_LCA_profile_plot_full <- ggplot(
  CSHQ_LCA_profile_df_full,
  aes(x = Subscale, y = Score, group = Group, color = Group)
) +
  geom_line(size = 1.2) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Sleep Profile by Latent Class",
    subtitle = "Mean Relevelled Subscale Scores by Group",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(size = 22, face = "bold"),
    plot.subtitle   = element_text(size = 18),
    legend.text     = element_text(size = 14)
  )



CSHQ_LCA_profile_plot_full_standardised <- ggplot(
  CSHQ_LCA_profile_df_full_standardised,
  aes(x = Subscale, y = Score, group = Group, color = Group)
) +
  geom_line(size = 1.2) +
  geom_point(size = 2.4) +
  theme_minimal(base_size = 18) +
  labs(
    title = "Sleep Profile by Latent Class ",
    subtitle = "Mean Relevelled Subscale Scores by Group (standardised)",
    x = "CSHQ Subscale",
    y = "Mean Score"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(
    legend.title    = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.title      = element_text(size = 22, face = "bold"),
    plot.subtitle   = element_text(size = 18),
    legend.text     = element_text(size = 14)
  )



# Comparing the old model with the new one by looking at modal classification on the new dataset


table(max.col(new_data_LCR_classification), full_minVI_cluster$cl)

#We should take a look at the individual subscales maybe for the smaller two groups in the old and the new model
#Also want to look at the correspondence between the clinical cohorts and the latent groups

#First looking at the correspondence between the clinical cohorts and the latent groups

table(full_minVI_cluster$cl, data_full$dataset)

chisq.test(x = full_minVI_cluster$cl, y = data_full$dataset, simulate.p.value = TRUE, B = 10000)

#No significant association found here anyway

#Looking at the heatmaps for the two smaller groups

#Should try to sort each group by total as well maybe
heatmap_order_groupA <- order(row_score_group_A) 
pheatmap(Y_full[which(full_minVI_cluster$cl == 1),][heatmap_order_groupA,], 
         cluster_cols = FALSE, cluster_rows = FALSE, color = viridis(3), 
         main = 'Heatmap of Latent Group 1 Responses', labels_row ='')

heatmap_order_groupB <- order(row_score_group_B)
pheatmap(Y_full[which(full_minVI_cluster$cl == 2),][heatmap_order_groupB,],
         cluster_cols = FALSE, cluster_rows = FALSE, color = viridis(3), 
         main = 'Heatmap of Latent Group 2 Responses', labels_row ='')

heatmap_order_groupC <- order(row_score_group_C)
pheatmap(Y_full[which(full_minVI_cluster$cl == 3),][heatmap_order_groupC,],
         cluster_cols = FALSE, cluster_rows = FALSE, color = viridis(3), 
         main = 'Heatmap of Latent Group 3 Responses', labels_row ='')

heatmap_order_groupD <- order(row_score_group_D)
pheatmap(Y_full[which(full_minVI_cluster$cl == 4),][heatmap_order_groupD,],
         cluster_cols = FALSE, cluster_rows = FALSE, color = viridis(3), 
         main = 'Heatmap of Latent Group 4 Responses', labels_row ='')



#Looking at all of the groups in a combined heatmap

Y_sorted_by_group <- rbind(
  Y_full[which(full_minVI_cluster$cl == 1),][heatmap_order_groupA,],
  Y_full[which(full_minVI_cluster$cl == 2),][heatmap_order_groupB,],
  Y_full[which(full_minVI_cluster$cl == 3),][heatmap_order_groupC,],
  Y_full[which(full_minVI_cluster$cl == 4),][heatmap_order_groupD,]
)

n_A <- sum(full_minVI_cluster$cl == 1)
n_B <- sum(full_minVI_cluster$cl == 2)
n_C <- sum(full_minVI_cluster$cl == 3)
n_D <- sum(full_minVI_cluster$cl == 4)

row_gaps <- rep(cumsum(c(n_A, n_B, n_C)), each = 4)

group_anno <- data.frame(
  Group = factor(rep(c("A", "B", "C", "D"), times = c(n_A, n_B, n_C, n_D)))
)
rownames(group_anno) <- rownames(Y_sorted_by_group)

pheatmap(Y_sorted_by_group,
         color           = viridis(3),
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         gaps_row        = row_gaps,
         annotation_row  = group_anno,
         show_rownames   = FALSE,
         main = 'Heatmap of Responses Sorted by Latent Group')


#Comparing with the solution with Olives dataset to see how the heatmaps compare

Olive_heatmap_order_groupA <- order(rowSums(Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 1),]))

Olive_heatmap_order_groupB <- order(rowSums(Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 2),]))

Olive_heatmap_order_groupC <- order(rowSums(Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 3),]))

Olive_heatmap_order_groupD <- order(rowSums(Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 4),]))

Olive_Y_sorted_by_group <- rbind(
  Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 1),][Olive_heatmap_order_groupA,],
  Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 2),][Olive_heatmap_order_groupB,],
  Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 3),][Olive_heatmap_order_groupC,],
  Y_CSHQ[which(CSHQ_varsel_minVI_cluster$cl == 4),][Olive_heatmap_order_groupD,]
)

rownames(Olive_Y_sorted_by_group) <- paste0("obs_", seq_len(nrow(Olive_Y_sorted_by_group)))
rownames(group_anno_Olive)        <- rownames(Olive_Y_sorted_by_group)

Olive_n_A <- sum(CSHQ_varsel_minVI_cluster$cl == 1)
Olive_n_B <- sum(CSHQ_varsel_minVI_cluster$cl == 2)
Olive_n_C <- sum(CSHQ_varsel_minVI_cluster$cl == 3)
Olive_n_D <- sum(CSHQ_varsel_minVI_cluster$cl == 4)

row_gaps_Olive <- rep(cumsum(c(Olive_n_A, Olive_n_B, Olive_n_C)), each = 4)

group_anno_Olive <- data.frame(
  Group = factor(rep(c("A", "B", "C", "D"), times = c(Olive_n_A, Olive_n_B, Olive_n_C, Olive_n_D)))
)
rownames(group_anno_Olive) <- rownames(Olive_Y_sorted_by_group)

pheatmap(Olive_Y_sorted_by_group,
         color           = viridis(3),
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         gaps_row        = row_gaps_Olive,
         annotation_row  = group_anno_Olive,
         show_rownames   = FALSE,
         main = 'Heatmap of Responses Sorted by Latent Group (Olive\'s Data)')




#Should try to make a heatmap of each group with the old one below it for comparison maybe













# Looking at an LCA model where we consider Olive's and Eleanor's data combined

# We'll do as before and fit an overfitted mixture on the combined dataset for imputation

# Then we carry out a number of collapsed sampler runs with imputed datasets to get uncertainty 

# for variable and group number selection across imputations


Olive_Y <- Y_CSHQ

Eleanor_Y <- Y_full


Y_pooled <- rbind(Olive_Y, Eleanor_Y)



set.seed(123)
overfitted_gibbs_LCA_pooled <- blca.gibbs(X = Y_pooled, G = 10, impute.missing = TRUE, verbose = TRUE, delta = 1/10)

overfitted_gibbs_LCA_pooled$classprob
overfitted_gibbs_LCA_pooled$Z
max.col(overfitted_gibbs_LCA_pooled$Z)
table(max.col(overfitted_gibbs_LCA_pooled$Z))



missing_draws_num_pooled <- apply(overfitted_gibbs_LCA_pooled$samples$missing, 2, as.numeric)

blca_collapsed_imputation_run_pooled <- run_blca_collapsed_on_imputations(Y = Y_pooled, 
                                                                   missing_draws = missing_draws_num,
                                                                   m = 10)


G_post_across_runs_pooled <- get_G_posterior_matrix(blca_collapsed_imputation_run_pooled$fits)

G_post_mean_across_runs_pooled <- colMeans(G_post_across_runs_pooled)


varsel_mat_pooled <- get_varprob_matrix_from_mcmc(blca_collapsed_imputation_run_pooled$fits, 
                                           var_names = colnames(Y_full), drop_failed = TRUE)

varsel_mat_mean_across_runs_pooled <- colMeans(varsel_mat_pooled)


#Computing point estimates via min VI
valid_runs_pooled <- list()

for (run in 1:10) {
  print(run)
  
  labels <- blca_collapsed_imputation_run_pooled$fits[[run]]$samples$labels
  
  if (is.null(labels) || !is.list(labels)) {
    message("Skipping run ", run, ": labels is NULL or not a list")
    next
  }
  
  run_mat <- tryCatch(
    do.call(rbind, labels),
    error = function(e) NULL
  )
  
  if (is.null(run_mat)) {
    message("Skipping run ", run, ": could not bind labels")
    next
  }
  
  valid_runs_pooled[[length(valid_runs_pooled) + 1]] <- run_mat
}

full_cluster_mat_pooled <- do.call(rbind, valid_runs_pooled)

full_psm_mat_pooled <- comp.psm(full_cluster_mat_pooled)


# For the moment, I'll just use as the initial clustering point estimate the posterior mode from the overfitted mixture

full_minVI_cluster_pooled <- minVI(psm = full_psm_mat_pooled, method = 'greedy', start.cl = max.col(overfitted_gibbs_LCA_pooled$Z), suppress.comment = FALSE)

full_minVI_cluster_pooled$cl
table(full_minVI_cluster_pooled$cl)

# Comparing the classification with the old models (each dataset separately)
# with the pooled model


#Olive's

table(full_minVI_cluster_pooled$cl[1:nrow(Olive_Y)], CSHQ_varsel_minVI_cluster$cl)

#Eleanor's

table(full_minVI_cluster_pooled$cl[(nrow(Olive_Y) + 1):nrow(Y_pooled)], full_minVI_cluster$cl)

ari_olive <- adjustedRandIndex(full_minVI_cluster_pooled$cl[1:nrow(Olive_Y)], 
                               CSHQ_varsel_minVI_cluster$cl)
ari_eleanor <- adjustedRandIndex(full_minVI_cluster_pooled$cl[(nrow(Olive_Y)+1):nrow(Y_pooled)], 
                                 full_minVI_cluster$cl)


barplot(G_post_mean_across_runs_pooled, 
                 names.arg = 4:(3+length(G_post_mean_across_runs_pooled)),
                 xlab = "Number of latent classes", ylab = "Posterior probability",
                main = 'posterior for G across collapsed runs')









#Running a 6-group LCA model to see how the item responses look

pooled_LCA_6group <- blca.gibbs(X = Y_pooled, G = 6, impute.missing = TRUE,verbose = TRUE)




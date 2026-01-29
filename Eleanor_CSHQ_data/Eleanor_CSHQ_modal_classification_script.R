#This script classifies the observations from the firefly control datasets, the Pt firefly datasets 
# and the starfish preterm datasets according to the LCA model fitted to the Olive's ASD dataset 


#Looking at the raw data for the FIREFLY project

firefly_CSHQ_raw_data_CN <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/Firefly_rawscores_CN.csv")

colnames(firefly_CSHQ_raw_data_CN) <- as.character(firefly_CSHQ_raw_data_CN[1, ])

firefly_CSHQ_raw_data_CN <- firefly_CSHQ_raw_data_CN[-1, ]

colnames(firefly_CSHQ_raw_data_CN)[1:5] <- c('ID', 'bedtime', 'waketime', 'night_time_slept (hh:mm)', 'day_naps_time_slept (hh:mm)')

rownames(firefly_CSHQ_raw_data_CN) <- 1:nrow(firefly_CSHQ_raw_data_CN)

firefly_CN_CSHQ_Y <- firefly_CSHQ_raw_data_CN[1:20,c(6:23, 25:39)]

firefly_CN_CSHQ_Y <- matrix(as.numeric(unlist(firefly_CN_CSHQ_Y)), ncol = ncol(firefly_CN_CSHQ_Y))

firefly_CN_CSHQ_Y[,c(7,8)] <- firefly_CN_CSHQ_Y[,c(7,8)] + 1

#We have a number of issues around missingness for this dataset

rowSums(is.na(firefly_CN_CSHQ_Y))

#Observations 12 and 15 have large amounts of missingness currently (25 and 9 cols respectively) We'll drop these for the moment

firefly_CN_CSHQ_Y <- firefly_CN_CSHQ_Y[-c(12,15),]

#We still have a few missing values here, a lot of which are for the same variables. 

#Considering a classification with the fitted model:

firefly_CN_CSHQ_classification_4group <- LCR_posterior_membership_prob(theta = CSHQ_LCR_varsel$itemprob, beta = t(matrix(CSHQ_LCR_varsel$beta_estimate[1,], nrow = 4, ncol = 2)), X = matrix(0, nrow = nrow(firefly_CN_CSHQ_Y), ncol = 1), Y = firefly_CN_CSHQ_Y)





#Looking at the Pt data now

firefly_CSHQ_raw_data_Pt <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/Firefly_rawscores_Pt.csv")

colnames(firefly_CSHQ_raw_data_Pt) <- as.character(firefly_CSHQ_raw_data_Pt[1, ])

firefly_CSHQ_raw_data_Pt <- firefly_CSHQ_raw_data_Pt[-1, ]

colnames(firefly_CSHQ_raw_data_Pt)[1:5] <- c('ID', 'bedtime', 'waketime', 'night_time_slept (hh:mm)', 'day_naps_time_slept (hh:mm)')

rownames(firefly_CSHQ_raw_data_Pt) <- 1:nrow(firefly_CSHQ_raw_data_Pt)

firefly_Pt_CSHQ_Y <- firefly_CSHQ_raw_data_Pt[,c(6:23, 25:39)]

#Rearranging by subscales
firefly_Pt_CSHQ_Y <- firefly_Pt_CSHQ_Y[,c(1,2,9,10,3,13,4,5,11,14,12,15,16,17,18,19:25,26:28,6,29:33,7,8)]

firefly_Pt_CSHQ_Y <- matrix(as.numeric(unlist(firefly_Pt_CSHQ_Y)), ncol = ncol(firefly_Pt_CSHQ_Y))

firefly_Pt_CSHQ_Y[,c(32,33)] <- firefly_Pt_CSHQ_Y[,c(32,33)] + 1

colnames(firefly_Pt_CSHQ_Y) <- c('BR1','BR2','BR3','BR5','SOD1','SD1','SD2','SD3',
                                 'SA1','SA2','SA3','SA4','NW1','NW2','NW3','P1',
                                 'P2','P3','P4','P5','P6','P7','SDB1','SDB2','SDB3',
                                 'DS1','DS2','DS3','DS4','DS5','DS6','DS7','DS8')

rowSums(is.na(firefly_Pt_CSHQ_Y))

#Observations 1 and 3 have 25 missing variables, and 34 has 9. Other than that, there are not more than 1 or 2 for the rest of the observations
#Can remove rows 1,3,34 since these are missing a substantial amount of information

firefly_Pt_CSHQ_Y <- firefly_Pt_CSHQ_Y[-c(1,3,34),]

#We have an entry of 4 for observation 31, item DS8 - the question is do we change this to a 3 or an NA? For the moment I'll change this to a 4
firefly_Pt_CSHQ_Y[31,33] <- 3


#For the moment, just considering classification with the fitted model: 

firefly_Pt_CSHQ_classification_4group <- LCR_posterior_membership_prob(theta = CSHQ_LCR_varsel$itemprob, beta = t(matrix(CSHQ_LCR_varsel$beta_estimate[1,], nrow = 4, ncol = 2)), X = matrix(0, nrow = nrow(firefly_Pt_CSHQ_Y), ncol = 1), Y = firefly_Pt_CSHQ_Y)






#Looking at the classification for the starfish preterm followup data

data_starfish_preterm_followup <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_data_Eleanor_starfish_preterm_followup.csv")
starfish_preterm_followup_CSHQ_subscale_data <- data_starfish_preterm_followup[, c(10:15, 18, 19, 22:46)]

starfish_preterm_followup_CSHQ_Y <- starfish_preterm_followup_CSHQ_subscale_data[,c(1,2,9,10,3,13,4,5,11,14,12,15,16,17,18,22,19,20,21,23,25,24,26,27,28,6,29,30,31,32,33,7,8)]

#starfish_preterm_followup_CSHQ_subscale_data <- starfish_preterm_followup_CSHQ_subscale_data[complete.cases(starfish_preterm_followup_CSHQ_subscale_data),]

#So we now have a dataset of 22 complete CSHQ responses - we can classify by posterior probability of group membership, computed directly from model parameters.

#Observation 17 has values of 4 for some of the items, not sure why. We can either assume that this should be a 3 or NA.
#For the moment we'll consider it a 3

#starfish_preterm_followup_CSHQ_Y <- starfish_preterm_followup_CSHQ_Y[-17,]
starfish_preterm_followup_CSHQ_Y[17, which(starfish_preterm_followup_CSHQ_Y[17,] == 4)] <- 3


starfish_preterm_followup_CSHQ_classification_4group <- LCR_posterior_membership_prob(theta = CSHQ_LCR_varsel$itemprob, beta = t(matrix(CSHQ_LCR_varsel$beta_estimate[1,], nrow = 4, ncol = 2)), X = matrix(0, nrow = nrow(starfish_preterm_followup_CSHQ_Y), ncol = 1), Y = starfish_preterm_followup_CSHQ_Y)








#Looking at the classification of the 3 datasets pooled together by modal classification

firefly_pooled_classification <- rbind(firefly_CN_CSHQ_classification_4group, firefly_Pt_CSHQ_classification_4group)

pooled_classification <- rbind(firefly_pooled_classification, starfish_preterm_followup_CSHQ_classification_4group)

modal_classification_firefly_pooled <- max.col(firefly_pooled_classification)

modal_classification_pooled_dataset <- max.col(pooled_classification)


chisq.test(x = table(modal_classification_firefly_pooled), p = CSHQ_LCR_varsel$pi)

chisq.test(x = table(modal_classification_pooled_dataset), p = CSHQ_LCR_varsel$pi)




source('C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/LCR_posterior_membership_prob.R')
source('C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_classification_from_subscale_totals.R')

#Modifying some of the item probabilities to account for having response 3 for some of the items where only 1 and 2 were observed
# CSHQ_LCR_varsel$itemprob[[19]] <- cbind(CSHQ_LCR_varsel$itemprob[[19]], 0)
# CSHQ_LCR_varsel$itemprob[[21]] <- cbind(CSHQ_LCR_varsel$itemprob[[21]], 0)
# CSHQ_LCR_varsel$itemprob[[24]] <- cbind(CSHQ_LCR_varsel$itemprob[[24]], 0)
# CSHQ_LCR_varsel$itemprob[[25]] <- cbind(CSHQ_LCR_varsel$itemprob[[25]], 0)



data_firefly_controls <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_data_Eleanor_firefly_controls.csv")
View(data_firefly_controls)
firefly_controls_CHSQ_subscale_data <- data_firefly_controls[, 8:16]
#There seems to be a number of issues with this data

#1. There are 9 columns in total, however one of these seems to be 'minutes a night waking usually lasts', and the daytime sleepiness subscale column seems to contain the CSHQ total scores (I think?) 
#for all of the observations, except 1 of them with a score of 9.

#2. subscale 1 (bedtime resistance) has 6 items originally, but since 2 of the items are also contained in sleep anxiety, these are dropped. Here they have not been dropped.

#3. For some reason, subscale 8 (daytime sleepiness) says subscale totals range from 6-22, which can't be correct. Firstly, there are 8 items in this subscale, so the minimum must be 8.
# Even if there were 6 items, each item has a max possible score of 3, so the max score would then be 18, not 22. it couldn't be 22 in any case since 3 does not divide 22. 
# Anyway, the data from this column seems to be missing.

#These issues might make it tricky to attempt classification using the fitted LCR model. I will leave this for the moment and come back to it.


data_firefly_NE <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_data_Eleanor_firefly_NE.csv")
firefly_NE_CHSQ_subscale_data <- data_firefly_NE[, c(8:12, 14:16)]

#This dataset has some similar issues around the subscale total ranges for subscale 1 and subscale 8. I'll need to come back to this as well. 
#I think maybe if we find out the method used to calculate the totals, we can modify the classification functions to deal with this?




data_champion_NE <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_data_Eleanor_champion_NE.csv")
champion_NE_CHSQ_subscale_data <- data_champion_NE[, 3:10]


#This dataset has some small problems, there aresome of the subscale totals that are less than the minimum possible score, 
#e.g. some are 0, and for Parasomnias, sleep disordered breathing and daytime sleepiness, there are a few scores less than the minimum possible value.


#The champion CP scores only contain the total CSHQ scores, so Im not sure how useful that would be in this context

data_starfish_preterm_followup <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_data_Eleanor_starfish_preterm_followup.csv")
starfish_preterm_followup_CSHQ_subscale_data <- data_starfish_preterm_followup[, c(10:15, 18, 19, 22:46)]

starfish_preterm_followup_CSHQ_subscale_data <- starfish_preterm_followup_CSHQ_subscale_data[,c(1,2,9,10,3,13,4,5,11,14,12,15,16,17,18,22,19,20,21,23,25,24,26,27,28,6,29,30,31,32,33,7,8)]

starfish_preterm_followup_CSHQ_subscale_data <- starfish_preterm_followup_CSHQ_subscale_data[complete.cases(starfish_preterm_followup_CSHQ_subscale_data),]

#So we now have a dataset of 22 complete CSHQ responses - we can classify by posterior probability of group membership, computed directly from model parameters.

#Observation 17 has values of 4 for some of the items, not sure why, but I'll drop this row for the moment.

starfish_preterm_followup_CSHQ_subscale_data <- starfish_preterm_followup_CSHQ_subscale_data[-17,]

starfish_preterm_followup_CSHQ_classification_4group <- LCR_posterior_membership_prob(theta = CSHQ_LCR_varsel$itemprob, beta = t(matrix(CSHQ_LCR_varsel$beta_estimate[1,], nrow = 4, ncol = 2)), X = matrix(0, nrow = 21, ncol = 1), Y = starfish_preterm_followup_CSHQ_subscale_data)

#Looking at classification by modal assignment via posterior classification probability

modal_assignment_starfish_preterm_followup_CSHQ_LCR_varsel <- max.col(starfish_preterm_followup_CSHQ_classification_4group)




data_serenity <- read.csv("C:/Users/matth/Documents/AIM CP Project/LCA/BLCR_and_var_sel/Eleanor_CSHQ_data/CSHQ_data_Eleanor_serenity.csv")

#There seems to be issues again her with how the subscale totals are computed...


serenity_subscale_data <- data_serenity[,5:12]

serenity_subscale_data <- serenity_subscale_data[c(8, 11:15, 17:22, 24:26, 28:30),]

colnames(serenity_subscale_data) <- c('BR','SOD','SD','SA','NW','P','SDB','DS')

serenity_subscale_data <- apply(serenity_subscale_data, 2, function(x) as.numeric(as.character(x)))


serenity_subscale_classification_4group_LCA <- classify_from_subscale_matrix_CSHQ(subscale_totals_mat = as.matrix(serenity_subscale_data), 
                                                                                  pi = CSHQ_LCR_varsel$pi, 
                                                                                  item_probs = CSHQ_LCR_varsel$itemprob)








library(pheatmap)
library(viridis)
library(car)      
library(FSA)      
library(mice)
library(ggplot2)
library(GGally)   
library(Hmisc)    
library(dplyr)
library(stringr)



# paths 
data_dir <- "./Eleanor_CSHQ_data"
out_dir  <- "./Eleanor_CSHQ_data"

# item order & subscale definition
ord <- c("BR1","BR2","BR3","BR5",
         "SOD1",
         "SD1","SD2","SD3",
         "SA1","SA2","SA3","SA4",
         "NW1","NW2","NW3",
         "P1","P2","P3","P4","P5","P6","P7",
         "SDB1","SDB2","SDB3",
         "DS1","DS2","DS3","DS4","DS5","DS6","DS7","DS8")

subscale_items <- list(
  BR  = c("BR1","BR2","BR3","BR5"),
  SOD = "SOD1",
  SD  = c("SD1","SD2","SD3"),
  SA  = c("SA1","SA2","SA3","SA4"),
  NW  = c("NW1","NW2","NW3"),
  P   = c("P1","P2","P3","P4","P5","P6","P7"),
  SDB = c("SDB1","SDB2","SDB3"),
  DS  = c("DS1","DS2","DS3","DS4","DS5","DS6","DS7","DS8")
)
subscale_names <- names(subscale_items)








# column names shared by every raw item-level file (same order in all four)
cshq_colnames <- c("ID","BR1","BR2","SOD1","SD2","SD3","DS1","DS7","DS8",
                   "BR3","BR5","SA1","SA3","SD1","SA2","SA4","NW1",
                   "NW2","NW3","P2","P3","P4","P1","P5","P7","P6",
                   "SDB1","SDB2","SDB3","DS2","DS3","DS4","DS5","DS6")


#clean_id <- function(x) gsub("\\s+", "", toupper(trimws(as.character(x))))














clean_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("\\(.*?\\)", "", x)     
  x <- gsub("[[:space:]_-]+", "", x)
  x[x %in% c("", "NA", "NANA", "NULL")] <- NA_character_
  x
}


# ---- 2. cohort-specific recoding --------------------------------------
# This is to match some of the ID coding up between the datasets
harmonise_id <- function(id, cohort) {
  id <- clean_id(id)
  ch <- as.character(cohort)
  
  # champion_NE: CSHQ file uses bare numbers, cytokine + age files use NE-prefix
  i <- ch == "champion_NE" & grepl("^[0-9]+$", id) & !is.na(id)
  id[i] <- paste0("NE", id[i])
  
  # serenity CP: THU is a transposition of TUH (TUH is correct)
  i <- ch %in% c("serenity", "serenity_cp") & !is.na(id)
  id[i] <- sub("^THU", "TUH", id[i])
  
  # serenity controls: age file uses CJ01..CJ14, cytokine file uses C1..C14
  i <- ch == "serenity_control" & grepl("^CJ[0-9]+$", id) & !is.na(id)
  id[i] <- sub("^CJ0*", "C", id[i])
  
  id
}






# ---- 4. report genuine within-file duplicates -------------------------

show_dups <- function(d, nm, value_cols) {
  k <- paste(d$cohort, d$ID_clean, sep = " | ")
  dk <- unique(k[duplicated(k)])
  if (!length(dk)) { message(nm, ": no duplicates"); return(invisible(NULL)) }
  message(nm, ": ", length(dk), " duplicated key(s)")
  print(d[k %in% dk, c("cohort", "ID_clean", value_cols)])
  invisible(dk)
}



# ---- 5. collapse duplicates (run only after inspecting step 4) --------

mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

collapse_dups <- function(d, value_cols) {
  d %>%
    group_by(cohort, ID_clean) %>%
    summarise(across(all_of(value_cols), mean_na), .groups = "drop") %>%
    as.data.frame()
}

# scores_all  <- collapse_dups(scores_all, c(subscale_names, "CSHQ_total"))
# age_df      <- collapse_dups(age_df,     "age_months")
# cyto_all    <- collapse_dups(cyto_all,   marker_raw)



cohort_meta <- data.frame(
  cohort = c("firefly_control", "firefly_NE",
             "champion_control", "champion_NE", "champion_CP",
             "serenity_control", "serenity_cp",
             "starfish_control", "starfish_preterm"),
  study  = c("firefly", "firefly",
             "champion", "champion", "champion",
             "serenity", "serenity",
             "starfish", "starfish"),
  group  = c("control", "NE",
             "control", "NE", "CP",
             "control", "CP",
             "control", "preterm"),
  stringsAsFactors = FALSE
)


# ---- 7. post-join check
link_check <- function(data = master) {
  data %>%
    mutate(cshq = !is.na(CSHQ_total),
           age  = !is.na(age_months),
           cyt  = !is.na(EPO) | !is.na(EPO_LPS)) %>%
    group_by(cohort) %>%
    summarise(n = n(),
              cshq = sum(cshq), age = sum(age), cyt = sum(cyt),
              all_three = sum(cshq & age & cyt), .groups = "drop") %>%
    as.data.frame()
}

























load_cshq_cohort <- function(path, col_idx, dataset, row_idx = NULL) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  if (is.null(row_idx)) row_idx <- seq_len(nrow(raw))
  d <- raw[row_idx, col_idx]
  colnames(d) <- cshq_colnames
  d[, -1] <- lapply(d[, -1], function(x) as.numeric(trimws(as.character(x))))
  d[d == ""] <- NA
  d[, c("DS7", "DS8")] <- d[, c("DS7", "DS8")] + 1   # scoring offset for these 2 items
  d <- d[, c("ID", ord)]
  d$dataset <- dataset
  d
}

# Subscale totals from an item-level frame. Items not present (e.g. NW1 once
# dropped) are ignored, so NW automatically becomes NW2 + NW3.
subscale_totals <- function(items_df) {
  out <- sapply(subscale_items, function(cols) {
    cols <- intersect(cols, colnames(items_df))
    rowSums(items_df[, cols, drop = FALSE])
  })
  as.data.frame(out)
}

# Min/max proportion above a CSHQ cut-off given missingness (items score 1-3).
prop_above_threshold <- function(Y, thresh = 41) {
  score   <- rowSums(Y, na.rm = TRUE)
  missing <- rowSums(is.na(Y))
  c(min = mean((score + missing * 1) > thresh),
    max = mean((score + missing * 3) > thresh))
}


#load item-level cohorts
firefly_control <- load_cshq_cohort(file.path(data_dir, "Firefly_rawscores_CN.csv"),
                                    col_idx = c(1, 6:23, 25:39),
                                    dataset = "firefly_control", row_idx = 2:21)
firefly_NE <- load_cshq_cohort(file.path(data_dir, "Firefly_rawscores_Pt.csv"),
                                    col_idx = c(1, 6:23, 25:39),
                                    dataset = "firefly_NE", row_idx = 2:35)
serenity <- load_cshq_cohort(file.path(data_dir, "SERENITY_raw_data.csv"),
                                    col_idx = c(1, 6:23, 25:39),
                                    dataset = "serenity", row_idx = 3:21)
starfish <- load_cshq_cohort(file.path(data_dir, "starfish_preterm_followup_raw_CSHQ.csv"),
                                    col_idx = c(1, 10:15, 18, 19, 22:46),
                                    dataset = "starfish_preterm")

# per-cohort fixes

firefly_NE$DS8[which(firefly_NE$DS8 == 4)] <- 3    # stray 4 (items are 1-3) -> 3
starfish <- starfish[starfish$ID != "S17", ]       # S17 had out-of-range responses

# combine 
data_full <- rbind(firefly_control, firefly_NE, starfish, serenity)
Y_full <- data_full[, ord]

keep <- rowSums(is.na(Y_full)) != ncol(Y_full)     # drop entirely-NA rows
data_full <- data_full[keep, , drop = FALSE]
Y_full <- Y_full[keep, , drop = FALSE]
rownames(data_full) <- rownames(Y_full) <- NULL


# ---- age distributions -----------------------------------------------
# These totals files are only read for ages.
firefly_control_tot <- read.csv(file.path(data_dir, "CSHQ_data_Eleanor_firefly_controls.csv"))
firefly_NE_tot      <- read.csv(file.path(data_dir, "CSHQ_data_Eleanor_firefly_NE.csv"))
serenity_tot        <- read.csv(file.path(data_dir, "CSHQ_data_Eleanor_serenity.csv"))
starfish_raw        <- read.csv(file.path(data_dir, "starfish_preterm_followup_raw_CSHQ.csv"))
champion_CP_raw <- read.csv(file.path(data_dir, "CSHQ_data_Eleanor_champion_CP.csv"))

## CHECK: serenity age units unconfirmed. serenity_age_scale = 1 treats the
##        column as months; set to 0.230137 if it is actually in weeks.
serenity_age_scale <- 1

ages <- list(
  "Firefly controls" = list(x = as.numeric(firefly_control_tot$X.1[3:22]), col = "lightblue"),
  "Firefly NE"       = list(x = as.numeric(firefly_NE_tot$X.1[3:31]),      col = "lightgreen"),
  "Starfish"         = list(x = starfish_raw$AGE..years. * 12,             col = "orange"),
  "Serenity CP"         = list(x = as.numeric(serenity_tot$X.1[3:30]) * serenity_age_scale,
                               col = "purple"),
  "Champion CP"         = list(x = as.numeric(gsub(" YRS", "", champion_CP_raw$X.1[3:29])) * 12,
                               col = "red")
)

draw_age_hists <- function(ages) {
  xr <- range(unlist(lapply(ages, `[[`, "x")), na.rm = TRUE)
  par(mfrow = c(2, 2))
  for (nm in names(ages))
    hist(ages[[nm]]$x, main = paste("Age:", nm), xlab = "Age (months)",
         ylab = "Density", col = ages[[nm]]$col, xlim = xr, probability = TRUE)
}

# Draw once per device. (The original nested pdf() then png(), which sent all
# plotting to the PNG and left the PDF blank.)
png(file.path(out_dir, "age_hist.png"), width = 1000, height = 800); draw_age_hists(ages); dev.off()
pdf(file.path(out_dir, "age_hist.pdf"), width = 10, height = 8);     draw_age_hists(ages); dev.off()


#Creating age dfs for use later on in the regression models

firefly_control_age <- data.frame(
  cohort = "firefly_control",
  ID = firefly_control$ID,
  age_months = as.numeric(firefly_control_tot$X.1[3:22])
)

# We have some issues with inconsistent coding to attach age with the firefly NE group
# firefly_NE_age <- data.frame(
#   cohort = "firefly_NE",
#   ID = as.character(firefly_NE$ID),
#   age_months = as.numeric(firefly_NE_tot$X.1[3:31])
# )

starfish_age <- data.frame(
  cohort = "starfish_preterm",
  ID = toupper(as.character(starfish_raw$STARFISH.ID)),
  age_months = as.numeric(starfish_raw$AGE..years.) * 12
)


serenity_cp_age_csv_with_ID <- read.csv(file.path(data_dir, "SERENITY_raw_data.csv"))
serenity_cp_age_csv_with_ID$Nimbus.ID <- gsub("^THU", "TUH", serenity_cp_age_csv_with_ID$Nimbus.ID)
colnames(serenity_cp_age_csv_with_ID)[1] <- 'ID'
serenity_cp_tot_for_age_merge <- serenity_tot[3:30,c(1,3)]
colnames(serenity_cp_tot_for_age_merge) <- c('ID', 'age')
serenity_cp_age_csv_with_ID <- inner_join(serenity_cp_age_csv_with_ID, serenity_cp_tot_for_age_merge, by = 'ID')
serenity_cp_age <- data.frame(
  cohort = "serenity_cp",
  ID = as.character(serenity_cp_age_csv_with_ID$ID),
  age_months = as.numeric(serenity_cp_age_csv_with_ID$age) * serenity_age_scale
)

champion_cp_age <- data.frame(
  cohort = "champion_CP",
  ID = as.character(champion_CP_raw$X[3:21]),
  age_months = as.numeric(
    gsub(" YRS", "", champion_CP_raw$X.1[3:21])
  ) * 12
)

champion_ne_age_frame <- read.csv(file.path(data_dir, "champion_NE_age_data.csv"), header = FALSE)[,-2]
colnames(champion_ne_age_frame) <- c('ID', 'age')
champion_ne_age_frame$ID <- as.integer(gsub(".*?(\\d+).*", "\\1", champion_ne_age_frame$ID))
champion_ne_age_frame$age <- as.numeric(gsub("[^0-9.]", "", champion_ne_age_frame$age))
champion_ne_age <- data.frame(
  cohort = "champion_NE",
  ID = paste0("NE", as.character(champion_ne_age_frame$ID)),   
  age_months = as.numeric(champion_ne_age_frame$age) * 12
)

champion_control_age_frame <- read.csv(file.path(data_dir, "champion_control_age_data.csv"), header = FALSE)[,-c(2,4)]
colnames(champion_control_age_frame) <- c('ID', 'age')
champion_control_age <- data.frame(
  cohort = "champion_control",
  ID = as.character(champion_control_age_frame$ID),
  age_months = as.numeric(champion_control_age_frame$age)*12
)

serenity_control_age_frame <- read.csv(file.path(data_dir, "serenity_control_age_data.csv"))[,c(1,2)]
colnames(serenity_control_age_frame) <- c('ID', 'age')
serenity_control_age <- data.frame(
  cohort = "serenity_control",
  ID = as.character(serenity_control_age_frame$ID),
  age_months = as.numeric(serenity_control_age_frame$age)*12
)

starfish_control_age_frame <- read.csv(file.path(data_dir, "starfish_control_age_data.csv"))
colnames(starfish_control_age_frame) <- c('ID','age')
starfish_control_age_frame <- starfish_control_age_frame %>%
  mutate(
    age_months = as.numeric(str_extract(age, "\\d+(?=y)")) * 12 +
      as.numeric(str_extract(age, "(?<=y)\\d+(?=m)"))
  )
starfish_control_age <- data.frame(
  cohort = "starfish_control",
  ID = as.character(starfish_control_age_frame$ID),
  age_months = as.numeric(starfish_control_age_frame$age_months)
)


age_df <- bind_rows(
  firefly_control_age,
  # firefly_NE_age,
  starfish_age,
  serenity_cp_age,
  champion_cp_age,
  champion_ne_age,
  champion_control_age,
  serenity_control_age,
  starfish_control_age
)

age_df$ID_clean <- gsub("\\s+", "", toupper(trimws(age_df$ID)))






# ---- missingness + response heatmaps ---------------------------------
gaps_col <- c(4, 5, 8, 12, 15, 22, 25, 33)   # subscale boundaries within `ord`
gaps_row <- c(20, 54, 75) # cohort boundaries

# missingness (TRUE = missing). Simplified to is.na(): equivalent intent to the
# original misty::na.indicator()[,34:66] slice, but dimension-safe, so misty is
# no longer needed.
png(file.path(out_dir, "missingness_heatmap.png"), width = 1000, height = 800)
pheatmap(is.na(Y_full) * 1,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = c("yellow", "red"),          # observed = yellow, missing = red
         gaps_row = gaps_row, gaps_col = gaps_col,
         labels_row = rep("", nrow(Y_full)), labels_col = colnames(Y_full),
         main = "Missing Indicator Heatmap")
dev.off()

png(file.path(out_dir, "CSHQ_response_heatmap.png"), width = 1000, height = 800)
pheatmap(Y_full, cluster_rows = FALSE, cluster_cols = FALSE,
         color = viridis(3), gaps_row = gaps_row, gaps_col = gaps_col,
         labels_row = rep("", nrow(Y_full)), labels_col = colnames(Y_full),
         main = "CSHQ Response Heatmap")
dev.off()


# ---- cohort mean item-score heatmap ----------------------------------
cohort_item_mean <- t(sapply(split(as.data.frame(Y_full), data_full$dataset),
                             function(d) colMeans(d, na.rm = TRUE)))
cohort_item_mean <- cohort_item_mean[c("firefly_control", "firefly_NE",
                                       "starfish_preterm", "serenity"), ]
cohort_item_mean["firefly_control", "NW1"] <- NA   # control cohort almost all missing NW1

png(file.path(out_dir, "CSHQ_cohort_mean_response_heatmap.png"),
    width = 14, height = 3, units = "in", res = 600)
pheatmap(cohort_item_mean, display_numbers = round(cohort_item_mean, 2),
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("green", "white", "red"))(100),
         main = "Mean CSHQ item scores by cohort",
         fontsize_number = 9, fontcolor_number = "gray30",
         angle_col = 45, border_color = "grey60",
         cellwidth = 25, cellheight = 35,
         labels_row = c("firefly controls", "firefly NE", "starfish", "serenity"))
dev.off()


# ---- proportion above clinical cut-off (41) --------------------------
keep_clean      <- rowSums(is.na(Y_full)) <= 16    # drop severe missingness
data_full_clean <- data_full[keep_clean, , drop = FALSE]
Y_full_clean    <- Y_full[keep_clean, , drop = FALSE]

prop_above_threshold(Y_full_clean)                                   # overall
sapply(split(Y_full_clean, data_full_clean$dataset),                 # per cohort
       prop_above_threshold)


# ---- champion cohorts (subscale / total only) ------------------------
# champion_NE: per-subscale totals provided.
champion_NE_raw <- read.csv(file.path(data_dir, "CSHQ_data_Eleanor_champion_NE.csv"))
champ_ne <- champion_NE_raw[3:55, c(1, 3:10)]
colnames(champ_ne) <- c("ID", "BR", "SOD", "SD", "SA", "NW", "P", "SDB", "DS")
champ_ne <- champ_ne[!is.na(champ_ne$BR), ]          # drop empty rows
champ_ne$DS <- champ_ne$DS + 2                        # DS scoring offset (2 items)
champ_ne$SD[champ_ne$SD == 0] <- NA                   # stray 0 -> treat as missing
champ_ne <- complete(mice(champ_ne, method = "pmm", seed = 123, printFlag = FALSE))
# Remove contributions the other cohorts don't have, assuming equal item weights:
champ_ne$NW <- champ_ne$NW - champ_ne$NW / 3          # drop NW1-equivalent (NW has 3 items)
champ_ne$BR <- champ_ne$BR - 2 * champ_ne$BR / 6      # drop 2 items shared with SA (BR has 6)

# champion_CP: only a CSHQ total provided.
champion_CP_raw <- read.csv(file.path(data_dir, "CSHQ_data_Eleanor_champion_CP.csv"))
champ_cp <- champion_CP_raw[3:21, c(1, 4)]
colnames(champ_cp) <- c("ID", "total")
champ_cp$total <- as.numeric(champ_cp$total) + 2       # scoring offset (2 items)
champ_cp$total <- champ_cp$total - champ_cp$total / 33  # drop NW1-equivalent (33 items total)

# champion_control: only a CSHQ total provided (for most of the observations)
champion_control_raw <- read.csv(file.path(data_dir, "champion_control_CSHQ_data.csv"))
champ_control <- champion_control_raw[, c(1, 5)]
colnames(champ_control) <- c("ID", "total")
champ_control$total <- as.numeric(champ_control$total) + 2       # scoring offset (2 items)
champ_control$total <- champ_control$total - champ_control$total / 33  # drop NW1-equivalent (33 items total)
champ_control$ID <- paste0("CON", champ_control$ID)





# ---- single imputation -> master per-participant scores table --------
# Impute once on the clean item set (NW1 dropped to match the other cohorts),
# using cohort as an auxiliary variable, then reuse these scores everywhere.
items_no_nw1 <- setdiff(ord, "NW1")
imp       <- mice(cbind(Y_full_clean[, items_no_nw1], dataset = data_full_clean$dataset),
                  m = 5, method = "pmm", seed = 123, printFlag = FALSE)
items_imp <- complete(imp)[, items_no_nw1]

# item-level cohorts
sc_items <- data.frame(
  ID     = as.character(data_full_clean$ID),
  cohort = data_full_clean$dataset,
  subscale_totals(items_imp),
  CSHQ_total = rowSums(items_imp),
  stringsAsFactors = FALSE
)

# champion_NE (subscales provided -> CSHQ total is their sum)
sc_champ_ne <- data.frame(
  ID = as.character(champ_ne$ID), cohort = "champion_NE",
  champ_ne[, subscale_names],
  CSHQ_total = rowSums(champ_ne[, subscale_names]),
  stringsAsFactors = FALSE
)

# champion_CP (total only -> subscales NA)
sc_champ_cp <- data.frame(
  ID = as.character(champ_cp$ID), cohort = "champion_CP",
  setNames(as.list(rep(NA_real_, length(subscale_names))), subscale_names),
  CSHQ_total = champ_cp$total,
  stringsAsFactors = FALSE
)

# champion_control (total only -> subscales NA)
sc_champ_control <- data.frame(
  ID = as.character(champ_control$ID), cohort = "champion_control",
  setNames(as.list(rep(NA_real_, length(subscale_names))), subscale_names),
  CSHQ_total = champ_control$total,
  stringsAsFactors = FALSE
)

scores_all <- bind_rows(sc_items, sc_champ_ne, sc_champ_cp, sc_champ_control)

scores_all$cohort[scores_all$cohort == "serenity"] <- "serenity_cp"
scores_all$ID <- harmonise_id(scores_all$ID, scores_all$cohort)

scores_all$study    <- unname(study_of[scores_all$cohort])
scores_all$ID_clean <- clean_id(scores_all$ID)

# coarse `study` grouping (used to join cytokines safely - see cytokine section)
study_of <- c(firefly_control = "firefly", firefly_NE = "firefly",
              starfish_preterm = "starfish", serenity_cp = "serenity",
              serenity_control = "serenity",champion_NE = "champion", 
              champion_CP = "champion", champion_control = "champion")
scores_all$study <- unname(study_of[scores_all$cohort])
#Converting the THUXX codes to TUHXX for consistency
scores_all$ID[scores_all$cohort == 'serenity_cp'] <- sub("^THU", "TUH", scores_all$ID[scores_all$cohort == 'serenity_cp'])




# ---- cohort comparisons ----------------------------------------------
sub_df  <- scores_all[!is.na(scores_all$BR), ]   # subscale tests (excludes champion_CP)
cshq_df <- scores_all                            # CSHQ-total test (all cohorts)





# 1. variance check (Levene) across subscales
levene_p <- sapply(subscale_names, function(s)
  leveneTest(sub_df[[s]] ~ factor(sub_df$cohort))[1, "Pr(>F)"])

# 2. ANOVA residual normality was checked (Shapiro/QQ) and rejected, so
#    Kruskal-Wallis is used throughout. Quick optional re-check:
# for (s in c(subscale_names, "CSHQ_total")) {
#   r <- residuals(aov(reformulate("cohort", s),
#                      data = if (s == "CSHQ_total") cshq_df else sub_df))
#   print(c(subscale = s, shapiro_p = signif(shapiro.test(r)$p.value, 3)))
# }

# 3. Kruskal-Wallis
kw_subscale <- sapply(subscale_names, function(s)
  kruskal.test(sub_df[[s]] ~ factor(sub_df$cohort))$p.value)
kw_cshq <- kruskal.test(CSHQ_total ~ factor(cohort), data = cshq_df)$p.value

# 4. FDR across the 8 subscales (CSHQ total tested separately)
kw_subscale_fdr <- p.adjust(kw_subscale, method = "fdr")
round(kw_subscale_fdr, 4)
round(kw_cshq, 4)

# 5. Dunn post-hoc for subscales significant after FDR
#    (SOD, SD, NW, SDB in the original run - update to match kw_subscale_fdr)
dunn_sig <- c("SOD", "SD", "NW", "SDB")
dunn_res <- lapply(setNames(dunn_sig, dunn_sig), function(s)
  dunnTest(sub_df[[s]] ~ factor(sub_df$cohort), method = "bh")$res)
dunn_res


# ---- cohort mean subscale (per-item) heatmap -------------------------
# mean subscale total per cohort / number of items summed -> per-item score
n_items   <- sapply(subscale_items, function(cols) length(intersect(cols, items_no_nw1)))
sub_means <- t(sapply(
  split(scores_all[!is.na(scores_all$BR), subscale_names],
        scores_all$cohort[!is.na(scores_all$BR)]),
  colMeans))
sub_means <- sweep(sub_means, 2, n_items, "/")
row_order <- intersect(c("firefly_control", "firefly_NE", "starfish_preterm",
                         "serenity", "champion_NE"), rownames(sub_means))
sub_means <- sub_means[row_order, ]

png(file.path(out_dir, "CSHQ_cohort_mean_subtotal_heatmap.png"),
    width = 5, height = 4, units = "in", res = 600)
pheatmap(sub_means, display_numbers = round(sub_means, 2),
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("green", "white", "red"))(100),
         main = "Mean CSHQ per-item score by cohort",
         fontsize_number = 9, fontcolor_number = "gray30",
         angle_col = 45, border_color = "grey60",
         cellwidth = 25, cellheight = 35,
         labels_row = c("firefly controls", "firefly NE", "starfish",
                        "serenity", "champion NE")[seq_along(row_order)])
dev.off()


# ---- plots: CSHQ total + significant subscales -----------------------
cohort_levels <- c("firefly_control", "firefly_NE", "champion_NE",
                   "starfish_preterm", "serenity", "champion_CP",
                   "champion_control")
cohort_labels <- c(firefly_control = "Firefly controls", firefly_NE = "Firefly NE",
                   champion_NE = "Champion NE", starfish_preterm = "Starfish (preterm)",
                   serenity = "Serenity", champion_CP = "Champion CP", champion_control = "Champion controls")
relabel_cohort <- function(x) {
  present <- cohort_levels[cohort_levels %in% unique(x)]
  factor(x, levels = present, labels = cohort_labels[present])
}
ang_x <- theme(axis.text.x = element_text(angle = 25, hjust = 1))

# CSHQ total
cshq_plot_df <- cshq_df
cshq_plot_df$cohort <- relabel_cohort(cshq_plot_df$cohort)
p_cshq <- ggplot(cshq_plot_df, aes(cohort, CSHQ_total, fill = cohort)) +
  geom_violin(alpha = 0.30, colour = NA, width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.08, height = 0, size = 1.1, alpha = 0.5) +
  geom_hline(yintercept = 41, linetype = "dashed", colour = "red") +
  annotate("text", x = Inf, y = 41, label = "Clinical cut-off (41)",
           hjust = 1.05, vjust = -0.6, size = 3, colour = "red") +
  scale_fill_viridis_d(guide = "none") +
  labs(title = "CSHQ total score by cohort",
       subtitle = paste0("Kruskal\u2013Wallis p = ", signif(kw_cshq, 2)),
       x = NULL, y = "CSHQ total score") +
  theme_minimal(base_size = 12) + ang_x
ggsave(file.path(out_dir, "cshq_total_by_cohort.png"), p_cshq,
       width = 8.5, height = 5.5, dpi = 600)

# significant subscales (long format)
subscale_labels <- c(SOD = "Sleep onset delay", SD = "Sleep duration",
                     NW = "Night waking", SDB = "Sleep-disordered breathing")
long <- do.call(rbind, lapply(dunn_sig, function(s)
  data.frame(cohort = sub_df$cohort, subscale = subscale_labels[[s]],
             score = sub_df[[s]], stringsAsFactors = FALSE)))
long$subscale <- factor(long$subscale, levels = unname(subscale_labels))
long$cohort   <- relabel_cohort(long$cohort)

p_sub <- ggplot(long, aes(cohort, score, fill = cohort)) +
  geom_violin(alpha = 0.30, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(width = 0.10, height = 0, size = 0.7, alpha = 0.35) +
  facet_wrap(~ subscale, scales = "free_y") +
  scale_fill_viridis_d(guide = "none") +
  labs(title = "CSHQ subscale totals by cohort", x = NULL, y = "Subscale total") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "subscale_violin_plot.png"), p_sub,
       width = 7.5, height = 5, dpi = 600)
















# ==========
# CYTOKINES
# ==========




firefly_cyto_raw <- read.csv(file.path(data_dir, "FIREFLY_cytokines.csv"),
                             stringsAsFactors = FALSE)
firefly_cyto_raw_vehicle <- firefly_cyto_raw[2:52, ] 
firefly_cyto_raw_LPS <- firefly_cyto_raw[58:108, ] 
firefly_cyto <- data.frame(
  cohort = c(rep("firefly_NE",31), rep('firefly_control', 20)),
  ID    = as.character(firefly_cyto_raw_vehicle[[2]]),
  EPO   = as.numeric(firefly_cyto_raw_vehicle[[3]]),
  GMCSF = as.numeric(firefly_cyto_raw_vehicle[[4]]),
  VEGF  = as.numeric(firefly_cyto_raw_vehicle[[5]]),
  IL8   = as.numeric(firefly_cyto_raw_vehicle[[6]]),
  EPO_LPS = as.numeric(firefly_cyto_raw_LPS[[3]]),
  GMCSF_LPS = as.numeric(firefly_cyto_raw_LPS[[4]]),
  VEGF_LPS = as.numeric(firefly_cyto_raw_LPS[[5]]),
  IL8_LPS = as.numeric(firefly_cyto_raw_LPS[[6]]),
  stringsAsFactors = FALSE
)
#We have DIV/0! values for the GMCSF, will impute with the LLOD
firefly_cyto$GMCSF[is.na(firefly_cyto$GMCSF)] <- 0.24696
firefly_cyto$GMCSF_LPS[is.na(firefly_cyto$GMCSF_LPS)] <- 0.24696




# Rest of the cytokines - raw data

champion_cp_cyto_raw <- read.csv(file.path(data_dir, "champion_CP_cytokines.csv"), stringsAsFactors = FALSE)
champion_cp_cyto_raw <- champion_cp_cyto_raw[1:13,]
champion_cp_cyto <- data.frame(
  cohort = 'champion_CP',
  ID = as.character(champion_cp_cyto_raw$sample.name),
  EPO = as.numeric(champion_cp_cyto_raw$EPO),
  GMCSF = as.numeric(champion_cp_cyto_raw$GM.CSF),
  VEGF = as.numeric(champion_cp_cyto_raw$VEGF),
  IL8 = as.numeric(champion_cp_cyto_raw$IL.8),
  EPO_LPS = as.numeric(champion_cp_cyto_raw$EPO.1),
  GMCSF_LPS = as.numeric(champion_cp_cyto_raw$GM.CSF.1),
  VEGF_LPS = as.numeric(champion_cp_cyto_raw$VEGF.1),
  IL8_LPS = as.numeric(champion_cp_cyto_raw$IL.8.1),
  stringsAsFactors = FALSE
)



champion_ne_cyto_raw <- read.csv(file.path(data_dir, "champion_NE_cytokines.csv"), stringsAsFactors = FALSE)


champion_ne_cyto <- data.frame(
  cohort = 'champion_NE',
  ID = as.character((champion_ne_cyto_raw$Baseline)),
  EPO = as.numeric(champion_ne_cyto_raw$EPO),
  GMCSF = as.numeric(champion_ne_cyto_raw$GM.CSF),
  VEGF = as.numeric(champion_ne_cyto_raw$VEGF),
  IL8 = as.numeric(champion_ne_cyto_raw$IL.8),
  EPO_LPS = as.numeric(champion_ne_cyto_raw$EPO.1),
  GMCSF_LPS = as.numeric(champion_ne_cyto_raw$GM.CSF.1),
  VEGF_LPS = as.numeric(champion_ne_cyto_raw$VEGF.1),
  IL8_LPS = as.numeric(champion_ne_cyto_raw$IL.8.1),
  stringsAsFactors = FALSE
)


champion_control_cyto_raw <- read.csv(file.path(data_dir, "champion_control_cytokines.csv"), stringsAsFactors = FALSE)

champion_control_cyto <- data.frame(
  cohort = 'champion_control',
  ID = paste0('CON', as.character(champion_control_cyto_raw$sample.name)),
  EPO = as.numeric(champion_control_cyto_raw$EPO),
  GMCSF = as.numeric(champion_control_cyto_raw$GM.CSF),
  VEGF = as.numeric(champion_control_cyto_raw$VEGF),
  IL8 = as.numeric(champion_control_cyto_raw$IL.8),
  EPO_LPS = as.numeric(champion_control_cyto_raw$EPO.1),
  GMCSF_LPS = as.numeric(champion_control_cyto_raw$GM.CSF.1),
  VEGF_LPS = as.numeric(champion_control_cyto_raw$VEGF.1),
  IL8_LPS = as.numeric(champion_control_cyto_raw$IL.8.1),
  stringsAsFactors = FALSE
)


serenity_cp_cyto_raw <- read.csv(file.path(data_dir, "serenity_CP_cytokines.csv"), stringsAsFactors = FALSE)
#Query around ID labelling - not consistent with the other serenity data

#We can relabel the age using the correspondence obtained from Johana
serenity_cp_id_correspondence <- c(
  p1 = 'TUH03',
  p2 = 'TUH06',
  p3 = 'TUH07',
  p4 = 'TUH09',
  p5 = 'TUH11',
  p6 = 'TUH14',
  p7 = 'TUH16',
  p8 = 'TUH17',
  p9 = 'TUH18',
  p10 = 'TUH20',
  p11 = 'TUH21',
  p12 = 'TUH22',
  p13 = 'TUH24',
  p14 = 'TUH26',
  p15 = 'TUH27'
)


serenity_cp_cyto <- data.frame(
  cohort = 'serenity_cp',
  ID = c(unname(serenity_cp_id_correspondence[serenity_cp_cyto_raw$Sample[2:16]]), toupper(as.character(serenity_cp_cyto_raw$Sample[17:26]))),
  EPO = as.numeric(serenity_cp_cyto_raw$EPO[2:26]),
  GMCSF = as.numeric(serenity_cp_cyto_raw$GM.CSF[2:26]),
  VEGF = as.numeric(serenity_cp_cyto_raw$VEGF[2:26]),
  IL8 = as.numeric(serenity_cp_cyto_raw$IL.8[2:26]),
  EPO_LPS = as.numeric(serenity_cp_cyto_raw$EPO.1[2:26]),
  GMCSF_LPS = as.numeric(serenity_cp_cyto_raw$GM.CSF.1[2:26]),
  VEGF_LPS = as.numeric(serenity_cp_cyto_raw$VEGF.1[2:26]),
  IL8_LPS = as.numeric(serenity_cp_cyto_raw$IL.8.1[2:26]),
  stringsAsFactors = FALSE
)
# We have one value of DIV/0! in the GMCSF_LPS column, so we replace this with the LLOD
serenity_cp_cyto$GMCSF_LPS[is.na(serenity_cp_cyto$GMCSF_LPS)] <- 0.24696



serenity_control_cyto_raw <- read.csv(file.path(data_dir, "serenity_control_cytokines.csv"), stringsAsFactors = FALSE)
#No sleep data currently available for the serenity control group 

serenity_control_cyto <- data.frame(
  cohort = 'serenity_control',
  ID = toupper(as.character(serenity_control_cyto_raw$BASELINE[2:27])),
  EPO = as.numeric(serenity_control_cyto_raw$EPO[2:27]),
  GMCSF = as.numeric(serenity_control_cyto_raw$GM.CSF..Human.[2:27]),
  VEGF = as.numeric(serenity_control_cyto_raw$VEGF[2:27]),
  IL8 = as.numeric(serenity_control_cyto_raw$IL.8[2:27]),
  EPO_LPS = as.numeric(serenity_control_cyto_raw$EPO.1[2:27]),
  GMCSF_LPS = as.numeric(serenity_control_cyto_raw$GM.CSF..Human..1[2:27]),
  VEGF_LPS = as.numeric(serenity_control_cyto_raw$VEGF.1[2:27]),
  IL8_LPS = as.numeric(serenity_control_cyto_raw$IL.8.1[2:27]),
  stringsAsFactors = FALSE
)

# We have DIV/0! for a number of the EPO rows
# We have DIV/0! for a number of the VEGF rows
# We have DIV/0! for a number of the IL8 rows
# We have DIV/0! for a number of the GMCSF_LPS rows

serenity_control_cyto$EPO[is.na(serenity_control_cyto$EPO)] <- 2.430539
serenity_control_cyto$VEGF[is.na(serenity_control_cyto$VEGF)] <- 4.984724
serenity_control_cyto$IL8[is.na(serenity_control_cyto$IL8)] <- 0.198944
serenity_control_cyto$GMCSF_LPS[is.na(serenity_control_cyto$GMCSF_LPS)] <- 0.24696




starfish_preterm_cyto_raw <- read.csv(file.path(data_dir, "starfish_preterm_cytokines.csv"), stringsAsFactors = FALSE)

starfish_preterm_cyto <- data.frame(
  cohort = 'starfish_preterm',
  ID = toupper(as.character(starfish_preterm_cyto_raw$X[2:24])),
  EPO = as.numeric(starfish_preterm_cyto_raw$X.11[2:24]),
  GMCSF = as.numeric(starfish_preterm_cyto_raw$X.23[2:24]),
  VEGF = as.numeric(starfish_preterm_cyto_raw$X.143[2:24]),
  IL8 = as.numeric(starfish_preterm_cyto_raw$X.107[2:24]),
  EPO_LPS = as.numeric(starfish_preterm_cyto_raw$X.8[2:24]),
  GMCSF_LPS = as.numeric(starfish_preterm_cyto_raw$X.20[2:24]),
  VEGF_LPS = as.numeric(starfish_preterm_cyto_raw$X.140[2:24]),
  IL8_LPS = as.numeric(starfish_preterm_cyto_raw$X.104[2:24]),
  stringsAsFactors = FALSE
)

#imputing the LLOD for DIV/0! in starfish preterm
starfish_preterm_cyto$GMCSF[is.na(starfish_preterm_cyto$GMCSF)] <- 0.24696
starfish_preterm_cyto$GMCSF_LPS[is.na(starfish_preterm_cyto$GMCSF_LPS)] <- 0.24696



starfish_control_cyto_raw <- read.csv(file.path(data_dir, "starfish_control_cytokines.csv"), stringsAsFactors = FALSE)
#No sleep data currently available for the starfish control group 


# GMCSF have a number of DIV0 entries - there are missing ones as well though so watch out for this
# same for the LPS version of GMCSF

starfish_control_cyto_raw$X.23[starfish_control_cyto_raw$X.23 == '#DIV/0!'] <- 0.24696
starfish_control_cyto_raw$X.20[starfish_control_cyto_raw$X.20 == '#DIV/0!'] <- 0.24696

starfish_control_cyto <- data.frame(
  cohort = 'starfish_control',
  ID = toupper(as.character(starfish_control_cyto_raw$X[2:18])),
  EPO = as.numeric(starfish_control_cyto_raw$X.11[2:18]),
  GMCSF = as.numeric(starfish_control_cyto_raw$X.23[2:18]),
  VEGF = as.numeric(starfish_control_cyto_raw$X.143[2:18]),
  IL8 = as.numeric(starfish_control_cyto_raw$X.107[2:18]),
  EPO_LPS = as.numeric(starfish_control_cyto_raw$X.8[2:18]),
  GMCSF_LPS = as.numeric(starfish_control_cyto_raw$X.20[2:18]),
  VEGF_LPS = as.numeric(starfish_control_cyto_raw$X.140[2:18]),
  IL8_LPS = as.numeric(starfish_control_cyto_raw$X.104[2:18]),
  stringsAsFactors = FALSE
)










# Stack every cohort table. bind_rows() fills any missing marker with NA.
cyto_all <- bind_rows(
  firefly_cyto,
  champion_cp_cyto,
  champion_ne_cyto,
  champion_control_cyto,
  starfish_preterm_cyto,
  serenity_cp_cyto,
  serenity_control_cyto,
  starfish_control_cyto
)


marker_cols <- setdiff(names(cyto_all), c("cohort", "ID"))


cyto_all <- cyto_all %>%
  mutate(across(all_of(marker_cols), ~ as.numeric(scale(.)), .names = "{.col}_z"))










































#OLD CODE HERE:

# # ---- join to CSHQ scores ---------------------------------------------
# 
scores_all$ID_clean <- clean_id(scores_all$ID)
age_df$ID_clean     <- clean_id(age_df$ID)
cyto_all$ID_clean   <- clean_id(cyto_all$ID)
#Adding in 'NE' to the champion NE ids for consistency with the cytokine data
idx <- scores_all$cohort == "champion_NE"

scores_all$ID[idx] <-
  ifelse(
    grepl("^NE\\s*", scores_all$ID[idx]),
    scores_all$ID[idx],
    paste0("NE ", scores_all$ID[idx])
  )


#Ensuring ID values are standard across data frames - capitalising, removing spaces

cyto_all$ID_clean <- gsub("\\s+", "", toupper(trimws(cyto_all$ID)))
scores_all$ID_clean <- gsub("\\s+", "", toupper(trimws(scores_all$ID)))

#Adding in the ages that we currently have
scores_all <- full_join(
  scores_all,
  age_df[, c("ID_clean","cohort","age_months")],
  by = c("ID_clean","cohort")
)
scores_all$age_months_z <- scale(scores_all$age_months)
scores_all$cohort[scores_all$cohort == 'serenity'] <- 'serenity_cp'





#Looking at the number of IDs that we have across the two datasets
common <- intersect(cyto_all$ID_clean, scores_all$ID_clean)
length(common)

only_cyto <- setdiff(cyto_all$ID_clean, scores_all$ID_clean)
only_scores <- setdiff(scores_all$ID_clean, cyto_all$ID_clean)

#Again, we have a Query about the 112T, 85T correspondence, but we'll continue for the moment

data_cyto <- full_join(scores_all, cyto_all, by = c("ID_clean", 'cohort'))
#There are only a few (4 in total) ID values common between the cyto data and the sleep data with the champion CP group

#For the moment, we focus only on the vehicle markers and we drop any of the rows with missing values in these markers
marker_cols_standardised <- c('EPO_z', 'GMCSF_z', 'VEGF_z', 'IL8_z')
#data_cyto_clean <- data_cyto[complete.cases(data_cyto[, c(marker_cols_standardised, 'cohort', 'CSHQ_total')]), ]


#NEED TO CONSIDER DUPLICATES HERE


# All rows involved in a duplicate ID_clean (first occurrence + repeats)
dupe_ids <- data_cyto$ID_clean[duplicated(data_cyto$ID_clean) | 
                                       duplicated(data_cyto$ID_clean, fromLast = TRUE)]
dupe_ids <- unique(dupe_ids)

# View all the actual rows for those IDs, sorted so duplicates sit together
dupes_df <- data_cyto[data_cyto$ID_clean %in% dupe_ids, ]
dupes_df <- dupes_df[order(dupes_df$ID_clean), ]
dupes_df

# Dealing with these duplicates
# For the ones with an NA id, we can just substitute UNLABELLED_X in
blank <- is.na(data_cyto$ID_clean) | data_cyto$ID_clean == ""
data_cyto$ID_clean[blank] <- paste0("UNLABELLED_", seq_len(sum(blank)))

#The duplicates for NE31 and NE50 don't make too much sense (different ages and cytokines) - we'll drop for now
drop_ids <- c("NE31", "NE50")
data_cyto <- data_cyto[!(data_cyto$cohort == "champion_NE" & data_cyto$ID_clean %in% drop_ids), ]





















#Pairs plots from before - commenting this out for the moment...


# 
# 
# ggpairs(data_cyto_clean[, c(marker_cols_standardised, "cohort")], columns = seq_along(marker_cols_standardised), aes(colour = cohort))
# ggsave(file.path(out_dir, "cytokine_pairs_plot.png"), width = 10, height = 10)
# 
# 
# # Running a linear regression for inflammation markers predicting CSHQ totals, adjusting for cohort membership
# 
# 
# CSHQ_lm_4_markers <- lm(CSHQ_total ~ cohort + EPO_z + GMCSF_z + VEGF_z + IL8_z, data = data_cyto_clean)
# 
# 
# 
# # EPO vs CSHQ total (stats read straight from the fitted model so they stay correct)
# epo_co <- coef(summary(CSHQ_lm_4_markers))["EPO_z", ]
# p_epo <- ggplot(data_cyto_clean, aes(EPO_z, CSHQ_total)) +
#   geom_smooth(method = "lm", se = TRUE, colour = "black") +
#   geom_point(aes(colour = cohort), size = 2, alpha = 0.85) +
#   scale_colour_brewer(palette = "Set1", name = "Cohort") +
#   labs(title = "EPO vs CSHQ total score",
#        subtitle = sprintf("beta = %.2f, p = %.3f", epo_co["Estimate"], epo_co["Pr(>|t|)"]),
#        x = "EPO (standardised)", y = "CSHQ total score") +
#   theme_minimal(base_size = 12)
# ggsave(file.path(out_dir, "epo_vs_cshq_total_more_cohorts.png"), p_epo, width = 7.5, height = 5, dpi = 600)
# 
# 
# 
# 
# #Looking at the difference between the stimulated cytokine response and the baseline response
# 
# #EPO
# m_EPO <- max(data_cyto$EPO, data_cyto$EPO_LPS, na.rm = TRUE)
# ggplot(data_cyto, aes(x = EPO, y = EPO_LPS, color = cohort)) +
#   geom_point() +
#   labs(
#     x = "EPO Vehicle",
#     y = "EPO LPS",
#     title = "Plot of EPO Vehicle vs. LPS"
#   ) +
#   scale_x_continuous(limits = c(0, m_EPO)) +
#   scale_y_continuous(limits = c(0, m_EPO)) +
#   theme_minimal()
# ggsave(file.path(out_dir, "LPS_vehicle_EPO.png"), width = 7.5, height = 5, dpi = 600)
# 
# 
# 
# #GMCSF
# m_GMCSF <- max(data_cyto$GMCSF, data_cyto$GMCSF_LPS, na.rm = TRUE)
# ggplot(data_cyto, aes(x = GMCSF, y = GMCSF_LPS, color = cohort)) +
#   geom_point() +
#   labs(
#     x = "GM-CSF Vehicle",
#     y = "GM-CSF LPS",
#     title = "Plot of GM-CSF Vehicle vs. LPS"
#   ) +
#   scale_x_continuous(limits = c(0, m_GMCSF)) +
#   scale_y_continuous(limits = c(0, m_GMCSF)) +
#   theme_minimal()
# ggsave(file.path(out_dir, "LPS_vehicle_GMCSF.png"), width = 7.5, height = 5, dpi = 600)
# 
# 
# 
# #IL8
# m_IL8 <- max(data_cyto$IL8, data_cyto$IL8_LPS, na.rm = TRUE)
# ggplot(data_cyto, aes(x = IL8, y = IL8_LPS, color = cohort)) +
#   geom_point() +
#   labs(
#     x = "IL8 Vehicle",
#     y = "IL8 LPS",
#     title = "Plot of IL8 Vehicle vs. LPS"
#   ) +
#   scale_x_continuous(limits = c(0, m_IL8)) +
#   scale_y_continuous(limits = c(0, m_IL8)) +
#   theme_minimal()
# ggsave(file.path(out_dir, "LPS_vehicle_IL8.png"), width = 7.5, height = 5, dpi = 600)
# 
# 
# #VEGF
# 
# m_VEGF <- max(data_cyto$VEGF, data_cyto$VEGF_LPS, na.rm = TRUE)
# ggplot(data_cyto, aes(x = VEGF, y = VEGF_LPS, color = cohort)) +
#   geom_point() +
#   labs(
#     x = "VEGF Vehicle",
#     y = "VEGF LPS",
#     title = "Plot of VEGF Vehicle vs. LPS"
#   ) +
#   scale_x_continuous(limits = c(0, m_VEGF)) +
#   scale_y_continuous(limits = c(0, m_VEGF)) +
#   theme_minimal()
# ggsave(file.path(out_dir, "LPS_vehicle_VEGF.png"), width = 7.5, height = 5, dpi = 600)
# 
# 
# 
# 
# #We want to now consider the analysis split into groups of: controls, those at high risk of CP and those who have CP.
# 
# #Initially we'll look at controls vs. non-controls, and consider both LPS and vehicle in the regression
# 
# #Control groups are the cohorts of firefly_control and champion control, however we currently do not have CSHQ responses for the champion control group
# 
# data_cyto$CG <- ifelse(data_cyto$cohort == 'firefly_control', 1, 0)
# data_cyto_clean$CG <- ifelse(data_cyto_clean$cohort == 'firefly_control', 1, 0)
# 
# data_cyto$CG <- factor(data_cyto$CG,
#                        levels = c(1, 0),
#                        labels = c("control", "non-control"))
# 
# data_cyto_clean$CG <- factor(data_cyto_clean$CG,
#                              levels = c(1, 0),
#                              labels = c("control", "non-control"))
# 
# #Running a linear regression, including the control/non-control variable, and the standardised markers for both vehicle and LPS
# 
# #we'll drop any rows with missingness in the LPS columns now
# data_cyto_clean_LPS <- data_cyto_clean[complete.cases(data_cyto_clean[, c('EPO_LPS_z', 'GMCSF_LPS_z', 'IL8_LPS_z', 'VEGF_LPS_z', 'cohort', 'CSHQ_total')]), ]
# 
# 
# lm_4_markers_vehicle_LPS_control_vs_no_control <- lm(CSHQ_total ~ CG + EPO_z + GMCSF_z + VEGF_z + IL8_z + EPO_LPS_z + GMCSF_LPS_z + VEGF_LPS_z + IL8_LPS_z, 
#                                                      data = data_cyto_clean_LPS)
# 
# 
# 
# 
# 
# #Looking to include age in the analysis now
# data_cyto_clean_age <- data_cyto_clean[complete.cases(data_cyto_clean$age_months_z), ]
# data_cyto_clean_age_and_LPS  <- data_cyto_clean_LPS[ complete.cases(data_cyto_clean_LPS$age_months_z),]
# lm_4_markers_vehicle_LPS_control_vs_no_control_with_age <- lm(CSHQ_total ~ CG + EPO_z + GMCSF_z + VEGF_z + IL8_z + EPO_LPS_z + GMCSF_LPS_z + VEGF_LPS_z + IL8_LPS_z + age_months_z,
#                                                               data = data_cyto_clean_age_and_LPS) 
# 
# 
# #Generating a pairs plot with all of the different variables and the CSHQ totals, then we will look at it with a bit more granularity
# 
# ggpairs(data_cyto_clean_LPS[, c('EPO_z', 'GMCSF_z', 'IL8_z', 'VEGF_z','EPO_LPS_z', 'GMCSF_LPS_z', 'IL8_LPS_z', 'VEGF_LPS_z', 'CSHQ_total', 'CG')], columns = seq_along(c('EPO_z', 'GMCSF_z', 'IL8_z', 'VEGF_z','EPO_LPS_z', 'GMCSF_LPS_z', 'IL8_LPS_z', 'VEGF_LPS_z', 'CSHQ_total')), aes(colour = CG))
# ggsave(file.path(out_dir, "marker_pairs_plot.png"), width = 500, height = 500, units = 'px')
# ggsave(file.path(out_dir, "marker_pairs_plot.png"), width = 15, height = 15, units = "in", dpi = 600)















#We're going to run a model for each of the markers (along with their corresponding LPS) for predicting sleep score - excluding age for the moment

#We'll split into 3 groups - controls, those who were at high risk of CP and didn't get it and those who have CP

# Controls: firefly controls, champion controls

# High risk of CP: firefly NE, starfish preterm, champion NE (minus those who got CP)

# CP: Serenity (however we are missing serenity CSHQ since the ID labelling doesn't match - not usable currently), champion CP

# Along with this, there are specific members from firefly NE, starfish preterm and champion NE that move into CP.

# These are: S12 from starfish, and NE00224, NE00205 from firefly

data_cyto_with_riskgroup <- data_cyto

# CN = Control, HR = High Risk, CP = Cerebral Palsy
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'firefly_control'] <- 'CN' 
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'champion_control'] <- 'CN'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'serenity_control'] <- 'CN'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'firefly_NE'] <- 'HR'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'starfish_preterm'] <- 'HR'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'champion_NE'] <- 'HR'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'champion_CP'] <- 'CP'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$cohort == 'serenity_cp'] <- 'CP'

#Those who moved to CP from HR
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$ID_clean == 'S12'] <- 'CP'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$ID_clean == '56'] <- 'CP'
data_cyto_with_riskgroup$riskgroup[data_cyto_with_riskgroup$ID_clean == '30'] <- 'CP'


data_cyto_with_riskgroup$riskgroup <- as.factor(data_cyto_with_riskgroup$riskgroup)

data_cyto_with_riskgroup$riskgroup_CP <- ifelse(data_cyto_with_riskgroup$riskgroup == 'CP', 1, 0)

data_cyto_with_riskgroup$riskgroup_HR <- ifelse(data_cyto_with_riskgroup$riskgroup == 'HR' | data_cyto_with_riskgroup$riskgroup == 'CP', 1, 0)


#EPO

lm_EPO_vehicle_LPS_3groups <- lm(CSHQ_total ~ EPO_z + EPO_LPS_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_EPO_vehicle_LPS_3groups_with_age <- lm(CSHQ_total ~ EPO_z + EPO_LPS_z + riskgroup_HR + riskgroup_CP + age_months_z, data = data_cyto_with_riskgroup)


#GMCSF

lm_GMCSF_vehicle_LPS_3groups <- lm(CSHQ_total ~ GMCSF_z + GMCSF_LPS_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_GMCSF_vehicle_LPS_3groups_with_age <- lm(CSHQ_total ~ GMCSF_z + GMCSF_LPS_z + riskgroup_HR + riskgroup_CP + age_months_z, data = data_cyto_with_riskgroup)

#VEGF

lm_VEGF_vehicle_LPS_3groups <- lm(CSHQ_total ~ VEGF_z + VEGF_LPS_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_VEGF_vehicle_LPS_3groups_with_age <- lm(CSHQ_total ~ VEGF_z + VEGF_LPS_z + riskgroup_HR + riskgroup_CP + age_months_z, data = data_cyto_with_riskgroup)


#IL8

lm_IL8_vehicle_LPS_3groups <- lm(CSHQ_total ~ IL8_z + IL8_LPS_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_IL8_vehicle_LPS_3groups_with_age <- lm(CSHQ_total ~ IL8_z + IL8_LPS_z + riskgroup_HR + riskgroup_CP + age_months_z, data = data_cyto_with_riskgroup)




#Looking at regressions for the LPS markers by cohorts, accounting for age



lm_EPO_LPS_3groups <- lm(EPO_LPS_z ~ age_months_z + riskgroup_CP + riskgroup_HR, data = data_cyto_with_riskgroup)
lm_GMCSF_LPS_3groups <- lm(GMCSF_LPS_z ~ age_months_z + riskgroup_CP + riskgroup_HR, data = data_cyto_with_riskgroup)
lm_VEGF_LPS_3groups <- lm(VEGF_LPS_z ~ age_months_z + riskgroup_CP + riskgroup_HR, data = data_cyto_with_riskgroup)
lm_IL8_LPS_3groups <- lm(IL8_LPS_z ~ age_months_z + riskgroup_CP + riskgroup_HR, data = data_cyto_with_riskgroup)










# Looking at the (standardised) difference between LPS and vehicle for the regression

data_cyto_with_riskgroup$EPO_diff <- data_cyto_with_riskgroup$EPO_LPS - data_cyto_with_riskgroup$EPO
data_cyto_with_riskgroup$GMCSF_diff <- data_cyto_with_riskgroup$GMCSF_LPS - data_cyto_with_riskgroup$GMCSF
data_cyto_with_riskgroup$IL8_diff <- data_cyto_with_riskgroup$IL8_LPS - data_cyto_with_riskgroup$IL8
data_cyto_with_riskgroup$VEGF_diff <- data_cyto_with_riskgroup$VEGF_LPS - data_cyto_with_riskgroup$VEGF

data_cyto_with_riskgroup$EPO_diff_z <- scale(data_cyto_with_riskgroup$EPO_diff)
data_cyto_with_riskgroup$GMCSF_diff_z <- scale(data_cyto_with_riskgroup$GMCSF_diff)
data_cyto_with_riskgroup$IL8_diff_z <- scale(data_cyto_with_riskgroup$IL8_diff)
data_cyto_with_riskgroup$VEGF_diff_z <- scale(data_cyto_with_riskgroup$VEGF_diff)

#WITH AGE

lm_EPO_diff_3groups_with_age <- lm(EPO_diff_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_GMCSF_diff_3groups_with_age <- lm(GMCSF_diff_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_VEGF_diff_3groups_with_age <- lm(IL8_diff_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
# We have a significant (p = 0.004) coefficient for the high risk group for VEGF (beta = 0.5)
lm_IL8_diff_3groups_with_age <- lm(VEGF_diff_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)


#WITHOUT AGE

lm_EPO_diff_3groups <- lm(EPO_diff_z ~ riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_GMCSF_diff_3groups <- lm(GMCSF_diff_z ~ riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_VEGF_diff_3groups <- lm(IL8_diff_z ~ riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
# We have a significant (p < 0.001) coefficient for the high risk group for VEGF (beta = 0.66) and p = 0.01 for CP with beta = -0.4
lm_IL8_diff_3groups <- lm(VEGF_diff_z ~ riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
# We have a significant (p = 0.001) coefficient for the high risk group for IL8 (beta = 0.46) and p = 0.01 for CP with beta = -0.43



# Looking at the (standardised) difference between log LPS and log vehicle for the regression

data_cyto_with_riskgroup$EPO_diff_log2 <- log2(data_cyto_with_riskgroup$EPO_LPS) - log2(data_cyto_with_riskgroup$EPO)
data_cyto_with_riskgroup$GMCSF_diff_log2 <- log2(data_cyto_with_riskgroup$GMCSF_LPS) - log2(data_cyto_with_riskgroup$GMCSF)
data_cyto_with_riskgroup$IL8_diff_log2 <- log2(data_cyto_with_riskgroup$IL8_LPS) - log2(data_cyto_with_riskgroup$IL8)
data_cyto_with_riskgroup$VEGF_diff_log2 <- log2(data_cyto_with_riskgroup$VEGF_LPS) - log2(data_cyto_with_riskgroup$VEGF)

data_cyto_with_riskgroup$EPO_diff_log2_z <- scale(data_cyto_with_riskgroup$EPO_diff_log2)
data_cyto_with_riskgroup$GMCSF_diff_log2_z <- scale(data_cyto_with_riskgroup$GMCSF_diff_log2)
data_cyto_with_riskgroup$IL8_diff_log2_z <- scale(data_cyto_with_riskgroup$IL8_diff_log2)
data_cyto_with_riskgroup$VEGF_diff_log2_z <- scale(data_cyto_with_riskgroup$VEGF_diff_log2)

#WITH AGE

lm_EPO_diff_log_3groups_with_age <- lm(EPO_diff_log2_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_GMCSF_diff_log_3groups_with_age <- lm(GMCSF_diff_log2_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
# high risk group is a significant (p<0.001) predictor here for GMCSF, with beta = 0.78
lm_VEGF_diff_log_3groups_with_age <- lm(IL8_diff_log2_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_IL8_diff_log_3groups_with_age <- lm(VEGF_diff_log2_z ~ age_months_z + riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)


#WITHOUT AGE

lm_EPO_diff_log_3groups <- lm(EPO_diff_log2_z ~  riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_GMCSF_diff_log_3groups <- lm(GMCSF_diff_log2_z ~  riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_VEGF_diff_log_3groups <- lm(IL8_diff_log2_z ~ riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)
lm_IL8_diff_log_3groups <- lm(VEGF_diff_log2_z ~  riskgroup_HR + riskgroup_CP, data = data_cyto_with_riskgroup)

#ALL OF THESE SEEM TO GIVE SOME SIGNIFICANT RESULTS






# Visualisations: pairs plot with groups taken to be the (disjoint) high risk and CP groups along with controls

# First looking at baseline (vehicle)
ggpairs(data_cyto_with_riskgroup, columns = 48:51, aes(color = riskgroup))

# Now the LPS

# Now the difference

# Now the log-ratio 




# data_cyto$EPO_diff <- data_cyto$EPO_LPS - data_cyto$EPO
# data_cyto$GMCSF_diff <- data_cyto$GMCSF_LPS - data_cyto$GMCSF
# data_cyto$IL8_diff <- data_cyto$IL8_LPS - data_cyto$IL8
# data_cyto$VEGF_diff <- data_cyto$VEGF_LPS - data_cyto$VEGF
# 
# data_cyto$EPO_diff_z <- scale(data_cyto$EPO_diff)
# data_cyto$GMCSF_diff_z <- scale(data_cyto$GMCSF_diff)
# data_cyto$IL8_diff_z <- scale(data_cyto$IL8_diff)
# data_cyto$VEGF_diff_z <- scale(data_cyto$VEGF_diff)
# 
# data_cyto_clean_diff <- data_cyto[complete.cases(data_cyto[, c('EPO_diff_z', 'GMCSF_diff_z', 'IL8_diff_z', 'VEGF_diff_z', 'cohort', 'CSHQ_total')]), ]
# 
# 
# CSHQ_lm_4_markers_diff <- lm(CSHQ_total ~ cohort + EPO_diff_z + GMCSF_diff_z + IL8_diff_z + VEGF_diff_z, data = data_cyto_clean_diff)
















#The previous analysis when we only had the firefly cohort

# # ---- firefly cytokine analyses ---------------------------------------
# firefly_cyto_df <- subset(data_cyto, study == "firefly")
# firefly_cyto_df <- firefly_cyto_df[complete.cases(firefly_cyto_df[, marker_cols]), ]
# ## CHECK: original dropped a few firefly rows by position; here rows are dropped
# ##        only when a marker is missing - confirm the same participants leave.
# 
# # pairwise marker relationships, coloured by cohort
# ggpairs(firefly_cyto_df[, c(marker_cols, "cohort")],
#         columns = seq_along(marker_cols), aes(colour = cohort))
# ggsave(file.path(out_dir, "firefly_cytokine_pairs_plot.png"), width = 10, height = 10)
# 
# # correlation matrix (Pearson)
# res <- rcorr(as.matrix(firefly_cyto_df[, marker_cols]))
# res$r; res$P
# 
# # VEGF vs IL8 looked positively correlated -> confirm with Spearman, overall + by cohort
# cor.test(firefly_cyto_df$VEGF, firefly_cyto_df$IL8, method = "spearman")
# by(firefly_cyto_df, firefly_cyto_df$cohort,
#    function(d) cor.test(d$VEGF, d$IL8, method = "spearman"))
# 
# # subscale + CSHQ totals ~ cohort + standardised markers
# marker_z   <- paste0(marker_cols, "_z")
# resp_names <- c(subscale_names, "CSHQ_total")
# cyto_models <- lapply(setNames(resp_names, resp_names),
#                       function(y) lm(reformulate(c("cohort", marker_z), y),
#                                      data = firefly_cyto_df))
# 
# # FDR per marker across the 8 subscales
# sub_pvals <- sapply(subscale_names, function(s)
#   coef(summary(cyto_models[[s]]))[marker_z, "Pr(>|t|)"])   # rows = markers, cols = subscales
# sub_pvals_fdr <- t(apply(sub_pvals, 1, p.adjust, method = "fdr"))
# round(sub_pvals_fdr, 4)
# 
# # EPO vs CSHQ total (stats read straight from the fitted model so they stay correct)
# epo_co <- coef(summary(cyto_models[["CSHQ_total"]]))["EPO_z", ]
# p_epo <- ggplot(firefly_cyto_df, aes(EPO_z, CSHQ_total)) +
#   geom_smooth(method = "lm", se = TRUE, colour = "black") +
#   geom_point(aes(colour = cohort), size = 2, alpha = 0.85) +
#   scale_colour_brewer(palette = "Set1", name = "Cohort") +
#   labs(title = "EPO vs CSHQ total score",
#        subtitle = sprintf("beta = %.2f, p = %.3f", epo_co["Estimate"], epo_co["Pr(>|t|)"]),
#        x = "EPO (standardised)", y = "CSHQ total score") +
#   theme_minimal(base_size = 12)
# ggsave(file.path(out_dir, "epo_vs_cshq_total.png"), p_epo, width = 7.5, height = 5, dpi = 600)
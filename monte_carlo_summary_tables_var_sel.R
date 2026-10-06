source('./sim_study_monte_carlo_functions.R')
source('./sim_study_monte_carlo_summary_analysis.R')

library(xtable)

n_replicates <- 300
sample_sizes <- c(150, 300, 500)


sim1_results_N150 <- readRDS('./sim1_results_N150.rds')
sim1_results_N300 <- readRDS('./sim1_results_N300.rds')
sim1_results_N500 <- readRDS('./sim1_results_N500.rds')

sim2_results_N150 <- readRDS('./sim2_results_N150.rds')
sim2_results_N300 <- readRDS('./sim2_results_N300.rds')
sim2_results_N500 <- readRDS('./sim2_results_N500.rds')

sim3_results_N150 <- readRDS('./sim3_results_N150.rds')
sim3_results_N300 <- readRDS('./sim3_results_N300.rds')
sim3_results_N500 <- readRDS('./sim3_results_N500.rds')

sim1_pred_sel_N150 <- pred_var_sel_summarise(fit = sim1_results_N150, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_pred_sel_N300 <- pred_var_sel_summarise(fit = sim1_results_N300, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_pred_sel_N500 <- pred_var_sel_summarise(fit = sim1_results_N500, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)


sim1_item_sel_N150 <- item_var_sel_summarise(fit = sim1_results_N150, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_item_sel_N300 <- item_var_sel_summarise(fit = sim1_results_N300, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim1_item_sel_N500 <- item_var_sel_summarise(fit = sim1_results_N500, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)

sim2_pred_sel_N150 <- pred_var_sel_summarise(fit = sim2_results_N150, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_pred_sel_N300 <- pred_var_sel_summarise(fit = sim2_results_N300, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_pred_sel_N500 <- pred_var_sel_summarise(fit = sim2_results_N500, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE), 
                                             threshold = 0.5)


sim2_item_sel_N150 <- item_var_sel_summarise(fit = sim2_results_N150, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_item_sel_N300 <- item_var_sel_summarise(fit = sim2_results_N300, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)
sim2_item_sel_N500 <- item_var_sel_summarise(fit = sim2_results_N500, 
                                             true_item_active = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 
                                             threshold = 0.5)


sim3_pred_sel_N150 <- pred_var_sel_summarise(fit = sim3_results_N150, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE), 
                                             threshold = 0.5)
sim3_pred_sel_N300 <- pred_var_sel_summarise(fit = sim3_results_N300, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE), 
                                             threshold = 0.5)
sim3_pred_sel_N500 <- pred_var_sel_summarise(fit = sim3_results_N500, 
                                             true_pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE), 
                                             threshold = 0.5)


sim3_item_sel_N150 <- item_var_sel_summarise(fit = sim3_results_N150, 
                                             true_item_active = c(rep(TRUE, 20), rep(FALSE, 20)), 
                                             threshold = 0.5)
sim3_item_sel_N300 <- item_var_sel_summarise(fit = sim3_results_N300, 
                                             true_item_active = c(rep(TRUE, 20), rep(FALSE, 20)), 
                                             threshold = 0.5)
sim3_item_sel_N500 <- item_var_sel_summarise(fit = sim3_results_N500, 
                                             true_item_active = c(rep(TRUE, 20), rep(FALSE, 20)), 
                                             threshold = 0.5)




get_sel <- function(x) {
  pr <- x$per_replicate
  
  cand <- setdiff(names(x), "per_replicate")
  has_pip <- vapply(cand, function(nm) is.data.frame(x[[nm]]) && "mean_pip" %in% names(x[[nm]]),
                    logical(1))
  if (!any(has_pip)) stop("No per-variable element with a mean_pip column. Names found: ",
                          paste(names(x), collapse = ", "))
  pv <- x[[cand[has_pip][1]]]
  
  list(tpr   = pr$TPR,
       fpr   = pr$FPR,
       exact = mean(pr$FN == 0 & pr$FP == 0),
       pip   = as.numeric(pv$mean_pip))
}
fmt_ms <- function(x, d = 3) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("--")
  sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x), sd(x))
}

fmt_num <- function(x, d = 3) if (is.na(x)) "--" else sprintf(paste0("%.", d, "f"), x)


sel_block <- function(title, sel_by_N, active) {
  header <- data.frame(a = title, tpr = "", fpr = "", pa = "", pi = "",
                       stringsAsFactors = FALSE)
  body <- do.call(rbind, lapply(names(sel_by_N), function(n) {
    s <- get_sel(sel_by_N[[n]])
    data.frame(a   = sprintf("$N = %s$", n),
               tpr = fmt_ms(s$tpr),
               fpr = if (any(!active)) fmt_ms(s$fpr) else "--",
               pa  = fmt_num(mean(s$pip[active])),
               pi  = if (any(!active)) fmt_num(mean(s$pip[!active])) else "--",
               stringsAsFactors = FALSE)
  }))
  rbind(header, body)
}

make_sel_main_table <- function(pred_by_N, item_by_N, pred_active, item_active,
                                caption, label, file = "") {
  tab <- rbind(sel_block("\\textit{Predictors}", pred_by_N, pred_active),
               sel_block("\\textit{Items}",      item_by_N, item_active))
  colnames(tab) <- c("", "TPR", "FPR", "Mean PIP (active)", "Mean PIP (inactive)")
  n_per <- length(pred_by_N) + 1
  print(xtable(tab, align = c("l", "l", "c", "c", "c", "c"),
               caption = caption, label = label),
        file = file, include.rownames = FALSE, booktabs = TRUE,
        table.placement = "H", caption.placement = "bottom",
        sanitize.text.function = identity,
        sanitize.colnames.function = identity,
        hline.after = c(-1, 0, n_per, 2 * n_per))
}

# Full per-variable mean PIP table (supplement)
make_pip_table <- function(sel_by_N, active, var_label, caption, label,
                           file = "", longtable = FALSE) {
  pips <- sapply(sel_by_N, function(x) get_sel(x)$pip)   # variables x N
  tab <- data.frame(Var = as.character(seq_len(nrow(pips))),
                    Active = ifelse(active, "Yes", "No"),
                    apply(pips, 2, function(v) sprintf("%.3f", v)),
                    stringsAsFactors = FALSE)
  colnames(tab) <- c(var_label, "Active", sprintf("$N = %s$", names(sel_by_N)))
  print(xtable(tab, align = c("l", "l", "c", rep("c", ncol(pips))),
               caption = caption, label = label),
        file = file, include.rownames = FALSE, booktabs = TRUE,
        floating = !longtable,
        tabular.environment = if (longtable) "longtable" else "tabular",
        table.placement = if (longtable) NULL else "H",
        caption.placement = "bottom",
        sanitize.text.function = identity,
        sanitize.colnames.function = identity)
}

# ---- Run for all three simulations ---------------------------------------

dir.create("./sim_study_tables", showWarnings = FALSE)

sel_specs <- list(
  sim1 = list(
    pred = list("150" = sim1_pred_sel_N150, "300" = sim1_pred_sel_N300, "500" = sim1_pred_sel_N500),
    item = list("150" = sim1_item_sel_N150, "300" = sim1_item_sel_N300, "500" = sim1_item_sel_N500),
    pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    item_active = c(rep(TRUE, 4), rep(FALSE, 4)),
    longtable = FALSE, num = 1),
  sim2 = list(
    pred = list("150" = sim2_pred_sel_N150, "300" = sim2_pred_sel_N300, "500" = sim2_pred_sel_N500),
    item = list("150" = sim2_item_sel_N150, "300" = sim2_item_sel_N300, "500" = sim2_item_sel_N500),
    pred_active = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    item_active = c(rep(TRUE, 8), rep(FALSE, 5)),
    longtable = FALSE, num = 2),
  sim3 = list(
    pred = list("150" = sim3_pred_sel_N150, "300" = sim3_pred_sel_N300, "500" = sim3_pred_sel_N500),
    item = list("150" = sim3_item_sel_N150, "300" = sim3_item_sel_N300, "500" = sim3_item_sel_N500),
    pred_active = c(rep(TRUE, 4), FALSE),
    item_active = c(rep(TRUE, 20), rep(FALSE, 20)),
    longtable = TRUE, num = 3)
)

for (s in names(sel_specs)) {
  sp <- sel_specs[[s]]
  
  make_sel_main_table(
    sp$pred, sp$item, sp$pred_active, sp$item_active,
    caption = sprintf("Variable selection performance for simulation %d. TPR and FPR are means (SD) across replicates; PIP is the posterior inclusion probability averaged over replicates and over variables in each group.", sp$num),
    label   = sprintf("table:%s_var_sel_summary", s),
    file    = sprintf("./sim_study_tables/%s_var_sel_summary.tex", s))
  
  make_pip_table(
    sp$item, sp$item_active, "Item",
    caption = sprintf("Mean posterior inclusion probability of each item across replicates, simulation %d.", sp$num),
    label   = sprintf("table:%s_item_pip", s),
    file    = sprintf("./sim_study_tables/%s_item_pip.tex", s),
    longtable = sp$longtable)
  
  make_pip_table(
    sp$pred, sp$pred_active, "Predictor",
    caption = sprintf("Mean posterior inclusion probability of each predictor across replicates, simulation %d.", sp$num),
    label   = sprintf("table:%s_pred_pip", s),
    file    = sprintf("./sim_study_tables/%s_pred_pip.tex", s))
}








































# Per-variable table (variable, true_active, mean_pip, mean_ind) from a summary object

sel_summary <- function(pr) {
  fpr <- pr$FPR[!is.na(pr$FPR)]
  c(TPR   = mean(pr$TPR),
    FPR   = if (length(fpr)) mean(fpr) else NA)
}

get_pv <- function(x) {
  cand <- setdiff(names(x), "per_replicate")
  ok <- vapply(cand, function(nm) is.data.frame(x[[nm]]) && "mean_pip" %in% names(x[[nm]]),
               logical(1))
  x[[cand[ok][1]]]
}

cols <- gray(c(0.9, 0.6, 0.3))   


pip_barplot <- function(sim, type, metric = "mean_pip", ylab = "Mean PIP",
                        main = NULL) {
  sp  <- sel_specs[[sim]]
  sel <- sp[[type]]
  act <- if (type == "pred") sp$pred_active else sp$item_active
  
  if (is.null(main))
    main <- sprintf("Simulation %d: Posterior Inclusion Probabilities (%s)",
                    sp$num, if (type == "pred") "Predictors" else "Items")
  
  pips <- sapply(sel, function(x) get_pv(x)[[metric]])    # variables x N
  mid  <- barplot(t(pips), beside = TRUE, col = cols, ylim = c(0, 1),
                  names.arg = seq_len(nrow(pips)),
                  xlab = "Variable", ylab = ylab, main = main)
  abline(h = 0.5, lty = 2)
  if (any(!act))
    abline(v = (max(mid[, max(which(act))]) + min(mid[, min(which(!act))])) / 2, lty = 3)
  legend("topright", legend = paste("N =", names(sel)), fill = cols)
}

rates_barplot <- function(sim, type, main = NULL) {
  sp  <- sel_specs[[sim]]
  sel <- sp[[type]]
  
  if (is.null(main))
    main <- sprintf("Simulation %d: Variable Selection Performance (%s)",
                    sp$num, if (type == "pred") "Predictors" else "Items")
  
  mat <- sapply(sel, function(x) sel_summary(x$per_replicate))   
  barplot(t(mat), beside = TRUE, col = cols, ylim = c(0, 1), main = main)
  legend("topright", legend = paste("N =", names(sel)), fill = cols)
}
# ---- Use: run one line at a time ---------------------------------------------
pip_barplot("sim1", "pred")
pip_barplot("sim1", "item")
rates_barplot("sim1", "pred")
rates_barplot("sim1", "item")

pip_barplot("sim2", "pred")
pip_barplot("sim2", "item")
rates_barplot("sim2", "pred")
rates_barplot("sim2", "item")

pip_barplot("sim3", "pred")
pip_barplot("sim3", "item")
rates_barplot("sim3", "pred")
rates_barplot("sim3", "item")

for (s in c("sim1", "sim2", "sim3")) {
  for (type in c("pred", "item")) {
    # wider for Sim 3 items (40 items)
    w <- if (s == "sim3" && type == "item") 12 else 7
    
    pdf(sprintf("./sim_study_plots/%s_%s_pip.pdf", s, type), width = w, height = 4.5)
    pip_barplot(s, type)
    dev.off()
    
    pdf(sprintf("./sim_study_plots/%s_%s_rates.pdf", s, type), width = 7, height = 4.5)
    rates_barplot(s, type)
    dev.off()
  }
}

# pip_barplot("sim2", "item")
# pip_barplot("sim3", "item")

# selection frequency instead of PIP:
# pip_barplot("sim1", "item", metric = "mean_ind", ylab = "Selection frequency")

# to save one:
# pdf("sim1_item_pip.pdf", width = 6, height = 4); pip_barplot("sim1", "item"); dev.off()
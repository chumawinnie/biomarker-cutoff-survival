# =====================================================================
# biomarker_cutoff.R  (reusable functions)
# Prognostic biomarker cutoff & survival subgroup analysis
# Data: curatedOvarianData (TCGA ovarian) | Biomarker: configurable
# =====================================================================
set.seed(42)

.need <- c("survival", "survminer", "maxstat", "pROC")
.miss <- .need[!.need %in% rownames(installed.packages())]
if (length(.miss)) install.packages(.miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(Biobase); library(curatedOvarianData)
  library(survival); library(survminer); library(maxstat); library(pROC)
})

# ---------------------------- functions ------------------------------

load_cohort <- function(gene, eset_name = "TCGA_eset") {
  data(list = eset_name, package = "curatedOvarianData")
  eset <- get(eset_name)
  if (!gene %in% featureNames(eset))
    stop(sprintf("'%s' not in %s; try CD8A, GZMB, PDCD1", gene, eset_name))
  pd <- pData(eset)
  dat <- data.frame(
    sample_id = sampleNames(eset),
    biomarker = as.numeric(exprs(eset)[gene, ]),
    age       = as.numeric(pd$age_at_initial_pathologic_diagnosis),
    time      = pd$days_to_death / 365.25,
    event     = as.integer(pd$vital_status == "deceased"),
    stringsAsFactors = FALSE
  )
  dat[is.finite(dat$biomarker) & is.finite(dat$time) &
      dat$time > 0 & dat$event %in% c(0, 1), ]
}

determine_cutpoint <- function(dat, marker = "biomarker", time = "time",
                               event = "event", minprop = 0.10, maxprop = 0.90) {
  f  <- as.formula(sprintf("Surv(%s, %s) ~ %s", time, event, marker))
  ms <- maxstat.test(f, data = dat, smethod = "LogRank", pmethod = "Lau94",
                     minprop = minprop, maxprop = maxprop)
  list(cutoff = as.numeric(ms$estimate), corrected_p = ms$p.value, maxstat = ms)
}

roc_cutpoint <- function(dat, marker = "biomarker", time = "time",
                         event = "event", horizon = 3) {
  elig <- dat[dat[[time]] >= horizon | dat[[event]] == 1, ]
  y    <- as.integer(elig[[time]] <= horizon & elig[[event]] == 1)
  ro   <- roc(y, elig[[marker]], quiet = TRUE)
  best <- coords(ro, "best", best.method = "youden", transpose = FALSE)
  list(threshold = best$threshold[1], auc = as.numeric(auc(ro)), roc = ro)
}

assign_groups <- function(dat, cutoff, marker = "biomarker") {
  dat$group <- factor(ifelse(dat[[marker]] > cutoff, "High", "Low"),
                      levels = c("Low", "High"))
  dat
}

run_survival <- function(dat, covariates = "age", time = "time", event = "event") {
  base <- sprintf("Surv(%s, %s) ~ group", time, event)
  rhs  <- paste(c("group", covariates), collapse = " + ")
  list(km      = survfit(as.formula(base), data = dat),
       cox_uni = coxph(as.formula(base), data = dat),
       cox_adj = coxph(as.formula(sprintf("Surv(%s, %s) ~ %s", time, event, rhs)),
                       data = dat))
}

validate_cutpoint <- function(dat, marker = "biomarker", time = "time",
                              event = "event", train_frac = 0.7, n_boot = 500,
                              minprop = 0.10, maxprop = 0.90) {
  idx   <- sample(nrow(dat), floor(train_frac * nrow(dat)))
  train <- dat[idx, ]; test <- dat[-idx, ]
  cut_tr <- determine_cutpoint(train, marker, time, event, minprop, maxprop)$cutoff
  cx     <- coxph(Surv(time, event) ~ group, data = assign_groups(test, cut_tr, marker))
  boot <- replicate(n_boot, {
    b <- dat[sample(nrow(dat), replace = TRUE), ]
    tryCatch(determine_cutpoint(b, marker, time, event, minprop, maxprop)$cutoff,
             error = function(e) NA_real_)
  })
  boot <- boot[is.finite(boot)]
  list(train_cutoff = cut_tr, test_hr = unname(exp(coef(cx))[1]),
       test_p = summary(cx)$coefficients[1, "Pr(>|z|)"],
       boot = boot, boot_ci = quantile(boot, c(.025, .5, .975)))
}

save_outputs <- function(dat, rc, surv, val, gene, outdir = "results") {
  dir.create(outdir, showWarnings = FALSE)
  km <- ggsurvplot(surv$km, data = dat, pval = TRUE, risk.table = TRUE,
                   conf.int = TRUE, legend.labs = c("Low", "High"),
                   palette = c("#2E6B9E", "#C0392B"),
                   xlab = "Years", ylab = "Overall survival",
                   title = sprintf("%s in HGSOC (TCGA)", gene))
  png(file.path(outdir, "kaplan_meier.png"), 720, 760); print(km); dev.off()
  png(file.path(outdir, "roc.png"), 700, 600)
  plot(rc$roc, print.auc = TRUE, print.thres = "best",
       main = "3-year mortality (secondary)"); dev.off()
  png(file.path(outdir, "bootstrap_cutpoint.png"), 720, 500)
  hist(val$boot, breaks = 30, col = "#2E6B9E", border = "white",
       main = "Bootstrap distribution of optimal cutpoint", xlab = gene)
  abline(v = val$boot_ci[c(1, 3)], lty = 2, lwd = 2); dev.off()
  sink(file.path(outdir, "cox_summary.txt"))
  cat("== Univariable ==\n");   print(summary(surv$cox_uni))
  cat("\n== Age-adjusted ==\n"); print(summary(surv$cox_adj))
  sink()
  invisible(outdir)
}

# ------------------------------- run ---------------------------------
GENE <- "CXCL13"

dat <- load_cohort(GENE)
cat(sprintf("Loaded %d patients | %d deaths | median FU %.1f y\n",
            nrow(dat), sum(dat$event), median(dat$time)))

cp <- determine_cutpoint(dat)
cat(sprintf("Cutpoint (maxstat): %.3f | corrected p = %s\n",
            cp$cutoff, signif(cp$corrected_p, 3)))

rc <- roc_cutpoint(dat)
cat(sprintf("ROC/Youden cutoff: %.3f | AUC = %.3f\n", rc$threshold, rc$auc))

dat  <- assign_groups(dat, cp$cutoff)
surv <- run_survival(dat, covariates = "age")
cat(sprintf("Univariable HR (High vs Low): %.2f\n", exp(coef(surv$cox_uni))[1]))

val <- validate_cutpoint(dat)
cat(sprintf("Train cutpoint: %.3f | TEST HR: %.2f (p = %s)\n",
            val$train_cutoff, val$test_hr, signif(val$test_p, 3)))
cat(sprintf("Bootstrap cutpoint: %.3f [%.3f, %.3f]\n",
            val$boot_ci[2], val$boot_ci[1], val$boot_ci[3]))

outdir <- save_outputs(dat, rc, surv, val, GENE)
cat("\nDone. Figures + tables in:", normalizePath(outdir), "\n")

library(GEOquery)

# Load data
# Download them (creates a GSE310998/ folder in working directory)
getGEOSuppFiles("GSE310998")

# Extract the tar archive
untar("GSE310998/GSE310998_RAW.tar", exdir = "GSE310998/raw/")

files <- list.files("GSE310998/raw/", full.names = TRUE)

# Read all files into a list
counts_list <- lapply(files, function(f) {
  read.table(gzfile(f), header = TRUE, row.names = 1, sep = "\t")
})

# Combine by column
counts <- do.call(cbind, counts_list)

# Checks counts
dim(counts)
head(counts)

gset <- getGEO("GSE310998", GSEMatrix = TRUE)[[1]]
pd <- pData(gset)

# Log transform to normalize counts data
log.counts <- log2(counts + 1)   # +1 avoids log(0) = -Inf

# Extract Gender data
pd$gender <- sub(".*gender:\\s*([MFmf]).*", "\\1", pd$characteristics_ch1.1)

# Encode gender as 0/1
pd$gender_binary <- ifelse(pd$gender == "M", 1, 0)

# Transpose so rows = samples, cols = genes
expr <- t(log.counts)
expr_non <- t(counts)

# Reassign rownames of pd to title
rownames(pd) <- pd$title

# Make sure sample order matches
expr <- expr[rownames(pd), ]
expr_non <- expr_non[rownames(pd), ]

# Run logistic regression for every gene. Function for both data sets
run_logistic_regression <- function(expr_input, pd) {
  
  genes   <- c()
  betas   <- c()
  ses     <- c()
  zs      <- c()
  pvalues <- c()
  
  for (gene in colnames(expr_input)) {
    tryCatch({
      df  <- data.frame(gender = pd$gender_binary, expression = expr_input[, gene])
      fit <- glm(gender ~ expression, data = df, family = binomial)
      coef_summary <- summary(fit)$coefficients
      
      if (!"expression" %in% rownames(coef_summary)) next
      
      genes   <- c(genes,   gene)
      betas   <- c(betas,   as.numeric(coef_summary["expression", "Estimate"]))
      ses     <- c(ses,     as.numeric(coef_summary["expression", "Std. Error"]))
      zs      <- c(zs,      as.numeric(coef_summary["expression", "z value"]))
      pvalues <- c(pvalues, as.numeric(coef_summary["expression", "Pr(>|z|)"]))
      
    }, error = function(e) NULL)
  }
  
  results_df       <- data.frame(gene = genes, beta = betas, se = ses,
                                 z = zs, p_value = pvalues,
                                 stringsAsFactors = FALSE)
  results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
  results_df       <- results_df[order(results_df$p_adj), ]
  
  return(results_df)
}

results.norm.df <- run_logistic_regression(expr, pd)
results.non.df  <- run_logistic_regression(t(counts)[rownames(pd), ], pd)

# Verify ordering
head(results.norm.df$p_value, 10)   # should be small values
tail(results.norm.df$p_value, 10)   # should be large values

# If not sorted, re-sort explicitly
results.norm.df <- results.norm.df[order(results.norm.df$p_value), ]
rownames(results.norm.df) <- NULL

# Check top 10
print(results.norm.df[1:10, ])

# Histogram of only the top 500
hist(results.norm.df$p_value[1:500], breaks = 50,
     main = "P-value Distribution: Top 500 Genes",
     xlab = "p-value")

# ///Multiple hypothesis testing with leave-n-out

# Filter to top 500 most significant genes from ranked results
keep.genes.norm <- head(results.norm.df$gene, 500)
keep.genes.non <- head(results.non.df$gene, 500)

cat("Genes kept for resampling:", length(keep.genes.norm), "\n")
cat("Genes kept for resampling:", length(keep.genes.non), "\n")

expr.filt.norm <- expr[, keep.genes.norm]
expr.filt.non <- expr[, keep.genes.non]

n_iter <- 100
all_pvalues <- vector("list", n_iter)

set.seed(42)  # set seed for reproducibility

# Function for running raw and normalized data frames
run_loo_resampling <- function(expr_input, pd, n_iter = 100) {
  
  set.seed(42)
  leave_out_indices <- sample(1:nrow(expr_input), n_iter, replace = TRUE)
  
  all_pvalues <- lapply(leave_out_indices, function(leave_out) {
    expr_sub <- expr_input[-leave_out, , drop = FALSE]
    pd_sub   <- pd[-leave_out, , drop = FALSE]
    
    iter_results <- lapply(colnames(expr_sub), function(gene) {
      tryCatch({
        df  <- data.frame(
          gender     = pd_sub$gender_binary,
          expression = as.numeric(expr_sub[, gene, drop = TRUE])
        )
        fit <- glm(gender ~ expression, data = df, family = binomial)
        coef_summary <- summary(fit)$coefficients
        if (!"expression" %in% rownames(coef_summary)) return(NULL)
        data.frame(gene = gene, p_value = coef_summary["expression", "Pr(>|z|)"])
      }, error = function(e) NULL)
    })
    
    iter_results <- Filter(Negate(is.null), iter_results)
    do.call(rbind, iter_results)
  })
  
  all_pvalues_df <- do.call(rbind, lapply(seq_along(all_pvalues), function(i) {
    df <- all_pvalues[[i]]
    df$iteration <- i
    df
  }))
  
  avg_pval   <- aggregate(p_value ~ gene, data = all_pvalues_df, FUN = mean)
  var_pval   <- aggregate(p_value ~ gene, data = all_pvalues_df, FUN = var)
  summary_df <- merge(avg_pval, var_pval, by = "gene", suffixes = c("_mean", "_var"))
  
  # Forces numeric to avoid a named vector issue
  summary_df <- summary_df[sort.int(summary_df$p_value_mean, index.return = TRUE)$ix, ]
  rownames(summary_df) <- NULL
  
  return(list(
    summary    = summary_df,
    iterations = all_pvalues_df
  ))
}

loo.results.norm <- run_loo_resampling(expr.filt.norm, pd)
loo.results.non  <- run_loo_resampling(expr.filt.non, pd)

stability.norm <- loo.results.norm$summary
stability.non  <- loo.results.non$summary

# Manual sort of results
stability.norm <- stability.norm[sort.int(stability.norm$p_value_mean, index.return = TRUE)$ix, ]
stability.non  <- stability.non[sort.int(stability.non$p_value_mean, index.return = TRUE)$ix, ]

top10.norm <- head(stability.norm, 10)
top10.non  <- head(stability.non, 10)

#Reports the top 10 genes from each data set
top10.norm
top10.non

# Visualization
# Label top 10 most significant stable genes
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))  # wider margins

# Normalized
top10.norm <- head(stability.norm, 10)
plot(
  x    = stability.norm$p_value_mean,
  y    = stability.norm$p_value_var,
  pch  = 20,
  col  = ifelse(stability.norm$p_adj < 0.05, "red", "grey"),
  xlab = "Mean Raw P-value",
  ylab = "Variance of P-value",
  main = "Gene Stability: Normalized",
  xlim = c(0, 1),
  ylim = c(0, max(stability.norm$p_value_var) * 1.2)  # add headroom for labels
)
abline(v = 0.05, col = "blue", lty = 2)
text(top10.norm$p_value_mean, top10.norm$p_value_var, labels = top10.norm$gene, 
     pos = 4, cex = 0.6, col = "darkred")  # pos=4 puts labels to the right

# Non-normalized
top10.non <- head(stability.non, 10)
plot(
  x    = stability.non$p_value_mean,
  y    = stability.non$p_value_var,
  pch  = 20,
  col  = ifelse(stability.non$p_adj < 0.05, "red", "grey"),
  xlab = "Mean Raw P-value",
  ylab = "Variance of P-value",
  main = "Gene Stability: Non-Normalized",
  xlim = c(0, 1),
  ylim = c(0, max(stability.non$p_value_var) * 1.2)
)
abline(v = 0.05, col = "blue", lty = 2)
text(top10.non$p_value_mean, top10.non$p_value_var, labels = top10.non$gene, 
     pos = 4, cex = 0.6, col = "darkred")

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))  # reset to defaults


# /// Question 2: Survival Analysis with heart transplant data

library(survival)

# txt file load from directory
heart <- read.table("Heart Transplant Data.txt", 
                    header = FALSE, 
                    skip = 10,          # skip the header lines
                    col.names = c("age", "status", "time"))

# Creates survival object. status = 1 for event (death), 0 for censored
surv_obj <- Surv(time = heart$time, event = heart$status)

# Fit Kaplan-Meier estimator
km_fit <- survfit(surv_obj ~ 1)

# Summary
summary(km_fit)

# Plot the survival curve
plot(km_fit, 
     xlab = "Time (days)", 
     ylab = "Survival Probability",
     main = "Kaplan-Meier Survival Curve\nHeart Transplant Patients")

print(km_fit)  # shows median survival time etc.

# Fit Cox model
cox_model <- coxph(surv_obj ~ age, data = heart)

# Summary (hazard ratios, p-values, etc.)
summary(cox_model)

# Check proportional hazards assumption
cox.zph(cox_model)

# Plot hazard ratio effect of age
plot(survfit(cox_model), 
     xlab = "Time (days)", 
     ylab = "Survival Probability",
     col = c("blue", "red"),  # you can stratify if needed
     main = "Cox Model Survival Curves by Age")






















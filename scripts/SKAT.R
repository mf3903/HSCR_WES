library(dplyr)
library(SKAT)

# ---- Argument Parsing ----
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) stop("Usage: Rscript SKAT.R <file_in> <file_out> <file_covariate> <file_pheno>")

file_in <- args[1]
file_out <- args[2]
file_covariate <- args[3]
file_pheno <- args[4]

# ---- Data Import ----
gt <- read.table(file_in, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
Z <- gt
X <- read.table(file_covariate, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
pheno <- read.table(file_pheno, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

y.b <- pheno$disease
X <- data.matrix(X)
Z <- data.matrix(Z)

# ---- SKAT Null Model ----
obj <- SKAT_Null_Model(y.b ~ X, out_type = "D")

# ---- SKAT Tests ----
out1 <- SKATBinary_Robust(Z, obj, kernel = "linear.weighted", impute.method = "bestguess", method = "SKATO", weights.beta = c(1, 25))
out2 <- SKATBinary_Robust(Z, obj, kernel = "linear.weighted", impute.method = "bestguess", method = "SKATO", weights.beta = c(0.5, 0.5))

# ---- Output Results ----
out_df <- data.frame(out1 = out1$p.value, out2 = out2$p.value)
write.table(out_df, file_out, row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)




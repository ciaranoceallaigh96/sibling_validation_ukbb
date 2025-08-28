library(data.table)

file_path  <- args[1] #GWAS Summary statistics
prefix <- args[2] #save out
log_path <- paste0(prefix, "_genomic_control.log")

gwas_data <- read.csv(file_path, header=FALSE, sep="")

lambda_gc <- median(qchisq(1-gwas_data$P,1))/qchisq(0.5,1)
print(paste("Genomic Inflation Factor (Lambda GC):", round(lambda_gc, 3)))

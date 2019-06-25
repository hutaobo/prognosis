library(EnsDb.Hsapiens.v86)
all_gene <- genes(EnsDb.Hsapiens.v86)
all_gene_id <- all_gene$gene_id[all_gene$gene_biotype == "lincRNA"]

library(SummarizedExperiment)
study_type <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")
gene_id <- as.character(gene_id)

p_val <- c()
mean_diff <- c()
for (study in study_type) {
  load(paste0("~/R/GDCdata/", study, "-Gene-HTSeq_FPKM_UQ.RData"), verbose = TRUE)
  
  normal_gr <- data$definition == "Solid Tissue Normal"
  tumor_gr <- data$definition == "Primary solid Tumor"
  patient_id <- stringr::str_sub(colnames(data), 1, 12)
  
  normal_data <- data[, normal_gr]
  tumor_data <- data[, tumor_gr]
  matched_tumor_data <- tumor_data[, match(patient_id[normal_gr], patient_id[tumor_gr])]
  
  normal_exp <- log2(assay(normal_data) + 1)
  matched_tumor_exp <- log2(assay(matched_tumor_data) + 1)
  all_p_val <- sapply(1:nrow(normal_exp), function(x) t.test(normal_exp[x, ], matched_tumor_exp[x, ], paired = TRUE)$p.value)
  
  mean_normal_exp <- rowMeans(normal_exp, na.rm = TRUE)
  mean_matched_tumor_exp <- rowMeans(matched_tumor_exp, na.rm = TRUE)
  all_mean_diff <- mean_normal_exp - mean_matched_tumor_exp
  
  p_val <- cbind(p_val, all_p_val)
  mean_diff <- cbind(mean_diff, all_mean_diff)
}
rownames(p_val) <- rownames(mean_diff)

is_elevated <- mean_diff < 0 & p_val < 5e-8
is_elevated_sum <- rowSums(is_elevated)
table(is_elevated_sum)
is_elevated_lncR <- is_elevated[rownames(is_elevated) %in% all_gene_id, ]
is_elevated_lncR_sum <- rowSums(is_elevated_lncR)
table(is_elevated_lncR_sum)

is_reduced <- mean_diff > 0 & p_val < 5e-8
is_reduced_sum <- rowSums(is_reduced)
table(is_reduced_sum)
is_reduced_lncR <- is_reduced[rownames(is_reduced) %in% all_gene_id, ]
is_reduced_lncR_sum <- rowSums(is_reduced_lncR)
table(is_reduced_lncR_sum)

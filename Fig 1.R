# PRDM2高低表达对预后的影响

diff_surv_2geneIndex <- function(gene_id_1, gene_id_2, pdf_path, width = 9, height = 7, ThreshTop = 0.67, ThreshDown = 0.33) {
  library(SummarizedExperiment)
  study_type <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PRAD", "TCGA-STAD", "TCGA-THCA", "TCGA-UCEC")
  gene_id <- as.character(gene_id)
  
  plot_list <- list()
  for (study in study_type) {
    load(paste0("/Volumes/BioInfo/R/GDCdata/", study, "-Gene-HTSeq_FPKM_UQ.RData"), verbose = TRUE)
    data <- data[, data$definition == "Primary solid Tumor"]
    gene_id_values_1 <- assay(data[gene_id_1, ])
    gene_id_values_2 <- assay(data[gene_id_2, ])
    gene_id_values <- gene_id_values_1 / gene_id_values_2
    gene_id_values_top <- as.numeric(quantile(as.numeric(gene_id_values), ThreshTop)[1])
    gene_id_values_down <- as.numeric(quantile(as.numeric(gene_id_values), ThreshDown)[1])
    data$dex <- ifelse(gene_id_values > gene_id_values_top, "high", ifelse(gene_id_values <= gene_id_values_down, "low", NA))
    data <- data[, !is.na(data$dex)]
    
    ttime <- ifelse(data$vital_status == "alive", data$days_to_last_follow_up, data$days_to_death)
    status <- ifelse(data$vital_status == "dead", TRUE, FALSE)
    Pheno <- data.frame(ttime = ttime, status = status, dex = data$dex, stringsAsFactors = FALSE)
    Pheno <- rbind(Pheno[Pheno$dex == "high", ], Pheno[Pheno$dex == "low", ])
    Pheno$dex <- factor(Pheno$dex, levels = c("high", "low"))
    
    library(survival)
    surv <- survfit(Surv(ttime, status) ~ dex, data = Pheno)
    
    tabSurv_pvalue <- tryCatch({
      tabSurv <- survdiff(Surv(ttime, status) ~ dex, data = Pheno)
      tabSurv_chis <- unlist(tabSurv)$chisq
      tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
    }, error = function(e) {
      return(Inf)
    })
    
    library(GGally)
    p_surv <- ggsurv(surv, order.legend = FALSE) +
      guides(fill = guide_legend("Type")) +
      labs(title = study) +
      theme(legend.position = c(0.8, 0.8)) +
      theme(legend.title = element_blank()) +
      theme(plot.title = element_text(hjust = 0.5)) +
      annotate("text", x = 0, y = 1, label = paste0("                    p = ", format(tabSurv_pvalue, scientific = TRUE, digits = 2)))
    
    plot_list[[study]] <- p_surv
  }
  
  library(ggtree)
  pdf(pdf_path, width, height)
  ggtree::multiplot(plotlist = plot_list, ncol = 4)
  dev.off()
}

# C1orf112 & AKT3
diff_surv_2geneIndex('ENSG00000000460', 'ENSG00000117020', "/Volumes/BioInfo/RProjects/项目/prognosis/result/Fig 1.pdf")
# NFRKB & AKT3
diff_surv_2geneIndex('ENSG00000170322', 'ENSG00000117020', "/Volumes/BioInfo/RProjects/项目/prognosis/result/Fig 2.pdf")
# ERCC1 & AKT3
diff_surv_2geneIndex('ENSG00000012061', 'ENSG00000117020', "/Volumes/BioInfo/RProjects/项目/prognosis/result/Fig 3.pdf")

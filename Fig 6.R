library(clusterProfiler)
library(org.Hs.eg.db)

genelist <- c('STIL', 'KIF4A', 'KIF23', 'PRC1', 'KIF14', 'KIF2C', 'CCNA2', 'KIF11', 'AURKB', 'MELK')
ids <- bitr(genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ego <- enrichGO(OrgDb = "org.Hs.eg.db", gene = ids$ENTREZID, ont = "BP", pvalueCutoff = 0.01, readable = TRUE)
pdf('~/RProjects/项目/prognosis/result/Fig 6.pdf', width = 10, height = 5)
dotplot(ego, showCategory = 10)
dev.off()

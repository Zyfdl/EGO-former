library(graphite)
ls("package:graphite")
databases <- pathwayDatabases()
kpaths <- pathways("hsapiens", "kegg")

insulin_pathway_converted <-kpaths
insulin_pathway_converted <- convertIdentifiers(kpaths, "symbol")

pathway_edges_list <- list()

# 使用循环遍历每个通路
for (i in seq_along(insulin_pathway_converted)) {
  pathway_edges_list[[i]] <- edges(insulin_pathway_converted[[i]], which = "proteins")
}

result_df <- do.call(rbind, pathway_edges_list)
library(dplyr)
result_df <- distinct(result_df)
result_df<-result_df[result_df$src_type=='SYMBOL',]
result_df<-result_df[result_df$dest_type=='SYMBOL',]
save(result_df,file = "hsa_kegg_all.Rdata")
table(result_df$type)
filtered_df <- result_df[grepl("activation|inhibition|expression|repression", result_df$type), ]
filtered_df<-filtered_df[filtered_df$direction=='directed',]

library(dplyr)
filtered_df <- distinct(filtered_df)
table(filtered_df$type)
mmu_symbol<-filtered_df[filtered_df$src_type=='SYMBOL',]
mmu_sym_sym<-mmu_symbol[mmu_symbol$dest_type=='SYMBOL',]
hsa_regulation <- mmu_sym_sym[, c("src", "dest")]
library(dplyr)
hsa_regulation <- distinct(hsa_regulation)
write.csv(hsa_regulation, "hsa_gene_regulation_kegg.csv", row.names = FALSE)
library(org.Hs.eg.db)
hsa_regulation$src <- mapIds(org.Hs.eg.db, keys = hsa_regulation$src, column = "ENTREZID", keytype = "SYMBOL")
hsa_regulation$dest <- mapIds(org.Hs.eg.db, keys = hsa_regulation$dest, column = "ENTREZID", keytype = "SYMBOL")
library(dplyr)
hsa_regulation <- distinct(hsa_regulation)
write.csv(hsa_regulation, "hsa_gene_regulation_kegg_NCBI.csv", row.names = FALSE)
mmu_reverse<-mmu_sym_sym
names(mmu_reverse)<-c("dest_type","dest","src_type","src","direction","type")
mmu_del <- mmu_sym_sym %>%
  anti_join(mmu_reverse, by = c("src","dest"))
save(mmu_del,file = "hsa_kegg_del_reverse.Rdata")
load("hsa_kegg_del_reverse.Rdata")
inhibition<-mmu_del[grepl("inhibition|repression", mmu_del$type), ]
activation<-mmu_del[grepl("activation|expression", mmu_del$type), ]
only_act <- activation %>%
  anti_join(inhibition, by = c("src","dest"))
only_inh<- inhibition %>%
  anti_join(activation, by = c("src","dest"))
only_inh$type<-'1'
only_act$type<-'0'

regulation_data<-rbind(only_act,only_inh)
regulation_data <- regulation_data[, c("src", "dest","type")]
write.table(regulation_data, file = "hsa_kegg_act_inh_symbol.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
regulation_data$src <- mapIds(org.Hs.eg.db, keys = regulation_data$src, column = "ENTREZID", keytype = "SYMBOL")
regulation_data$dest <- mapIds(org.Hs.eg.db, keys = regulation_data$dest, column = "ENTREZID", keytype = "SYMBOL")
write.table(regulation_data, file = "hsa_kegg_act_inh_NCBI.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


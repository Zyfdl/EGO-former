
library(graphite)
names(reactome)[1:3]
reactome_pathways <- pathways(species = "hsapiens", database = "reactome")
insulin_pathway_converted <- convertIdentifiers(reactome_pathways, "symbol")
pathway_edges <- edges(insulin_pathway_converted[[1]], which = "proteins")

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
save(result_df,file = "hsa_reactome_all.Rdata")
table(result_df$type)
filtered_df <- result_df[grepl("ACTIVATION|INHIBITION", result_df$type), ]
filtered_df<-filtered_df[filtered_df$direction=='directed',]
library(dplyr)
filtered_df <- distinct(filtered_df)
mmu_out <- filtered_df[grepl("Out", filtered_df$type), ]
mmu_in  <- filtered_df[grepl("In", filtered_df$type), ]
names(mmu_in)<-c("dest_type","dest","src_type","src","direction","type")
mmu_in <- dplyr::select(mmu_in,src_type,src,dest_type,dest,direction,type)
mmu_out_d <- mmu_out %>%
  anti_join(mmu_in, by = c("src","dest"))
new_data<-rbind(mmu_out_d,mmu_in)
reverse<-new_data
names(reverse)<-c("dest_type","dest","src_type","src","direction","type")
reverse <- dplyr::select(reverse,src_type,src,dest_type,dest,direction,type)
del_reverse<- new_data %>%
  anti_join(reverse, by = c("src","dest"))
table(del_reverse$type)
save(del_reverse,file = "hsa_reactome_del_reverse.Rdata")
new_data <- dplyr::select(new_data,src,dest)
write.csv(new_data, "hsa_gene_regulation_reactome_symbol.csv", row.names = FALSE)
library(org.Hs.eg.db)
new_data$src <- mapIds(org.Hs.eg.db, keys = new_data$src, column = "ENTREZID", keytype = "SYMBOL")
new_data$dest <- mapIds(org.Hs.eg.db, keys = new_data$dest, column = "ENTREZID", keytype = "SYMBOL")
library(dplyr)
new_data <- distinct(new_data)
write.csv(new_data, "hsa_gene_regulation_reactome_NCBI.csv", row.names = FALSE)

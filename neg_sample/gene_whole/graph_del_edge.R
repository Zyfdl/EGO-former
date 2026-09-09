load(file = "use_all_negative_gene_whole.Rdata")
library(readr)
pre_pair <- read_csv("../pre_pair.csv")
names(pre_pair)<-c("gene2","gene1")
pre_pair <- data.frame(sapply(pre_pair, as.character))
library(dplyr)
del_possible <- del_possible %>%
  anti_join(pre_pair, by = c("gene1","gene2"))
library(readr)
true_pair <- read_csv("../true_pair.csv")
names(true_pair)<-c("gene2","gene1")
true_pair <- data.frame(sapply(true_pair, as.character))
del_possible <- del_possible %>%
  anti_join(true_pair, by = c("gene1","gene2"))

load(file ="../htftarget.Rdata" )
names(htftarget)<-c("gene2","gene1")
htftarget <- data.frame(sapply(htftarget, as.character))
del_possible <- del_possible %>%
  anti_join(htftarget, by = c("gene1","gene2"))


library(readr)
hsa_kegg <- read_csv("D:/rworkspace/Reactome/hsa/hsa_gene_regulation_kegg_NCBI.csv")
names(hsa_kegg)<-c("gene2","gene1")
hsa_kegg <- data.frame(sapply(hsa_kegg, as.character))
del_possible <- del_possible %>%
  anti_join(hsa_kegg, by = c("gene1","gene2"))

library(readr)
hsa_reactome <- read_csv("D:/rworkspace/Reactome/hsa/hsa_gene_regulation_reactome_NCBI.csv")
names(hsa_reactome)<-c("gene2","gene1")
hsa_reactome <- data.frame(sapply(hsa_reactome, as.character))
del_possible <- del_possible %>%
  anti_join(hsa_reactome, by = c("gene1","gene2"))


load(file = "possible_all.Rdata")
names(possible_all)<-c("gene2","gene1")
possible_all <- data.frame(sapply(possible_all, as.character))
del_possible <- del_possible %>%
  anti_join(possible_all, by = c("gene1","gene2"))

reverse_data<-del_possible
names(reverse_data)<-c("gene2","gene1")
possible_one <- del_possible %>%
  anti_join(reverse_data, by = c("gene1","gene2"))

possible_two <- del_possible %>%
  semi_join(reverse_data, by = c("gene1","gene2"))
# 创建一个新列，将每行的a和b值排序，并合并为一个字符串
possible_two$sorted_pair <- apply(possible_two, 1, function(x) paste(sort(c(x['gene1'], x['gene2'])), collapse = "-"))

# 去除排序后的对中重复的行
df_unique <- possible_two[!duplicated(possible_two$sorted_pair), ]

# 删除辅助的sorted_pair列
df_unique$sorted_pair <- NULL

# 查看结果
print(df_unique)
del_possible<-rbind(possible_one,df_unique)
save(del_possible,file = "graph_for_irr.Rdata")
target_irr<-unique(del_possible$gene2)
source_irr<-unique(del_possible$gene1)

combined_data <- del_possible[sample(nrow(del_possible)), ]
irr_one <- combined_data %>%
  group_by(gene1) %>%
  filter(if (n() >= 17) row_number() <=17  else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(irr_one$gene2)
source_irr<-unique(irr_one$gene1)

irr_two <- combined_data %>%
  group_by(gene2) %>%
  filter(if (n() >= 2) row_number() <=2 else TRUE) %>%  # 条件过滤
  ungroup()
irr_data<-rbind(irr_one,irr_two)
library(dplyr)
irr_data <- distinct(irr_data)

irr_data <- irr_data[sample(nrow(irr_data)), ]
irr_data <- irr_data[1:14586, ]
target_irr<-unique(irr_data$gene2)
source_irr<-unique(irr_data$gene1)


library(readr)
hsa_index <- read_csv("hsa_transcript_regulation_edge_index.csv")
hsa_reg<-hsa_index[hsa_index$trans_type=='regulation',]
hsa_irr<-irr_data
hsa_irr$trans_type<-'irrelation'
hsa_irr$data_type<-hsa_reg$data_type
load(file ="all_gene.Rdata" )
all_gene$index <- seq(0, by = 1, length.out = nrow(all_gene))
names(all_gene)<-c("gene1","source_type","source_index")
library(dplyr)
merged_df <- left_join(hsa_irr, all_gene, by = "gene1")
names(all_gene)<-c("gene2","target_type","target_index")
all_map <- left_join(merged_df, all_gene, by = "gene2")
hsa_irr<-all_map[, c("source_index", "target_index","trans_type","data_type")]
hsa_index<-rbind(hsa_reg,hsa_irr)
write.csv(hsa_index, "hsa_transcript_regulation_edge_index.csv", row.names = FALSE)
reverse_reg<-hsa_reg
names(reverse_reg)<-c("target_index", "source_index","trans_type","data_type")
library(dplyr)
semi_reg <- hsa_reg %>%
  semi_join(reverse_reg, by = c("source_index", "target_index"))

library(org.Hs.eg.db)
all_gene$go_term <- mapIds(org.Hs.eg.db, keys = all_gene$target_gene, column = "GO", keytype = "ENTREZID", multiVals = "list")
table(all_gene$go_term)
library(plyr)
go_terms_df <- data.frame(
  gene = rep(all_gene$target_gene, sapply(all_gene$go_term, length)),
  go_term = unlist(all_gene$go_term)
)
unique(go_terms_df$go_term)
go_terms_df$go_term <- ifelse(is.na(go_terms_df$go_term), "GO:", go_terms_df$go_term)

go_terms_df <- na.omit(go_terms_df)
library(dplyr)
go_terms_df <- distinct(go_terms_df)
write.csv(go_terms_df, "gene_go_terms.csv", row.names = FALSE)
library(readr)
gene_go_terms <- read_csv("gene_go_terms.csv")
 
library(readr)
possible_low <- read_csv("human_possible_low.csv")
possible_low<-possible_low[, c("source_name", "Gene symbol")]
library(org.Hs.eg.db)
possible_low$source_name <- mapIds(org.Hs.eg.db, keys = possible_low$source_name, column = "ENTREZID", keytype = "UNIPROT")
possible_low$`Gene symbol` <- mapIds(org.Hs.eg.db, keys = possible_low$`Gene symbol`, column = "ENTREZID", keytype = "SYMBOL")
# 假设 df 是你的数据框
possible_low <- na.omit(possible_low)
names(possible_low)<-c("source_name","target_name")
save(possible_low,file = "possible_low.Rdata")


library(readr)
possible_mid <- read_csv("human_possible_mid.csv")
possible_mid<-possible_mid[, c("source_name", "Gene symbol")]
library(org.Hs.eg.db)
possible_mid$source_name <- mapIds(org.Hs.eg.db, keys = possible_mid$source_name, column = "ENTREZID", keytype = "UNIPROT")
possible_mid$`Gene symbol` <- mapIds(org.Hs.eg.db, keys = possible_mid$`Gene symbol`, column = "ENTREZID", keytype = "SYMBOL")
# 假设 df 是你的数据框
possible_mid <- na.omit(possible_mid)
names(possible_mid)<-c("source_name","target_name")
save(possible_mid,file = "possible_mid.Rdata")


library(readr)
possible_high <- read_csv("human_possible_high.csv")
possible_high<-possible_high[, c("source_name", "Gene symbol")]
library(org.Hs.eg.db)
possible_high$source_name <- mapIds(org.Hs.eg.db, keys = possible_high$source_name, column = "ENTREZID", keytype = "UNIPROT")
possible_high$`Gene symbol` <- mapIds(org.Hs.eg.db, keys = possible_high$`Gene symbol`, column = "ENTREZID", keytype = "SYMBOL")
# 假设 df 是你的数据框
possible_high <- na.omit(possible_high)
names(possible_high)<-c("source_name","target_name")
save(possible_high,file = "possible_high.Rdata")


load(file ="possible_high.Rdata" )
load(file ="possible_low.Rdata" )
load(file ="possible_mid.Rdata" )
possible_gtrd<-rbind(possible_high,possible_low,possible_mid)
library(dplyr)
possible_gtrd <- distinct(possible_gtrd)
possible_high <- distinct(possible_high)
save(possible_gtrd,file = "possible_gtrd.Rdata")


load(file = "../del_possible.Rdata")
library(readr)
hsa_kegg <- read_csv("D:/rworkspace/Reactome/hsa/hsa_gene_regulation_kegg_NCBI.csv")
names(hsa_kegg)<-c("gene1","gene2")
hsa_kegg <- data.frame(sapply(hsa_kegg, as.character))
del_possible <- del_possible %>%
  anti_join(hsa_kegg, by = c("gene1","gene2"))
library(readr)
hsa_reactome <- read_csv("D:/rworkspace/Reactome/hsa/hsa_gene_regulation_reactome_NCBI.csv")
names(hsa_reactome)<-c("gene1","gene2")
hsa_reactome <- data.frame(sapply(hsa_reactome, as.character))
del_possible <- del_possible %>%
  anti_join(hsa_reactome, by = c("gene1","gene2"))
load(file ="possible_gtrd.Rdata" )
names(possible_gtrd)<-c("gene1","gene2")
possible_gtrd <- data.frame(sapply(possible_gtrd, as.character))
del_possible <- del_possible %>%
  anti_join(possible_gtrd, by = c("gene1","gene2"))
library(readr)
biao_zhun <- read_csv("D:/rworkspace/trans_pair/human/hsa_transcript_regulation_part_pre_test.csv")
biao_zhun<-biao_zhun[, c("source_gene", "target_gene")]
names(biao_zhun)<-c("gene1","gene2")
biao_zhun <- data.frame(sapply(biao_zhun, as.character))
del_possible<- del_possible %>%
  semi_join(biao_zhun, by = c("gene1"))
del_possible<- del_possible %>%
  semi_join(biao_zhun, by = c("gene2"))
target_irr<-unique(biao_zhun$gene2)
source_irr<-unique(biao_zhun$gene1)
target_irr<-unique(del_possible$gene2)
source_irr<-unique(del_possible$gene1)
save(del_possible,file = "use_all_negative.Rdata")
combined_data <- del_possible[sample(nrow(del_possible)), ]
irr_one <- combined_data %>%
  group_by(gene1) %>%
  filter(if (n() >= 12) row_number() <=12  else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(irr_one$gene2)
source_irr<-unique(irr_one$gene1)

irr_two <- combined_data %>%
  group_by(gene2) %>%
  filter(if (n() >= 3) row_number() <=3 else TRUE) %>%  # 条件过滤
  ungroup()
irr_data<-rbind(irr_one,irr_two)
library(dplyr)
irr_data <- distinct(irr_data)
target_irr<-unique(irr_data$gene2)
source_irr<-unique(irr_data$gene1)
irr_data <- irr_data[sample(nrow(irr_data)), ]
irr_data <- irr_data[1:14586, ]
library(readr)
true_pair <- read_csv("D:/rworkspace/trans_pair/human/true_pair.csv")
names(true_pair)<-c("gene1","gene2")
target_irr<-unique(true_pair$gene2)
source_irr<-unique(true_pair$gene1)
true_pair <- distinct(true_pair)
true_pair$type<-'0'
irr_data$type<-'1'
all_data<-rbind(true_pair,irr_data)
write.table(all_data, file = "hsa_transcript_regulation_depart_gtrd.txt", sep = "\t", row.names = FALSE, col.names = FALSE)



library(readr)
hsa_data<- read_csv("hsa_transcript_regulation_depart_gtrd.csv")
positive<-hsa_data[hsa_data$trans_type=='regulation',]
negative<-hsa_data[hsa_data$trans_type=='irrelation',]
positive <- positive[sample(nrow(positive)), ]
negative <- negative[sample(nrow(negative)), ]
new_column <- c(rep('test', 2917), rep('val', 2917), rep('train', 8752))
positive$data_type<-new_column
negative$data_type<-new_column
com_data<-rbind(positive,negative)
write.csv(com_data, "hsa_transcript_regulation_depart_gtrd.csv", row.names = FALSE)


library(readr)
possible_all <- read_csv("human_possible_all.csv")
possible_all<-possible_all[, c("source_name", "Gene symbol")]
library(org.Hs.eg.db)
possible_all$source_name <- mapIds(org.Hs.eg.db, keys = possible_all$source_name, column = "ENTREZID", keytype = "UNIPROT")
possible_all$`Gene symbol` <- mapIds(org.Hs.eg.db, keys = possible_all$`Gene symbol`, column = "ENTREZID", keytype = "SYMBOL")
# 假设 df 是你的数据框
possible_all <- na.omit(possible_all)
names(possible_all)<-c("source_name","target_name")
save(possible_all,file = "possible_all.Rdata")

load(file = "use_all_negative.Rdata")
load(file = "possible_all.Rdata")
names(possible_all)<-c("gene1","gene2")
possible_all <- data.frame(sapply(possible_all, as.character))
del_possible <- del_possible %>%
  anti_join(possible_all, by = c("gene1","gene2"))
target_irr<-unique(del_possible$gene2)
source_irr<-unique(del_possible$gene1)
library(readr)
biao_zhun <- read_csv("hsa_transcript_regulation_depart_gtrd.csv")
biao_zhun<-biao_zhun[, c("source_gene", "target_gene")]
names(biao_zhun)<-c("gene1","gene2")
biao_zhun <- data.frame(sapply(biao_zhun, as.character))
target_irr<-unique(biao_zhun$gene2)
source_irr<-unique(biao_zhun$gene1)
load("D:/rworkspace/protein/biogrid.Rdata")
names(biogrid)<-c("gene1","gene2")
library(dplyr)
del_possible <- del_possible %>%
  anti_join(biogrid, by = c("gene1","gene2"))
save(del_possible,file = "use_all_negative_gene_whole.Rdata")
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
target_irr<-unique(irr_data$gene2)
source_irr<-unique(irr_data$gene1)
irr_data <- irr_data[sample(nrow(irr_data)), ]
irr_data <- irr_data[1:14586, ]
irr_data <- irr_data[sample(nrow(irr_data)), ]
hsa_data <- read_csv("hsa_transcript_regulation_depart_gene_whole.csv")
hsa_reg<-hsa_data[hsa_data$trans_type=='regulation',]
hsa_reg$data<-NULL
hsa_reg$data_type<-NULL
load(file ="all_gene.Rdata" )
all_gene$index <- seq(0, by = 1, length.out = nrow(all_gene))
names(irr_data)<-c("source_gene","target_gene")
irr_data$trans_type<-'irrelation'
hsa_data<-rbind(hsa_reg,irr_data)
names(all_gene)<-c("source_gene","source_type","source_index")
hsa_data <- data.frame(sapply(hsa_data, as.character))
library(dplyr)
merged_df <- left_join(hsa_data, all_gene, by = "source_gene")
names(all_gene)<-c("target_gene","target_type","target_index")
all_map <- left_join(merged_df, all_gene, by = "target_gene")
all_map<-all_map[, c("source_gene","target_gene","source_index", "target_index","trans_type")]
write.csv(all_map, "hsa_edge_index.csv", row.names = FALSE)


library(readr)
all_data <- read_csv("hsa_edge_index.csv")
new_column <- c(rep('test', 2917), rep('val', 2917), rep('train', 8752),rep('test', 2917), rep('val', 2917), rep('train', 8752))
all_data$data_type<-new_column
write.table(all_data, file = "hsa_bar_image.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

##插曲
library(readr)
all_data <- read_delim("hsa_bar_image.txt", 
                            delim = "\t", escape_double = FALSE, 
                            col_names = FALSE, trim_ws = TRUE)
new_column<-c(rep('one', 2917), rep('two', 2917), rep('three', 2918),rep('four', 2917),rep('five', 2917),rep('one', 2917), rep('two', 2917), rep('three', 2918),rep('four', 2917),rep('five', 2917))
all_data$X6<-NULL
all_data$X6<-new_column
write.table(all_data, file = "hsa_bar_image_5k.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
library(readr)
has_negative <- read_csv("has_edge_index_all_negative.csv")
has_negative$trans_type<-'irrelation'
has_negative <- has_negative[1:14586, ]
unique(has_negative$target_gene)
new_column<-c(rep('one', 2917), rep('two', 2917), rep('three', 2918),rep('four', 2917),rep('five', 2917))
has_negative$data_type<-new_column
names(has_negative)<-c("X1","X2","X3","X4","X5","X6")
hsa_positive<-all_data[all_data$X5=='regulation',]
hsa_random<-rbind(hsa_positive,has_negative)
write.table(hsa_random, file = "hsa_random_image_5k.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


all_data$source_index<-NULL
all_data$target_index<-NULL
all_data$trans_type[ all_data$trans_type == "regulation" ] <- '0'
all_data$trans_type[ all_data$trans_type == "irrelation" ] <- '1'
library(dplyr)
all_data <- distinct(all_data)
write.table(all_data, file = "hsa_tflink_image.txt", sep = "\t", row.names = FALSE, col.names = FALSE)













library(readr)
true_pair <- read_csv("D:/rworkspace/trans_pair/human/true_pair.csv")
names(true_pair)<-c("gene1","gene2")
target_irr<-unique(true_pair$gene2)
source_irr<-unique(true_pair$gene1)
true_pair <- distinct(true_pair)
true_pair$type<-'0'
irr_data$type<-'1'
all_data<-rbind(true_pair,irr_data)
write.table(all_data, file = "hsa_transcript_regulation_depart_gene_whole.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


library(readr)
hsa_data <- read_csv("hsa_transcript_regulation_depart_gene_whole.csv")
positive<-hsa_data[hsa_data$trans_type=='regulation',]
negative<-hsa_data[hsa_data$trans_type=='irrelation',]
positive <- positive[sample(nrow(positive)), ]
negative <- negative[sample(nrow(negative)), ]
new_column <- c(rep('test', 2917), rep('val', 2917), rep('train', 8752))
positive$data_type<-new_column
negative$data_type<-new_column
com_data<-rbind(positive,negative)
write.csv(com_data, "hsa_transcript_regulation_depart_gene_whole.csv", row.names = FALSE)

all_data <- dplyr::select(com_data,source_gene,target_gene,trans_type,data_type)
all_data$trans_type[ all_data$trans_type == "regulation" ] <- '0'
all_data$trans_type[ all_data$trans_type == "irrelation" ] <- '1'
library(dplyr)
all_data <- distinct(all_data)
write.table(all_data, file = "hsa_image_tf_link.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


library(readr)
hsa_data <- read_csv("hsa_transcript_regulation_depart_gene_whole.csv")
trans_factor<-unique(hsa_data$source_gene)
trans_factor<-as.data.frame(trans_factor)
trans_factor$type<-NULL

library(readr)
all_TF <- read_delim("Homo_sapiens_TF.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
all_TF<-unique(all_TF$Entrez_ID)
all_TF<-as.data.frame(all_TF)
names(all_TF)<-c("trans_factor")
load(file ="possible_gtrd.Rdata" )
gtrd_tf<-unique(possible_gtrd$source_name)
gtrd_tf<-as.data.frame(gtrd_tf)
names(gtrd_tf)<-c("trans_factor")
all_TF<-rbind(gtrd_tf,all_TF,trans_factor)
library(dplyr)
all_TF <- distinct(all_TF)
type(all_TF$trans_factor)
target_gene<-unique(hsa_data$target_gene)
target_gene<-as.data.frame(target_gene)
names(target_gene)<-c("trans_factor")
target_gene <- data.frame(sapply(target_gene, as.character))
add_TF<- all_TF %>%
  semi_join(target_gene, by = c("trans_factor"))
other_type<-target_gene %>%
  anti_join(all_TF, by = c("trans_factor"))
trans_factor<-rbind(trans_factor,add_TF)
library(dplyr)
trans_factor <- distinct(trans_factor)
trans_factor$type<-'0'
other_type$type<-'1'
all_gene<-rbind(trans_factor,other_type)
trans_factor<-unique(hsa_data$source_gene)
trans_factor<-as.data.frame(trans_factor)
com_pare<-rbind(trans_factor,target_gene)
library(dplyr)
com_pare <- distinct(com_pare)
names(all_gene)<-c("name","type")
save(all_gene,file = "all_gene.Rdata")




load(file ="all_gene.Rdata" )
all_gene$index <- seq(0, by = 1, length.out = nrow(all_gene))
library(readr)
hsa_data <- read_csv("hsa_transcript_regulation_depart_gene_whole.csv")
names(all_gene)<-c("source_gene","source_type","source_index")
hsa_data <- data.frame(sapply(hsa_data, as.character))
library(dplyr)
merged_df <- left_join(hsa_data, all_gene, by = "source_gene")
names(all_gene)<-c("target_gene","target_type","target_index")
all_map <- left_join(merged_df, all_gene, by = "target_gene")
positive<-all_map[all_map$trans_type=='regulation',]
negative<-all_map[all_map$trans_type=='irrelation',]
pos_train<-positive[positive$data_type=='train',]
all_map<-all_map[, c("source_gene","target_gene","source_index", "target_index","trans_type","data_type")]
write.table(all_map, file = "add_GO_term.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

write.csv(all_map, "hsa_transcript_regulation_edge_index_cnn.csv", row.names = FALSE)
names(all_gene)<-c("gene_id","gene_type","gene_index")
pattern <- rep(c("True", "True", "True", "False", "False"), length.out = nrow(all_gene))

# 将这个向量作为新列添加到数据框中
all_gene$train_mask <- pattern

val_vetor <- rep(c("False", "False", "False", "True", "False"), length.out = nrow(all_gene))

# 将这个向量作为新列添加到数据框中
all_gene$val_mask <- val_vetor
test_vetor <- rep(c("False", "False", "False", "False", "True"), length.out = nrow(all_gene))

# 将这个向量作为新列添加到数据框中
all_gene$test_mask <- test_vetor
write.csv(all_gene, "hsa_transcript_regulation_gene_node.csv", row.names = FALSE)
use_gene<-all_gene[, c("gene_id")]
use_gene<-as.data.frame(use_gene)
write.table(use_gene, file = "hsa_gene_node_id.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


library(readr)
gene_node <- read_csv("hsa_transcript_regulation_gene_node.csv")
tf_gene<- gene_node[1:973,]
other_gene<- gene_node[974:4398,]
tf_gene$train_mask<-'True'
random_indices <- sample(seq_len(nrow(tf_gene)), size = 195)
tf_gene$train_mask[random_indices] <- 'False'
table(tf_gene$train_mask)
tf_gene$val_mask<-NULL
tf_gene$test_mask <- tf_gene$train_mask
tf_gene$test_mask <- ifelse(tf_gene$test_mask == 'False', 'F', tf_gene$test_mask)
tf_gene$test_mask <- ifelse(tf_gene$test_mask == 'True', 'False', tf_gene$test_mask)
tf_gene$test_mask <- ifelse(tf_gene$test_mask == 'F', 'True', tf_gene$test_mask)


other_gene$train_mask<-'True'
random_indices <- sample(seq_len(nrow(other_gene)), size = 195)
other_gene$train_mask[random_indices] <- 'False'
table(other_gene$train_mask)
other_gene$val_mask<-NULL
other_gene$test_mask <- other_gene$train_mask
other_gene$test_mask <- ifelse(other_gene$test_mask == 'False', 'F', other_gene$test_mask)
other_gene$test_mask <- ifelse(other_gene$test_mask == 'True', 'False', other_gene$test_mask)
other_gene$test_mask <- ifelse(other_gene$test_mask == 'F', 'True', other_gene$test_mask)
all_gene<- rbind(tf_gene,other_gene)
write.csv(all_gene, "hsa_transcript_regulation_gene_node.csv", row.names = FALSE)


load(file = "use_all_negative_gene_whole.Rdata")

load(file ="all_gene.Rdata" )
all_gene$index <- seq(0, by = 1, length.out = nrow(all_gene))



library(readr)
hsa_data <- read_csv("hsa_transcript_regulation_depart_gene_whole.csv")
names(del_possible)<-c("source_gene","target_gene")
hsa_data <- data.frame(sapply(hsa_data, as.character))
library(dplyr)
del_possible<-del_possible %>%
  anti_join(hsa_data, by = c("source_gene","target_gene"))
has_irr<-hsa_data[hsa_data$trans_type=='irrelation',]
irr_train<-has_irr[has_irr$data_type=='train',]
irr_train<-irr_train[, c("source_gene", "target_gene")]
del_possible<-rbind(irr_train,del_possible)




del_possible<-del_possible[sample(nrow(del_possible)), ]
names(all_gene)<-c("source_gene","source_type","source_index")
names(del_possible)<-c("source_gene","target_gene")
library(dplyr)
merged_df <- left_join(del_possible, all_gene, by = "source_gene")
names(all_gene)<-c("target_gene","target_type","target_index")
all_map <- left_join(merged_df, all_gene, by = "target_gene")
all_map<-all_map[, c("source_gene","target_gene","source_index", "target_index")]
#write.csv(all_map, "negative_edge_all_possible.csv", row.names = FALSE)
write.csv(all_map, "has_edge_index_all_negative.csv", row.names = FALSE)



library(readr)
hsa_index <- read_csv("hsa_transcript_regulation_edge_index.csv")
has_irr<-hsa_index[hsa_index$trans_type=='irrelation',]
has_reg<-hsa_index[hsa_index$trans_type=='regulation',]
names(has_irr)<-c("target_index", "source_index","trans_type","data_type")
mutual_gene<- has_reg %>%
  semi_join(has_irr, by = c("source_index","target_index"))
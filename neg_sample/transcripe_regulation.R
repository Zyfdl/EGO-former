library(readr)
true_data<- read_delim("TFLink_Homo_sapiens_interactions_SS_simpleFormat_v1.0.tsv", 
                                                                    delim = "\t", escape_double = FALSE, 
                                                                    trim_ws = TRUE)
problems(true_data)
table(true_data$`Small-scale.evidence`)
true_pair <- true_data[, c("NCBI.GeneID.TF", "NCBI.GeneID.Target")]
write.csv(true_pair, "true_pair.csv", row.names = FALSE)

library(readr)
pre_data <- read_delim("TFLink_Homo_sapiens_interactions_LS_simpleFormat_v1.0.tsv", 
                                                                    delim = "\t", escape_double = FALSE, 
                                                                    trim_ws = TRUE)
table(pre_data$`Small-scale.evidence`)
pre_pair <- pre_data[, c("NCBI.GeneID.TF", "NCBI.GeneID.Target")]
write.csv(pre_pair, "pre_pair.csv", row.names = FALSE)


library(readr)
pre_pair <- read_csv("pre_pair.csv")
library(stringr)
library(dplyr)

df_pair <- pre_pair[!grepl("-", pre_pair$NCBI.GeneID.TF), ]
df_pair<- df_pair[!grepl("-", df_pair$NCBI.GeneID.Target), ]
rm(pre_pair)
df_TF <- df_pair[grepl(";", df_pair$NCBI.GeneID.TF), ]
df_target <- df_pair[grepl(";", df_pair$NCBI.GeneID.Target), ]
# 筛选出含有";"的行
df_pair <- df_pair[!grepl(";", df_pair$NCBI.GeneID.TF), ]
df_pair <- df_pair[!grepl(";", df_pair$NCBI.GeneID.Target), ]

library(tidyverse)

long_TF <- df_TF %>%
  separate_rows(NCBI.GeneID.TF, sep = ";")
long_TF <- long_TF %>%
  separate_rows(NCBI.GeneID.Target, sep = ";")
long_target <- df_target %>%
  separate_rows(NCBI.GeneID.TF, sep = ";")
long_target <- long_target %>%
  separate_rows(NCBI.GeneID.Target, sep = ";")

all_df<-rbind(df_pair,long_TF,long_target)
write.csv(all_df, "pre_pair.csv", row.names = FALSE)


library(readr)
true_pair <- read_csv("true_pair.csv")
df_pair <- true_pair[!grepl("-", true_pair$NCBI.GeneID.TF), ]
df_pair<- df_pair[!grepl("-", df_pair$NCBI.GeneID.Target), ]
df_TF <- df_pair[grepl(";", df_pair$NCBI.GeneID.TF), ]
df_target <- df_pair[grepl(";", df_pair$NCBI.GeneID.Target), ]
# 筛选出含有";"的行
df_pair <- df_pair[!grepl(";", df_pair$NCBI.GeneID.TF), ]
df_pair <- df_pair[!grepl(";", df_pair$NCBI.GeneID.Target), ]
library(tidyverse)
long_target <- df_target %>%
  separate_rows(NCBI.GeneID.TF, sep = ";")
long_target <- long_target %>%
  separate_rows(NCBI.GeneID.Target, sep = ";")
all_df<-rbind(df_pair,long_target)
library(dplyr)
all_df <- distinct(all_df)
write.csv(all_df, "true_pair.csv", row.names = FALSE)
target_irr<-unique(all_df$NCBI.GeneID.Target)
source_irr<-unique(all_df$NCBI.GeneID.TF)
trans_factor<-as.data.frame(source_irr)
dest_type<-as.data.frame(target_irr)
save(trans_factor,file = "trans_factor.Rdata")
save(dest_type,file = "dest_type.Rdata")


library(readr)
pre_pair <- read_csv("pre_pair.csv")
load(file ="trans_factor.Rdata" )
load(file ="dest_type.Rdata")
target_irr<-unique(pre_pair$NCBI.GeneID.Target)
source_irr<-unique(pre_pair$NCBI.GeneID.TF)
source_name<-as.data.frame(source_irr)
target_name<-as.data.frame(target_irr)
source_name<-rbind(source_name,trans_factor)
library(dplyr)
source_name <- distinct(source_name)
pre_pair<-distinct(pre_pair)
target_name <- rbind(dest_type,target_name)
library(dplyr)
target_name <- distinct(target_name)
names(source_name)<-c("name")
names(target_name)<-c("name")
target_name$name<-as.character(target_name$name)
source_name$name<-as.character(source_name$name)
target_name <- target_name %>%
  anti_join(source_name, by = c("name"))
target_r <- target_name[sample(nrow(target_name)), ]
source_r <- source_name[sample(nrow(source_name)), ]
source_name<-as.data.frame(source_r)
target_name<-as.data.frame(target_r)
names(source_name)<-c("name")
names(target_name)<-c("name")
source_one<-target_name[1:5000, ]
source_one<-as.data.frame(source_one)
names(source_one)<-c("name")
source_name<-rbind(source_one,source_name)
target_one<-target_name[5001:13000, ]
target_one<-as.data.frame(target_one)
save(source_name,file = "source_one.Rdata")
save(target_one,file = "target_one.Rdata")



load(file ="source_one.Rdata" )
load(file ="target_one.Rdata")
combined_data <- expand.grid(gene1=source_name$name, gene2=target_one$target_one)
combined_data <- data.frame(sapply(combined_data, as.character))
combined_data <- subset(combined_data, gene1 != gene2)
library(readr)
true_pair <- read_csv("true_pair.csv")
names(true_pair)<-c("gene1","gene2")
true_pair <- data.frame(sapply(true_pair, as.character))
library(dplyr)
combined_data <- combined_data %>%
  anti_join(true_pair, by = c("gene1","gene2"))
pre_pair <- read_csv("pre_pair.csv")
names(pre_pair)<-c("gene1","gene2")
pre_pair <- data.frame(sapply(pre_pair, as.character))
library(dplyr)
combined_data <- combined_data %>%
  anti_join(pre_pair, by = c("gene1","gene2"))
combined_data <- combined_data[sample(nrow(combined_data)), ]
irr_one <- combined_data %>%
  group_by(gene1) %>%
  filter(if (n() >= 4) row_number() <=4  else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(irr_one$gene2)
source_irr<-unique(irr_one$gene1)

irr_two <- irr_one %>%
  group_by(gene2) %>%
  filter(if (n() >= 4) row_number() <=4 else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(irr_one$gene2)
source_irr<-unique(irr_one$gene1)

true_pair$type<-'0'
irr_two$type<-'1'
all_data<-rbind(true_pair,irr_two)
write.table(all_data, file = "hsa_transcripte_regulation.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


library(readr)
hsa_data <- read_csv("hsa_transcript_regulation.csv")
positive<-hsa_data[hsa_data$trans_type=='positive',]
negative<-hsa_data[hsa_data$trans_type=='negative',]
target_irr<-unique(negative$target_gene)
source_irr<-unique(negative$source_gene)
other_one <- negative %>%
  group_by(source_gene) %>%
  filter(if (n() >= 3) row_number() <=3  else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(other_one$target_gene)
source_irr<-unique(other_one$source_gene)
other_two <- other_one %>%
  group_by(target_gene) %>%
  filter(if (n() >= 4) row_number() <=4  else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(other_two$target_gene)
source_irr<-unique(other_two$source_gene)
other_two <- other_two[sample(nrow(other_two)), ]
negative <- other_two[1:14586, ]
positive <- positive[sample(nrow(positive)), ]
new_column <- c(rep('test', 2917), rep('val', 2917), rep('train', 8752))
positive$data_type<-new_column
negative$data_type<-new_column
com_data<-rbind(positive,negative)
write.csv(com_data, "hsa_transcript_regulation_for_train.csv", row.names = FALSE)




library(readr)
htftarget <- read_delim("TF-Target-information.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
target_irr<-unique(htftarget$target)
source_irr<-unique(htftarget$TF)
library(org.Hs.eg.db)
htftarget$TF <- mapIds(org.Hs.eg.db, keys = htftarget$TF, column = "ENTREZID", keytype = "SYMBOL")
htftarget$target <- mapIds(org.Hs.eg.db, keys = htftarget$target, column = "ENTREZID", keytype = "SYMBOL")
# 假设 df 是你的数据框
htftarget <- na.omit(htftarget)
htftarget$tissue<-NULL
save(htftarget,file = "htftarget.Rdata")


library(readr)
pre_pair <- read_csv("pre_pair.csv")
load(file ="trans_factor.Rdata" )
load(file ="dest_type.Rdata")
combined_data <- expand.grid(gene1=trans_factor$source_irr, gene2=dest_type$target_irr)
combined_data <- data.frame(sapply(combined_data, as.character))
combined_data <- subset(combined_data, gene1 != gene2)
library(readr)
pre_pair <- read_csv("pre_pair.csv")
names(pre_pair)<-c("gene1","gene2")
pre_pair <- data.frame(sapply(pre_pair, as.character))
library(dplyr)
irr_for_pre <- combined_data %>%
  semi_join(pre_pair, by = c("gene1","gene2"))
target_irr<-unique(irr_for_pre$gene2)
source_irr<-unique(irr_for_pre$gene1)
del_pre <- combined_data %>%
  anti_join(pre_pair, by = c("gene1","gene2"))
library(readr)
true_pair <- read_csv("true_pair.csv")
names(true_pair)<-c("gene1","gene2")
true_pair <- data.frame(sapply(true_pair, as.character))
del_yes <- del_pre %>%
  anti_join(true_pair, by = c("gene1","gene2"))
load(file ="htftarget.Rdata" )
names(htftarget)<-c("gene1","gene2")
htftarget <- data.frame(sapply(htftarget, as.character))
del_possible <- del_yes %>%
  anti_join(htftarget, by = c("gene1","gene2"))
irr_for_pre_htf<- del_yes %>%
  semi_join(htftarget, by = c("gene1","gene2"))
save(irr_for_pre_htf,file = "htftarget_for_pre.Rdata")
save(irr_for_pre,file = "tflink_for_pre.Rdata")
save(del_possible,file = "del_possible.Rdata")
del_true_pair <- combined_data %>%
  anti_join(true_pair, by = c("gene1","gene2"))
save(del_true_pair,file = "del_true_pair.Rdata")


load(file = "del_possible.Rdata")
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
combined_data <- del_possible[sample(nrow(del_possible)), ]
irr_one <- combined_data %>%
  group_by(gene1) %>%
  filter(if (n() >= 14) row_number() <=14  else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(irr_one$gene2)
source_irr<-unique(irr_one$gene1)

irr_two <- combined_data %>%
  group_by(gene2) %>%
  filter(if (n() >= 4) row_number() <=4 else TRUE) %>%  # 条件过滤
  ungroup()
target_irr<-unique(irr_two$gene2)
source_irr<-unique(irr_two$gene1)

irr_data<-rbind(irr_one,irr_two)
library(dplyr)
irr_data <- distinct(irr_data)
target_irr<-unique(irr_data$gene2)
source_irr<-unique(irr_data$gene1)
library(readr)
true_pair <- read_csv("true_pair.csv")
names(true_pair)<-c("gene1","gene2")
target_irr<-unique(true_pair$gene2)
source_irr<-unique(true_pair$gene1)
true_pair <- distinct(true_pair)
true_pair$type<-'0'
irr_data$type<-'1'
all_data<-rbind(true_pair,irr_data)
write.table(all_data, file = "hsa_transcripte_regulation_part.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

library(readr)
true_pair <- read_csv("true_pair.csv")
TF_num<-table(true_pair$NCBI.GeneID.TF)
TF_num<-as.data.frame(TF_num)
# 假设 df 是你的数据框，column_name 是你想要按照其值进行降序排序的列名
df_sorted <- TF_num[order(TF_num$Freq, decreasing = TRUE), ]
row.names(df_sorted) <- NULL
zoshi<-table(df_sorted$Freq)
zoshi<-as.data.frame(zoshi)
zoshi<-zoshi[6:96, ]
# 假设你的数据框叫做df，a列和b列是你要绘制的横坐标和纵坐标
plot(zoshi$Var1, zoshi$Freq, type = "l", main = "Line Plot", xlab = "转录因子调控靶基因数量", ylab = "次数")

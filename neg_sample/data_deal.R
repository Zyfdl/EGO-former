library(readr)
hsa_data <- read_csv("hsa_transcript_regulation_part.csv")
positive<-hsa_data[hsa_data$trans_type=='regulation',]
negative<-hsa_data[hsa_data$trans_type=='irrelation',]
target_irr<-unique(positive$target_gene)
source_irr<-unique(positive$source_gene)
negative_one <- negative %>%
  semi_join(positive, by = c("source_gene"))
negative_two <- negative_one %>%
  semi_join(positive, by = c("target_gene"))

target_irr<-unique(negative_two$target_gene)
source_irr<-unique(negative_two$source_gene)

irr_one <- negative_two %>%
  group_by(source_gene) %>%
  filter(if (n() >= 2) row_number() <=2  else TRUE) %>%  # 条件过滤
  ungroup()

irr_two <- negative_two %>%
  group_by(target_gene) %>%
  filter(if (n() >= 4) row_number() <=4 else TRUE) %>%  # 条件过滤
  ungroup()
source_irr<-unique(irr_two$source_gene)
irr_data<-rbind(irr_two,irr_one)
library(dplyr)
irr_data <- distinct(irr_data)
irr_data <- irr_data[1:14586, ]
target_irr<-unique(irr_data$target_gene)
source_irr<-unique(irr_data$source_gene)
positive <- positive[sample(nrow(positive)), ]
irr_data <- irr_data[sample(nrow(irr_data)), ]
new_column <- c(rep('test', 2917), rep('val', 2917), rep('train', 8752))
positive$data_type<-new_column
irr_data$data_type<-new_column
com_data<-rbind(positive,irr_data)
write.csv(com_data, "hsa_transcript_regulation_part_pre_test.csv", row.names = FALSE)
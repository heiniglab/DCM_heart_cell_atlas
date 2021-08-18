# Packages
library(ggplot2)

# Code
table_ <- openxlsx::read.xlsx("/Users/Eric/Desktop/aggregation_all_patients.xlsx")
head(table_)
table_$patient <- gsub("b'(.*)'", "\\1", table_$patient)
table_$Gene_pred <- gsub("b'(.*)'", "\\1", table_$Gene_pred)
table_$true.label <- gsub("b'(.*)'", "\\1", table_$true.label)

table_[which(table_$patient=="H25"),]


table_$Gene_pred <- gsub("mutation negative", "PVneg", table_$Gene_pred)
table_$true.label <- gsub("mutation negative", "PVneg", table_$true.label)


table_subset <- table_[which(table_$Gene_pred==table_$true.label),]
table_$patient <- factor(table_$patient, levels=table_subset[order(table_subset[,"proba"]),"patient"])

table_$Gene_pred <- factor(table_$Gene_pred, levels=c("LMNA", 'TTN', 'RBM20', 'PKP2', 'PVneg'))
table_[,'true.label'] <- factor(table_[,'true.label'], levels=c("LMNA", 'TTN', 'RBM20', 'PKP2', 'PVneg'))

sum(table_[which(table_$patient=="H25"),"proba"])
g <- ggplot(table_, aes(x=patient, y=proba, fill=Gene_pred)) + geom_bar(stat="identity", colour='black') + facet_grid(.~true.label, scale="free_x") + 
  ylab("Aggregated Probability") + xlab("Patient") + theme(axis.text.x = element_text(color = "black", size = 13, angle = 90, hjust = 1, vjust = 1, face = "bold"),
                                                           #axis.text.y = element_blank(),  
                                                           axis.text.y = element_text(color = "black", size = 13, angle = 0, hjust = 0, vjust = 0, face = "bold"),  
                                                           axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
                                                           axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "bold"),
                                                           legend.text = element_text(size=15, face = "bold"),
                                                           legend.title = element_text(size=15, face = "bold"),
                                                           strip.text.x = element_text(size=20, face = "bold"),
                                                           legend.position = "none"
  ) + scale_fill_manual(values = c('#FC8D59', '#FEE08B', '#E6F598', '#99D594', '#3288BD'),
                        name="Primary Genetic\nDiagnosis") +
  scale_y_continuous(limits=c(0,1.01), breaks=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), expand = c(0.01, 0.01))
g
ggsave("/Users/Eric/Desktop/MainFig_6D_II.svg", g, width = 10, height = 7)
ggsave("/Users/Eric/Desktop/MainFig_6D_II.png", g, width = 10, height = 7)





table_MY <- openxlsx::read.xlsx("/Users/Eric/Desktop/results.xlsx", sheet=1)
table_EC <- openxlsx::read.xlsx("/Users/Eric/Desktop/results.xlsx", sheet=2)
table_FB <- openxlsx::read.xlsx("/Users/Eric/Desktop/results.xlsx", sheet=3)
table_CM <- openxlsx::read.xlsx("/Users/Eric/Desktop/results.xlsx", sheet=4)
table_MY[,"cell_type"] <- "MY"
table_EC[,"cell_type"] <- "EC"
table_FB[,"cell_type"] <- "FB"
table_CM[,"cell_type"] <- "CM"
table_2 <- rbind(table_MY, table_EC, table_FB, table_CM)

head(table_2)
tapply(table_2$quantaty.predicted, table_2$patient, sum)==table_2[-which(duplicated(table_2$patient)),"quantaty.true"]

table_2$Gene_pred <- gsub("mutation negative", "PVneg", table_2$Gene_pred)
table_2$true.label <- gsub("mutation negative", "PVneg", table_2$true.label)
table_2$Gene_pred <- gsub("b'(.*)'", "\\1", table_2$Gene_pred)
table_2$true.label <- gsub("b'(.*)'", "\\1", table_2$true.label)
table_2$patient <- gsub("b'(.*)'", "\\1", table_2$patient)
table_2[,'Gene_pred'] <- factor(table_2[,'Gene_pred'], levels=c("LMNA", 'TTN', 'RBM20', 'PKP2', 'PVneg'))
table_2[,'true.label'] <- factor(table_2[,'true.label'], levels=c("LMNA", 'TTN', 'RBM20', 'PKP2', 'PVneg'))
table_2[,'cell_type'] <- factor(table_2[,'cell_type'], levels=c("CM", 'FB', 'EC', 'MY'))

library(ggplot2)
g <- ggplot(table_2, aes(x=patient, y=share, fill=Gene_pred))+ geom_bar(stat="identity") + 
  facet_grid(cell_type~true.label, scales = "free_x")+ theme(axis.text.x = element_text(color = "black", size = 13, angle = 90, hjust = 1, vjust = 1, face = "bold"),
                                                             #axis.text.y = element_blank(),  
                                                             axis.text.y = element_text(color = "black", size = 13, angle = 0, hjust = 0, vjust = 0, face = "bold"),  
                                                             axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
                                                             axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "bold"),
                                                             legend.text = element_text(size=15, face = "bold"),
                                                             legend.title = element_text(size=15, face = "bold"),
                                                             strip.text.x = element_text(size=13, face = "bold"),
                                                             strip.text.y = element_text(size=13, face = "bold")
  ) + scale_fill_manual(values = c('#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b'),
                        name="Primary Genetic\nDiagnosis") + xlab("Patient") + ylab("% Probability")

head(table_2)
table_2_sub <- table_2[which(table_2$true.label==table_2$Gene_pred),]


g <- ggplot(table_2_sub, aes(x=true.label, y=share, fill=Gene_pred))+ geom_boxplot(outlier.shape=NA) + geom_jitter(size=1, width = 0.2) + 
  facet_wrap(cell_type~., ncol=2)+ theme(axis.text.x = element_text(color = "black", size = 13, angle = 0, face = "bold"),
                                                             #axis.text.y = element_blank(),  
                                                             axis.text.y = element_text(color = "black", size = 13, angle = 0, hjust = 0, vjust = 0, face = "bold"),  
                                                             axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
                                                             axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "bold"),
                                                             legend.text = element_text(size=15, face = "bold"),
                                                             legend.title = element_text(size=15, face = "bold"),
                                                             strip.text.x = element_text(size=13, face = "bold"),
                                                             strip.text.y = element_text(size=13, face = "bold"),
                                                             legend.position = "none"
  ) + scale_fill_manual(values = c('#FC8D59', '#FEE08B', '#E6F598', '#99D594', '#3288BD'),
                        name="Primary Genetic\nDiagnosis") + xlab("") + ylab("% Probability")
g
ggsave("/Users/Eric/Desktop/MainFig_6D_I.svg", g, width = 10, height = 7)
ggsave("/Users/Eric/Desktop/MainFig_6D_I.png", g, width = 10, height = 7)

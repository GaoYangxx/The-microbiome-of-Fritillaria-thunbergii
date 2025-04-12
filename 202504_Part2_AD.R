#图2c α多样性
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
install.packages(c("dplyr", "tidyverse", "data.table", "hillR",
                   "viridis", "hrbrthemes", "paletteer", "mgcv"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "breakaway", "microbiome"))
install.packages("remotes")
remotes::install_github("kstagaman/phyloseqCompanion")
library(dplyr)
library(tidyverse)
library(breakaway)
library(data.table)
library(phyloseq)
library(microViz)
library(viridis)
library(hrbrthemes)
library(paletteer)
library(microbiome)
library(mgcv)
library(hillR)
library(phyloseqCompanion)
#细菌阿尔法
bacteria_ASV <- read.csv(file = "D:/study/master/meiji/bacteria_ASV.csv",
                         sep=",",header=TRUE,check.names = FALSE)
# 定义文件的绝对路径（使用正斜杠）
metadata <- "D:/study/master/metadata.tsv"
# 导入元数据文件
metadata <- import_qiime_sample_data(metadata)
# 查看元数据
head(metadata)
# 读取 ASV 表（特征表）
b_ASV <- bacteria_ASV[, c("ASV",metadata$Sample.ID)]
b_ASV <- b_ASV %>% column_to_rownames(var = "ASV")
head(b_ASV)
# 读取 taxonomy（分类信息）
b_tax <- bacteria_ASV[, 2:9]
b_tax <- b_tax %>% column_to_rownames(var = "ASV")
head(b_tax)
# 合并为 phyloseq 对象
b_ASV <- otu_table(b_ASV, taxa_are_rows = TRUE)
b_tax <- tax_table(as.matrix(b_tax))
b_merged <- merge_phyloseq(b_ASV, b_tax, metadata)
#去除非细菌
b_merged <- subset_taxa(b_merged, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
b_merged <- subset_taxa(b_merged, (Kingdom!="d__Archaea") | is.na(Kingdom))
b_merged <- subset_taxa(b_merged, (Order!="o__Chloroplast") )
b_merged <- subset_taxa(b_merged, (Family!="f__Mitochondria"))
b_merged <- subset_taxa(b_merged, (Family!="NA"))
# 提取 ASV 表
b_ASV_df <- as.matrix(otu_table(b_merged, taxa_are_rows=T))
b_ASV_df <- b_ASV_df[rowSums(b_ASV_df[])>0,]
# 提取分类表
b_tax <- tax_table(b_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(b_merged))
b_diversity_nomis <- phyloseq::sample_data(estimate_richness(b_merged,measures=c("Observed","Shannon")))#计算样本的丰富度指数
b_alphadiv_nomis <- merge_phyloseq(b_merged, b_diversity_nomis)# 合并物种丰度数据和多样性数据
b_alphadt_nomis_df<- sample_data(b_alphadiv_nomis)# 提取样本数据
#计算丰富度
b_OR <- b_alphadt_nomis_df[, c("Sample.ID", "Group", "Observed")]#提取3列
b_ASVrichness <- b_OR %>% #每个组（Group）的平均物种丰富度和标准差
  group_by(Group) %>% 
  summarise(average=mean(Observed), std=sd(Observed))
b_ASVrichness
b_average_richness <- b_OR %>% #整个数据集的平均物种丰富度和标准差
  summarise(average=mean(Observed), std=sd(Observed))
b_average_richness
b_median_richness <- b_OR %>% #整个数据集的中位数和四分位数（25%、50%、75%分位数
  summarise(median=median(Observed), x = quantile(Observed, c(0.25, 0.5, 0.75)))
b_median_richness
#画图
habitat_labeller <- c("JRG" = "Jurong\nRhizosphere\nSoil", "JJG" = "Jingjiang\nRhizosphere\nSoil","TZG" = "Tongzhou\nRhizosphere\nSoil","PAG" = "Panan\nRhizosphere\nSoil","JRN" = "Jurong\nBulb", "JJN" = "Jingjiang\nBulb","TZN" = "Tongzhou\nBulb","PAN" = "Panan\nBulb")
b_OR$Group <- factor(b_OR$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))# 按照你想要的顺序设置 Group 因子的 levels
plot_OR_bacteria <- ggplot(b_OR,aes(x= Group,y=Observed, color= Group)) + 
  geom_violin(width=1.4, alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) +  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16), # 坐标轴标题（x、y）
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14), # 坐标轴刻度标签（数字）
    legend.position = "none"
  ) +
  labs(y = "Bacterial ASV richness", x = "") +
  scale_x_discrete(labels = habitat_labeller) +
  scale_colour_manual(values= c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_OR_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2c_alpha_bacteria.png", plot = plot_OR_bacteria, width = 8, height = 8, dpi = 600, bg = "transparent")
#saveRDS(plot_OR_bacteria,file="D:/study/master/Main_Figure_tables/Figure_2/2c_alpha_bacteria.rds")#二进制文件.rds为R 对象
#Shannon多样性
b_Shannon <- b_alphadt_nomis_df[, c("Sample.ID", "Group", "Observed", "Shannon")]
b_shannon <- b_Shannon %>% 
  group_by(Group) %>% 
  summarise(average=mean(Shannon), std=sd(Shannon))
b_shannon
b_median_shannon <- b_Shannon %>% 
  reframe(median=median(Shannon), x = quantile(Shannon, c(0.25, 0.5, 0.75)))
b_median_shannon
b_Shannon$Group <- factor(b_Shannon$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))# 按照你想要的顺序设置 Group 因子的 levels
#画图
plot_Shannon_bacteria <- ggplot(b_Shannon,aes(x= Group,y=Shannon, color= Group)) + 
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) +  # 设置 alpha 值 fShannon 透明度
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(labels = habitat_labeller)+
  labs(x = "")+
  scale_colour_manual(values= c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_Shannon_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2c_Shannon_bacteria.png", plot = plot_Shannon_bacteria, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_Shannon_bacteria, file = " D:/study/master/Main_Figure_tables/Figure_2/2c_Shannon_bacteria.rds")
#计算 Hill 数量指数（Hill number），当 q = 1 时，Hill 数等价于 Shannon 指数的指数形式，即 exp(Shannon)
b_hill_q1_nomis<-as.data.frame(hill_taxa(t(b_ASV_df), q = 1))
# 数据框转换成样本元数据（sample_data）
b_hill_q1_indices_nomis <- phyloseq::sample_data(b_hill_q1_nomis) 
#合并，样本元数据增加额外的列 hill_q1 数据
b_merged_hill_nomis <- phyloseq::merge_phyloseq(b_merged, b_hill_q1_indices_nomis) 
b_meta_diversity_nomis_hill <- sample.data.frame(b_merged_hill_nomis)#提取样本元数据
b_median_shannon_hill<- b_meta_diversity_nomis_hill %>% #计算中位数和四分位数，列名乱
  summarise(median=median(hill_taxa.t.b_ASV_df...q...1.), x = quantile(hill_taxa.t.b_ASV_df...q...1..1, c(0.25, 0.5, 0.75))) 
b_median_shannon_hill
#用布拉指数来衡量均匀度
b_bulla_estimate <- phyloseq::sample_data(microbiome::evenness(b_merged, index="all"))
b_alphadiv_nomis <- merge_phyloseq(b_alphadiv_nomis, b_bulla_estimate)
b_alphadt_nomis_df<- sample_data(b_alphadiv_nomis)
hist(b_alphadt_nomis_df$bulla)
b_evenness <- b_alphadt_nomis_df %>% 
  group_by(Group) %>% 
  summarise(average=mean(bulla), std=sd(bulla))
b_evenness
b_evenness_median <- b_alphadt_nomis_df %>% 
  reframe(median=median(bulla), x = quantile(bulla, c(0.25, 0.5, 0.75)))
b_evenness_median
b_evenness_total <- b_alphadt_nomis_df %>% 
  summarise(average=mean(bulla), std=sd(bulla))
#真菌阿尔法
fungi_ASV <- read.csv(file = "D:/study/master/meiji/fungi_ASV.csv",
                      sep=",",header=TRUE,check.names = FALSE)
# 定义文件的绝对路径（使用正斜杠）
metadata <- "D:/study/master/metadata.tsv"
# 导入元数据文件
metadata <- import_qiime_sample_data(metadata)
# 查看元数据
head(metadata)
# 读取 ASV 表（特征表）
f_ASV <- fungi_ASV[, c("ASV",metadata$Sample.ID)]
f_ASV <- f_ASV %>% column_to_rownames(var = "ASV")
head(f_ASV)
# 读取 taxonomy（分类信息）
f_tax <- fungi_ASV[, 2:9]
f_tax <- f_tax %>% column_to_rownames(var = "ASV")
head(f_tax)
# 合并为 phyloseq 对象
f_ASV <- otu_table(f_ASV, taxa_are_rows = TRUE)
f_tax <- tax_table(as.matrix(f_tax))
f_merged <- merge_phyloseq(f_ASV, f_tax, metadata)
#去除非细菌
f_merged <- subset_taxa(f_merged, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
f_merged <- subset_taxa(f_merged, (Kingdom!="d__Archaea") | is.na(Kingdom))
f_merged <- subset_taxa(f_merged, (Order!="o__Chloroplast") )
f_merged <- subset_taxa(f_merged, (Family!="f__Mitochondria"))
f_merged <- subset_taxa(f_merged, (Family!="NA"))
# 提取 ASV 表
f_ASV_df <- as.matrix(otu_table(f_merged, taxa_are_rows=T))
f_ASV_df <- f_ASV_df[rowSums(f_ASV_df[])>0,]
# 提取分类表
f_tax <- tax_table(f_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(f_merged))
f_diversity_nomis <- phyloseq::sample_data(estimate_richness(f_merged,measures=c("Observed","Shannon")))#计算样本的丰富度指数
f_alphadiv_nomis <- merge_phyloseq(f_merged, f_diversity_nomis)# 合并物种丰度数据和多样性数据
f_alphadt_nomis_df<- sample_data(f_alphadiv_nomis)# 提取样本数据
#计算丰富度
f_OR <- f_alphadt_nomis_df[, c("Sample.ID", "Group", "Observed")]#提取3列
f_ASVrichness <- f_OR %>% #每个组（Group）的平均物种丰富度和标准差
  group_by(Group) %>% 
  summarise(average=mean(Observed), std=sd(Observed))
f_ASVrichness
f_average_richness <- f_OR %>% #整个数据集的平均物种丰富度和标准差
  summarise(average=mean(Observed), std=sd(Observed))
f_average_richness
f_median_richness <- f_OR %>% #整个数据集的中位数和四分位数（25%、50%、75%分位数
  summarise(median=median(Observed), x = quantile(Observed, c(0.25, 0.5, 0.75)))
f_median_richness
#画图
habitat_labeller <- c("JRG" = "Jurong\nRhizosphere\nSoil", "JJG" = "Jingjiang\nRhizosphere\nSoil","TZG" = "Tongzhou\nRhizosphere\nSoil","PAG" = "Panan\nRhizosphere\nSoil","JRN" = "Jurong\nBulb", "JJN" = "Jingjiang\nBulb","TZN" = "Tongzhou\nBulb","PAN" = "Panan\nBulb")
f_OR$Group <- factor(f_OR$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))# 按照你想要的顺序设置 Group 因子的 levels
plot_OR_fungi <- ggplot(f_OR,aes(x= Group,y=Observed, color= Group)) + 
  geom_violin(width=1.4, alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) +  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16), # 坐标轴标题（x、y）
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14), # 坐标轴刻度标签（数字）
    legend.position = "none"
  ) +
  labs(y = "Fungal ASV richness", x = "") +
  scale_x_discrete(labels = habitat_labeller) +
  scale_colour_manual(values= c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_OR_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2c_alpha_fungi.png", plot = plot_OR_fungi, width = 8, height = 8, dpi = 600, bg = "transparent")
#saveRDS(plot_OR_fungi,file="D:/study/master/Main_Figure_tables/Figure_2/2c_alpha_fungi.rds")#二进制文件.rds为R 对象
#Shannon多样性
f_Shannon <- f_alphadt_nomis_df[, c("Sample.ID", "Group", "Observed", "Shannon")]
f_shannon <- f_Shannon %>% 
  group_by(Group) %>% 
  summarise(average=mean(Shannon), std=sd(Shannon))
f_shannon
f_median_shannon <- f_Shannon %>% 
  reframe(median=median(Shannon), x = quantile(Shannon, c(0.25, 0.5, 0.75)))
f_median_shannon
f_Shannon$Group <- factor(f_Shannon$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))# 按照你想要的顺序设置 Group 因子的 levels
#画图
plot_Shannon_fungi <- ggplot(f_Shannon,aes(x= Group,y=Shannon, color= Group)) + 
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1, outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) +  # 设置 alpha 值 fShannon 透明度
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1),
        legend.position = "none") +
  scale_x_discrete(labels = habitat_labeller)+
  labs(x = "")+
  scale_colour_manual(values= c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_Shannon_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2c_Shannon_fungi.png", plot = plot_Shannon_fungi, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_Shannon_fungi, file = " D:/study/master/Main_Figure_tables/Figure_2/2c_Shannon_fungi.rds")
#计算 Hill 数量指数（Hill number），当 q = 1 时，Hill 数等价于 Shannon 指数的指数形式，即 exp(Shannon)
f_hill_q1_nomis<-as.data.frame(hill_taxa(t(f_ASV_df), q = 1))
# 数据框转换成样本元数据（sample_data）
f_hill_q1_indices_nomis <- phyloseq::sample_data(f_hill_q1_nomis) 
#合并，样本元数据增加额外的列 hill_q1 数据
f_merged_hill_nomis <- phyloseq::merge_phyloseq(f_merged, f_hill_q1_indices_nomis) 
f_meta_diversity_nomis_hill <- sample.data.frame(f_merged_hill_nomis)#提取样本元数据
f_median_shannon_hill<- f_meta_diversity_nomis_hill %>% #计算中位数和四分位数，列名乱
  summarise(median=median(hill_taxa.t.f_ASV_df...q...1.), x = quantile(hill_taxa.t.f_ASV_df...q...1..1, c(0.25, 0.5, 0.75))) 
f_median_shannon_hill
#用布拉指数来衡量均匀度
f_bulla_estimate <- phyloseq::sample_data(microbiome::evenness(f_merged, index="all"))
f_alphadiv_nomis <- merge_phyloseq(f_alphadiv_nomis, f_bulla_estimate)
f_alphadt_nomis_df<- sample_data(f_alphadiv_nomis)
hist(f_alphadt_nomis_df$bulla)
f_evenness <- f_alphadt_nomis_df %>% 
  group_by(Group) %>% 
  summarise(average=mean(bulla), std=sd(bulla))
f_evenness
f_evenness_median <- f_alphadt_nomis_df %>% 
  reframe(median=median(bulla), x = quantile(bulla, c(0.25, 0.5, 0.75)))
f_evenness_median
f_evenness_total <- f_alphadt_nomis_df %>% 
  summarise(average=mean(bulla), std=sd(bulla))
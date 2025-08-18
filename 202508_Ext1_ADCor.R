#图1 b α多样性的多样性、均匀度和PD及其相关性
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
install.packages(c("dplyr", "tidyverse", "data.table", "hillR",
                   "viridis", "hrbrthemes", "paletteer", "mgcv", "readxl", "picante", "corrplot"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "breakaway", "microbiome"))
install.packages("remotes")
remotes::install_github("kstagaman/phyloseqCompanion")
devtools::install_github("jbisanz/qiime2R")
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
library(qiime2R)
library(readxl)
library(picante)
library(corrplot)
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
# 筛选：每个 ASV 在所有样本中的总和 ≥ 4
#b_ASV <- b_ASV[rowSums(b_ASV) >= 4, ]
# 读取 taxonomy（分类信息）
b_tax <- bacteria_ASV[, 2:9]
b_tax <- b_tax %>% column_to_rownames(var = "ASV")
head(b_tax) 
b_tax <- data.frame(b_tax)
#b_tax_df <- b_tax %>%
#  filter(!(stringr::str_detect(Genus, "g__uncultured") | stringr::str_detect(Genus, "g__unclassified") | stringr::str_detect(Genus, "g__norank")))
#b_ASV <- b_ASV[rownames(b_ASV) %in% rownames(b_tax_df), ]
#分步生成和读取系统发育树
b_tree <- read_qza("D:/study/master/meiji/b_rooted-tree.qza")
b_phylo_tree <- b_tree$data
b_md5<-read_excel("D:/study/master/meiji/b_ASV_md5.xlsx")#md5值
b_rename_vector <- setNames(b_md5$`ASV ID`, b_md5$md5)
b_phylo_tree$tip.label <- b_rename_vector[b_phylo_tree$tip.label] #替换tip.label
# 合并为 phyloseq 对象
b_ASV <- otu_table(b_ASV, taxa_are_rows = TRUE)
b_tax <- tax_table(as.matrix(b_tax))
b_merged <- merge_phyloseq(b_ASV, b_tax, metadata, phy_tree(b_phylo_tree))
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
habitat_labeller <- c("JRG" = "Jurong\nRhizosphere\nSoil", "JJG" = "Jingjiang\nRhizosphere\nSoil","TZG" = "Tongzhou\nRhizosphere\nSoil","PAG" = "Panan\nRhizosphere\nSoil","JRN" = "Jurong\nBulb", "JJN" = "Jingjiang\nBulb","TZN" = "Tongzhou\nBulb","PAN" = "Panan\nBulb")
plot_Shannon_bacteria <-
  ggplot(b_Shannon,aes(x= Group,y=Shannon, color= Group)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1,
               outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) + # 设置 alpha 值 fShannon 透明度
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    panel.border = element_blank(),     # 移除面板边框
    axis.line.x = element_line(color = "black"), # 设置 X 轴线颜色为黑色
    axis.line.y = element_line(color = "black"), # 设置 Y 轴线颜色为黑色
    axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 14), # 调整 X 轴刻度标签样式和大小
    legend.position = "none", # 移除图例
    axis.title.x = element_text(size = 16), # X 轴标题大小
    axis.title.y = element_text(size = 16), # Y 轴标题大小
    axis.text.y = element_text(size = 14) # Y 轴刻度标签大小
  ) +
  scale_x_discrete(labels = habitat_labeller)+
  labs(y = "Bacterial Shannon diversity", x = "")+ # 修改 Y 轴标签为 Shannon 多样性
  scale_colour_manual(values=
                        c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                          "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_Shannon_bacteria
# ggsave("D:/study/master/Extended_Data_tables/Figure_1/1b_Shannon_bacteria.png", plot = plot_Shannon_bacteria, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_Shannon_bacteria, file = " D:/study/master/Main_Figure_tables/Figure_1/1b_Shannon_bacteria.rds")
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
b_alphadt_nomis_df$Group <- factor(b_alphadt_nomis_df$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))# 按照你想要的顺序设置 Group 因子的 levels
plot_bulla_bacteria <-
  ggplot(b_alphadt_nomis_df, aes(x = Group, y = bulla, color = Group)) +
  geom_violin(width = 1.4, alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black", alpha = 1,
               outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.9) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    panel.border = element_blank(),     # 移除面板边框
    axis.line.x = element_line(color = "black"), # 设置 X 轴线颜色为黑色
    axis.line.y = element_line(color = "black"), # 设置 Y 轴线颜色为黑色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), # 调整 X 轴刻度标签样式和大小
    legend.position = "none", # 移除图例
    axis.title.x = element_text(size = 16), # X 轴标题大小
    axis.title.y = element_text(size = 16), # Y 轴标题大小
    axis.text.y = element_text(size = 14) # Y 轴刻度标签大小
  ) +
  scale_x_discrete(labels = habitat_labeller) + # 应用 habitat_labeller 来设置 X 轴标签
  labs(y = "Bacterial Bulla's evenness", x = "") + # 修改 Y 轴标签为布拉均匀度
  scale_colour_manual(values = c("#8C57A2FF", "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#DC9445FF", "#bee183"))
plot_bulla_bacteria
# ggsave("D:/study/master/Extended_Data_tables/Figure_1/1b_bulla_bacteria.png", plot = plot_bulla_bacteria, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_bulla_bacteria, file = " D:/study/master/Main_Figure_tables/Figure_1/1b_bulla_bacteria.rds")
# 计算 Faith's PD
# 从 phyloseq 对象中提取 OTU 表和树
b_otu_mat <- as.data.frame(otu_table(b_merged))
b_tree_obj <- phy_tree(b_merged)
b_pd_results <- picante::pd(t(b_otu_mat), b_tree_obj, include.root = FALSE)
# 将 PD 结果添加到样本元数据中
b_alphadt_nomis_df <- sample_data(b_alphadiv_nomis) %>%
  as_tibble() %>%
  left_join(as_tibble(b_pd_results, rownames = "Sample.ID") %>% dplyr::select(Sample.ID, PD), by = "Sample.ID")
# PD 多样性数据框
b_PD <- b_alphadt_nomis_df[, c("Sample.ID", "Group", "PD")]
# 计算 PD 的平均值和标准差
b_pd_summary <- b_PD %>%
  group_by(Group) %>%
  summarise(average = mean(PD, na.rm = TRUE), std = sd(PD, na.rm = TRUE)) # 添加 na.rm = TRUE 处理缺失值
print(b_pd_summary)
# 计算 PD 的中位数和四分位数
b_median_pd <- b_PD %>%
  reframe(median = median(PD, na.rm = TRUE), x = quantile(PD, c(0.25, 0.5, 0.75), na.rm = TRUE))
print(b_median_pd)
# 按照你想要的顺序设置 Group 因子的 levels
b_PD$Group <- factor(b_PD$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))
#绘制 PD 图 
plot_PD_bacteria <-
  ggplot(b_PD, aes(x = Group, y = PD, color = Group)) +
  geom_violin(width = 1.4, alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black", alpha = 1,
               outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.9) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14)
  ) +
  scale_x_discrete(labels = habitat_labeller) +
  labs(y = "Bacterial Faith's PD", x = "") + # 修改 Y 轴标签
  scale_colour_manual(values = c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                 "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_PD_bacteria
# 保存 PD 图
# ggsave("D:/study/master/Extended_Data_tables/Figure_1/1b_PD_bacteria.png",
#        plot = plot_PD_bacteria, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_PD_bacteria, file = "D:/study/master/Main_Figure_tables/Figure_1/1b_PD_bacteria.rds")
# 选择用于相关性分析的 alpha 多样性指数列
b_alpha_diversity_corr <- b_alphadt_nomis_df %>%
  dplyr::select(Observed, Shannon, PD, bulla) %>%
  na.omit()
#计算 Spearman 相关矩阵 
b_correlation_matrix_alpha <- cor(b_alpha_diversity_corr, method = "spearman")
# 重命名列名（和行名）
colnames(b_correlation_matrix_alpha) <- c("Richness", "Shannon", "PD", "Evenness")
rownames(b_correlation_matrix_alpha) <- c("Richness", "Shannon", "PD", "Evenness")
# 定义你想要的最终顺序
desired_order <- c("Richness", "Shannon", "Evenness", "PD")
# 重新排列矩阵的行和列，使其按照 desired_order 的顺序
b_correlation_matrix_alpha_ordered <-
  b_correlation_matrix_alpha[desired_order, desired_order]
# 定义新的颜色渐变
my_colors <- colorRampPalette(c("#F7CAD0", "#F8E19B", "#90E0C2", "#00B4D8", "#7A4D9B"))(200)
# 绘制相关性热图
#png(file = "D:/study/master/Extended_Data_tables/Figure_1/1c_alpha_correlation_bacteria.png", width = 10, height = 8, units = "in", res = 300)
corrplot(b_correlation_matrix_alpha_ordered, # 使用重新排序后的矩阵
         method = "color",       # 显示颜色方块
         addCoef.col = "white",  # 添加相关系数并设置为黑色
         number.cex = 1.4, 
         tl.col = "black",       # 标签颜色
         tl.srt = 45,            # 标签旋转角度
         tl.cex = 1.4,
         cl.pos = "r",           # 色阶条位置（右侧）
         cl.ratio = 0.2,         # 色阶条宽度
         cl.cex = 1.4,    
         col = my_colors, # 应用新的颜色渐变
         #col.lim = NULL,       # 设置颜色条的范围
         #is.corr = FALSE,        # 关键更改：视为非相关矩阵 
         title = "Correlations among bacterial\nAlpha diversity indices", # 图表标题
         cex.main = 1.8,   
         mar = c(0,0,3,0)        # 调整边距，使标题显示完整
)
#dev.off()
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
#分步生成和读取系统发育树
f_tree <- read_qza("D:/study/master/meiji/f_rooted-tree.qza")
f_phylo_tree <- f_tree$data
f_md5<-read_excel("D:/study/master/meiji/f_ASV_md5.xlsx")#md5值
f_rename_vector <- setNames(f_md5$`ASV ID`, f_md5$md5)
f_phylo_tree$tip.label <- f_rename_vector[f_phylo_tree$tip.label] #替换tip.label
# 合并为 phyloseq 对象
f_ASV <- otu_table(f_ASV, taxa_are_rows = TRUE)
f_tax <- tax_table(as.matrix(f_tax))
f_merged <- merge_phyloseq(f_ASV, f_tax, metadata, phy_tree(f_phylo_tree))
#去除非真菌
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
f_tax <- data.frame(f_tax)
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
habitat_labeller <- c("JRG" = "Jurong\nRhizosphere\nSoil", "JJG" = "Jingjiang\nRhizosphere\nSoil","TZG" = "Tongzhou\nRhizosphere\nSoil","PAG" = "Panan\nRhizosphere\nSoil","JRN" = "Jurong\nBulb", "JJN" = "Jingjiang\nBulb","TZN" = "Tongzhou\nBulb","PAN" = "Panan\nBulb")
plot_Shannon_fungi <-
  ggplot(f_Shannon,aes(x= Group,y=Shannon, color= Group)) +
  geom_violin(width=1.4,alpha=0.5) +
  geom_boxplot(width=0.1, color="black", alpha=1,
               outlier.shape=NA) +
  geom_jitter(position=position_jitter(0.2), alpha=0.9) + # 设置 alpha 值 fShannon 透明度
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    panel.border = element_blank(),     # 移除面板边框
    axis.line.x = element_line(color = "black"), # 设置 X 轴线颜色为黑色
    axis.line.y = element_line(color = "black"), # 设置 Y 轴线颜色为黑色
    axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 14), # 调整 X 轴刻度标签样式和大小
    legend.position = "none", # 移除图例
    axis.title.x = element_text(size = 16), # X 轴标题大小
    axis.title.y = element_text(size = 16), # Y 轴标题大小
    axis.text.y = element_text(size = 14) # Y 轴刻度标签大小
  ) +
  scale_x_discrete(labels = habitat_labeller)+
  labs(y = "Fungal Shannon diversity", x = "")+ # 修改 Y 轴标签为 Shannon 多样性
  scale_colour_manual(values=
                        c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                          "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_Shannon_fungi
# ggsave("D:/study/master/Extended_Data_tables/Figure_1/1b_Shannon_fungi.png", plot = plot_Shannon_fungi, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_Shannon_fungi, file = " D:/study/master/Main_Figure_tables/Figure_1/1b_Shannon_fungi.rds")
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
f_alphadt_nomis_df$Group <- factor(f_alphadt_nomis_df$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))# 按照你想要的顺序设置 Group 因子的 levels
plot_bulla_fungi <-
  ggplot(f_alphadt_nomis_df, aes(x = Group, y = bulla, color = Group)) +
  geom_violin(width = 1.4, alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black", alpha = 1,
               outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.9) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), # 移除主网格线
    panel.grid.minor = element_blank(), # 移除次网格线
    panel.border = element_blank(),     # 移除面板边框
    axis.line.x = element_line(color = "black"), # 设置 X 轴线颜色为黑色
    axis.line.y = element_line(color = "black"), # 设置 Y 轴线颜色为黑色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), # 调整 X 轴刻度标签样式和大小
    legend.position = "none", # 移除图例
    axis.title.x = element_text(size = 16), # X 轴标题大小
    axis.title.y = element_text(size = 16), # Y 轴标题大小
    axis.text.y = element_text(size = 14) # Y 轴刻度标签大小
  ) +
  scale_x_discrete(labels = habitat_labeller) + # 应用 habitat_labeller 来设置 X 轴标签
  labs(y = "Fungal Bulla's evenness", x = "") + # 修改 Y 轴标签为布拉均匀度
  scale_colour_manual(values = c("#8C57A2FF", "#3EBCB6", "#82581FFF", "#2F509EFF",
                                 "#E5614CFF", "#97A1A7FF", "#DC9445FF", "#bee183"))
plot_bulla_fungi
# ggsave("D:/study/master/Extended_Data_tables/Figure_1/1b_bulla_fungi.png", plot = plot_bulla_fungi, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_bulla_fungi, file = " D:/study/master/Main_Figure_tables/Figure_1/1b_bulla_fungi.rds")
# 计算 Faith's PD
# 从 phyloseq 对象中提取 OTU 表和树
f_otu_mat <- as.data.frame(otu_table(f_merged))
f_tree_obj <- phy_tree(f_merged)
f_pd_results <- picante::pd(t(f_otu_mat), f_tree_obj, include.root = FALSE)
# 将 PD 结果添加到样本元数据中
f_alphadt_nomis_df <- sample_data(f_alphadiv_nomis) %>%
  as_tibble() %>%
  left_join(as_tibble(f_pd_results, rownames = "Sample.ID") %>% dplyr::select(Sample.ID, PD), by = "Sample.ID")
# PD 多样性数据框
f_PD <- f_alphadt_nomis_df[, c("Sample.ID", "Group", "PD")]
# 计算 PD 的平均值和标准差
f_pd_summary <- f_PD %>%
  group_by(Group) %>%
  summarise(average = mean(PD, na.rm = TRUE), std = sd(PD, na.rm = TRUE)) # 添加 na.rm = TRUE 处理缺失值
print(f_pd_summary)
# 计算 PD 的中位数和四分位数
f_median_pd <- f_PD %>%
  reframe(median = median(PD, na.rm = TRUE), x = quantile(PD, c(0.25, 0.5, 0.75), na.rm = TRUE))
print(f_median_pd)
# 按照你想要的顺序设置 Group 因子的 levels
f_PD$Group <- factor(f_PD$Group, levels = c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN"))
#绘制 PD 图 
plot_PD_fungi <-
  ggplot(f_PD, aes(x = Group, y = PD, color = Group)) +
  geom_violin(width = 1.4, alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black", alpha = 1,
               outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.9) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14)
  ) +
  scale_x_discrete(labels = habitat_labeller) +
  labs(y = "Fungal Faith's PD", x = "") + # 修改 Y 轴标签
  scale_colour_manual(values = c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                 "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183"))
plot_PD_fungi
# 保存 PD 图
# ggsave("D:/study/master/Extended_Data_tables/Figure_1/1b_PD_fungi.png",
#        plot = plot_PD_fungi, width = 8, height = 8, dpi = 600, bg = "transparent")
# saveRDS(plot_PD_fungi, file = "D:/study/master/Main_Figure_tables/Figure_1/1b_PD_fungi.rds")
# 选择用于相关性分析的 alpha 多样性指数列
f_alpha_diversity_corr <- f_alphadt_nomis_df %>%
  dplyr::select(Observed, Shannon, PD, bulla) %>%
  na.omit()
#计算 Spearman 相关矩阵 
f_correlation_matrix_alpha <- cor(f_alpha_diversity_corr, method = "spearman")
# 重命名列名（和行名）
colnames(f_correlation_matrix_alpha) <- c("Richness", "Shannon", "PD", "Evenness")
rownames(f_correlation_matrix_alpha) <- c("Richness", "Shannon", "PD", "Evenness")
# 定义你想要的最终顺序
desired_order <- c("Richness", "Shannon", "Evenness", "PD")
# 重新排列矩阵的行和列，使其按照 desired_order 的顺序
f_correlation_matrix_alpha_ordered <-
  f_correlation_matrix_alpha[desired_order, desired_order]
# 定义新的颜色渐变
my_colors <- colorRampPalette(c("#F7CAD0", "#F8E19B", "#90E0C2", "#00B4D8", "#7A4D9B"))(200)
# 绘制相关性热图
#png(file = "D:/study/master/Extended_Data_tables/Figure_1/1c_alpha_correlation_fungi.png", width = 10, height = 8, units = "in", res = 300)
corrplot(f_correlation_matrix_alpha_ordered, # 使用重新排序后的矩阵
         method = "color",       # 显示颜色方块
         addCoef.col = "white",  # 添加相关系数并设置为黑色
         number.cex = 1.4, 
         tl.col = "black",       # 标签颜色
         tl.srt = 45,            # 标签旋转角度
         tl.cex = 1.4,
         cl.pos = "r",           # 色阶条位置（右侧）
         cl.ratio = 0.2,         # 色阶条宽度
         cl.cex = 1.4,    
         col = my_colors, # 应用新的颜色渐变
         #col.lim = NULL,       # 设置颜色条的范围
         #is.corr = FALSE,        # 关键更改：视为非相关矩阵 
         title = "Correlations among fungal\nAlpha diversity indices", # 图表标题
         cex.main = 1.8,   
         mar = c(0,0,3,0)        # 调整边距，使标题显示完整
)
#dev.off()

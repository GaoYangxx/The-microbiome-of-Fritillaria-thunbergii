#图2b 伽马多样性
install.packages(c("ggplot2", "vegan", "knitr", "dplyr", "tibble", "iNEXT", "grid", "cowplot"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "biomformat"))
# 加载必要的 R 包
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(iNEXT)
library(tibble)
library(grid) 
library(cowplot)
#细菌伽马
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
row_sums <- rowSums(b_ASV)#去除过少ASV
rows_to_remove <- row_sums <= 10
b_ASV <- b_ASV[!rows_to_remove, ]
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
b_ASV <- otu_table(b_merged, taxa_are_rows=TRUE)
# 提取分类表
b_tax <- tax_table(b_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(b_merged))
groups <- c("JRG", "JJG", "TZG", "PAG","JRN", "JJN", "TZN", "PAN") 
b_sample_lists <- list()
b_asv_lists <- list()
b_vec_lists <- list()
sample_data(b_merged)$Group <- as.factor(sample_data(b_merged)$Group)
sample_data(b_merged)$Sample.ID <- as.factor(sample_data(b_merged)$Sample.ID)
sample_data(b_merged)$Origin <- as.factor(sample_data(b_merged)$Origin)
sample_data(b_merged)$Niche <- as.factor(sample_data(b_merged)$Niche)
for (group in groups) {
  b_subset_result <- subset_samples(b_merged, Group == group)#每组的phyloseq对象
  b_merged_result <- merge_samples(b_subset_result,"Sample.ID")#以样本列重新分组
  # 每组样本asv数量
  b_asv_table_group <- otu_table(b_merged_result, taxa_are_rows=T)
  b_asv_table_group_t <- t(b_asv_table_group)
  b_asv_lists[[group]] <- b_asv_table_group_t
  #每组asv出现样本数
  b_asv_table_df <- b_asv_table_group_t
  b_sumrow_group <- unname(rowSums(b_asv_table_df>0))
  b_sort_group<- sort(b_sumrow_group, decreasing=T)
  b_vec_group <- b_sort_group[b_sort_group >0]
  b_vec_lists[[group]] <- b_vec_group
}
# 准备 iNEXT 输入数据
b_list_exped_all <- list(jrg=c(ncol(b_asv_lists$JRG),b_vec_lists$JRG),jjg=c(ncol(b_asv_lists$JJG),b_vec_lists$JJG), tzg=c(ncol(b_asv_lists$TZG),b_vec_lists$TZG), pag=c(ncol(b_asv_lists$PAG),b_vec_lists$PAG), jrn=c(ncol(b_asv_lists$JRN),b_vec_lists$JRN),jjn=c(ncol(b_asv_lists$JJN),b_vec_lists$JJN), tzn=c(ncol(b_asv_lists$TZN),b_vec_lists$TZN), pan=c(ncol(b_asv_lists$PAN),b_vec_lists$PAN))
# 执行 iNEXT 分析丰度，基于 q 阶的希尔数（Hill numbers）进行插值与外推分析
b_out_all_exped <- iNEXT(b_list_exped_all, q=0, datatype="incidence_freq", se=TRUE, conf=0.95, nboot=99)
# 将 iNEXT 结果转换为数据框
b_df <- fortify(b_out_all_exped, type = 1)
# 分割出观测数据和预测数据
b_df.point <- b_df[which(b_df$Method == "Observed"),]
b_df.line <- b_df[which(b_df$Method != "Observed"),]
b_df.line$Method <- factor(b_df.line$Method, c("Rarefaction", "Extrapolation"))
b_df.asympote <- data.frame(y = c(24,8),
                            Asymptote = c("jrg", "jjg", "tzg", "pag","jrn", "jjn", "tzn", "pan"))
# 绘图
plot_bacteria <-ggplot(b_df, aes(x = x, y = y, colour = Assemblage)) + 
  geom_line(aes(linetype = Method), lwd = 1.5, data = b_df.line) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = Assemblage, colour = NULL), alpha = 0.2) +
  labs(x = "Number of samples", y = "Bacterial ASV richness") +
  scale_fill_manual(values = c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                               "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")) +
  scale_color_manual(values = c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")) +
  scale_linetype_discrete(name = NULL) +
  scale_x_continuous(breaks = c(3, 6, 9, 12)) +
  guides(
    colour = "none", 
    fill = "none",
    linetype = guide_legend(override.aes = list(size = 2))  # 控制图例线条粗细
  ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = c(0.01, 1),
        legend.justification = c(0, 1),
        legend.key.width = unit(2, "cm"), # 控制横线长度
        axis.title.x = element_text(size = 16), # 坐标轴标题（x、y）
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14), # 坐标轴刻度标签（数字）
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14), # 图例文字
        legend.title = element_text(size = 14) )
plot_bacteria
# 自定义标签和颜色（与你主图保持一致）
legend_df <- data.frame(
  x = rep(c(1, 3), each = 4),
  y = rep(4:1, 2),
  label = c("Jurong Rhizosphere Soil", "Jingjiang Rhizosphere Soil", "Tongzhou Rhizosphere Soil", "Panan Rhizosphere Soil",
            "Jurong Bulb", "Jingjiang Bulb", "Tongzhou Bulb", "Panan Bulb"),
  color = c("#8C57A2", "#3EBCB6", "#82581F", "#2F509E", "#E5614C", "#97A1A7", "#DC9445", "#BEE183")
)
# 图例图层（圆点 + 标签）
legend_plot <- ggplot(legend_df, aes(x, y)) +
  geom_point(aes(color = color), size = 4.5, shape = 16, show.legend = FALSE) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.1, color = "black", size = 4.5) +
  scale_color_identity() +  # 直接使用提供的颜色
  theme_void() +
  xlim(1, 4) + ylim(0, 6)
# 合并主图和自定义图例（右下角 inset）
final_plot_bacteria <- ggdraw() +
  draw_plot(plot_bacteria, 0, 0, 1, 1) +
  draw_plot(legend_plot, 0.46, 0.1, 0.55, 0.22)  # 控制位置和大小（x, y, width, height）
#print(final_plot_bacteria)
#ggsave("D:/study/master/Main_Figure_tables/Figure_2/2a_gamma_bacteria.png", plot = plot_bacteria, width = 8, height = 6, dpi = 600, bg = "transparent")
# 计算频率比
b_inext_freq_results <- b_out_all_exped$AsyEst  #提取iNEXT 结果中的渐近估计表格，包括观测值与估计值
b_inext_freq_results$prop <- b_inext_freq_results$Observed / b_inext_freq_results$Estimator#计算频率比，观测值 / 估计值
b_inext_freq_results <- b_inext_freq_results[b_inext_freq_results$Diversity == 'Species richness',]# 只保留“物种丰富度”指标
# 计算中位数和四分位数
b_median_GD_freq <- b_inext_freq_results %>% 
  summarise(med = median(prop), 
            lower_quartile = quantile(prop, 0.25),
            median = quantile(prop, 0.5),
            upper_quartile = quantile(prop, 0.75))
#真菌伽马
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
#row_sums <- rowSums(f_ASV)#去除过少ASV
#rows_to_remove <- row_sums <= 1
#f_ASV <- f_ASV[!rows_to_remove, ]
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
f_ASV <- otu_table(f_merged, taxa_are_rows=TRUE)
# 提取分类表
f_tax <- tax_table(f_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(f_merged))
groups <- c("JRG", "JJG", "TZG", "PAG","JRN", "JJN", "TZN", "PAN") 
f_sample_lists <- list()
f_asv_lists <- list()
f_vec_lists <- list()
sample_data(f_merged)$Group <- as.factor(sample_data(f_merged)$Group)
sample_data(f_merged)$Sample.ID <- as.factor(sample_data(f_merged)$Sample.ID)
sample_data(f_merged)$Origin <- as.factor(sample_data(f_merged)$Origin)
sample_data(f_merged)$Niche <- as.factor(sample_data(f_merged)$Niche)
for (group in groups) {
  f_subset_result <- subset_samples(f_merged, Group == group)#每组的phyloseq对象
  f_merged_result <- merge_samples(f_subset_result,"Sample.ID")#以样本列重新分组
  # 每组样本asv数量
  f_asv_table_group <- otu_table(f_merged_result, taxa_are_rows=T)
  f_asv_table_group_t <- t(f_asv_table_group)
  f_asv_lists[[group]] <- f_asv_table_group_t
  #每组asv出现样本数
  f_asv_table_df <- f_asv_table_group_t
  f_sumrow_group <- unname(rowSums(f_asv_table_df>0))
  f_sort_group<- sort(f_sumrow_group, decreasing=T)
  f_vec_group <- f_sort_group[f_sort_group >0]
  f_vec_lists[[group]] <- f_vec_group
}
# 准备 iNEXT 输入数据
f_list_exped_all <- list(jrg=c(ncol(f_asv_lists$JRG),f_vec_lists$JRG),jjg=c(ncol(f_asv_lists$JJG),f_vec_lists$JJG), tzg=c(ncol(f_asv_lists$TZG),f_vec_lists$TZG), pag=c(ncol(f_asv_lists$PAG),f_vec_lists$PAG), jrn=c(ncol(f_asv_lists$JRN),f_vec_lists$JRN),jjn=c(ncol(f_asv_lists$JJN),f_vec_lists$JJN), tzn=c(ncol(f_asv_lists$TZN),f_vec_lists$TZN), pan=c(ncol(f_asv_lists$PAN),f_vec_lists$PAN))
# 执行 iNEXT 分析丰度，基于 q 阶的希尔数（Hill numbers）进行插值与外推分析
f_out_all_exped <- iNEXT(f_list_exped_all, q=0, datatype="incidence_freq", se=TRUE, conf=0.95, nboot=99)
# 将 iNEXT 结果转换为数据框
f_df <- fortify(f_out_all_exped, type = 1)
# 分割出观测数据和预测数据
f_df.point <- f_df[which(f_df$Method == "Observed"),]
f_df.line <- f_df[which(f_df$Method != "Observed"),]
f_df.line$Method <- factor(f_df.line$Method, c("Rarefaction", "Extrapolation"))
f_df.asympote <- data.frame(y = c(24,8),
                            Asymptote = c("jrg", "jjg", "tzg", "pag","jrn", "jjn", "tzn", "pan"))
# 绘图
plot_fungi <-ggplot(f_df, aes(x = x, y = y, colour = Assemblage)) + 
  geom_line(aes(linetype = Method), lwd = 1.5, data = f_df.line) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = Assemblage, colour = NULL), alpha = 0.2) +
  labs(x = "Number of samples", y = "Fungal ASV richness") +
  scale_fill_manual(values = c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                               "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")) +
  scale_color_manual(values = c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                                "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")) +
  scale_linetype_discrete(name = NULL) +
  scale_x_continuous(breaks = c(3, 6, 9, 12)) +
  guides(
    colour = "none", 
    fill = "none",
    linetype = guide_legend(override.aes = list(size = 2))  # 控制图例线条粗细
  ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = c(0.01, 1),
        legend.justification = c(0, 1),
        legend.key.width = unit(2, "cm"), # 控制横线长度
        axis.title.x = element_text(size = 16), # 坐标轴标题（x、y）
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14), # 坐标轴刻度标签（数字）
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14), # 图例文字
        legend.title = element_text(size = 14) )
plot_fungi
#ggsave("D:/study/master/Main_Figure_tables/Figure_2/2a_gamma_fungi.png", plot = plot_fungi, width = 8, height = 6, dpi = 600, bg = "transparent")
# 自定义图例数据：4列2行
legend_df <- data.frame(
  x = rep(seq(1,4), each = 2),  # 设置为4列，增大列间距，by = 1.5表示间隔1.5
  y = rep(2:1, 4),  # 2 行
  label = c("Jurong Rhizosphere Soil", "Jingjiang Rhizosphere Soil","Jurong Bulb","Jingjiang Bulb", "Tongzhou Rhizosphere Soil", "Panan Rhizosphere Soil","Tongzhou Bulb",
            "Panan Bulb"),
  color = c("#82581F","#8C57A2","#2F509E", "#3EBCB6", "#DC9445","#E5614C", "#BEE183", "#97A1A7")
)
# 图例图层（圆点 + 标签）
legend_plot <- ggplot(legend_df, aes(x, y)) +
  geom_point(aes(color = color), size = 4.5, shape = 16, show.legend = FALSE) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.05, color = "black", size = 4.5) +
  scale_color_identity() +  # 直接使用提供的颜色
  theme_void() +
  xlim(1, 5) + ylim(0.5, 2.5)  # 调整范围以适应 4 列 2 行布局
legend_plot
#ggsave("D:/study/master/Main_Figure_tables/Figure_2/2a_gamma_legend.png", plot = legend_plot, width = 12, height = 0.8, dpi = 2000, bg = "transparent")
# 合并主图和自定义图例（右下角 inset）
#final_plot_fungi <- ggdraw() +
#draw_plot(plot_fungi, 0, 0, 1, 1) +
#draw_plot(legend_plot, 0.46, 0.1, 0.55, 0.22)  # 控制位置和大小（x, y, width, height）
#print(final_plot_fungi)
# 计算频率比
f_inext_freq_results <- f_out_all_exped$AsyEst  
f_inext_freq_results$prop <- f_inext_freq_results$Observed / f_inext_freq_results$Estimator
f_inext_freq_results <- f_inext_freq_results[f_inext_freq_results$Diversity == 'Species richness',]
# 计算中位数和四分位数
f_median_GD_freq <- f_inext_freq_results %>% 
  summarise(med = median(prop), 
            lower_quartile = quantile(prop, 0.25),
            median = quantile(prop, 0.5),
            upper_quartile = quantile(prop, 0.75))

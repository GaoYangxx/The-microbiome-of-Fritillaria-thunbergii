#图3b 甜甜圈图
install.packages(c("RColorBrewer", "ggplot2", "tidyverse", "reshape2", "ggpubr", " dplyr", "Matrix"))
install.packages("phyloseqCompanion")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("phyloseq", "indicspecies"))
library("RColorBrewer")
library("ggplot2")
library(tidyverse)
library(reshape2)
library(ggpubr)
library(phyloseq)
library(indicspecies)
library(phyloseqCompanion)
library(dplyr)
library(Matrix)
library(scales) # 用于格式化百分比标签
library(ggrepel) # 用于处理标签重叠
#细菌甜甜圈
#读取数据
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
# 如果 b_ASV 和 b_tax 的行名是 ASV ID，需要先将其转为列
b_ASV_df <- b_ASV %>% rownames_to_column(var = "ASV")
b_tax_df <- b_tax %>% rownames_to_column(var = "ASV")
# 将 b_ASV 与 b_tax 中的 'Species' 列连接
b_species <- b_ASV_df %>%
  left_join(select(b_tax_df, ASV, Species), by = "ASV")
# 按照 'Species' 列进行分组，并对每个样本的丰度求和
b_species <- b_species %>%
  select(all_of(c("Species", metadata$Sample.ID))) %>%
  group_by(Species) %>%
  summarise(across(all_of(metadata$Sample.ID), sum)) %>%
  ungroup() # 解除分组
# 排序
b_species <- b_species %>%
  arrange(Species)
b_species_df <- b_species # 复制操作，R中通常直接赋值即可，除非您需要独立副本
print(paste("Length of b_species_df (number of species):", nrow(b_species_df)))
# 按照 'Species' 分组，并保留每个物种的第一个分类信息（排除ASV列）
b_tax_species <- b_tax_df %>%
  group_by(Species) %>%
  slice(1) %>% # 获取每个Species组的第一个行
  ungroup() %>%
  select(-ASV) %>% # 排除ASV列
  arrange(Species) # 按照Species排序，与Python中的sort_index()对应
# 进行替换操作和去除 'p__' 前缀
b_tax_species <- b_tax_species %>%
  mutate(
    Phylum = case_when(
      Phylum == "p__unclassified_k__norank_d__Bacteria" ~ "Unclassified",
      Phylum == "p__SAR324_cladeMarine_group_B" ~ "SAR324",
      TRUE ~ Phylum # 保持其他值不变
    ),
    # 去除 'p__' 前缀
    Phylum = str_replace(Phylum, "^p__", "")
  )
# 创建物种学名到 'species_数字' 格式的映射并应用
b_unique_species_names <- b_species %>% pull(Species) %>% unique() %>% sort()
# 创建一个从物种学名到 'species_数字' 格式的映射数据框
b_species_to_numeric_id <- tibble(
  original_species = b_unique_species_names,
  numeric_id = paste0("species_", 1:length(b_unique_species_names))
)
# 转换 b_species 的 'Species' 列
b_species <- b_species %>%
  left_join(b_species_to_numeric_id, by = c("Species" = "original_species")) %>%
  select(-Species) %>% # 移除原始的Species列
  rename(Species = numeric_id) %>% # 将numeric_id重命名为Species
  select(Species, everything()) # 确保Species列在第一位
# 转换 b_tax_species 的 'Species' 列
b_tax_species <- b_tax_species %>%
  left_join(b_species_to_numeric_id, by = c("Species" = "original_species")) %>%
  select(-Species) %>% # 移除原始的Species列
  rename(Species = numeric_id) %>% # 将numeric_id重命名为Species
  select(Species, everything()) # 确保Species列在第一位
b_species_df <- b_species
print(paste("Length of b_species_df after species ID mapping:", nrow(b_species_df)))
metadata <- data.frame(metadata)
#核心asv（在大于2个居群生态位中相对丰度大于0.1%）
# 计算相对丰度
b_species_matrix <- b_species_df %>%
  column_to_rownames(var = "Species") %>% # 将 'Species' 列转换为行名
  as.matrix() # 转换为矩阵以便进行高效的列操作
# 计算相对丰度：将每个样本（列）的丰度除以该样本的总丰度
b_species_rel_matrix <- t(t(b_species_matrix) / colSums(b_species_matrix))
# 将结果转换回 tibble，并将行名（物种名称）转回列
b_species_rel <- b_species_rel_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "species") %>% # 将行名转回 'species' 列
  as_tibble() # 转换为 tibble
# 将数据重塑为长格式并合并元数据
b_rel_long <- b_species_rel %>%
  pivot_longer(
    cols = all_of(metadata$Sample.ID), # 仅选择样本ID列进行重塑
    names_to = "Sample",       # 新的列，存放样本名称
    values_to = "value"        # 新的列，存放丰度值
  )
# 创建从元数据中提取 Origin 和 Niche 的映射数据框
sample_origins_map_df <- metadata %>% select(`Sample.ID`, Origin) %>% rename(Sample = `Sample.ID`)
sample_niches_map_df <- metadata %>% select(`Sample.ID`, Niche) %>% rename(Sample = `Sample.ID`)
# 通过 left_join 将 Origin 和 Niche 信息添加到长格式表中
b_rel_long <- b_rel_long %>%
  left_join(sample_origins_map_df, by = "Sample") %>%
  left_join(sample_niches_map_df, by = "Sample")
# 计算物种整体相对丰度并进行筛选 
b_species_raw_matrix_for_sum <- b_species_df %>%
  column_to_rownames(var = "Species") %>%
  as.matrix()
b_species_total_counts <- rowSums(b_species_raw_matrix_for_sum)
# 计算所有物种的总丰度
b_total_counts_across_all_species <- sum(b_species_total_counts)
# 计算每个物种的整体相对丰度（相对于所有物种的总和）
b_species_relative_abundance <- b_species_total_counts / b_total_counts_across_all_species
# 定义整体相对丰度阈值 (0.1% 转换为小数)
overall_abundance_threshold <- 0.001
# 筛选出满足整体丰度阈值的物种名称列表
b_species_above_overall_threshold <- names(b_species_relative_abundance[b_species_relative_abundance >= overall_abundance_threshold])
# 根据整体丰度阈值过滤长格式丰度表
b_rel_long_filtered_by_overall_abundance <- b_rel_long %>%
  filter(species %in% b_species_above_overall_threshold)
#计算物种在每个 Origin 和 Niche 中的出现次数
# 计算每个物种在每个 Origin 中的平均丰度
b_species_origin_mean_filtered <- b_rel_long_filtered_by_overall_abundance %>%
  group_by(species, Origin) %>%
  summarise(mean_abundance_in_origin = mean(value), .groups = 'drop') # .groups = 'drop' 用于解除分组
# 识别在特定 Origin 中平均丰度大于 0 的物种（即有出现）
b_species_in_origin <- b_species_origin_mean_filtered %>%
  filter(mean_abundance_in_origin > 0)
# 统计每个物种出现的 Origin 数量
b_species_origin_count <- b_species_in_origin %>%
  group_by(species) %>%
  summarise(origin_count = n(), .groups = 'drop')
# 计算每个物种在每个 Niche 中的平均丰度
b_species_niche_mean_filtered <- b_rel_long_filtered_by_overall_abundance %>%
  group_by(species, Niche) %>%
  summarise(mean_abundance_in_niche = mean(value), .groups = 'drop')
# 识别在特定 Niche 中平均丰度大于 0 的物种（即有出现）
b_species_in_niche <- b_species_niche_mean_filtered %>%
  filter(mean_abundance_in_niche > 0)
# 统计每个物种出现的 Niche 数量
b_species_niche_count <- b_species_in_niche %>%
  group_by(species) %>%
  summarise(niche_count = n(), .groups = 'drop')
# 合并 Origin 和 Niche 的出现次数统计
b_species_counts_combined <- b_species_origin_count %>%
  inner_join(b_species_niche_count, by = "species")
# 定义出现次数阈值
origin_occurrence_threshold <- 3
niche_occurrence_threshold <- 2
# 筛选满足所有条件（整体丰度，以及在 Origin 和 Niche 中的出现次数）的物种，得到核心物种的名称列表
b_core_species_names <- b_species_counts_combined %>%
  filter(origin_count >= origin_occurrence_threshold &
           niche_count >= niche_occurrence_threshold) %>%
  pull(species) # 提取物种名称作为向量
# 从相对丰度数据中提取核心物种的数据
b_core_abundance <- b_species_rel %>%
  filter(species %in% b_core_species_names)
# 从原始计数数据中提取核心物种的数据
b_core_species_raw_counts <- b_species_df %>%
  filter(Species %in% b_core_species_names)
# 打印核心物种的数量
print(paste("找到的核心物种数量:", nrow(b_core_species_raw_counts)))
#最特殊
# 计算所有物种整体相对丰度的中位数并筛选稀有物种 
b_median_relative_abundance <- median(b_species_relative_abundance)
print(paste("所有物种整体相对丰度的中位数:", b_median_relative_abundance))
# 定义稀有类群的丰度阈值
b_rare_threshold <- 0.00001
# 筛选出稀有物种
b_rare_species_names <- names(b_species_relative_abundance[b_species_relative_abundance < b_rare_threshold])
print(paste("稀有物种的数量:", length(b_rare_species_names)))
# 将 b_species_df 转换为长格式并添加 Group 信息
# 转置 b_species_df 并转换为长格式
b_speciesdfmelt <- b_species_df %>%
  pivot_longer(
    cols = -Species, # 排除 Species 列，将其余列转换为长格式
    names_to = "Sample",
    values_to = "value"
  ) %>%
  rename(species = Species) %>% # 重命名 Species 列为 species 以匹配 Python 代码中的命名
  filter(value > 0) # 过滤掉丰度为 0 的行
# 创建从 Sample ID 到 Group 的映射数据框
sample_group_map_df <- metadata %>% select(`Sample.ID`, Group) %>% rename(Sample = `Sample.ID`)
# 添加 Group 信息到长格式表中
b_control_species <- b_speciesdfmelt %>%
  left_join(sample_group_map_df, by = "Sample") %>%
  filter(!is.na(Group)) # 移除 Group 为 NA 的行 (对应 Python 的 dropna(subset=['Group']))
# 统计物种出现在多少个 Group 中 
# 提取 Group 和 species 列并去重
b_controlspecies_unique <- b_control_species %>%
  select(Group, species) %>%
  distinct() # 相当于 drop_duplicates()
# 统计每个 species 出现在多少个 Group 中
b_controlspecies_group_count <- b_controlspecies_unique %>%
  group_by(species) %>%
  summarise(group_count = n(), .groups = 'drop')
# 只保留那些只在一个 Group 中出现的物种
b_unique_species_single_group_names <- b_controlspecies_group_count %>%
  filter(group_count == 1) %>%
  pull(species) # 提取物种名称向量
print(paste("只在一个 Group 中出现的物种数量:", length(b_unique_species_single_group_names)))
# 找到同时满足两个条件的“最特殊稀有物种
b_rare_unique_species_names <- intersect(b_rare_species_names, b_unique_species_single_group_names)
print(paste("最特殊稀有物种（Rare AND Uniques to single group）的数量:", length(b_rare_unique_species_names)))
# 提取这些“最特殊稀有物种”的原始丰度数据
b_rare_unique_abundance <- b_species_df %>%
  filter(Species %in% b_rare_unique_species_names)
b_rare_unique_species_raw_counts <- b_species_df %>%
  filter(Species %in% b_rare_unique_species_names)
#统一读取
b_species_df <- data.frame(b_species)
rownames(b_species_df) <- b_species_df$Species
b_species_df <- b_species_df[,-1]
b_core_species_raw_counts <- data.frame(b_core_species_raw_counts)
rownames(b_core_species_raw_counts) <- b_core_species_raw_counts$Species
b_core_species_raw_counts <- b_core_species_raw_counts[,-1]
b_rare_unique_species_raw_counts <- data.frame(b_rare_unique_species_raw_counts)
rownames(b_rare_unique_species_raw_counts) <- b_rare_unique_species_raw_counts$Species
b_rare_unique_species_raw_counts <- b_rare_unique_species_raw_counts[,-1]
b_dat	<- b_species_df
b_cores <- b_core_species_raw_counts 
b_uniq <- b_rare_unique_species_raw_counts
b_generalists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/b_generalists_raw_counts.csv", header = TRUE, row.names = 1)
b_specialists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/b_specialists_raw_counts.csv", header = TRUE, row.names = 1)
b_species <- rownames(b_dat)
b_species <- as.data.frame(b_species)
colnames(b_species) <- "species"
b_cores$species<-rownames(b_cores)
rownames(b_cores)<-NULL
b_cores <- b_cores %>%
  select(species)%>%# 只保留 species 这一列
  mutate(type = "core")#增加一列 type，所有行都标记为 "core"
b_uniq$species<-rownames(b_uniq)
rownames(b_uniq)<-NULL
b_uniq <- b_uniq %>%
  select(species)%>%
  mutate(type = "uniq")
colnames(b_uniq) <- c("species","type")
b_generalists$species<-rownames(b_generalists)
rownames(b_generalists)<-NULL
b_generalists <- b_generalists %>%
  select(species)%>%
  mutate(type = "generalists")
colnames(b_generalists) <- c("species", "type")
b_specialists$species<-rownames(b_specialists)
rownames(b_specialists)<-NULL
b_specialists <- b_specialists %>%
  select(species)%>%
  mutate(type = "specialists")
colnames(b_specialists) <- c("species","type")
b_dat_numb <- rbind(b_cores, b_uniq, b_generalists, b_specialists)#合并
b_species_other <- b_species %>%#筛选出那些没有出现在四类标签里的species
  filter(!(species %in% b_dat_numb$species))%>%
  mutate(type = "Other")
#计算每类有多少种species
b_datToPlot <- b_dat_numb %>%
  arrange(species, type) %>%#先对 species 和 type 排序，确保拼接顺序一致
  group_by(species)%>%# 每个species为一组
  summarise(type = str_c(type, collapse="_"))%>%#将同一个species的多个 type 合并成一个字符串，比如 "core_indicator"
  ungroup()%>%#取消分组状态
  bind_rows(b_species_other)%>%#合并"Other" 的species
  group_by(type)%>%# 按 type 统计每一类species的种类数量
  summarise(sum = n())
b_datToPlot$type <- factor(b_datToPlot$type, levels = c("core_generalists","core",  "core_specialists",  "specialists", "specialists_uniq", "uniq","generalists_uniq","generalists", "Other"))#设置分类顺序与颜色
colors <- c(
  "core_generalists" = "#4AC0FF",
  "core" = "#5499C7",
  "core_specialists" = "#2F609F",
  "specialists" = "#45B39D",
  "specialists_uniq" = "#006400",
  "uniq" = "#D4AC0D",
  "generalists_uniq" = "#30308B",
  "generalists" = "#EC7063",
  "Other" = "#737373"
)
b_datToPlot <- b_datToPlot[order(b_datToPlot$type),]#排序
b_datToPlot$fraction = b_datToPlot$sum / sum(b_datToPlot$sum)#计算每一类所占百分比
b_datToPlot$ymax = cumsum(b_datToPlot$fraction)#计算累积百分比
b_datToPlot$ymin = c(0, head(b_datToPlot$ymax, n=-1))# 计算累积最小值
b_datToPlot$labelPosition <- (b_datToPlot$ymax + b_datToPlot$ymin) / 2#计算标签的中心点
# 绘制第一张图 (b_p1): 绝对分布
# 准备 b_p1 的标签数据
b_datToPlot$label <- paste0(b_datToPlot$sum)
# 绘制 b_p1
b_p1 <- ggplot(b_datToPlot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = type)) +
  geom_rect() +
  # 使用 geom_text_repel 替代 geom_text
  geom_text_repel(aes(x = 4.5, y = labelPosition, label = label),
                  size = 4,
                  # 调整推力，确保标签不会重叠
                  force = 0,
                  # 调整 x轴方向的推力
                  nudge_x = 0,
                  # 连接线样式
                  segment.color = "gray60",
                  segment.size = 0) +
  scale_fill_manual(
    values = colors,
    labels = c("Cores and generalists", "Cores",
               "Cores and specialists", "Specialists",
               " Uniques and specialists", "Uniques",
               " Uniques and generalists", "Generalists", "Other")
  ) +
  coord_polar(theta = "y") +
  xlim(c(1, 6)) + # 为标签留出更多空间
  theme_void() +
  theme(
    legend.title = element_blank(), # legend.position = "none", #取消图例
    legend.text = element_text(size = 14)
  ) +
  labs(title = "Absolute distribution of bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = -25)))
#+guides(fill = guide_legend(nrow = 1, byrow = TRUE))
b_p1
b_dat_sum <- rowSums(b_dat[1:48]) # 对每个 species 求所有样本丰度之和
b_dat_sum <- as.data.frame(b_dat_sum)
b_dat_sum$species <- rownames(b_dat_sum) # 加上 species 名称
#计算每类的相对丰度
b_datToPlotCov <- b_dat_numb %>%
  arrange(species, type) %>%
  group_by(species)%>%
  summarise(type = str_c(type, collapse="_"))%>% # 合并同一个species多个type
  ungroup()%>%
  bind_rows(b_species_other)%>% 
  left_join(b_dat_sum)%>%# 添加丰度信息（b_dat_sum 是每个 species 的总丰度）
  filter(!(b_dat_sum == 0))%>% # 移除丰度为0的species
  mutate(perc = b_dat_sum / sum(b_dat_sum))%>%#计算相对丰度
  group_by(type)%>%
  summarise(fraction = sum(perc)) # 统计每种类型的丰度总占比
sum(b_datToPlotCov$fraction)
b_datToPlotCov$type <- factor(b_datToPlotCov$type, levels = c("core_generalists","core", "core_specialists", "specialists", "specialists_uniq", "uniq","generalists_uniq","generalists", "Other"))
b_datToPlotCov <- b_datToPlotCov[order(b_datToPlotCov$type),]
# b_datToPlotCov$fraction = b_datToPlotCov$sum / sum(b_datToPlotCov$sum)
b_datToPlotCov$ymax = cumsum(b_datToPlotCov$fraction)
b_datToPlotCov$ymin = c(0, head(b_datToPlotCov$ymax, n=-1))
b_datToPlotCov$labelPosition <- (b_datToPlotCov$ymax + b_datToPlotCov$ymin) / 2
# 绘制第二张图 (b_p2): 相对分布
# 准备 b_p2 的标签数据
b_datToPlotCov$label <- paste0(round(b_datToPlotCov$fraction * 100, 2), "%")
# 绘制 b_p2
b_p2 <- ggplot(b_datToPlotCov, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = type)) +
  geom_rect() +
  # 使用 geom_text_repel 替代 geom_text
  geom_text_repel(aes(x = 4.5, y = labelPosition, label = label),
                  size = 4,
                  # 调整推力，确保标签不会重叠
                  force = 0.5,
                  # 调整 x轴方向的推力
                  nudge_x = 1,
                  # 连接线样式
                  segment.color = "gray60",
                  segment.size = 0.5) +
  scale_fill_manual(
    values = colors,
    labels = c("Cores and generalists", "Cores",
               "Cores and specialists", "Specialists",
               " Uniques and specialists", "Uniques",
               " Uniques and generalists", "Generalists", "Other")
  ) +
  coord_polar(theta = "y") +
  xlim(c(1, 6)) + # 为标签留出更多空间
  theme_void() +
  theme(
    legend.title = element_blank(), # legend.position = "none", #取消图例
    legend.text = element_text(size = 14)
  ) +
  labs(title = "Relative distribution of bacteria") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = -25)))
#+guides(fill = guide_legend(nrow = 1, byrow = TRUE))
b_p2
# 组合并显示图表 (p3_bacteria)
p3_bacteria <- ggarrange(
  b_p1, b_p2,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)
# 显示最终组合图
p3_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3a_donut_bacteria.png", plot = p3_bacteria, width = 12, height = 6, dpi = 600, bg = "transparent")
#最特殊在各生态位相对丰度，各生态位相对丰度的平均值为总相对丰度
b_uniq2 <- b_rare_unique_species_raw_counts
b_species2 <- b_species_df
G_samples <- rownames(metadata[metadata$Niche == "G", ])
N_samples <- rownames(metadata[metadata$Niche == "N", ])
# G组特异species的相对丰度
b_G_rel_ab <- sum(b_uniq2[, G_samples]) / sum(b_species2[, G_samples])
b_G_rel_ab
# N组特异species的相对丰度
b_N_rel_ab <- sum(b_uniq2[, N_samples]) / sum(b_species2[, N_samples])
b_N_rel_ab
#真菌甜甜圈
#读取数据
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
# 如果 f_ASV 和 f_tax 的行名是 ASV ID，需要先将其转为列
f_ASV_df <- f_ASV %>% rownames_to_column(var = "ASV")
f_tax_df <- f_tax %>% rownames_to_column(var = "ASV")
# 将 f_ASV 与 f_tax 中的 'Species' 列连接
f_species <- f_ASV_df %>%
  left_join(select(f_tax_df, ASV, Species), by = "ASV")
# 按照 'Species' 列进行分组，并对每个样本的丰度求和
f_species <- f_species %>%
  select(all_of(c("Species", metadata$Sample.ID))) %>%
  group_by(Species) %>%
  summarise(across(all_of(metadata$Sample.ID), sum)) %>%
  ungroup() # 解除分组
# 排序
f_species <- f_species %>%
  arrange(Species)
f_species_df <- f_species # 复制操作，R中通常直接赋值即可，除非您需要独立副本
print(paste("Length of f_species_df (number of species):", nrow(f_species_df)))
# 按照 'Species' 分组，并保留每个物种的第一个分类信息（排除ASV列）
f_tax_species <- f_tax_df %>%
  group_by(Species) %>%
  slice(1) %>% # 获取每个Species组的第一个行
  ungroup() %>%
  select(-ASV) %>% # 排除ASV列
  arrange(Species) # 按照Species排序，与Python中的sort_index()对应
# 进行替换操作和去除 'p__' 前缀
f_tax_species <- f_tax_species %>%
  mutate(
    Phylum = case_when(
      Phylum == "p__unclassified_k__norank_d__Fungi" ~ "Unclassified",
      Phylum == "p__SAR324_cladeMarine_group_B" ~ "SAR324",
      TRUE ~ Phylum # 保持其他值不变
    ),
    # 去除 'p__' 前缀
    Phylum = str_replace(Phylum, "^p__", "")
  )
# 创建物种学名到 'species_数字' 格式的映射并应用
f_unique_species_names <- f_species %>% pull(Species) %>% unique() %>% sort()
# 创建一个从物种学名到 'species_数字' 格式的映射数据框
f_species_to_numeric_id <- tibble(
  original_species = f_unique_species_names,
  numeric_id = paste0("species_", 1:length(f_unique_species_names))
)
# 转换 f_species 的 'Species' 列
f_species <- f_species %>%
  left_join(f_species_to_numeric_id, by = c("Species" = "original_species")) %>%
  select(-Species) %>% # 移除原始的Species列
  rename(Species = numeric_id) %>% # 将numeric_id重命名为Species
  select(Species, everything()) # 确保Species列在第一位
# 转换 f_tax_species 的 'Species' 列
f_tax_species <- f_tax_species %>%
  left_join(f_species_to_numeric_id, by = c("Species" = "original_species")) %>%
  select(-Species) %>% # 移除原始的Species列
  rename(Species = numeric_id) %>% # 将numeric_id重命名为Species
  select(Species, everything()) # 确保Species列在第一位
f_species_df <- f_species
print(paste("Length of f_species_df after species ID mapping:", nrow(f_species_df)))
metadata <- data.frame(metadata)
#核心asv（在大于2个居群生态位中相对丰度大于0.1%）
# 计算相对丰度
f_species_matrix <- f_species_df %>%
  column_to_rownames(var = "Species") %>% # 将 'Species' 列转换为行名
  as.matrix() # 转换为矩阵以便进行高效的列操作
# 计算相对丰度：将每个样本（列）的丰度除以该样本的总丰度
f_species_rel_matrix <- t(t(f_species_matrix) / colSums(f_species_matrix))
# 将结果转换回 tibble，并将行名（物种名称）转回列
f_species_rel <- f_species_rel_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "species") %>% # 将行名转回 'species' 列
  as_tibble() # 转换为 tibble
# 将数据重塑为长格式并合并元数据
f_rel_long <- f_species_rel %>%
  pivot_longer(
    cols = all_of(metadata$Sample.ID), # 仅选择样本ID列进行重塑
    names_to = "Sample",       # 新的列，存放样本名称
    values_to = "value"        # 新的列，存放丰度值
  )
# 创建从元数据中提取 Origin 和 Niche 的映射数据框
sample_origins_map_df <- metadata %>% select(`Sample.ID`, Origin) %>% rename(Sample = `Sample.ID`)
sample_niches_map_df <- metadata %>% select(`Sample.ID`, Niche) %>% rename(Sample = `Sample.ID`)
# 通过 left_join 将 Origin 和 Niche 信息添加到长格式表中
f_rel_long <- f_rel_long %>%
  left_join(sample_origins_map_df, by = "Sample") %>%
  left_join(sample_niches_map_df, by = "Sample")
# 计算物种整体相对丰度并进行筛选 
f_species_raw_matrix_for_sum <- f_species_df %>%
  column_to_rownames(var = "Species") %>%
  as.matrix()
f_species_total_counts <- rowSums(f_species_raw_matrix_for_sum)
# 计算所有物种的总丰度
f_total_counts_across_all_species <- sum(f_species_total_counts)
# 计算每个物种的整体相对丰度（相对于所有物种的总和）
f_species_relative_abundance <- f_species_total_counts / f_total_counts_across_all_species
# 定义整体相对丰度阈值 (0.1% 转换为小数)
overall_abundance_threshold <- 0.001
# 筛选出满足整体丰度阈值的物种名称列表
f_species_above_overall_threshold <- names(f_species_relative_abundance[f_species_relative_abundance >= overall_abundance_threshold])
# 根据整体丰度阈值过滤长格式丰度表
f_rel_long_filtered_by_overall_abundance <- f_rel_long %>%
  filter(species %in% f_species_above_overall_threshold)
#计算物种在每个 Origin 和 Niche 中的出现次数
# 计算每个物种在每个 Origin 中的平均丰度
f_species_origin_mean_filtered <- f_rel_long_filtered_by_overall_abundance %>%
  group_by(species, Origin) %>%
  summarise(mean_abundance_in_origin = mean(value), .groups = 'drop') # .groups = 'drop' 用于解除分组
# 识别在特定 Origin 中平均丰度大于 0 的物种（即有出现）
f_species_in_origin <- f_species_origin_mean_filtered %>%
  filter(mean_abundance_in_origin > 0)
# 统计每个物种出现的 Origin 数量
f_species_origin_count <- f_species_in_origin %>%
  group_by(species) %>%
  summarise(origin_count = n(), .groups = 'drop')
# 计算每个物种在每个 Niche 中的平均丰度
f_species_niche_mean_filtered <- f_rel_long_filtered_by_overall_abundance %>%
  group_by(species, Niche) %>%
  summarise(mean_abundance_in_niche = mean(value), .groups = 'drop')
# 识别在特定 Niche 中平均丰度大于 0 的物种（即有出现）
f_species_in_niche <- f_species_niche_mean_filtered %>%
  filter(mean_abundance_in_niche > 0)
# 统计每个物种出现的 Niche 数量
f_species_niche_count <- f_species_in_niche %>%
  group_by(species) %>%
  summarise(niche_count = n(), .groups = 'drop')
# 合并 Origin 和 Niche 的出现次数统计
f_species_counts_combined <- f_species_origin_count %>%
  inner_join(f_species_niche_count, by = "species")
# 定义出现次数阈值
origin_occurrence_threshold <- 3
niche_occurrence_threshold <- 2
# 筛选满足所有条件（整体丰度，以及在 Origin 和 Niche 中的出现次数）的物种，得到核心物种的名称列表
f_core_species_names <- f_species_counts_combined %>%
  filter(origin_count >= origin_occurrence_threshold &
           niche_count >= niche_occurrence_threshold) %>%
  pull(species) # 提取物种名称作为向量
# 从相对丰度数据中提取核心物种的数据
f_core_abundance <- f_species_rel %>%
  filter(species %in% f_core_species_names)
# 从原始计数数据中提取核心物种的数据
f_core_species_raw_counts <- f_species_df %>%
  filter(Species %in% f_core_species_names)
# 打印核心物种的数量
print(paste("找到的核心物种数量:", nrow(f_core_species_raw_counts)))
#最特殊
# 计算所有物种整体相对丰度的中位数并筛选稀有物种 
f_median_relative_abundance <- median(f_species_relative_abundance)
print(paste("所有物种整体相对丰度的中位数:", f_median_relative_abundance))
# 定义稀有类群的丰度阈值
f_rare_threshold <- 0.000008
# 筛选出稀有物种
f_rare_species_names <- names(f_species_relative_abundance[f_species_relative_abundance < f_rare_threshold])
print(paste("稀有物种的数量:", length(f_rare_species_names)))
# 将 f_species_df 转换为长格式并添加 Group 信息
# 转置 f_species_df 并转换为长格式
f_speciesdfmelt <- f_species_df %>%
  pivot_longer(
    cols = -Species, # 排除 Species 列，将其余列转换为长格式
    names_to = "Sample",
    values_to = "value"
  ) %>%
  rename(species = Species) %>% # 重命名 Species 列为 species 以匹配 Python 代码中的命名
  filter(value > 0) # 过滤掉丰度为 0 的行
# 创建从 Sample ID 到 Group 的映射数据框
sample_group_map_df <- metadata %>% select(`Sample.ID`, Group) %>% rename(Sample = `Sample.ID`)
# 添加 Group 信息到长格式表中
f_control_species <- f_speciesdfmelt %>%
  left_join(sample_group_map_df, by = "Sample") %>%
  filter(!is.na(Group)) # 移除 Group 为 NA 的行 (对应 Python 的 dropna(subset=['Group']))
# 统计物种出现在多少个 Group 中 
# 提取 Group 和 species 列并去重
f_controlspecies_unique <- f_control_species %>%
  select(Group, species) %>%
  distinct() # 相当于 drop_duplicates()
# 统计每个 species 出现在多少个 Group 中
f_controlspecies_group_count <- f_controlspecies_unique %>%
  group_by(species) %>%
  summarise(group_count = n(), .groups = 'drop')
# 只保留那些只在一个 Group 中出现的物种
f_unique_species_single_group_names <- f_controlspecies_group_count %>%
  filter(group_count == 1) %>%
  pull(species) # 提取物种名称向量
print(paste("只在一个 Group 中出现的物种数量:", length(f_unique_species_single_group_names)))
# 找到同时满足两个条件的“最特殊稀有物种
f_rare_unique_species_names <- intersect(f_rare_species_names, f_unique_species_single_group_names)
print(paste("最特殊稀有物种（Rare AND Uniques to single group）的数量:", length(f_rare_unique_species_names)))
# 提取这些“最特殊稀有物种”的原始丰度数据
f_rare_unique_abundance <- f_species_df %>%
  filter(Species %in% f_rare_unique_species_names)
f_rare_unique_species_raw_counts <- f_species_df %>%
  filter(Species %in% f_rare_unique_species_names)
#统一读取
f_species_df <- data.frame(f_species)
rownames(f_species_df) <- f_species_df$Species
f_species_df <- f_species_df[,-1]
f_core_species_raw_counts <- data.frame(f_core_species_raw_counts)
rownames(f_core_species_raw_counts) <- f_core_species_raw_counts$Species
f_core_species_raw_counts <- f_core_species_raw_counts[,-1]
f_rare_unique_species_raw_counts <- data.frame(f_rare_unique_species_raw_counts)
rownames(f_rare_unique_species_raw_counts) <- f_rare_unique_species_raw_counts$Species
f_rare_unique_species_raw_counts <- f_rare_unique_species_raw_counts[,-1]
f_dat	<- f_species_df
f_cores <- f_core_species_raw_counts 
f_uniq <- f_rare_unique_species_raw_counts
f_generalists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/f_generalists_raw_counts.csv", header = TRUE, row.names = 1)
f_specialists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/f_specialists_raw_counts.csv", header = TRUE, row.names = 1)
f_species <- rownames(f_dat)
f_species <- as.data.frame(f_species)
colnames(f_species) <- "species"
f_cores$species<-rownames(f_cores)
rownames(f_cores)<-NULL
f_cores <- f_cores %>%
  select(species)%>%# 只保留 species 这一列
  mutate(type = "core")#增加一列 type，所有行都标记为 "core"
f_uniq$species<-rownames(f_uniq)
rownames(f_uniq)<-NULL
f_uniq <- f_uniq %>%
  select(species)%>%
  mutate(type = "uniq")
colnames(f_uniq) <- c("species","type")
f_generalists$species<-rownames(f_generalists)
rownames(f_generalists)<-NULL
f_generalists <- f_generalists %>%
  select(species)%>%
  mutate(type = "generalists")
colnames(f_generalists) <- c("species", "type")
f_specialists$species<-rownames(f_specialists)
rownames(f_specialists)<-NULL
f_specialists <- f_specialists %>%
  select(species)%>%
  mutate(type = "specialists")
colnames(f_specialists) <- c("species","type")
f_dat_numb <- rbind(f_cores, f_uniq, f_generalists, f_specialists)#合并
f_species_other <- f_species %>%#筛选出那些没有出现在四类标签里的species
  filter(!(species %in% f_dat_numb$species))%>%
  mutate(type = "Other")
#计算每类有多少种species
f_datToPlot <- f_dat_numb %>%
  arrange(species, type) %>%#先对 species 和 type 排序，确保拼接顺序一致
  group_by(species)%>%# 每个species为一组
  summarise(type = str_c(type, collapse="_"))%>%#将同一个species的多个 type 合并成一个字符串，比如 "core_indicator"
  ungroup()%>%#取消分组状态
  bind_rows(f_species_other)%>%#合并"Other" 的species
  group_by(type)%>%# 按 type 统计每一类species的种类数量
  summarise(sum = n())
f_datToPlot$type <- factor(f_datToPlot$type, levels = c("core_generalists","core",  "core_specialists",  "specialists", "specialists_uniq", "uniq","generalists_uniq","generalists", "Other"))#设置分类顺序与颜色
colors <- c(
  "core_generalists" = "#4AC0FF",
  "core" = "#5499C7",
  "core_specialists" = "#2F609F",
  "specialists" = "#45B39D",
  "specialists_uniq" = "#006400",
  "uniq" = "#D4AC0D",
  "generalists_uniq" = "#30308B",
  "generalists" = "#EC7063",
  "Other" = "#737373"
)
f_datToPlot <- f_datToPlot[order(f_datToPlot$type),]#排序
f_datToPlot$fraction = f_datToPlot$sum / sum(f_datToPlot$sum)#计算每一类所占百分比
f_datToPlot$ymax = cumsum(f_datToPlot$fraction)#计算累积百分比
f_datToPlot$ymin = c(0, head(f_datToPlot$ymax, n=-1))# 计算累积最小值
f_datToPlot$labelPosition <- (f_datToPlot$ymax + f_datToPlot$ymin) / 2#计算标签的中心点
# 绘制第一张图 (f_p1): 绝对分布
# 准备 f_p1 的标签数据
f_datToPlot$label <- paste0(f_datToPlot$sum)
# 绘制 f_p1
f_p1 <- ggplot(f_datToPlot, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = type)) +
  geom_rect() +
  # 使用 geom_text_repel 替代 geom_text
  geom_text_repel(aes(x = 4.5, y = labelPosition, label = label),
                  size = 4,
                  # 调整推力，确保标签不会重叠
                  force = 0,
                  # 调整 x轴方向的推力
                  nudge_x = 0,
                  # 连接线样式
                  segment.color = "gray60",
                  segment.size = 0) +
  scale_fill_manual(
    values = colors,
    labels = c("Cores and generalists", "Cores",
               "Cores and specialists", "Specialists",
               " Uniques and specialists", "Uniques",
               " Uniques and generalists", "Generalists", "Other")
  ) +
  coord_polar(theta = "y") +
  xlim(c(1, 6)) + # 为标签留出更多空间
  theme_void() +
  theme(
    legend.position = "none", #取消图例
    legend.text = element_text(size = 14)
  ) +
  labs(title = "Absolute distribution of fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = -25))) 
f_p1
f_dat_sum <- rowSums(f_dat[1:48]) # 对每个 species 求所有样本丰度之和
f_dat_sum <- as.data.frame(f_dat_sum)
f_dat_sum$species <- rownames(f_dat_sum) # 加上 species 名称
#计算每类的相对丰度
f_datToPlotCov <- f_dat_numb %>%
  arrange(species, type) %>%
  group_by(species)%>%
  summarise(type = str_c(type, collapse="_"))%>% # 合并同一个species多个type
  ungroup()%>%
  bind_rows(f_species_other)%>% 
  left_join(f_dat_sum)%>%# 添加丰度信息（f_dat_sum 是每个 species 的总丰度）
  filter(!(f_dat_sum == 0))%>% # 移除丰度为0的species
  mutate(perc = f_dat_sum / sum(f_dat_sum))%>%#计算相对丰度
  group_by(type)%>%
  summarise(fraction = sum(perc)) # 统计每种类型的丰度总占比
sum(f_datToPlotCov$fraction)
f_datToPlotCov$type <- factor(f_datToPlotCov$type, levels = c("core_generalists","core", "core_specialists", "specialists", "specialists_uniq", "uniq","generalists_uniq","generalists", "Other"))
f_datToPlotCov <- f_datToPlotCov[order(f_datToPlotCov$type),]
# f_datToPlotCov$fraction = f_datToPlotCov$sum / sum(f_datToPlotCov$sum)
f_datToPlotCov$ymax = cumsum(f_datToPlotCov$fraction)
f_datToPlotCov$ymin = c(0, head(f_datToPlotCov$ymax, n=-1))
f_datToPlotCov$labelPosition <- (f_datToPlotCov$ymax + f_datToPlotCov$ymin) / 2
# 绘制第二张图 (f_p2): 相对分布
# 准备 f_p2 的标签数据
f_datToPlotCov$label <- paste0(round(f_datToPlotCov$fraction * 100, 2), "%")
# 绘制 f_p2
f_p2 <- ggplot(f_datToPlotCov, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = type)) +
  geom_rect() +
  # 使用 geom_text_repel 替代 geom_text
  geom_text_repel(aes(x = 4.5, y = labelPosition, label = label),
                  size = 4,
                  # 调整推力，确保标签不会重叠
                  force = 0.5,
                  # 调整 x轴方向的推力
                  nudge_x = 1,
                  # 连接线样式
                  segment.color = "gray60",
                  segment.size = 0.5) +
  scale_fill_manual(
    values = colors,
    labels = c("Cores and generalists", "Cores",
               "Cores and specialists", "Specialists",
               " Uniques and specialists", "Uniques",
               " Uniques and generalists", "Generalists", "Other")
  ) +
  coord_polar(theta = "y") +
  xlim(c(1, 6)) + # 为标签留出更多空间
  theme_void() +
  theme(
    legend.position = "none", #取消图例
    legend.text = element_text(size = 14)
  ) +
  labs(title = "Relative distribution of fungi") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = -25)))
f_p2
# 组合并显示图表 (p3_fungi)
p3_fungi <- ggarrange(
  f_p1, f_p2,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)
# 显示最终组合图
p3_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3a_donut_fungi.png", plot = p3_fungi, width = 12, height = 6, dpi = 600, bg = "transparent")
#最特殊在各生态位相对丰度，各生态位相对丰度的平均值为总相对丰度
f_uniq2 <- f_rare_unique_species_raw_counts
f_species2 <- f_species_df
G_samples <- rownames(metadata[metadata$Niche == "G", ])
N_samples <- rownames(metadata[metadata$Niche == "N", ])
# G组特异species的相对丰度
f_G_rel_ab <- sum(f_uniq2[, G_samples]) / sum(f_species2[, G_samples])
f_G_rel_ab
# N组特异species的相对丰度
f_N_rel_ab <- sum(f_uniq2[, N_samples]) / sum(f_species2[, N_samples])
f_N_rel_ab

#图5跨界网络
install.packages("BiocManager")
library(BiocManager)
install("remotes")
install("tidyverse")
install("tidyfst")
install("igraph")
install("sna")
install("phyloseq")
install("ggalluvial")
install("ggraph")
install("WGCNA")
install("ggnewscale")
install("pulsar")
install("patchwork")
remotes::install_github("taowenmicro/EasyStat")
remotes::install_github("taowenmicro/ggClusterNet")
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
library(RColorBrewer)
#细菌真菌跨界
#创建文件夹
Envnetplot1<- paste("D:/study/master/Main_Figure_tables/Figure_5/16S_ITS_network",sep = "")
dir.create(Envnetplot1)
#读取细菌数据
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
b_tax_species <- b_tax %>% rownames_to_column(var = "ASV")
# 将 b_ASV 与 b_tax 中的 'Species' 列连接
b_species <- b_ASV_df %>%
  left_join(select(b_tax_species, ASV, Species), by = "ASV")
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
b_tax_species <- b_tax_species %>%
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
  # 不删除原始的Species列
  # 将 numeric_id 列重命名为 ID
  rename(ID = numeric_id) %>%
  # 将 ID 列设置为行名
  column_to_rownames(var = "ID")
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
  rownames_to_column(var = "Species") %>% # 将行名转回 'species' 列
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
  filter(Species %in% b_species_above_overall_threshold)
#计算物种在每个 Origin 和 Niche 中的出现次数
# 计算每个物种在每个 Origin 中的平均丰度
b_species_origin_mean_filtered <- b_rel_long_filtered_by_overall_abundance %>%
  group_by(Species, Origin) %>%
  summarise(mean_abundance_in_origin = mean(value), .groups = 'drop') # .groups = 'drop' 用于解除分组
# 识别在特定 Origin 中平均丰度大于 0 的物种（即有出现）
b_species_in_origin <- b_species_origin_mean_filtered %>%
  filter(mean_abundance_in_origin > 0)
# 统计每个物种出现的 Origin 数量
b_species_origin_count <- b_species_in_origin %>%
  group_by(Species) %>%
  summarise(origin_count = n(), .groups = 'drop')
# 计算每个物种在每个 Niche 中的平均丰度
b_species_niche_mean_filtered <- b_rel_long_filtered_by_overall_abundance %>%
  group_by(Species, Niche) %>%
  summarise(mean_abundance_in_niche = mean(value), .groups = 'drop')
# 识别在特定 Niche 中平均丰度大于 0 的物种（即有出现）
b_species_in_niche <- b_species_niche_mean_filtered %>%
  filter(mean_abundance_in_niche > 0)
# 统计每个物种出现的 Niche 数量
b_species_niche_count <- b_species_in_niche %>%
  group_by(Species) %>%
  summarise(niche_count = n(), .groups = 'drop')
# 合并 Origin 和 Niche 的出现次数统计
b_species_counts_combined <- b_species_origin_count %>%
  inner_join(b_species_niche_count, by = "Species")
# 定义出现次数阈值
origin_occurrence_threshold <- 3
niche_occurrence_threshold <- 2
# 筛选满足所有条件（整体丰度，以及在 Origin 和 Niche 中的出现次数）的物种，得到核心物种的名称列表
b_core_species_names <- b_species_counts_combined %>%
  filter(origin_count >= origin_occurrence_threshold &
           niche_count >= niche_occurrence_threshold) %>%
  pull(Species) # 提取物种名称作为向量
# 从相对丰度数据中提取核心物种的数据
b_core_abundance <- b_species_rel %>%
  filter(Species %in% b_core_species_names)
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
  rename(Species = Species) %>% # 重命名 Species 列为 species 以匹配 Python 代码中的命名
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
  select(Group, Species) %>%
  distinct() # 相当于 drop_duplicates()
# 统计每个 species 出现在多少个 Group 中
b_controlspecies_group_count <- b_controlspecies_unique %>%
  group_by(Species) %>%
  summarise(group_count = n(), .groups = 'drop')
# 只保留那些只在一个 Group 中出现的物种
b_unique_species_single_group_names <- b_controlspecies_group_count %>%
  filter(group_count == 1) %>%
  pull(Species) # 提取物种名称向量
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
b_species_rel <- data.frame(b_species_rel)
rownames(b_species_rel) <- b_species_rel$Species
b_species_rel <- b_species_rel[,-1]
b_core_abundance <- data.frame(b_core_abundance)
rownames(b_core_abundance) <- b_core_abundance$Species
b_core_abundance <- b_core_abundance[,-1]
b_rare_unique_abundance <- data.frame(b_rare_unique_abundance)
rownames(b_rare_unique_abundance) <- b_rare_unique_abundance$Species
b_rare_unique_abundance <- b_rare_unique_abundance[,-1]
b_generalists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/b_generalists_raw_counts.csv", header = TRUE, row.names = 1)
b_specialists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/b_specialists_raw_counts.csv", header = TRUE, row.names = 1)
# 提取 b_generalists 的行名
b_generalist_rownames <- rownames(b_generalists)
# 从 b_species_rel 中筛选出具有相同行名的行
b_generalist_rel <- b_species_rel[rownames(b_species_rel) %in% b_generalist_rownames, , drop = FALSE]
# 提取 b_specialists 的行名
b_specialist_rownames <- rownames(b_specialists)
# 从 b_species_rel 中筛选出具有相同行名的行
b_specialist_rel <- b_species_rel[rownames(b_species_rel) %in% b_specialist_rownames, , drop = FALSE]
# 合并为 phyloseq 对象
b_species_rel <- otu_table(b_species_rel, taxa_are_rows = TRUE)
b_tax_species <- tax_table(as.matrix(b_tax_species))
metadata <- sample_data(metadata)
b_merged_rel <- merge_phyloseq(b_species_rel, b_tax_species, metadata)
#找出相对丰度总和小于阈值为other门
b_merged_phylum_ra <- b_merged_rel
#聚合到门 (Phylum) 水平
b_merged_phylum <- tax_glom(b_merged_phylum_ra, taxrank = "Phylum", NArm = FALSE)
# 提取并处理数据框 
b_phylum_ra_df <- psmelt(b_merged_phylum)
# 确保 Phylum 列是字符类型
b_phylum_ra_df$Phylum <- as.character(b_phylum_ra_df$Phylum)
# 计算当前 Phylum 级别的平均相对丰度
b_phylum_mean_ra <- b_phylum_ra_df %>%
  group_by(Phylum) %>% 
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选出最终的“Top Phyla”列表
b_top_phyla_final <- b_phylum_mean_ra %>%
  filter(MeanAbundance > 0.001) %>%
  pull(Phylum)
#最终的 Phylum 分组：将所有不符合条件的门归入“Other”
b_tax_species <- data.frame(b_tax_species) %>%
  mutate(Phylum = ifelse(
    Phylum %in% b_top_phyla_final,
    Phylum, # 如果是 Top Phyla，则保留其名称
    "Other"             # 否则，归类为最终的“Other”
  ))
b_tax_species <- tax_table(as.matrix(b_tax_species))
b_merged_rel <- merge_phyloseq(b_species_rel, b_tax_species, metadata)
#读取真菌数据
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
f_tax_species <- f_tax %>% rownames_to_column(var = "ASV")
# 将 f_ASV 与 f_tax 中的 'Species' 列连接
f_species <- f_ASV_df %>%
  left_join(select(f_tax_species, ASV, Species), by = "ASV")
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
f_tax_species <- f_tax_species %>%
  group_by(Species) %>%
  slice(1) %>% # 获取每个Species组的第一个行
  ungroup() %>%
  select(-ASV) %>% # 排除ASV列
  arrange(Species) # 按照Species排序，与Python中的sort_index()对应
# 进行替换操作和去除 'p__' 前缀
f_tax_species <- f_tax_species %>%
  mutate(
    Phylum = case_when(
      Phylum == "p__unclassified_k__Fungi" ~ "Unclassified",
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
  # 不删除原始的Species列
  # 将 numeric_id 列重命名为 ID
  rename(ID = numeric_id) %>%
  # 将 ID 列设置为行名
  column_to_rownames(var = "ID")
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
  rownames_to_column(var = "Species") %>% # 将行名转回 'species' 列
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
  filter(Species %in% f_species_above_overall_threshold)
#计算物种在每个 Origin 和 Niche 中的出现次数
# 计算每个物种在每个 Origin 中的平均丰度
f_species_origin_mean_filtered <- f_rel_long_filtered_by_overall_abundance %>%
  group_by(Species, Origin) %>%
  summarise(mean_abundance_in_origin = mean(value), .groups = 'drop') # .groups = 'drop' 用于解除分组
# 识别在特定 Origin 中平均丰度大于 0 的物种（即有出现）
f_species_in_origin <- f_species_origin_mean_filtered %>%
  filter(mean_abundance_in_origin > 0)
# 统计每个物种出现的 Origin 数量
f_species_origin_count <- f_species_in_origin %>%
  group_by(Species) %>%
  summarise(origin_count = n(), .groups = 'drop')
# 计算每个物种在每个 Niche 中的平均丰度
f_species_niche_mean_filtered <- f_rel_long_filtered_by_overall_abundance %>%
  group_by(Species, Niche) %>%
  summarise(mean_abundance_in_niche = mean(value), .groups = 'drop')
# 识别在特定 Niche 中平均丰度大于 0 的物种（即有出现）
f_species_in_niche <- f_species_niche_mean_filtered %>%
  filter(mean_abundance_in_niche > 0)
# 统计每个物种出现的 Niche 数量
f_species_niche_count <- f_species_in_niche %>%
  group_by(Species) %>%
  summarise(niche_count = n(), .groups = 'drop')
# 合并 Origin 和 Niche 的出现次数统计
f_species_counts_combined <- f_species_origin_count %>%
  inner_join(f_species_niche_count, by = "Species")
# 定义出现次数阈值
origin_occurrence_threshold <- 3
niche_occurrence_threshold <- 2
# 筛选满足所有条件（整体丰度，以及在 Origin 和 Niche 中的出现次数）的物种，得到核心物种的名称列表
f_core_species_names <- f_species_counts_combined %>%
  filter(origin_count >= origin_occurrence_threshold &
           niche_count >= niche_occurrence_threshold) %>%
  pull(Species) # 提取物种名称作为向量
# 从相对丰度数据中提取核心物种的数据
f_core_abundance <- f_species_rel %>%
  filter(Species %in% f_core_species_names)
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
  rename(Species = Species) %>% # 重命名 Species 列为 species 以匹配 Python 代码中的命名
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
  select(Group, Species) %>%
  distinct() # 相当于 drop_duplicates()
# 统计每个 species 出现在多少个 Group 中
f_controlspecies_group_count <- f_controlspecies_unique %>%
  group_by(Species) %>%
  summarise(group_count = n(), .groups = 'drop')
# 只保留那些只在一个 Group 中出现的物种
f_unique_species_single_group_names <- f_controlspecies_group_count %>%
  filter(group_count == 1) %>%
  pull(Species) # 提取物种名称向量
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
f_species_rel <- data.frame(f_species_rel)
rownames(f_species_rel) <- f_species_rel$Species
f_species_rel <- f_species_rel[,-1]
f_core_abundance <- data.frame(f_core_abundance)
rownames(f_core_abundance) <- f_core_abundance$Species
f_core_abundance <- f_core_abundance[,-1]
f_rare_unique_abundance <- data.frame(f_rare_unique_abundance)
rownames(f_rare_unique_abundance) <- f_rare_unique_abundance$Species
f_rare_unique_abundance <- f_rare_unique_abundance[,-1]
f_generalists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/f_generalists_raw_counts.csv", header = TRUE, row.names = 1)
f_specialists <- read.csv("D:/study/master/Main_Figure_tables/Figure_3/f_specialists_raw_counts.csv", header = TRUE, row.names = 1)
# 提取 f_generalists 的行名
f_generalist_rownames <- rownames(f_generalists)
# 从 f_species_rel 中筛选出具有相同行名的行
f_generalist_rel <- f_species_rel[rownames(f_species_rel) %in% f_generalist_rownames, , drop = FALSE]
# 提取 f_specialists 的行名
f_specialist_rownames <- rownames(f_specialists)
# 从 f_species_rel 中筛选出具有相同行名的行
f_specialist_rel <- f_species_rel[rownames(f_species_rel) %in% f_specialist_rownames, , drop = FALSE]
# 合并为 phyloseq 对象
f_species_rel <- otu_table(f_species_rel, taxa_are_rows = TRUE)
f_tax_species <- tax_table(as.matrix(f_tax_species))
metadata <- sample_data(metadata)
f_merged_rel <- merge_phyloseq(f_species_rel, f_tax_species, metadata)
#找出相对丰度总和小于阈值为other门
f_merged_phylum_ra <- f_merged_rel
#聚合到门 (Phylum) 水平
f_merged_phylum <- tax_glom(f_merged_phylum_ra, taxrank = "Phylum", NArm = FALSE)
# 提取并处理数据框 
f_phylum_ra_df <- psmelt(f_merged_phylum)
# 确保 Phylum 列是字符类型
f_phylum_ra_df$Phylum <- as.character(f_phylum_ra_df$Phylum)
# 计算当前 Phylum 级别的平均相对丰度
f_phylum_mean_ra <- f_phylum_ra_df %>%
  group_by(Phylum) %>% 
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选出最终的“Top Phyla”列表
f_top_phyla_final <- f_phylum_mean_ra %>%
  filter(MeanAbundance > 0.001) %>%
  pull(Phylum)
#最终的 Phylum 分组：将所有不符合条件的门归入“Other”
f_tax_species <- data.frame(f_tax_species) %>%
  mutate(Phylum = ifelse(
    Phylum %in% f_top_phyla_final,
    Phylum, # 如果是 Top Phyla，则保留其名称
    "Other"             # 否则，归类为最终的“Other”
  ))
f_tax_species <- tax_table(as.matrix(f_tax_species))
f_merged_rel <- merge_phyloseq(f_species_rel, f_tax_species, metadata)
# 确保细菌 (b_merged_rel) 和真菌 (f_merged_rel) 的 phyloseq 对象拥有相同的样本元数据 (map file)。
# 这一步是关键，因为 merge16S_ITS 函数需要两个 phyloseq 对象具有一致的样本信息。
ps.merge <- ggClusterNet::merge16S_ITS(
  ps16s = b_merged_rel, # 细菌的 phyloseq 对象
  psITS = f_merged_rel, # 真菌的 phyloseq 对象
  N16s = 2500, # 可选参数：如果需要，可以指定16S数据保留的最高丰度OTU/ASV数量
  NITS = 2500  # 可选参数：如果需要，可以指定ITS数据保留的最高丰度OTU/ASV数量
)
# 打印合并后的 phyloseq 对象概览
ps.merge
# 定义一个包含12种基础颜色的向量，这些颜色来自 RColorBrewer 的 "Paired" 调色板
base_colors <- c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
                 "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
                 "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                 "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                 "#FFFF99","#33A02C")
# 确保 'Other' 是 Phylum 列中的最后一个因子水平
# 首先获取所有独特的 Phylum 名称
b_Phylum <- data.frame(b_tax_species)
b_Phylum <- sort(unique(b_Phylum$Phylum))
f_Phylum <- data.frame(f_tax_species)
f_Phylum <- sort(unique(f_Phylum$Phylum))
all_unique_phylums <- unique(c(b_Phylum, f_Phylum))
# 将 'Other' 从列表中移除
phylums_sorted_without_other <- setdiff(all_unique_phylums, "Other")
# 将 'Other' 放到排序后的 Phylum 列表的最后
ordered_phylum_levels <- c(phylums_sorted_without_other, "Other")
# 为 Phylum 类别创建自定义颜色映射
# 生成除 "Other" 之外的其他 Phylum 的颜色
generated_non_other_colors <- colorRampPalette(base_colors)(length(phylums_sorted_without_other))
# 创建一个命名颜色向量，将所有颜色与 Phylum 名称对应起来
# 将生成的颜色分配给非 "Other" 的 Phylum
phylum_color_map <- generated_non_other_colors
names(phylum_color_map) <- phylums_sorted_without_other
# 为 "Other" 分配特定颜色
phylum_color_map["Other"] <-"#cccccc"
corBionetwork.st = function(
    ps.st= ps.merge,# phyloseq对象
    g1 = "Group",# 分组1
    g2 = NULL,# 分组2
    g3 = NULL,# 分组3
    ord.g1 = NULL, # 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL, # 排序顺序
    order = NULL, # 出图每行代表的变量
    fill = "filed",
    size = "igraph.degree",
    method = "spearman",
    lab = NULL,
    label = TRUE,
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,
    layout_net = "model_maptree2",
    r.threshold=0.8,
    p.threshold=0.01,
    maxnode = 5,
    N= 500,
    scale = TRUE,
    env = NULL,
    bio = TRUE,
    minsize = 4,
    maxsize = 14,
    group.node = NULL,
    model.node = FALSE
){
  
  if (scale) {
    ps.st  = ps.st %>% ggClusterNet::scale_micro()
  }
  
  
  ps.all = ps.st
  map = sample_data(ps.all)
  
  # g2 = NULL
  if (is.null(g2)) {
    sp = ""
  } else if (is.null(ord.g2)){
    sp = map[,g2] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    sp = ord.g2
    print(sp)
  }
  # g3 = NULL
  if (is.null(g3)) {
    ti = ""
  } else if (is.null(ord.g3)){
    ti = map[,g3] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    ti = ord.g3
    print(ti)
  }
  
  if (is.null(ord.g1)) {
    group = map[,g1] %>% as.matrix() %>% as.vector() %>% unique()
    
  } else{
    group = ord.g1
    print(group)
  }
  
  
  
  #-构造两两组合全部情况
  for (i in 1:length(sp)) {
    dat = data.frame(g2 = sp[i],g3 = ti)
    
    if (i ==1) {
      dat.f = dat
    } else{
      dat.f = rbind(dat.f,dat)
    }
    
  }
  
  if (!is.null(env)) {
    colnames(env)[1] = "ID"
    env_sub <-  env[match(mapsub$ID,env$ID),]
    head(env_sub)
  }
  
  
  j = 1
  n = 1
  
  cor.all = list()
  # 存储不同分组的相关矩阵拮即可
  for (j in 1:nrow(dat.f)) {
    if (dat.f[j,1] == "") {
      ps.t = ps.all
    } else{
      ps.t = ps.all %>% subset_samples.wt(g2,dat.f[j,1])
    }
    
    if (dat.f[j,2] == "") {
      ps.f = ps.t
    } else{
      ps.f = ps.t  %>% subset_samples.wt(g3,dat.f[j,2])
    }
    
    
    for (n in 1:length(group)) {
      map = sample_data(ps.f)
      head(map)
      map$Group
      ps.g = ps.f  %>% subset_samples.wt(g1,group[n]) %>% filter_OTU_ps(N)
      big = TRUE
      if (bio) {
        if (!is.null(env)) {
          
          if (big) {
            
            result <- corBiostripeBig(data =  env_sub,
                                      group = envGroup,
                                      ps = ps.g,
                                      r.threshold = r.threshold,
                                      p.threshold = p.threshold,
                                      method = method)
          } else {
            result <- corBiostripe(data =  env_sub,
                                   group = envGroup,
                                   ps = ps.g,
                                   r.threshold = r.threshold,
                                   p.threshold = p.threshold,
                                   method = method)
          }
          
          
          #-- extract cor matrix
          occor.r = result[[1]]
          tax = as.data.frame((ggClusterNet::vegan_tax(ps.g)))
          
          if (length(tax$filed) != 0) {
            group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
          } else {
            group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
          }
          
          colnames(envGroup) <-c("SampleID","Group")
          netClu = rbind(envGroup,group2)
          colnames(netClu) <- c("ID","group")
        } else {
          
          if (big) {
            result <- corBiostripeBig(ps = ps.g,
                                      r.threshold = r.threshold, p.threshold = p.threshold, method = method)
          } else {
            result <- corBiostripe(ps = ps.g,
                                   r.threshold = r.threshold, p.threshold = p.threshold, method = method)
          }
          
          #-- extract cor matrix
          occor.r = result[[1]]
          tax = as.data.frame((ggClusterNet::vegan_tax(ps.g)))
          
          if (length(tax$filed) != 0) {
            group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
          } else {
            group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
          }
          
          netClu = group2
          colnames(netClu) <- c("ID","group")
          
        }
        
      }
      
      # head(netClu)
      
      cor = result[[1]]
      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[tem]] = cor
      # gru.all[[tem]] = cor
      
      if (!is.null(group.node)) {
        group.node = group.node %>% dplyr::filter(ID %in% row.names(cor))
        netClu = group.node
        
      }
      
      res = node.edge(
        cor = cor,
        select_layout = TRUE,
        clu_method=clu_method,
        layout_net = layout_net,
        group.node = group.node,
        model.node = model.node
      )
      
      nod = res[[1]]
      nod$group = tem
      # head(nod)
      # ggplot(nod) + geom_point(aes(X1,X2))
      
      edg = res[[2]]
      edg$group = tem
      
      netClu$group2 = tem
      head(nod)
      head(edg)
      
      
      if (j ==1 & n == 1) {
        node = nod
        edge = edg
        netClu2 = netClu
      } else{
        node = rbind(node,nod)
        edge = rbind(edg,edge)
        netClu2 = rbind(netClu,netClu2)
      }
      
    }
    
  }
  
  #   head(edge)
  #   head(node)
  #  edge$group %>% unique()
  #  netClu2$group2%>% unique()
  # netClu2$group = as.factor(netClu2$group)
  
  #-统计边的节点数量 node link
  tem = edge$group %>% table() %>% as.data.frame()
  colnames(tem) = c("group","links")
  i = 1
  id = edge$group %>% unique()
  aa = c()
  for (i in 1:length(id)) {
    aa[i] = edge %>% filter(group == id[i]) %>%
      select("OTU_2", "OTU_1") %>% as.matrix() %>%
      as.vector() %>% unique() %>% length()
  }
  tem2 = data.frame(group = id,nodes = aa)
  
  tem3 = tem %>% full_join(tem2,by = "group")
  tem3$label= paste(tem3$group,": (nodes: ",
                    tem3$nodes,"; links: ",tem3$links,")",sep = "")
  
  
  if (is.null(order)){
    order = ""
    
  }
  
  #-行--空间
  if (order == "space"|order == "g2") {
    row.id = g3
    row.num = length(group) * length(ti)
    a = c()
    for (i in 1:length(sp)) {
      for (j in 1:length(ti)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }
    
  } else if (order == "time"|order == "g3") {
    row.id = g2
    row.num = length(group) * length(sp)
    a = c()
    for (j in 1:length(ti)) {
      for (i in 1:length(sp)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }
  } else{
    a = NULL
    row.num = length(group)
  }
  
  if (!is.null(a)) {
    node$group = factor(node$group,levels = a)
    edge$group = factor(edge$group,levels = a)
    tem3 = tem3[match(a,tem3$group),]
  }else{
    node$group = factor(node$group)
    edge$group = factor(edge$group)
  }
  
  head(edge)
  head(node)
  
  # ggplot(node %>% filter(group == "..Group1")) +geom_point(aes(X1,X2))
  
  tax = ps.st %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("elements")
  node = node %>% left_join(tax,by = "elements")
  
  
  tem3$label = factor(tem3$label,levels = tem3$label)
  edge = edge %>% left_join(tem3,by = "group")
  head(edge)
  edge$label = factor(edge$label,levels = as.character(tem3$label))
  
  head(node)
  
  node = node %>% left_join(tem3,by = "group")
  node$label = factor(node$label,levels = as.character(tem3$label))
  
  
  
  net.dat = list(
    cortab = cor.all,
    node = node,
    edge = edge
  )
  
  
  # 确保 'Other' 是 Phylum 列中的最后一个因子水平
  node$Phylum <- factor(node$Phylum, levels = ordered_phylum_levels)#
  col1 = c("#b10026","#4AC0FF")
  names(col1) = c("+","-")
  p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.3,alpha = 0.5) +
    geom_point(aes(x = X1, y = X2,size = !!sym(size), fill = Phylum, shape = !!sym(fill)),data =  node,alpha = 0.5) +#fill = Ecotype
    #scale_colour_brewer(palette = "Set1") +
    scale_color_manual(values = col1, labels = c("Negative", "Positive")) +
    scale_shape_manual(values = c(21,25),labels = c("Bacteria", "Fungi")) + 
    scale_size(range = c(minsize, maxsize)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_manual(values = phylum_color_map) +#values = brewer.pal(3, "Pastel1")
    # labs( title = current_title) +
    #facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    theme_void()+
    guides(
      fill = guide_legend(title = "Phylum", override.aes = list(shape = 22, size = 5), order = 2), # 增加 size = 5# title = "Ecotype"
      shape = guide_legend(title = "Microbiomes", override.aes = list(size = 5), order = 1),       # 增加 size = 5
      color = guide_legend(title = "Edge", override.aes = list(size = 8), order = 3),             # 增加 size = 1 (边的粗细)
      size = "none"
    ) +
    theme(
      #plot.title = element_text(size = 60, hjust = 0.5
      legend.title = element_text(size = 24), # 调整图例标题的大小，例如 14
      legend.text = element_text(size = 22),   # 调整图例内容的大小，例如 12
      legend.spacing.y = unit(5.5, "cm")
    )
  
  p0  
  
  tem2 = node %>% dplyr::filter(elements %in% lab[[1]])
  if (label == TRUE ) {
    
    if (!is.null(lab)) {
      p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = tem2)
    } else {
      p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = node)
    }
    
  }
  
  return(list(network.plot = p0,network.data = net.dat,p1))
  
}
# 提取合并后的 phyloseq 对象的样本元数据
map = phyloseq::sample_data(ps.merge)
# 查看元数据的前几行
head(map)
# 添加一个新的样本分组列 "Group1"，并将其所有值设置为 "one"。
# 这通常用于创建一个单一的整体网络，而不是按特定分组进行子网络分析。
# 如果有其他分组变量，可以替换 "one" 或添加更多条件。
map$Group1 = "one"
# 将修改后的元数据重新赋值给 ps.merge 对象
phyloseq::sample_data(ps.merge) <- map
map = phyloseq::sample_data(ps.merge)
ps.merge %>% vegan_tax()
result <- corBionetwork.st(
  ps.st= ps.merge,
  g1 = "Group1",
  g2 = NULL,
  g3 = NULL,
  ord.g1 = NULL,
  ord.g2 = NULL,
  ord.g3 = NULL,
  order = NULL,
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_Gephi.3",
  r.threshold=0.8,
  p.threshold=0.01,
  maxnode = 5,
  N= 4000,scale = TRUE,env = NULL,
  bio = TRUE,minsize = 4,maxsize = 14)
summary(result[[2]])
# 提取 corBionetwork 返回的绘图对象 (ggplot2 对象)
p = result[[1]]
# 显示绘图 (如果在一个交互式环境中，这会显示图表)
p
# 构建输出PDF文件的完整路径和文件名
plotname1 = paste(Envnetplot1,"/network_all.pdf",sep = "")
# 保存网络图为PDF文件
ggsave(plotname1, p, width = 24, height = 19)
# 使用 model_maptree 函数对网络相关性矩阵进行社区检测和布局
tem <- model_maptree(
  cor = result$network.data$cortab$..one,          # 输入相关性矩阵
  method = "cluster_fast_greedy", # 社区检测方法
  seed = 12                     # 设置随机种子，确保结果可重复
)
# 提取社区检测后的节点模型信息 (通常包含每个节点的社区ID等)
node_model = tem[[2]]
# 查看节点模型的前几行
head(node_model)
# 构建输出节点模型信息CSV文件的完整路径和文件名
tablename <- paste(Envnetplot1,"/node_model_imformation",".csv",sep = "")
# 将节点模型信息写入CSV文件
write.csv(node_model,tablename)
# 构建输出节点绘图数据CSV文件的完整路径和文件名
tablename <- paste(Envnetplot1,"/nodeG_plot",".csv",sep = "")
# 将用于绘图的节点数据写入CSV文件 (通常包含坐标、颜色、大小等信息)
write.csv(result$network.data$node,tablename)
# 构建输出边绘图数据CSV文件的完整路径和文件名
tablename <- paste(Envnetplot1,"/edge_plot",".csv",sep = "")
# 将用于绘图的边数据写入CSV文件 (通常包含起点、终点、相关性、颜色等信息)
write.csv(result$network.data$edge,tablename)
# 构建输出相关性矩阵CSV文件的完整路径和文件名
tablename <- paste(Envnetplot1,"/cor_matrix",".csv",sep = "")
# 将计算出的相关性矩阵写入CSV文件
write.csv(result$network.data$cortab,tablename)
#指向包含 edge_plot.csv 的文件夹。
tablename_edge_plot <- paste(Envnetplot1,"/edge_plot",".csv",sep = "")
# 读取表格，计算边密度
edge_data <- read.csv(file = tablename_edge_plot, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# 提取节点数和实际边数
n <- edge_data$nodes[1]
L <- edge_data$links[1]
# 计算最大可能边数
if (n > 1) {
  L_max <- (n * (n - 1)) / 2
} else {
  L_max <- 0
}
# 计算边密度
if (L_max > 0) {
  edge_density_df <- L / L_max
} else {
  edge_density_df <- 0
}
# 将结果写入 CSV 文件
write.csv(edge_density_df, paste(Envnetplot1,"/edge_density",".csv",sep = ""), row.names = FALSE)
#分组
# 定义输入和输出文件夹路径
input_base_dir <- "D:/study/master/Main_Figure_tables/Figure_5/"
output_base_dir <- "D:/study/master/Main_Figure_tables/Figure_5/"
# 定义每个处理的子文件夹名称
treatment_folders <- c("G", "N","JR", "JJ", "TZ", "PA")
output_suffixes <- c("g", "n", "jr","jj", "tz", "pa")
title <- c( "Rhizosphere Soil", "Bulb","Jurong ","Jingjiang","Tongzhou","Pan'an")
# --- 循环处理每个文件夹 ---
for (i in 1:length(treatment_folders)) {
  current_input_folder <- paste0(input_base_dir, treatment_folders[i], "/")
  current_output_suffix <- output_suffixes[i] 
  current_title <- title[i]
  # 创建当前循环的输出文件夹
  Envnetplot2 <- paste0(output_base_dir, "16S_ITS_network_", current_output_suffix)
  dir.create(Envnetplot2, showWarnings = FALSE, recursive = TRUE) # recursive=TRUE 确保创建所有父目录
  # 读取和处理细菌数据
  b_file_path <- paste0(current_input_folder, "b_filtered_otus.csv")
  b_filtered_otus <- read.csv(file = b_file_path, sep = ",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  b_filtered_otus <- t(b_filtered_otus)
  colnames(b_filtered_otus) <- b_filtered_otus[1,]
  b_filtered_otus <- b_filtered_otus[-1,]
  b_filtered_otus <- as.matrix(b_filtered_otus)
  mode(b_filtered_otus) <- "numeric" 
  b_filtered_otus <- otu_table(b_filtered_otus, taxa_are_rows = TRUE)
  b_filtered_merged <- merge_phyloseq(b_filtered_otus, b_tax_species, metadata)
  # 读取和处理真菌数据 
  f_file_path <- paste0(current_input_folder, "f_filtered_otus.csv")
  f_filtered_otus <- read.csv(file = f_file_path, sep = ",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  f_filtered_otus <- t(f_filtered_otus)
  colnames(f_filtered_otus) <- f_filtered_otus[1,]
  f_filtered_otus <- f_filtered_otus[-1,]
  f_filtered_otus <- as.matrix(f_filtered_otus)
  mode(f_filtered_otus) <- "numeric" 
  f_filtered_otus <- otu_table(f_filtered_otus, taxa_are_rows = TRUE)
  f_filtered_merged <- merge_phyloseq(f_filtered_otus, f_tax_species, metadata)
  ps.merge <- ggClusterNet::merge16S_ITS(
    ps16s = b_filtered_merged, # 细菌的 phyloseq 对象
    psITS = f_filtered_merged, # 真菌的 phyloseq 对象
    N16s = 2500, # 可选参数：如果需要，可以指定16S数据保留的最高丰度OTU/ASV数量
    NITS = 2500  # 可选参数：如果需要，可以指定ITS数据保留的最高丰度OTU/ASV数量
  )
  ps.merge
  map = phyloseq::sample_data(ps.merge)
  head(map)
  map$Group1 = "one"
  phyloseq::sample_data(ps.merge) <- map
  map = phyloseq::sample_data(ps.merge)
  ps.merge %>% vegan_tax()
  result <- corBionetwork.st(
    ps.st= ps.merge,
    g1 = "Group1",
    g2 = NULL,
    g3 = NULL,
    ord.g1 = NULL,
    ord.g2 = NULL,
    ord.g3 = NULL,
    order = NULL,
    fill = "filed",
    size = "igraph.degree",method = "spearman",
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,layout_net = "model_Gephi.3",
    r.threshold=0.8,
    p.threshold=0.01,
    maxnode = 5,
    N= 4000,scale = TRUE,env = NULL,
    bio = TRUE,minsize = 4,maxsize = 14)
  summary(result[[2]])
  # 提取 corBionetwork 返回的绘图对象 (ggplot2 对象)
  p = result[[1]]
  # 显示绘图 (如果在一个交互式环境中，这会显示图表)
  p
  # 构建输出PDF文件的完整路径和文件名
  plotname1 = paste(Envnetplot2,"/network_all.pdf",sep = "")
  # 保存网络图为PDF文件
  ggsave(plotname1, p, width = 24, height = 19)
  # 使用 model_maptree 函数对网络相关性矩阵进行社区检测和布局
  tem <- model_maptree(
    cor = result$network.data$cortab$..one,          # 输入相关性矩阵
    method = "cluster_fast_greedy", # 社区检测方法
    seed = 12                     # 设置随机种子，确保结果可重复
  )
  # 提取社区检测后的节点模型信息 (通常包含每个节点的社区ID等)
  node_model = tem[[2]]
  # 查看节点模型的前几行
  head(node_model)
  # 构建输出节点模型信息CSV文件的完整路径和文件名
  tablename <- paste(Envnetplot2,"/node_model_imformation",".csv",sep = "")
  # 将节点模型信息写入CSV文件
  write.csv(node_model,tablename)
  # 构建输出节点绘图数据CSV文件的完整路径和文件名
  tablename <- paste(Envnetplot2,"/nodeG_plot",".csv",sep = "")
  # 将用于绘图的节点数据写入CSV文件 (通常包含坐标、颜色、大小等信息)
  write.csv(result$network.data$node,tablename)
  # 构建输出边绘图数据CSV文件的完整路径和文件名
  tablename <- paste(Envnetplot2,"/edge_plot",".csv",sep = "")
  # 将用于绘图的边数据写入CSV文件 (通常包含起点、终点、相关性、颜色等信息)
  write.csv(result$network.data$edge,tablename)
  # 构建输出相关性矩阵CSV文件的完整路径和文件名
  tablename <- paste(Envnetplot2,"/cor_matrix",".csv",sep = "")
  # 将计算出的相关性矩阵写入CSV文件
  write.csv(result$network.data$cortab,tablename)
}
# --- 循环：遍历每个输出文件夹，读取 edge_plot.csv 并计算边密度 ---
for (i in 1:length(treatment_folders)) {
  # 构建当前循环的输出文件夹路径 (这里是查找已生成的edge_plot.csv的路径)
  current_output_suffix <- output_suffixes[i]
  current_title <- title[i]
  Envnetplot2 <- paste0(output_base_dir, "16S_ITS_network_", current_output_suffix)
  # 构建 edge_plot.csv 文件的完整路径
  tablename_edge_plot <- paste(Envnetplot2, "/edge_plot", ".csv", sep = "")
  # 构建当前分组的边密度输出文件路径
  output_edge_density_file <- paste0(Envnetplot2, "/edge_density.csv")
  # --- 读取 edge_plot.csv 并计算边密度 ---
  if (file.exists(tablename_edge_plot)) {
    edge_data <- read.csv(file = tablename_edge_plot, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    # 确保 edge_data 不为空
    if (nrow(edge_data) > 0) {
      # 从第一行获取 links 和 nodes (因为它们在每行都一样)
      L <- edge_data$links[1] # 实际存在的边数
      n <- edge_data$nodes[1] # 节点总数
      # 计算最大可能存在的边数 (无向图)
      L_max <- (n * (n - 1)) / 2
      # 避免除以零的情况 (虽然理论上 n >= 2 才会形成边)
      if (L_max > 0) {
        current_edge_density <- L / L_max
      } else {
        current_edge_density <- 0 # 如果没有节点或只有一个节点，边密度为0
      }
      print(paste0("Group: ", current_title, ", Edge Density: ", round(current_edge_density, 5)))
      # --- 将边密度存储到当前分组的文件夹中 ---
      # 创建一个单行数据框，包含 'edge density' 列和对应的值
      edge_density_df <- data.frame(
        "edge density" = current_edge_density,
        stringsAsFactors = FALSE
      )
      colnames(edge_density_df) <- "edge density" 
      # 将结果写入 CSV 文件
      write.csv(edge_density_df, output_edge_density_file, row.names = FALSE)
      
    } else {
      print(paste0("Warning: edge_plot.csv for ", current_title, " is empty. Edge density cannot be calculated."))
      # 对于空文件的情况，也可以选择写入一个包含 NA 或 0 的文件
      edge_density_df <- data.frame(
        "edge density" = NA,
        stringsAsFactors = FALSE
      )
      write.csv(edge_density_df, output_edge_density_file, row.names = FALSE)
    }
  } else {
    print(paste0("Error: ", tablename_edge_plot, " not found. Please ensure your previous network analysis has been run and generated this file."))
    # 对于文件未找到的情况，也写入一个包含 NA 的文件
    edge_density_df <- data.frame(
      "edge density" = NA,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    write.csv(edge_density_df, output_edge_density_file, row.names = FALSE)
  }
}
# --- 循环结束后，不再统一查看所有边密度，因为已单独保存 ---
print("所有分组的边密度已计算并保存到各自的文件夹中。")
#生态型
corBionetwork.st_ecotype = function(
    ps.st= ps.merge,# phyloseq对象
    g1 = "Group",# 分组1
    g2 = NULL,# 分组2
    g3 = NULL,# 分组3
    ord.g1 = NULL, # 排序顺序
    ord.g2 = NULL, # 排序顺序
    ord.g3 = NULL, # 排序顺序
    order = NULL, # 出图每行代表的变量
    fill = "filed",
    size = "igraph.degree",
    method = "spearman",
    lab = NULL,
    label = TRUE,
    clu_method = "cluster_fast_greedy",
    select_layout = TRUE,
    layout_net = "model_maptree2",
    r.threshold=0.8,
    p.threshold=0.01,
    maxnode = 5,
    N= 500,
    scale = TRUE,
    env = NULL,
    bio = TRUE,
    minsize = 4,
    maxsize = 14,
    group.node = NULL,
    model.node = FALSE
){
  
  if (scale) {
    ps.st  = ps.st %>% ggClusterNet::scale_micro()
  }
  
  
  ps.all = ps.st
  map = sample_data(ps.all)
  
  # g2 = NULL
  if (is.null(g2)) {
    sp = ""
  } else if (is.null(ord.g2)){
    sp = map[,g2] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    sp = ord.g2
    print(sp)
  }
  # g3 = NULL
  if (is.null(g3)) {
    ti = ""
  } else if (is.null(ord.g3)){
    ti = map[,g3] %>% as.matrix() %>% as.vector() %>% unique()
  } else{
    ti = ord.g3
    print(ti)
  }
  
  if (is.null(ord.g1)) {
    group = map[,g1] %>% as.matrix() %>% as.vector() %>% unique()
    
  } else{
    group = ord.g1
    print(group)
  }
  
  
  
  #-构造两两组合全部情况
  for (i in 1:length(sp)) {
    dat = data.frame(g2 = sp[i],g3 = ti)
    
    if (i ==1) {
      dat.f = dat
    } else{
      dat.f = rbind(dat.f,dat)
    }
    
  }
  
  if (!is.null(env)) {
    colnames(env)[1] = "ID"
    env_sub <-  env[match(mapsub$ID,env$ID),]
    head(env_sub)
  }
  
  
  j = 1
  n = 1
  
  cor.all = list()
  # 存储不同分组的相关矩阵拮即可
  for (j in 1:nrow(dat.f)) {
    if (dat.f[j,1] == "") {
      ps.t = ps.all
    } else{
      ps.t = ps.all %>% subset_samples.wt(g2,dat.f[j,1])
    }
    
    if (dat.f[j,2] == "") {
      ps.f = ps.t
    } else{
      ps.f = ps.t  %>% subset_samples.wt(g3,dat.f[j,2])
    }
    
    
    for (n in 1:length(group)) {
      map = sample_data(ps.f)
      head(map)
      map$Group
      ps.g = ps.f  %>% subset_samples.wt(g1,group[n]) %>% filter_OTU_ps(N)
      big = TRUE
      if (bio) {
        if (!is.null(env)) {
          
          if (big) {
            
            result <- corBiostripeBig(data =  env_sub,
                                      group = envGroup,
                                      ps = ps.g,
                                      r.threshold = r.threshold,
                                      p.threshold = p.threshold,
                                      method = method)
          } else {
            result <- corBiostripe(data =  env_sub,
                                   group = envGroup,
                                   ps = ps.g,
                                   r.threshold = r.threshold,
                                   p.threshold = p.threshold,
                                   method = method)
          }
          
          
          #-- extract cor matrix
          occor.r = result[[1]]
          tax = as.data.frame((ggClusterNet::vegan_tax(ps.g)))
          
          if (length(tax$filed) != 0) {
            group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
          } else {
            group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
          }
          
          colnames(envGroup) <-c("SampleID","Group")
          netClu = rbind(envGroup,group2)
          colnames(netClu) <- c("ID","group")
        } else {
          
          if (big) {
            result <- corBiostripeBig(ps = ps.g,
                                      r.threshold = r.threshold, p.threshold = p.threshold, method = method)
          } else {
            result <- corBiostripe(ps = ps.g,
                                   r.threshold = r.threshold, p.threshold = p.threshold, method = method)
          }
          
          #-- extract cor matrix
          occor.r = result[[1]]
          tax = as.data.frame((ggClusterNet::vegan_tax(ps.g)))
          
          if (length(tax$filed) != 0) {
            group2 <- data.frame(SampleID = row.names(tax),Group = tax$filed)
          } else {
            group2 <- data.frame(SampleID = row.names(tax),Group = "OTU")
          }
          
          netClu = group2
          colnames(netClu) <- c("ID","group")
          
        }
        
      }
      
      # head(netClu)
      
      cor = result[[1]]
      tem = paste(dat.f[j,1],dat.f[j,2],group[n],sep = ".")
      cor.all[[tem]] = cor
      # gru.all[[tem]] = cor
      
      if (!is.null(group.node)) {
        group.node = group.node %>% dplyr::filter(ID %in% row.names(cor))
        netClu = group.node
        
      }
      
      res = node.edge(
        cor = cor,
        select_layout = TRUE,
        clu_method=clu_method,
        layout_net = layout_net,
        group.node = group.node,
        model.node = model.node
      )
      
      nod = res[[1]]
      nod$group = tem
      # head(nod)
      # ggplot(nod) + geom_point(aes(X1,X2))
      
      edg = res[[2]]
      edg$group = tem
      
      netClu$group2 = tem
      head(nod)
      head(edg)
      
      
      if (j ==1 & n == 1) {
        node = nod
        edge = edg
        netClu2 = netClu
      } else{
        node = rbind(node,nod)
        edge = rbind(edg,edge)
        netClu2 = rbind(netClu,netClu2)
      }
      
    }
    
  }
  
  #   head(edge)
  #   head(node)
  #  edge$group %>% unique()
  #  netClu2$group2%>% unique()
  # netClu2$group = as.factor(netClu2$group)
  
  #-统计边的节点数量 node link
  tem = edge$group %>% table() %>% as.data.frame()
  colnames(tem) = c("group","links")
  i = 1
  id = edge$group %>% unique()
  aa = c()
  for (i in 1:length(id)) {
    aa[i] = edge %>% filter(group == id[i]) %>%
      select("OTU_2", "OTU_1") %>% as.matrix() %>%
      as.vector() %>% unique() %>% length()
  }
  tem2 = data.frame(group = id,nodes = aa)
  
  tem3 = tem %>% full_join(tem2,by = "group")
  tem3$label= paste(tem3$group,": (nodes: ",
                    tem3$nodes,"; links: ",tem3$links,")",sep = "")
  
  
  if (is.null(order)){
    order = ""
    
  }
  
  #-行--空间
  if (order == "space"|order == "g2") {
    row.id = g3
    row.num = length(group) * length(ti)
    a = c()
    for (i in 1:length(sp)) {
      for (j in 1:length(ti)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }
    
  } else if (order == "time"|order == "g3") {
    row.id = g2
    row.num = length(group) * length(sp)
    a = c()
    for (j in 1:length(ti)) {
      for (i in 1:length(sp)) {
        tem = paste(sp[i],ti[j],group,sep = ".")
        a = c(a,tem)
      }
    }
  } else{
    a = NULL
    row.num = length(group)
  }
  
  if (!is.null(a)) {
    node$group = factor(node$group,levels = a)
    edge$group = factor(edge$group,levels = a)
    tem3 = tem3[match(a,tem3$group),]
  }else{
    node$group = factor(node$group)
    edge$group = factor(edge$group)
  }
  
  head(edge)
  head(node)
  
  # ggplot(node %>% filter(group == "..Group1")) +geom_point(aes(X1,X2))
  
  tax = ps.st %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("elements")
  node = node %>% left_join(tax,by = "elements")
  
  
  tem3$label = factor(tem3$label,levels = tem3$label)
  edge = edge %>% left_join(tem3,by = "group")
  head(edge)
  edge$label = factor(edge$label,levels = as.character(tem3$label))
  
  head(node)
  
  node = node %>% left_join(tem3,by = "group")
  node$label = factor(node$label,levels = as.character(tem3$label))
  
  
  
  net.dat = list(
    cortab = cor.all,
    node = node,
    edge = edge
  )
  
  
  # 确保 'Other' 是 Phylum 列中的最后一个因子水平
  node$Phylum <- factor(node$Phylum, levels = ordered_phylum_levels)#
  col1 = c("#b10026","#4AC0FF")
  names(col1) = c("+","-")
  p0 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.3,alpha = 0.5) +
    geom_point(aes(x = X1, y = X2,size = !!sym(size), fill = Ecotype, shape = !!sym(fill)),data =  node,alpha = 0.5) +#fill = Ecotype
    scale_colour_brewer(palette = "Set1") +
    scale_color_manual(values = col1, labels = c("Negative", "Positive")) +
    scale_shape_manual(values = c(21,25),labels = c("Bacteria", "Fungi")) + 
    scale_size(range = c(minsize, maxsize)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    scale_fill_manual(values = brewer.pal(3, "Pastel1"))+
    # labs( title = current_title) +
    #facet_wrap(.~ label,scales="free_y",ncol = row.num ) +
    theme_void()+
    guides(
      fill = guide_legend(title = "Phylum", override.aes = list(shape = 22, size = 5), order = 2), # 增加 size = 5# title = "Ecotype"
      shape = guide_legend(title = "Microbiomes", override.aes = list(size = 5), order = 1),       # 增加 size = 5
      color = guide_legend(title = "Edge", override.aes = list(size = 8), order = 3),             # 增加 size = 1 (边的粗细)
      size = "none"
    ) +
    theme(
      #plot.title = element_text(size = 60, hjust = 0.5
      legend.title = element_text(size = 24), # 调整图例标题的大小，例如 14
      legend.text = element_text(size = 22),   # 调整图例内容的大小，例如 12
      legend.spacing.y = unit(5.5, "cm")
    )
  
  p0  
  
  tem2 = node %>% dplyr::filter(elements %in% lab[[1]])
  if (label == TRUE ) {
    
    if (!is.null(lab)) {
      p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = tem2)
    } else {
      p1 <- p0 + ggrepel::geom_text_repel(aes(X1, X2,label= elements),size=4, data = node)
    }
    
  }
  
  return(list(network.plot = p0,network.data = net.dat,p1))
  
}
#核心、特殊、其他
b_tax_species_abundance <- data.frame(b_tax_species)
# 初始化所有物种为 'Others'
b_tax_species_abundance$Ecotype <- "Others"
# 获取 b_core_abundance 和 b_rare_unique_abundance 的行名（即 OTU ID）
b_core_otus <- row.names(b_core_abundance)
b_unique_otus <- row.names(b_rare_unique_abundance)
# 根据索引匹配，设置 'Cores'
# 使用 `match` 函数找到 b_tax_species_abundance 中与 core_otus 匹配的行
# 然后将这些行的Ecotype列设置为 'Cores'
b_tax_species_abundance[row.names(b_tax_species_abundance) %in% b_core_otus, "Ecotype"] <- "Cores"
# 根据索引匹配，设置 'Uniques'
b_tax_species_abundance[row.names(b_tax_species_abundance) %in% b_unique_otus, "Ecotype"] <- "Uniques"
b_tax_species_abundance <- tax_table(as.matrix(b_tax_species_abundance))
b_merged_abundance <- merge_phyloseq(b_species_rel, b_tax_species_abundance, metadata)
f_tax_species_abundance <- data.frame(f_tax_species)
# 初始化所有物种为 'Others'
f_tax_species_abundance$Ecotype <- "Others"
# 获取 f_core_abundance 和 f_rare_unique_abundance 的行名（即 OTU ID）
f_core_otus <- row.names(f_core_abundance)
f_unique_otus <- row.names(f_rare_unique_abundance)
# 根据索引匹配，设置 'Cores'
# 使用 `match` 函数找到 f_tax_species_abundance 中与 core_otus 匹配的行
# 然后将这些行的Ecotype列设置为 'Cores'
f_tax_species_abundance[row.names(f_tax_species_abundance) %in% f_core_otus, "Ecotype"] <- "Cores"
# 根据索引匹配，设置 'Uniques'
f_tax_species_abundance[row.names(f_tax_species_abundance) %in% f_unique_otus, "Ecotype"] <- "Uniques"
f_tax_species_abundance <- tax_table(as.matrix(f_tax_species_abundance))
f_merged_abundance <- merge_phyloseq(f_species_rel, f_tax_species_abundance, metadata)
#创建文件夹
Envnetplot3<- paste("D:/study/master/Main_Figure_tables/Figure_5/16S_ITS_network_abundance",sep = "")
dir.create(Envnetplot3)
# 确保细菌 (b_merged_rel) 和真菌 (f_merged_rel) 的 phyloseq 对象拥有相同的样本元数据 (map file)。
# 这一步是关键，因为 merge16S_ITS 函数需要两个 phyloseq 对象具有一致的样本信息。
ps.merge <- ggClusterNet::merge16S_ITS(
  ps16s = b_merged_abundance, # 细菌的 phyloseq 对象
  psITS = f_merged_abundance, # 真菌的 phyloseq 对象
  N16s = 2500, # 可选参数：如果需要，可以指定16S数据保留的最高丰度OTU/ASV数量
  NITS = 2500  # 可选参数：如果需要，可以指定ITS数据保留的最高丰度OTU/ASV数量
)
# 打印合并后的 phyloseq 对象概览
ps.merge
# 提取合并后的 phyloseq 对象的样本元数据
map = phyloseq::sample_data(ps.merge)
# 查看元数据的前几行
head(map)
# 添加一个新的样本分组列 "Group1"，并将其所有值设置为 "one"。
# 这通常用于创建一个单一的整体网络，而不是按特定分组进行子网络分析。
# 如果有其他分组变量，可以替换 "one" 或添加更多条件。
map$Group1 = "one"
# 将修改后的元数据重新赋值给 ps.merge 对象
phyloseq::sample_data(ps.merge) <- map
map = phyloseq::sample_data(ps.merge)
ps.merge %>% vegan_tax()
result <- corBionetwork.st_ecotype(
  ps.st= ps.merge,
  g1 = "Group1",
  g2 = NULL,
  g3 = NULL,
  ord.g1 = NULL,
  ord.g2 = NULL,
  ord.g3 = NULL,
  order = NULL,
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_Gephi.3",
  r.threshold=0.8,
  p.threshold=0.01,
  maxnode = 5,
  N= 4000,scale = TRUE,env = NULL,
  bio = TRUE,minsize = 4,maxsize = 14)
summary(result[[2]])
# 提取 corBionetwork 返回的绘图对象 (ggplot2 对象)
p = result[[1]]
# 显示绘图 (如果在一个交互式环境中，这会显示图表)
p
# 构建输出PDF文件的完整路径和文件名
plotname1 = paste(Envnetplot3,"/network_all.pdf",sep = "")
# 保存网络图为PDF文件
ggsave(plotname1, p, width = 21, height = 19)
# 使用 model_maptree 函数对网络相关性矩阵进行社区检测和布局
tem <- model_maptree(
  cor = result$network.data$cortab$..one,          # 输入相关性矩阵
  method = "cluster_fast_greedy", # 社区检测方法
  seed = 12                     # 设置随机种子，确保结果可重复
)
# 提取社区检测后的节点模型信息 (通常包含每个节点的社区ID等)
node_model = tem[[2]]
# 查看节点模型的前几行
head(node_model)
# 构建输出节点模型信息CSV文件的完整路径和文件名
tablename <- paste(Envnetplot3,"/node_model_imformation",".csv",sep = "")
# 将节点模型信息写入CSV文件
write.csv(node_model,tablename)
# 构建输出节点绘图数据CSV文件的完整路径和文件名
tablename <- paste(Envnetplot3,"/nodeG_plot",".csv",sep = "")
# 将用于绘图的节点数据写入CSV文件 (通常包含坐标、颜色、大小等信息)
write.csv(result$network.data$node,tablename)
# 构建输出边绘图数据CSV文件的完整路径和文件名
tablename <- paste(Envnetplot3,"/edge_plot",".csv",sep = "")
# 将用于绘图的边数据写入CSV文件 (通常包含起点、终点、相关性、颜色等信息)
write.csv(result$network.data$edge,tablename)
# 构建输出相关性矩阵CSV文件的完整路径和文件名
tablename <- paste(Envnetplot3,"/cor_matrix",".csv",sep = "")
# 将计算出的相关性矩阵写入CSV文件
write.csv(result$network.data$cortab,tablename)
#广适、狭适、其他
b_tax_species_breadth <- data.frame(b_tax_species)
# 初始化所有物种为 'Others'
b_tax_species_breadth$Ecotype <- "Others"
# 获取行名（即 OTU ID）
b_generalist_otus <- row.names(b_generalist_rel)
b_specialist_otus <- row.names(b_specialist_rel)
# 根据索引匹配，设置 'Cores'
# 使用 `match` 函数找到 b_tax_species_breadth 中与 b_generalist_otus 匹配的行
# 然后将这些行的Ecotype列设置为
b_tax_species_breadth[row.names(b_tax_species_breadth) %in% b_generalist_otus, "Ecotype"] <- "Generalists"
# 根据索引匹配，设置 
b_tax_species_breadth[row.names(b_tax_species_breadth) %in% b_specialist_otus, "Ecotype"] <- "Specialists"
b_tax_species_breadth <- tax_table(as.matrix(b_tax_species_breadth))
b_merged_breadth <- merge_phyloseq(b_species_rel, b_tax_species_breadth, metadata)
f_tax_species_breadth <- data.frame(f_tax_species)
# 初始化所有物种为 'Others'
f_tax_species_breadth$Ecotype <- "Others"
# 获取行名（即 OTU ID）
f_generalist_otus <- row.names(f_generalist_rel)
f_specialist_otus <- row.names(f_specialist_rel)
# 根据索引匹配，设置 
# 使用 `match` 函数找到 f_tax_species_breadth 中与 f_generalist_otus 匹配的行
# 然后将这些行的Ecotype列设置为 'Cores'
f_tax_species_breadth[row.names(f_tax_species_breadth) %in% f_generalist_otus, "Ecotype"] <- "Generalists"
# 根据索引匹配，设置
f_tax_species_breadth[row.names(f_tax_species_breadth) %in% f_specialist_otus, "Ecotype"] <- "Specialists"
f_tax_species_breadth <- tax_table(as.matrix(f_tax_species_breadth))
f_merged_breadth <- merge_phyloseq(f_species_rel, f_tax_species_breadth, metadata)
#创建文件夹
Envnetplot4<- paste("D:/study/master/Main_Figure_tables/Figure_5/16S_ITS_network_breadth",sep = "")
dir.create(Envnetplot4)
# 确保细菌 (b_merged_rel) 和真菌 (f_merged_rel) 的 phyloseq 对象拥有相同的样本元数据 (map file)。
# 这一步是关键，因为 merge16S_ITS 函数需要两个 phyloseq 对象具有一致的样本信息。
ps.merge <- ggClusterNet::merge16S_ITS(
  ps16s = b_merged_breadth, # 细菌的 phyloseq 对象
  psITS = f_merged_breadth, # 真菌的 phyloseq 对象
  N16s = 2500, # 可选参数：如果需要，可以指定16S数据保留的最高丰度OTU/ASV数量
  NITS = 2500  # 可选参数：如果需要，可以指定ITS数据保留的最高丰度OTU/ASV数量
)
# 打印合并后的 phyloseq 对象概览
ps.merge
# 提取合并后的 phyloseq 对象的样本元数据
map = phyloseq::sample_data(ps.merge)
# 查看元数据的前几行
head(map)
# 添加一个新的样本分组列 "Group1"，并将其所有值设置为 "one"。
# 这通常用于创建一个单一的整体网络，而不是按特定分组进行子网络分析。
# 如果有其他分组变量，可以替换 "one" 或添加更多条件。
map$Group1 = "one"
# 将修改后的元数据重新赋值给 ps.merge 对象
phyloseq::sample_data(ps.merge) <- map
map = phyloseq::sample_data(ps.merge)
ps.merge %>% vegan_tax()
result <- corBionetwork.st_ecotype(
  ps.st= ps.merge,
  g1 = "Group1",
  g2 = NULL,
  g3 = NULL,
  ord.g1 = NULL,
  ord.g2 = NULL,
  ord.g3 = NULL,
  order = NULL,
  fill = "filed",
  size = "igraph.degree",method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,layout_net = "model_Gephi.3",
  r.threshold=0.8,
  p.threshold=0.01,
  maxnode = 5,
  N= 4000,scale = TRUE,env = NULL,
  bio = TRUE,minsize = 4,maxsize = 14)
summary(result[[2]])
# 提取 corBionetwork 返回的绘图对象 (ggplot2 对象)
p = result[[1]]
# 显示绘图 (如果在一个交互式环境中，这会显示图表)
p
# 构建输出PDF文件的完整路径和文件名
plotname1 = paste(Envnetplot4,"/network_all.pdf",sep = "")
# 保存网络图为PDF文件
ggsave(plotname1, p, width = 21, height = 19)
# 使用 model_maptree 函数对网络相关性矩阵进行社区检测和布局
tem <- model_maptree(
  cor = result$network.data$cortab$..one,          # 输入相关性矩阵
  method = "cluster_fast_greedy", # 社区检测方法
  seed = 12                     # 设置随机种子，确保结果可重复
)
# 提取社区检测后的节点模型信息 (通常包含每个节点的社区ID等)
node_model = tem[[2]]
# 查看节点模型的前几行
head(node_model)
# 构建输出节点模型信息CSV文件的完整路径和文件名
tablename <- paste(Envnetplot4,"/node_model_imformation",".csv",sep = "")
# 将节点模型信息写入CSV文件
write.csv(node_model,tablename)
# 构建输出节点绘图数据CSV文件的完整路径和文件名
tablename <- paste(Envnetplot4,"/nodeG_plot",".csv",sep = "")
# 将用于绘图的节点数据写入CSV文件 (通常包含坐标、颜色、大小等信息)
write.csv(result$network.data$node,tablename)
# 构建输出边绘图数据CSV文件的完整路径和文件名
tablename <- paste(Envnetplot4,"/edge_plot",".csv",sep = "")
# 将用于绘图的边数据写入CSV文件 (通常包含起点、终点、相关性、颜色等信息)
write.csv(result$network.data$edge,tablename)
# 构建输出相关性矩阵CSV文件的完整路径和文件名
tablename <- paste(Envnetplot4,"/cor_matrix",".csv",sep = "")
# 将计算出的相关性矩阵写入CSV文件
write.csv(result$network.data$cortab,tablename)

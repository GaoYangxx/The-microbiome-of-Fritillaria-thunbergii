#图2a没分析的有筛选多于生物学重复的序列、
#图2a
# 安装 R 包
#https://github.com/grunwaldlab/metacoder/archive/refs/heads/master.zip
install.packages("remotes")
install.packages(c("ggplot2", "dplyr", "readr", "tibble", "vegan", "ape", "agricolae", "stringr"))
install.packages("BiocManager")
BiocManager::install("phyloseq")
library(remotes)
library(phyloseq)
library(metacoder)
library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
library(stringr)
remotes::install_github("mikemc/phyloseqCompanion") 
remotes::install_local("D:/study/micro/reference/metacoder-master.zip")
#细菌树状热图
bacteria_ASV <- read.csv(file = "D:/study/master/meiji/bacteria_ASV.csv",
                         sep=",",header=TRUE,check.names = FALSE)
# 定义文件的绝对路径（使用正斜杠）
metadata <- "D:/study/master/metadata.tsv"
# 导入元数据文件
metadata <- import_qiime_sample_data(metadata)
# 查看元数据
head(metadata)
# 读取 ASV 表（特征表）
#b_ASV <- bacteria_ASV[, c("ASV",metadata$Sample.ID)]#根际细菌
b_ASV <- bacteria_ASV[, c(9,13:60)]#根际和鳞茎细菌
colnames(b_ASV)
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
b_ASV <- otu_table(b_merged, taxa_are_rows=TRUE)
# 提取分类表
b_tax <- tax_table(b_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(b_merged))
#**转换 taxonomy 表为 metacoder 兼容格式**
b_tax_df<-as.data.frame(b_tax)
b_tax_df$ASV <- rownames(b_tax_df)
b_tax_df$lineage <- paste(      #参考parse_tax_data函数参数
  b_tax_df$Kingdom,          # 添加 Kingdom
  b_tax_df$Phylum,           # 添加 Phylum
  b_tax_df$Class,            # 添加 Class
  b_tax_df$Order,            # 添加 Order
  b_tax_df$Family,           # 添加 Family
  b_tax_df$Genus,            # 添加 Genus
  b_tax_df$Species,          # 添加 Species
  sep = ";"                 # 使用分号作为分隔符
)
b_tax_df$lineage <- str_replace_all(b_tax_df$lineage, ";uncultured_bacterium", "")
b_tax_df_clean <- b_tax_df %>%
  filter(!str_detect(Genus, "unclassified"))
b_tax_df_clean <- b_tax_df_clean %>%
  filter(!str_detect(Genus, "norank"))
#b_tax_df_clean <- b_tax_df_clean %>%
#add_count(Genus, name = "Genus_Count") %>%  # 添加计数列，列名为Genus_Count
#filter(Genus_Count >= 2) %>%                # 保留出现次数≥2的属
#select(-Genus_Count)                        # 移除计数列（可选）
b_ASV<-as.data.frame(b_ASV)
b_ASV$ASV <- rownames(b_ASV) 
b_otus <- merge(b_tax_df_clean[, c("ASV", "lineage")], b_ASV, by = "ASV", all.x = TRUE)
row_sums <- rowSums(b_otus[3:48])#去除过少ASV
rows_to_remove <- row_sums <= 3
b_otus <- b_otus[!rows_to_remove, ]
# **解析数据，使 metacoder 识别**
b_tree <- parse_tax_data(b_otus,
                         class_cols = "lineage", 
                         class_sep = ";", 
                         class_regex = "^(.+)__(.+)$", 
                         class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
# 移除无意义的分类名
b_tree <- metacoder::filter_taxa(b_tree, taxon_names != "")
b_tree <- metacoder::filter_taxa(b_tree, taxon_names != "Bacteria;;;;;")
#筛选属及以上
b_tree <- b_tree %>% 
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa=T)
# 计算丰度
b_tree$data$tax_abund <- calc_taxon_abund(b_tree, "tax_data")
# 计算出现次数
b_tree$data$tax_occ <- calc_n_samples(b_tree, "tax_abund")  
# **绘制分类热树（Heat Tree）**
set.seed(199)  # 保证图形的可复现性，可循环1到1000
b_heat_tree <- heat_tree(b_tree, 
                         #基础映射
                         node_color = n_obs,  # 颜色代表 ASV 频率
                         node_size = n_obs,   # 节点大小代表 ASV 频率
                         node_label = NA,#taxon_names,  # 节点标签为分类名称
                         edge_size= n_obs,
                         #颜色配置
                         edge_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                         node_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                         edge_color = n_samples,  # 边的颜色表示分类单元出现在多少个样本中
                         #图例控制
                         make_node_legend = FALSE,
                         make_edge_legend=FALSE,
                         edge_legend_title = NULL,  # 标题格式
                         edge_color_axis_label = NULL, # 隐藏左侧颜色条标签
                         edge_color_digits = 2,
                         #布局优化
                         initial_layout = "re", 
                         layout = "da")  # 选择合适的树形布局
b_heat_tree
# 导出图片
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2a_bacteria_heat_tree.png", plot = b_heat_tree, width = 10, height = 10, dpi = 300)
#真菌树状热图
fungi_ASV <- read.csv(file = "D:/study/master/meiji/fungi_ASV.csv",
                      sep=",",header=TRUE,check.names = FALSE)
# 定义文件的绝对路径（使用正斜杠）
metadata <- "D:/study/master/metadata.tsv"
# 导入元数据文件
metadata <- import_qiime_sample_data(metadata)
# 查看元数据
head(metadata)
# 读取 ASV 表（特征表）
#f_ASV <- fungi_ASV[, c("ASV",metadata$Sample.ID)]#根际真菌
f_ASV <- fungi_ASV[, c(9,13:60)]#根际和鳞茎真菌
colnames(f_ASV)
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
f_ASV <- otu_table(f_merged, taxa_are_rows=TRUE)
# 提取分类表
f_tax <- tax_table(f_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(f_merged))
#**转换 taxonomy 表为 metacoder 兼容格式**
f_tax_df<-as.data.frame(f_tax)
f_tax_df$ASV <- rownames(f_tax_df)
f_tax_df$lineage <- paste(      #参考parse_tax_data函数参数
  f_tax_df$Kingdom,          # 添加 Kingdom
  f_tax_df$Phylum,           # 添加 Phylum
  f_tax_df$Class,            # 添加 Class
  f_tax_df$Order,            # 添加 Order
  f_tax_df$Family,           # 添加 Family
  f_tax_df$Genus,            # 添加 Genus
  f_tax_df$Species,          # 添加 Species
  sep = ";"                 # 使用分号作为分隔符
)
f_tax_df$lineage <- str_replace_all(f_tax_df$lineage, ";uncultured_bacterium", "")
f_tax_df_clean <- f_tax_df %>%
  filter(!str_detect(Genus, "unclassified"))
f_tax_df_clean <- f_tax_df_clean %>%
  filter(!str_detect(Genus, "norank"))#clean有变化，ftree没变。
#f_tax_df_clean <- f_tax_df_clean %>%
#add_count(Genus, name = "Genus_Count") %>%  # 添加计数列，列名为Genus_Count
#filter(Genus_Count >= 2) %>%                # 保留出现次数≥2的属
#select(-Genus_Count)                        # 移除计数列（可选）
f_ASV<-as.data.frame(f_ASV)
f_ASV$ASV <- rownames(f_ASV) 
f_otus <- merge(f_tax_df_clean[, c("ASV", "lineage")], f_ASV, by = "ASV", all.x = TRUE)
#row_sums <- rowSums(f_otus[3:48])#去除过少ASV
#rows_to_remove <- row_sums <= 3
#f_otus <- f_otus[!rows_to_remove, ]
# **解析数据，使 metacoder 识别**
f_tree <- parse_tax_data(f_otus,
                         class_cols = "lineage", 
                         class_sep = ";", 
                         class_regex = "^(.+)__(.+)$", 
                         class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
# 移除无意义的分类名
f_tree <- metacoder::filter_taxa(f_tree, taxon_names != "")
f_tree <- metacoder::filter_taxa(f_tree, taxon_names != "Fungi;;;;;")
#筛选属及以上
f_tree <- f_tree %>% 
  metacoder::filter_taxa(taxon_ranks == "g", supertaxa=T)
# 计算丰度
f_tree$data$tax_abund <- calc_taxon_abund(f_tree, "tax_data")
# 计算出现次数
f_tree$data$tax_occ <- calc_n_samples(f_tree, "tax_abund")  
# **绘制分类热树（Heat Tree）**
set.seed(247)  # 保证图形的可复现性，可循环1到1000
f_heat_tree <- heat_tree(f_tree, 
                         #基础映射
                         node_color = n_obs,  # 颜色代表 ASV 频率
                         node_size = n_obs,   # 节点大小代表 ASV 频率
                         node_label = NA,#taxon_names,  # 节点标签为分类名称
                         edge_size= n_obs,
                         #颜色配置
                         edge_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                         node_color_range = c("#CC8394FF", "#AC563BFF", "#CDA97CFF","#7C8EC5FF","#2B3C51FF"),
                         edge_color = n_samples,  # 边的颜色表示分类单元出现在多少个样本中
                         #图例控制
                         make_node_legend = FALSE,
                         make_edge_legend=FALSE,
                         edge_legend_title = NULL,  # 标题格式
                         edge_color_axis_label = NULL, # 隐藏左侧颜色条标签
                         edge_color_digits = 2,
                         #布局优化
                         initial_layout = "re", 
                         layout = "da")  # 选择合适的树形布局
f_heat_tree
# 导出图片
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2a_fungi_heat_tree.png", plot = f_heat_tree, width = 10, height = 10, dpi = 300)

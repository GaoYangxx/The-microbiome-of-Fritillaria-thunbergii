#图5d
# 定义需要安装的包列表
packages_to_install <- c("picante", "progress", "ggplot2", "ggridges", "forcats", "tibble", "qiime2R", "readxl", "dplyr", "tidyr", "phyloseq", "stringr", "ggtext")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste("Package", pkg, "is already installed."))
  }
}
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("phyloseq")
library(picante)
library(progress)
library(ggplot2)
library(ggridges)
library(forcats)
library(tibble)
library(qiime2R)
library(readxl)
library(dplyr) 
library(tidyr)
library(phyloseq)
library(stringr)
library(ggtext)
#细菌系统发育距离
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
b_merged <- phyloseq::subset_taxa(b_merged, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
b_merged <- phyloseq::subset_taxa(b_merged, (Kingdom!="d__Archaea") | is.na(Kingdom))
b_merged <- phyloseq::subset_taxa(b_merged, (Order!="o__Chloroplast") )
b_merged <- phyloseq::subset_taxa(b_merged, (Family!="f__Mitochondria"))
b_merged <- phyloseq::subset_taxa(b_merged, (Family!="NA"))
# 提取 ASV 表
b_ASV <- as.matrix(otu_table(b_merged, taxa_are_rows=T))
b_ASV_df <- b_ASV[rowSums(b_ASV[])>0,]
b_ASVs <- data.frame(b_ASV)
# 读取系统发育树
b_tree <- phy_tree(b_merged)
# 读取分类信息表
b_tax <- data.frame(b_tax)
# 去除未知属
b_tax <- b_tax %>%
  filter(Class != "c__bacteriap25") %>%
  filter(!(Genus == "g__uncultured" | stringr::str_detect(Genus, "g__unclassified") | stringr::str_detect(Genus, "g__norank")))
# 查看Genus列的因子水平
levels(factor(b_tax$Genus))
# 统计每个Genus的频数
b_goi <- data.frame(table(b_tax$Genus))
# 按频数降序排序
b_goi <- b_goi[order(b_goi$Freq, decreasing=TRUE), ]
# 选择频数排在第3到32位的Genus
b_goi2 <- b_goi[3:32,]
# 排除某些特定的Genus
b_goi2 <- b_goi2[!b_goi2$Var1 %in% c(" g__OM190"," g__vadinHA49"," g__A21b"," g__0319-6G20", " g__TRA3-20"," g__AKYH767", " g__Anaeromyxobacter"," g__Blfdi19"," g__Candidatus_Udaeobacter"),]
# 选择频数排在第3到102位的Genus
b_goi3 <- b_goi[3:102,]
# 初始化进度条
b_pb <- progress_bar$new(total = nrow(b_goi3))
# 初始化结果存储列表
b_res1 <- c()
# 遍历b_goi3中的每个Genus
for(b_g in levels(factor(b_goi3$Var1))){
  # 筛选当前Genus的分类信息
  b_tax.g <- b_tax[b_tax$Genus== b_g,]
  # 筛选当前Genus的ASV丰度
  b_com.g <- b_ASVs[rownames(b_ASVs) %in% rownames(b_tax.g),]
  # 筛选在超过2个样本中存在的ASV
  b_com.g <- b_com.g[,colSums(b_com.g>0)>2]
  # 基于当前Genus的ASV修剪系统发育树
  b_tree.g <- prune.sample(t(b_com.g),b_tree)
  # 计算修剪后树的平均共表型距离
  b_coph.g <- mean(cophenetic(b_tree.g))
  # 存储结果：平均共表型距离和Genus名称
  b_res1 <- rbind(b_res1,c(b_coph.g, b_g))
  b_pb$tick()
}
# 初始化结果存储列表
b_res <- c()
# 遍历b_goi2中的每个Genus
for(b_g in levels(factor(b_goi2$Var1))){
  # 筛选当前Genus的分类信息
  b_tax.g <- b_tax[b_tax$Genus== b_g,]
  # 筛选当前Genus的ASV丰度
  b_com.g <- b_ASVs[rownames(b_ASVs) %in% rownames(b_tax.g),]
  # 筛选在超过2个样本中存在的ASV
  b_com.g <- b_com.g[,colSums(b_com.g>0)>2]
  # 基于当前Genus的ASV修剪系统发育树
  b_tree.g <- prune.sample(t(b_com.g),b_tree)
  # 遍历当前Genus的每个样本
  for(k in 1:ncol(b_com.g)){
    # 提取当前样本的丰度数据
    b_com.k <- b_com.g[,k, drop=FALSE]
    # 筛选当前样本中存在的ASV
    b_com.k <- b_com.k[b_com.k>0,,drop=FALSE]
    # 基于当前样本存在的ASV修剪系统发育树
    b_tree.k <- prune.sample(t(b_com.k),b_tree.g)
    # 计算修剪后树的平均共表型距离
    b_coph.k <- mean(cophenetic(b_tree.k))
    # 存储结果：Genus名称、平均共表型距离和样本名称
    b_res <- rbind(b_res, c(b_g, b_coph.k, colnames(b_com.k)))}
}
# 将结果转换为数据框
b_res <- data.frame(b_res)
# 设置列名
colnames(b_res) <- c("genus","phy.dist","sample")
b_res$phy.dist <- as.numeric(as.character(b_res$phy.dist))  
# 计算每个属的密度峰值 (使用新变量名)
# 这个自定义函数用于找到连续分布的众数（密度曲线的峰值）
get_density_mode <- function(x_values) {
  # 处理数据点过少无法计算密度的情况
  if (length(x_values) < 2) {
    return(NA)
  }
  # 计算密度
  density_obj <- density(x_values, na.rm = TRUE)
  # 找到密度最高的 x 值
  mode_x <- density_obj$x[which.max(density_obj$y)]
  return(mode_x)
}
# (可选) 计算并查看每个属的众数，以便核对
b_genus_modes <- b_res %>% 
  group_by(genus) %>%
  summarise(
    mode_phy_dist = get_density_mode(phy.dist)
  ) %>%
  ungroup() %>%
  arrange(mode_phy_dist) %>%
  mutate(order = row_number())
# 将 b_genus_modes 中的排序信息合并到 b_res
b_res_ordered <- b_res %>%
  left_join(b_genus_modes, by = "genus") %>%
  filter(mode_phy_dist < 0.5) %>%
  # 然后根据 mode_phy_dist 排序
  arrange(mode_phy_dist) %>%
  mutate(genus = str_replace(genus, "g__", "")) %>%# 移除 "g__" 前缀
  mutate(genus = ifelse(
    genus == "Burkholderia-Caballeronia-Paraburkholderia", # 如果 genus 是这个长名字
    "Burkholderia group", # 替换为 "Burkholderia group"
    genus # 否则保持当前 genus 的值
  ))
# 统计每个属的 ASV 数量
b_genus_asv_counts <- b_tax %>%
  group_by(Genus) %>%
  summarise(ASV_Count = n()) %>%
  ungroup()
# 重命名 Genus 列以方便后续合并，并确保名称与 b_res 中的 Genus 格式一致
# b_res 中的 genus 列在之前已经移除了 "g__" 并替换了长名称
# 所以这里我们需要在合并前，也对 genus_asv_counts 中的 Genus 进行相同的处理
b_genus_asv_counts <- b_genus_asv_counts %>%
  mutate(Genus = as.character(Genus)) %>%
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  mutate(Genus = ifelse(
    Genus == "Burkholderia-Caballeronia-Paraburkholderia",
    "Burkholderia group",
    Genus
  )) %>%
  rename(genus = Genus) # 重命名以便与 b_res_ordered 的列名匹配
#将 ASV 数量合并到 b_res_ordered 中
b_res_ordered <- b_res_ordered %>%
  left_join(b_genus_asv_counts, by = "genus")
# 创建新的 Y 轴标签列 (例如 'display_genus') 
# 格式化为 "属名 (ASV数量)"
b_res_ordered <- b_res_ordered %>%
  mutate(display_genus = paste0("*", genus, "*", " (", ASV_Count, ")"))
b_res_ordered$phy.dist <- as.numeric(as.character(b_res_ordered$phy.dist))
b_res_ordered$order <- as.numeric(as.character(b_res_ordered$order))
b_ordered_genera <- unique(b_res_ordered$genus)
# 定义一个包含12种不重复颜色的向量
specific_12_colors <- c(
  "#D99A95", "#FAC0BB", 
  "#F5B48B", "#EDC470", 
  "#E0D36B", "#D1D67B", 
  "#B7E4C7", "#87D3BC", 
  "#87CEEB", "#ADD8E6",
  "#957DAD", "#E0BBE4" 
)
# 创建一个命名向量，将颜色与属名精确关联起来
names(specific_12_colors) <- b_ordered_genera
# 使用 ggplot2 绘制山脊图，并应用这个精确的命名颜色调色板
phyl_distance_bacteria <- ggplot(b_res_ordered, aes(x = phy.dist,
                                                    y = fct_reorder(display_genus, order), # 保持按 order 排序
                                                    fill = genus)) +
  # 添加山脊图层
  geom_density_ridges() +
  # 使用自定义颜色调色板，精确映射
  scale_fill_manual(values = specific_12_colors) +
  # 使用山脊主题
  theme_ridges() +
  # 移除图例（因为你可能希望通过Y轴标签直接识别）
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0, size = 16), #调整横坐标标题位置和大小
        #修改 Y 轴文本外观和对齐
        axis.text.y = element_markdown(vjust = -0.5, # 设置垂直对齐方式为居中
                                       hjust = 0, # 设置为左对齐
                                       margin = margin(r = unit(-50, "pt")),# 设置右侧边距为0点，使其紧贴Y轴
                                       size = 14), 
        # 修改 X 轴文本外观和对齐
        axis.text.x = element_text(margin = margin(t = unit(-20, "pt")),#设置顶部边距为0点，使其紧贴X轴
                                   size = 14,
                                   vjust = 0 ) # 设置为底端对齐
  ) + 
  # 添加轴标签和标题
  labs(x = "Bacterial phylogenetic distance",
       y = "") +
  ggtitle("") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     labels = c("0.0", "0.2", "0.4", "0.6", "0.8"))
phyl_distance_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_5/5d_phyl_distance_bacteria.png", plot = phyl_distance_bacteria, width = 6, height = 6, dpi = 600, bg = "transparent")
# 对 b_res_ordered 数据按属进行分组
# 计算每个属的 phy.dist 的中位数和 IQR
# 格式化结果以匹配参考文献的样式
b_genus_phyl_summary <- b_res_ordered %>%
  group_by(genus) %>% # 按处理后的属名分组
  summarise(
    Median_Phy_Dist = median(phy.dist, na.rm = TRUE), # 计算中位数
    Q1_Phy_Dist = quantile(phy.dist, 0.25, na.rm = TRUE), # 计算第一四分位数
    Q3_Phy_Dist = quantile(phy.dist, 0.75, na.rm = TRUE)  # 计算第三四分位数
  ) %>%
  ungroup() %>%
  # 格式化输出为 "属名 (中位数; IQR: Q1–Q3)" 的形式
  # 使用 sprintf 来控制小数位数，这里保留三位小数
  mutate(
    Formatted_Summary = sprintf("%s (%.3f; IQR: %.3f–%.3f)",
                                genus,
                                Median_Phy_Dist,
                                Q1_Phy_Dist,
                                Q3_Phy_Dist)
  ) %>%
  # 按照中位数排序（可选，但通常有助于查看）
  arrange(Median_Phy_Dist)
# 打印结果
b_genus_phyl_summary
#真菌系统发育距离
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
f_merged <- phyloseq::subset_taxa(f_merged, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
f_merged <- phyloseq::subset_taxa(f_merged, (Kingdom!="d__Archaea") | is.na(Kingdom))
f_merged <- phyloseq::subset_taxa(f_merged, (Order!="o__Chloroplast") )
f_merged <- phyloseq::subset_taxa(f_merged, (Family!="f__Mitochondria"))
f_merged <- phyloseq::subset_taxa(f_merged, (Family!="NA"))
# 提取 ASV 表
f_ASV <- as.matrix(otu_table(f_merged, taxa_are_rows=T))
f_ASV_df <- f_ASV[rowSums(f_ASV[])>0,]
f_ASVs <- data.frame(f_ASV)
# 读取系统发育树
f_tree <- phy_tree(f_merged)
# 读取分类信息表
f_tax <- data.frame(f_tax)
# 去除未知属
f_tax <- f_tax %>%
  filter(!(Genus == "g__uncultured" |
             stringr::str_detect(Genus, "g__unclassified") |
             stringr::str_detect(Genus, "g__norank")))
# 查看Genus列的因子水平
levels(factor(f_tax$Genus))
# 统计每个Genus的频数
f_goi <- data.frame(table(f_tax$Genus))
# 按频数降序排序
f_goi <- f_goi[order(f_goi$Freq, decreasing=TRUE), ]
# 选择频数排在第3到32位的Genus
f_goi2 <- f_goi[1:32,]
# 选择频数排在第3到102位的Genus
f_goi3 <- f_goi[1:102,]
# 计算每个 Genus 的平均共表型距离
# 初始化进度条
f_pb1 <- progress_bar$new(total = nrow(f_goi3),
                          format = " [:bar] :percent :eta :current/:total (Genus level)")
# 初始化结果存储列表
f_res1 <- data.frame(
  Avg_Cophenetic_Dist = numeric(),
  Genus = character(),
  stringsAsFactors = FALSE
)
# 遍历f_goi3中的每个Genus
for(f_g in levels(factor(f_goi3$Var1))){
  # 筛选当前Genus的分类信息
  f_tax.g <- f_tax[f_tax$Genus == f_g, , drop = FALSE]
  # 筛选当前Genus的ASV丰度
  f_com.g <- f_ASVs[rownames(f_ASVs) %in% rownames(f_tax.g), , drop = FALSE]
  # 筛选在超过2个样本中存在的ASV
  # 检查行数以避免 rowSums(0行) 错误
  if (nrow(f_com.g) > 0) {
    f_com.g <- f_com.g[rowSums(f_com.g > 0) > 2, , drop = FALSE]
  } else {
    f_pb1$tick() # 推进进度条
    next # 没有ASV，跳过当前属
  }
  # 检查是否有足够的ASV来计算树和共表型距离 (至少2个ASV)
  if (nrow(f_com.g) < 2) {
    f_pb1$tick() # 推进进度条
    next # ASV太少，跳过当前属
  }
  # 基于当前Genus的ASV修剪系统发育树
  f_tree.g <- prune.sample(t(f_com.g),f_tree)
  # 检查修剪后的树是否有效 (非 NULL 且至少有2个叶子节点)
  if (is.null(f_tree.g) || Ntip(f_tree.g) < 2) {
    f_pb1$tick() # 推进进度条
    next # 树无效，跳过当前属
  }
  # 计算修剪后树的平均共表型距离
  f_coph.g <- mean(cophenetic(f_tree.g), na.rm = TRUE)
  # 存储结果
  new_row_f_res1 <- data.frame(
    Avg_Cophenetic_Dist = f_coph.g,
    Genus = f_g,
    stringsAsFactors = FALSE
  )
  f_res1 <- rbind(f_res1, new_row_f_res1)
  f_pb1$tick() # 推进进度条
}
head(f_res1)
# 计算每个样本在每个 Genus 内的平均共表型距离
# 初始化进度条
f_pb2 <- progress_bar$new(total = nrow(f_goi2),
                          format = " [:bar] :percent :eta :current/:total (Sample level)")
# 初始化结果存储列表
f_res <- data.frame(
  genus = character(),
  phy.dist = numeric(),
  sample = character(),
  stringsAsFactors = FALSE
)
# 遍历f_goi2中的每个Genus
for(f_g in levels(factor(f_goi2$Var1))){
  # 筛选当前Genus的分类信息
  f_tax.g <- f_tax[f_tax$Genus == f_g, , drop = FALSE]
  # 筛选当前Genus的ASV丰度
  f_com.g <- f_ASVs[rownames(f_ASVs) %in% rownames(f_tax.g), , drop = FALSE]
  # 修正: 筛选在超过2个样本中存在的ASV (如果ASV是行，样本是列)
  if (nrow(f_com.g) > 0) {
    f_com.g <- f_com.g[rowSums(f_com.g > 0) > 2, , drop = FALSE]
  } else {
    f_pb2$tick() # 推进进度条
    next # 没有ASV，跳过当前属
  }
  # 检查是否有足够的ASV来构建Genus-level的树 (至少2个ASV)
  if (nrow(f_com.g) < 2) {
    f_pb2$tick() # 推进进度条
    next # ASV太少，跳过当前属
  }
  # 基于当前Genus的ASV修剪系统发育树
  f_tree.g <- prune.sample(t(f_com.g), f_tree)
  
  # 检查修剪后的Genus树是否有效
  if (is.null(f_tree.g) || Ntip(f_tree.g) < 2) {
    f_pb2$tick() # 推进进度条
    next # Genus树无效，跳过当前属
  }
  # 遍历当前Genus的每个样本
  for(k in 1:ncol(f_com.g)){
    # 提取当前样本的丰度数据
    f_com.k <- f_com.g[, k, drop = FALSE]
    # 筛选当前样本中存在的ASV (丰度大于0)
    f_com.k <- f_com.k[f_com.k > 0, , drop = FALSE]
    # 检查当前样本是否有足够的ASV来计算距离 (至少2个ASV)
    if (nrow(f_com.k) < 2) {
      next # 样本中ASV太少，跳过当前样本
    }
    # 基于当前样本存在的ASV修剪系统发育树 (使用 f_tree.g 作为基础树)
    f_tree.k <- prune.sample(t(f_com.k), f_tree.g)
    # 检查修剪后的样本树是否有效
    if (is.null(f_tree.k) || Ntip(f_tree.k) < 2) {
      next # 样本树无效，跳过当前样本
    }
    # 计算修剪后树的平均共表型距离
    f_coph.k <- mean(cophenetic(f_tree.k), na.rm = TRUE)
    # 存储结果
    new_row_f_res <- data.frame(
      genus = f_g,
      phy.dist = f_coph.k,
      sample = colnames(f_com.k), # 确保这里只取到单个样本名
      stringsAsFactors = FALSE
    )
    f_res <- rbind(f_res, new_row_f_res)
  }
  f_pb2$tick() # 推进进度条
}
head(f_res)
summary(f_res$phy.dist)
# 将结果转换为数据框
f_res <- data.frame(f_res)
# 设置列名
colnames(f_res) <- c("genus","phy.dist","sample")
f_res$phy.dist <- as.numeric(as.character(f_res$phy.dist))  
# 计算每个属的密度峰值 (使用新变量名)
# 这个自定义函数用于找到连续分布的众数（密度曲线的峰值）
get_density_mode <- function(x_values) {
  # 处理数据点过少无法计算密度的情况
  if (length(x_values) < 2) {
    return(NA)
  }
  # 计算密度
  density_obj <- density(x_values, na.rm = TRUE)
  # 找到密度最高的 x 值
  mode_x <- density_obj$x[which.max(density_obj$y)]
  return(mode_x)
}
# (可选) 计算并查看每个属的众数，以便核对
f_genus_modes <- f_res %>% 
  group_by(genus) %>%
  summarise(
    mode_phy_dist = get_density_mode(phy.dist)
  ) %>%
  ungroup() %>%
  arrange(mode_phy_dist) %>%
  mutate(order = row_number())
# 将 f_genus_modes 中的排序信息合并到 f_res
f_res_ordered <- f_res %>%
  left_join(f_genus_modes, by = "genus") %>%
  filter(mode_phy_dist < 0.5) %>%
  # 然后根据 mode_phy_dist 排序
  arrange(mode_phy_dist) %>%
  mutate(genus = str_replace(genus, "g__", "")) %>%# 移除 "g__" 前缀
  mutate(genus = ifelse(
    genus == "Burkholderia-Caballeronia-Paraburkholderia", # 如果 genus 是这个长名字
    "Burkholderia group", # 替换为 "Burkholderia group"
    genus # 否则保持当前 genus 的值
  ))
# 统计每个属的 ASV 数量
f_genus_asv_counts <- f_tax %>%
  group_by(Genus) %>%
  summarise(ASV_Count = n()) %>%
  ungroup()
# 重命名 Genus 列以方便后续合并，并确保名称与 f_res 中的 Genus 格式一致
# f_res 中的 genus 列在之前已经移除了 "g__" 并替换了长名称
# 所以这里我们需要在合并前，也对 genus_asv_counts 中的 Genus 进行相同的处理
f_genus_asv_counts <- f_genus_asv_counts %>%
  mutate(Genus = as.character(Genus)) %>%
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  mutate(Genus = ifelse(
    Genus == "Burkholderia-Caballeronia-Paraburkholderia",
    "Burkholderia group",
    Genus
  )) %>%
  rename(genus = Genus) # 重命名以便与 f_res_ordered 的列名匹配
#将 ASV 数量合并到 f_res_ordered 中
f_res_ordered <- f_res_ordered %>%
  left_join(f_genus_asv_counts, by = "genus")
# 创建新的 Y 轴标签列 (例如 'display_genus') 
# 格式化为 "属名 (ASV数量)"
f_res_ordered <- f_res_ordered %>%
  mutate(display_genus = paste0("*", genus, "*", " (", ASV_Count, ")"))
f_res_ordered$phy.dist <- as.numeric(as.character(f_res_ordered$phy.dist))
f_res_ordered$order <- as.numeric(as.character(f_res_ordered$order))
f_ordered_genera <- unique(f_res_ordered$genus)
# 定义一个包含25种不重复颜色的向量
specific_25_colors <- c(
  "#D99A95", "#FAC0BB", 
  "#F5B48B", "#EDC470", 
  "#E0D36B", "#D1D67B", 
  "#B7E4C7", "#87D3BC", 
  "#87CEEB", "#ADD8E6",
  "#957DAD", "#E0BBE4" 
)
specific_25_colors <- colorRampPalette(specific_25_colors)
specific_25_colors <- specific_25_colors(length(f_ordered_genera))
# 创建一个命名向量，将颜色与属名精确关联起来
names(specific_25_colors) <- f_ordered_genera
# 使用 ggplot2 绘制山脊图，并应用这个精确的命名颜色调色板
phyl_distance_fungi <- ggplot(f_res_ordered, aes(x = phy.dist,
                                                 y = fct_reorder(display_genus, order), # 保持按 order 排序
                                                 fill = genus)) +
  # 添加山脊图层
  geom_density_ridges() +
  # 使用自定义颜色调色板，精确映射
  scale_fill_manual(values = specific_25_colors) +
  # 使用山脊主题
  theme_ridges() +
  # 移除图例（因为你可能希望通过Y轴标签直接识别）
  theme(legend.position = "none",
        axis.title.x = element_text(hjust = 0, size = 16), #调整横坐标标题位置和大小
        #修改 Y 轴文本外观和对齐
        axis.text.y = element_markdown(vjust = -0.1, # 设置垂直对齐方式为居中
                                       hjust = 0, # 设置为左对齐
                                       margin = margin(r = unit(-30, "pt")),# 设置右侧边距为0点，使其紧贴Y轴
                                       size = 14), 
        # 修改 X 轴文本外观和对齐
        axis.text.x = element_text(margin = margin(t = unit(-20, "pt")),#设置顶部边距为0点，使其紧贴X轴
                                   size = 14,
                                   vjust = -2 ) # 设置为底端对齐
  ) + 
  # 添加轴标签和标题
  labs(x = "Fungal phylogenetic distance",
       y = "") +
  ggtitle("") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8),
                     labels = c("0.0", "0.2", "0.4", "0.6", "0.8"))
phyl_distance_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_5/5d_phyl_distance_fungi.png", plot = phyl_distance_fungi, width = 6, height = 6, dpi = 600, bg = "transparent")
# 对 f_res_ordered 数据按属进行分组
# 计算每个属的 phy.dist 的中位数和 IQR
# 格式化结果以匹配参考文献的样式
f_genus_phyl_summary <- f_res_ordered %>%
  group_by(genus) %>% # 按处理后的属名分组
  summarise(
    Median_Phy_Dist = median(phy.dist, na.rm = TRUE), # 计算中位数
    Q1_Phy_Dist = quantile(phy.dist, 0.25, na.rm = TRUE), # 计算第一四分位数
    Q3_Phy_Dist = quantile(phy.dist, 0.75, na.rm = TRUE)  # 计算第三四分位数
  ) %>%
  ungroup() %>%
  # 格式化输出为 "属名 (中位数; IQR: Q1–Q3)" 的形式
  # 使用 sprintf 来控制小数位数，这里保留三位小数
  mutate(
    Formatted_Summary = sprintf("%s (%.3f; IQR: %.3f–%.3f)",
                                genus,
                                Median_Phy_Dist,
                                Q1_Phy_Dist,
                                Q3_Phy_Dist)
  ) %>%
  # 按照中位数排序（可选，但通常有助于查看）
  arrange(Median_Phy_Dist)
# 打印结果
f_genus_phyl_summary

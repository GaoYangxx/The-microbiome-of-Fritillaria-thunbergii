#图1a堆叠柱状图
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
cran_packages <- c(
  "dplyr", "tidyverse", "data.table", "hillR", "viridis",
  "hrbrthemes", "paletteer", "mgcv", "readxl", "picante", "corrplot",
  "ggplot2", "reshape2", "ggprism", "ggalluvial", "plyr", "ggpubr",
  "scales"
)
for (p in cran_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}
bioc_packages <- c("phyloseq", "breakaway", "microbiome")
BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("kstagaman/phyloseqCompanion", upgrade = "never")
remotes::install_github("jbisanz/qiime2R", upgrade = "never")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(data.table)
  library(readxl)
  library(scales)
  library(phyloseq)
  library(microViz)
  library(breakaway)
  library(microbiome)
  library(picante)
  library(corrplot)
  library(ggplot2)
  library(reshape2)
  library(ggprism)
  library(ggalluvial)
  library(ggpubr)
  library(viridis)
  library(hrbrthemes)
  library(paletteer)
  library(mgcv)
  library(hillR)
  library(phyloseqCompanion)
  library(qiime2R)
})
message("所有必要的R包已安装并加载。")
#细菌组成
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
b_tax <- data.frame(b_tax)
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
# 将数据转换为相对丰度
b_merged_phylum_ra <- transform_sample_counts(b_merged, function(x) x / sum(x))
#聚合到门 (Phylum) 水平
b_merged_phylum <- tax_glom(b_merged_phylum_ra, taxrank = "Phylum", NArm = FALSE)
# 提取并处理数据框 
b_phylum_ra_df <- psmelt(b_merged_phylum)
# 确保 Phylum 列是字符类型
b_phylum_ra_df$Phylum <- as.character(b_phylum_ra_df$Phylum)
# 定义一个函数，用于判断 Phylum 是否为“特殊”类型（未培养、未分类、无等级或NA）
is_special_phylum <- function(phylum_name) {
  is.na(phylum_name) | # 检查是否为NA
    stringr::str_detect(phylum_name, "uncultured") |
    stringr::str_detect(phylum_name, "unclassified") |
    stringr::str_detect(phylum_name, "norank")
}
# 如果 Phylum 属于“特殊”类型，则将其标记为“Other_Special_Phylum”。
b_phylum_pre_grouped_df <- b_phylum_ra_df %>%
  mutate(Phylum_Pre_Grouped = ifelse(is_special_phylum(Phylum), "Other_Special_Phylum", Phylum))
# 计算当前 Phylum 级别（已包含 "Other_Special_Phylum"）的平均相对丰度
b_phylum_mean_ra <- b_phylum_pre_grouped_df %>%
  group_by(Phylum_Pre_Grouped) %>% # 注意这里是对 Phylum_Pre_Grouped 列分组
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选出最终的“Top Phyla”列表
b_top_phyla_final <- b_phylum_mean_ra %>%
  filter(MeanAbundance > 0.001) %>%
  pull(Phylum_Pre_Grouped)
#最终的 Phylum 分组：将所有不符合条件的门归入“Other”
b_phylum_final_grouped_df <- b_phylum_pre_grouped_df %>%
  mutate(Phylum_Grouped = ifelse(
    Phylum_Pre_Grouped %in% b_top_phyla_final & Phylum_Pre_Grouped != "Other_Special_Phylum",
    Phylum_Pre_Grouped, # 如果是 Top Phyla 且不是特殊的“Other”标记，则保留其名称
    "Other"             # 否则，归类为最终的“Other”
  ))
# 新计算每个样本中（新）分组门的相对丰度
b_phylum_summarized_df <- b_phylum_final_grouped_df %>%
  group_by(Sample, Group, Phylum_Grouped) %>% # 确保 'Group' 列存在于 psmelt 结果中
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
b_phylum_summarized_df <- b_phylum_summarized_df %>%
  filter(!is.na(Group))
# 查看最终分组后的门类及其平均丰度
b_phylum_summarized_df %>%
  group_by(Phylum_Grouped) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(MeanAbundance)) %>%
  print()
# 计算每个门在所有样本中的平均相对丰度
b_phylum_mean_ra <- b_phylum_ra_df %>%
  group_by(Phylum) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选平均相对丰度大于 0的门 
b_top_phyla <- b_phylum_mean_ra %>%
  filter(MeanAbundance > 0) %>%
  pull(Phylum)
# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
b_phylum_summarized_df <- b_phylum_summarized_df %>%
  filter(!is.na(Group))
# 准备绘图数据 (计算每个组的平均相对丰度)
b_mean_phylum_abundance <- b_phylum_summarized_df %>%
  group_by(Group, Phylum_Grouped) %>% # 按组和处理后的门进行聚合
  summarise(MeanValue = mean(Abundance), .groups = 'drop') %>% # 计算平均丰度
  # 重命名列以匹配旧代码中的 'variable' 和 'value'
  rename(variable = Group, value = MeanValue, tax = Phylum_Grouped) %>%
  as.data.frame()
# 移除 tax 列中的 'p__' 前缀
b_mean_phylum_abundance$tax <- gsub("^p__", "", b_mean_phylum_abundance$tax)
# 首先获取所有唯一的门名称（包括 'Other'）
b_current_phylum <- unique(b_mean_phylum_abundance$tax)
# 将 'Other' 分离出来
b_other_category_phylum <- b_current_phylum[b_current_phylum == "Other" | b_current_phylum == "others"]
b_other_category_phylum <- b_other_category_phylum[1] # 确保只取一个 'Other' 或 'others'
# 排除 'Other'/'others' 后，对剩余的门进行字母排序
b_non_other_phylum <- b_current_phylum[b_current_phylum != "Other" & b_current_phylum != "others"]
b_sorted_phylum <- sort(b_non_other_phylum)
# 将排序后的门和 'Other' 组合起来，确保 'Other' 在最后
b_tax_order_phylum <- c(b_sorted_phylum)
if (length(b_other_category_phylum) > 0) {
  b_tax_order_phylum <- c(b_tax_order_phylum, b_other_category_phylum)
}
# 将 'tax' 列转换为因子，并应用排序
b_mean_phylum_abundance$tax <- factor(b_mean_phylum_abundance$tax, levels = b_tax_order_phylum)# 为绘图中的门和组设置顺序
# 确保 'Other' 类别在堆叠图和图例的最底部显示
if ("Other" %in% b_tax_order_phylum) {
  b_tax_order_phylum <- c(b_tax_order_phylum[b_tax_order_phylum != "Other"], "Other")
} else if ("others" %in% b_tax_order_phylum) { # 兼容你可能使用 'others' 的情况
  b_tax_order_phylum <- c(b_tax_order_phylum[b_tax_order_phylum != "others"], "others")
}
# 将 'tax' 和 'variable' 列转换为因子，并应用排序
b_mean_phylum_abundance$tax <- factor(b_mean_phylum_abundance$tax, levels = b_tax_order_phylum)
b_mean_phylum_abundance$variable <- factor(b_mean_phylum_abundance$variable, levels = unique(b_mean_phylum_abundance$variable))
# 定义自定义颜色板
# 获取绘图数据中实际存在的唯一 Phylum_Grouped 类别数量
b_unique_phyla_in_plot <- levels(b_mean_phylum_abundance$tax)
b_num_unique_phyla <- length(b_unique_phyla_in_plot)
# 你的自定义颜色列表
my_manual_colors_phylum <- c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
                             "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
                             "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                             "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                             "#cccccc") # 为 "Other" 添加一个灰色
# 检查颜色数量是否足够，如果不够则发出警告或进行调整
if (length(my_manual_colors_phylum) < b_num_unique_phyla) {
  warning("提供的颜色数量不足以覆盖所有门类别。某些颜色可能被重用或缺失。")
  # 可以在这里添加代码，自动生成更多颜色，例如使用 RColorBrewer
  # library(RColorBrewer)
  # my_manual_colors_phylum <- c(brewer.pal(n = b_num_unique_phyla - 1, name = "Set3"), "#cccccc") # 示例
} else if (length(my_manual_colors_phylum) > b_num_unique_phyla) {
  my_manual_colors_phylum <- my_manual_colors_phylum[1:b_num_unique_phyla] # 如果颜色太多则截取
}
# 确保 'Other' 类别被正确映射为灰色 (再次检查，以防裁剪后颜色被改变)
if ("Other" %in% b_unique_phyla_in_plot) {
  other_pos <- which(b_unique_phyla_in_plot == "Other")
  my_manual_colors_phylum[other_pos] <- "#cccccc"
}
# 将颜色名称映射到排序后的门类别，以确保正确的颜色分配
names(my_manual_colors_phylum) <- b_unique_phyla_in_plot
# 绘制堆叠柱状图 
b_mean_phylum_abundance$variable <- factor(b_mean_phylum_abundance$variable,
                                           levels = c("JRG", "JJG", "TZG", "PAG",
                                                      "JRN", "JJN", "TZN", "PAN"))
habitat_labeller <- c(
  "JRG" = "Jurong\nRhizosphere\nSoil",
  "JJG" = "Jingjiang\nRhizosphere\nSoil",
  "TZG" = "Tongzhou\nRhizosphere\nSoil",
  "PAG" = "Panan\nRhizosphere\nSoil",
  "JRN" = "Jurong\nBulb",
  "JJN" = "Jingjiang\nBulb",
  "TZN" = "Tongzhou\nBulb",
  "PAN" = "Panan\nBulb"
)
phylum_stack_bacteria <- ggplot(b_mean_phylum_abundance,
                                aes(x = variable, # 直接使用 'variable' (即 'Group')
                                    y = value,
                                    fill = factor(tax, levels = b_tax_order_phylum))) + # rev(b_tax_order_phylum) 用于图例顺序
  geom_bar(stat = "identity", position = "fill", width = 0.7) + # 'position = "fill"' 将Y轴缩放到0-1
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0,0)) + # Y轴标签为百分比
  coord_cartesian(ylim = c(0,1)) + # 强制 Y 轴范围为 0 到 1
  xlab("") + # X轴标签
  ylab("Bacterial relative abundance") + # Y轴标签
  theme_classic() + # 经典主题
  guides(fill = guide_legend(title = "Phylum")) + # 图例标题
  theme(legend.key.size = unit(0.4, "cm")) + # 图例键大小
  theme(text = element_text(family = "sans", size = 8)) + # 基础字体大小
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.background = element_blank(), # 移除面板背景
    panel.grid = element_blank(), # 移除网格线
    axis.text.y = element_text(size = 12, colour = "black", family = "sans", angle = 0),
    axis.text.x = element_text(size = 12, colour = "black", family = "sans", angle = 45, hjust = 1), # 旋转 X 轴标签
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # 标题居中加粗
    axis.line = element_line(size = 0.1, colour = "black") # 坐标轴线
  ) +
  scale_fill_manual(values = my_manual_colors_phylum) + # 应用自定义颜色
  scale_x_discrete(labels = habitat_labeller, # 应用自定义标签
                   guide = guide_axis(angle = 90))
phylum_stack_bacteria
# ggsave("D:/study/master/Extended_Data_tables/Figure_3/1a_phylum_stack_bacteria.png", plot = phylum_stack_bacteria, width = 10, height = 8, dpi = 600, bg = "transparent")
# 将数据转换为相对丰度
b_merged_genus_ra <- transform_sample_counts(b_merged, function(x) x / sum(x))
#聚合到属 (Genus) 水平
b_merged_genus <- tax_glom(b_merged_genus_ra, taxrank = "Genus", NArm = FALSE)
# 提取并处理数据框，筛选平均相对丰度 > 1% 的门 
b_genus_ra_df <- psmelt(b_merged_genus)
# 确保 Genus 列是字符类型
b_genus_ra_df$Genus <- as.character(b_genus_ra_df$Genus)
# 定义一个函数，用于判断 Genus 是否为“特殊”类型（未培养、未分类、无等级或NA）
is_special_genus <- function(genus_name) {
  is.na(genus_name) | # 检查是否为NA
    stringr::str_detect(genus_name, "uncultured") |
    stringr::str_detect(genus_name, "unclassified") |
    stringr::str_detect(genus_name, "norank")
}
# 如果 Genus 属于“特殊”类型，则将其标记为“Other_Special_Genus”。
b_genus_pre_grouped_df <- b_genus_ra_df %>%
  mutate(Genus_Pre_Grouped = ifelse(is_special_genus(Genus), "Other_Special_Genus", Genus))
# 计算当前 Genus 级别（已包含 "Other_Special_Genus"）的平均相对丰度
b_genus_mean_ra <- b_genus_pre_grouped_df %>%
  group_by(Genus_Pre_Grouped) %>% # 注意这里是对 Genus_Pre_Grouped 列分组
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选出最终的“Top Phyla”列表
b_top_phyla_final <- b_genus_mean_ra %>%
  filter(MeanAbundance > 0.005) %>%
  pull(Genus_Pre_Grouped)
#最终的 Genus 分组：将所有不符合条件的门归入“Other”
b_genus_final_grouped_df <- b_genus_pre_grouped_df %>%
  mutate(Genus_Grouped = ifelse(
    Genus_Pre_Grouped %in% b_top_phyla_final & Genus_Pre_Grouped != "Other_Special_Genus",
    Genus_Pre_Grouped, # 如果是 Top Phyla 且不是特殊的“Other”标记，则保留其名称
    "Other"             # 否则，归类为最终的“Other”
  ))
# 新计算每个样本中（新）分组门的相对丰度
b_genus_summarized_df <- b_genus_final_grouped_df %>%
  group_by(Sample, Group, Genus_Grouped) %>% # 确保 'Group' 列存在于 psmelt 结果中
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
b_genus_summarized_df <- b_genus_summarized_df %>%
  filter(!is.na(Group))
# 查看最终分组后的门类及其平均丰度
b_genus_summarized_df %>%
  group_by(Genus_Grouped) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(MeanAbundance)) %>%
  print()
# 计算每个门在所有样本中的平均相对丰度
b_genus_mean_ra <- b_genus_ra_df %>%
  group_by(Genus) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选平均相对丰度大于 0的门 
b_top_phyla <- b_genus_mean_ra %>%
  filter(MeanAbundance > 0) %>%
  pull(Genus)
# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
b_genus_summarized_df <- b_genus_summarized_df %>%
  filter(!is.na(Group))
# 准备绘图数据 (计算每个组的平均相对丰度)
b_mean_genus_abundance <- b_genus_summarized_df %>%
  group_by(Group, Genus_Grouped) %>% # 按组和处理后的门进行聚合
  summarise(MeanValue = mean(Abundance), .groups = 'drop') %>% # 计算平均丰度
  # 重命名列以匹配旧代码中的 'variable' 和 'value'
  rename(variable = Group, value = MeanValue, tax = Genus_Grouped) %>%
  as.data.frame()
# 移除 tax 列中的 'g__' 前缀
b_mean_genus_abundance$tax <- gsub("^g__", "", b_mean_genus_abundance$tax)
# 替换特定的长属名
b_mean_genus_abundance$tax <- recode(b_mean_genus_abundance$tax,
                                     "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" = "Rhizobium group",
                                     "Burkholderia-Caballeronia-Paraburkholderia" = "Burkholderia group",
                                     .default = as.character(b_mean_genus_abundance$tax))
# 首先获取所有唯一的门名称（包括 'Other'）
b_current_genus <- unique(b_mean_genus_abundance$tax)
# 将 'Other' 分离出来
b_other_category_genus <- b_current_genus[b_current_genus == "Other" | b_current_genus == "others"]
b_other_category_genus <- b_other_category_genus[1] # 确保只取一个 'Other' 或 'others'
# 排除 'Other'/'others' 后，对剩余的门进行字母排序
b_non_other_genus <- b_current_genus[b_current_genus != "Other" & b_current_genus != "others"]
b_sorted_genus <- sort(b_non_other_genus)
# 将排序后的门和 'Other' 组合起来，确保 'Other' 在最后
b_tax_order_genus <- c(b_sorted_genus)
if (length(b_other_category_genus) > 0) {
  b_tax_order_genus <- c(b_tax_order_genus, b_other_category_genus)
}
# 将 'tax' 列转换为因子，并应用排序
b_mean_genus_abundance$tax <- factor(b_mean_genus_abundance$tax, levels = b_tax_order_genus)# 为绘图中的门和组设置顺序
# 确保 'Other' 类别在堆叠图和图例的最底部显示
if ("Other" %in% b_tax_order_genus) {
  b_tax_order_genus <- c(b_tax_order_genus[b_tax_order_genus != "Other"], "Other")
} else if ("others" %in% tax_order) { # 兼容你可能使用 'others' 的情况
  b_tax_order_genus <- c(b_tax_order_genus[b_tax_order_genus != "others"], "others")
}
# 将 'tax' 和 'variable' 列转换为因子，并应用排序
b_mean_genus_abundance$tax <- factor(b_mean_genus_abundance$tax, levels = b_tax_order_genus)
b_mean_genus_abundance$variable <- factor(b_mean_genus_abundance$variable, levels = unique(b_mean_genus_abundance$variable))
# 定义自定义颜色板
# 获取绘图数据中实际存在的唯一 Genus_Grouped 类别数量
b_unique_genus_in_plot <- levels(b_mean_genus_abundance$tax)
b_num_unique_genus <- length(b_unique_genus_in_plot)
# 你的自定义颜色列表
my_manual_colors_genus <- c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
                            "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
                            "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                            "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                            "#F7CAD0","#F8E19B","#90E0C2","#00B4D8","#7A4D9B",
                            "#cccccc") # 为 "Other" 添加一个灰色
# 检查颜色数量是否足够，如果不够则发出警告或进行调整
if (length(my_manual_colors_genus) < b_num_unique_genus) {
  warning("提供的颜色数量不足以覆盖所有门类别。某些颜色可能被重用或缺失。")
} else if (length(my_manual_colors_genus) > b_num_unique_genus) {
  my_manual_colors_genus <- my_manual_colors_genus[1:b_num_unique_genus] # 如果颜色太多则截取
}
# 确保 'Other' 类别被正确映射为灰色 (再次检查，以防裁剪后颜色被改变)
if ("Other" %in% b_unique_genus_in_plot) {
  other_pos <- which(b_unique_genus_in_plot == "Other")
  my_manual_colors_genus[other_pos] <- "#cccccc"
}
# 将颜色名称映射到排序后的门类别，以确保正确的颜色分配
names(my_manual_colors_genus) <- b_unique_genus_in_plot
# 绘制堆叠柱状图 
b_mean_genus_abundance$variable <- factor(b_mean_genus_abundance$variable,
                                          levels = c("JRG", "JJG", "TZG", "PAG",
                                                     "JRN", "JJN", "TZN", "PAN"))
habitat_labeller <- c(
  "JRG" = "Jurong\nRhizosphere\nSoil",
  "JJG" = "Jingjiang\nRhizosphere\nSoil",
  "TZG" = "Tongzhou\nRhizosphere\nSoil",
  "PAG" = "Panan\nRhizosphere\nSoil",
  "JRN" = "Jurong\nBulb",
  "JJN" = "Jingjiang\nBulb",
  "TZN" = "Tongzhou\nBulb",
  "PAN" = "Panan\nBulb"
)
genus_stack_bacteria <- ggplot(b_mean_genus_abundance,
                               aes(x = variable, # 直接使用 'variable' (即 'Group')
                                   y = value,
                                   fill = factor(tax, levels = b_tax_order_genus))) + # rev(b_tax_order_genus) 用于图例顺序
  geom_bar(stat = "identity", position = "fill", width = 0.7) + # 'position = "fill"' 将Y轴缩放到0-1
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0,0)) + # Y轴标签为百分比
  coord_cartesian(ylim = c(0,1)) + # 强制 Y 轴范围为 0 到 1
  xlab("") + # X轴标签
  ylab("Bacterial relative abundance") + # Y轴标签
  theme_classic() + # 经典主题
  guides(fill = guide_legend(title = "Genus")) + # 图例标题
  theme(legend.key.size = unit(0.4, "cm")) + # 图例键大小
  theme(text = element_text(family = "sans", size = 8)) + # 基础字体大小
  theme(
    legend.text = element_text(size = 12 , face = "italic"),
    legend.title = element_text(size = 12),
    panel.background = element_blank(), # 移除面板背景
    panel.grid = element_blank(), # 移除网格线
    axis.text.y = element_text(size = 12, colour = "black", family = "sans", angle = 0),
    axis.text.x = element_text(size = 12, colour = "black", family = "sans", angle = 45, hjust = 1), # 旋转 X 轴标签
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # 标题居中加粗
    axis.line = element_line(size = 0.1, colour = "black") # 坐标轴线
  ) +
  scale_fill_manual(values = my_manual_colors_genus) + # 应用自定义颜色
  scale_x_discrete(labels = habitat_labeller, # 应用自定义标签
                   guide = guide_axis(angle = 90))
genus_stack_bacteria
# ggsave("D:/study/master/Extended_Data_tables/Figure_3/1a_genus_stack_bacteria.png", plot = genus_stack_bacteria, width = 10, height = 8, dpi = 600, bg = "transparent")
#真菌组成
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
f_tax <- data.frame(f_tax)
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
# 提取样本元数据
metadata <- as_tibble(sample_data(f_merged))
# 将数据转换为相对丰度
f_merged_phylum_ra <- transform_sample_counts(f_merged, function(x) x / sum(x))
#聚合到门 (Phylum) 水平
f_merged_phylum <- tax_glom(f_merged_phylum_ra, taxrank = "Phylum", NArm = FALSE)
# 提取并处理数据框 
f_phylum_ra_df <- psmelt(f_merged_phylum)
# 确保 Phylum 列是字符类型
f_phylum_ra_df$Phylum <- as.character(f_phylum_ra_df$Phylum)
# 定义一个函数，用于判断 Phylum 是否为“特殊”类型（未培养、未分类、无等级或NA）
is_special_phylum <- function(phylum_name) {
  is.na(phylum_name) | # 检查是否为NA
    stringr::str_detect(phylum_name, "uncultured") |
    stringr::str_detect(phylum_name, "unclassified") |
    stringr::str_detect(phylum_name, "norank")
}
# 如果 Phylum 属于“特殊”类型，则将其标记为“Other_Special_Phylum”。
f_phylum_pre_grouped_df <- f_phylum_ra_df %>%
  mutate(Phylum_Pre_Grouped = ifelse(is_special_phylum(Phylum), "Other_Special_Phylum", Phylum))
# 计算当前 Phylum 级别（已包含 "Other_Special_Phylum"）的平均相对丰度
f_phylum_mean_ra <- f_phylum_pre_grouped_df %>%
  group_by(Phylum_Pre_Grouped) %>% # 注意这里是对 Phylum_Pre_Grouped 列分组
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选出最终的“Top Phyla”列表
f_top_phyla_final <- f_phylum_mean_ra %>%
  filter(MeanAbundance > 0.001) %>%
  pull(Phylum_Pre_Grouped)
#最终的 Phylum 分组：将所有不符合条件的门归入“Other”
f_phylum_final_grouped_df <- f_phylum_pre_grouped_df %>%
  mutate(Phylum_Grouped = ifelse(
    Phylum_Pre_Grouped %in% f_top_phyla_final & Phylum_Pre_Grouped != "Other_Special_Phylum",
    Phylum_Pre_Grouped, # 如果是 Top Phyla 且不是特殊的“Other”标记，则保留其名称
    "Other"             # 否则，归类为最终的“Other”
  ))
# 新计算每个样本中（新）分组门的相对丰度
f_phylum_summarized_df <- f_phylum_final_grouped_df %>%
  group_by(Sample, Group, Phylum_Grouped) %>% # 确保 'Group' 列存在于 psmelt 结果中
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
f_phylum_summarized_df <- f_phylum_summarized_df %>%
  filter(!is.na(Group))
# 查看最终分组后的门类及其平均丰度
f_phylum_summarized_df %>%
  group_by(Phylum_Grouped) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(MeanAbundance)) %>%
  print()
# 计算每个门在所有样本中的平均相对丰度
f_phylum_mean_ra <- f_phylum_ra_df %>%
  group_by(Phylum) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选平均相对丰度大于 0的门 
f_top_phyla <- f_phylum_mean_ra %>%
  filter(MeanAbundance > 0) %>%
  pull(Phylum)
# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
f_phylum_summarized_df <- f_phylum_summarized_df %>%
  filter(!is.na(Group))
# 准备绘图数据 (计算每个组的平均相对丰度)
f_mean_phylum_abundance <- f_phylum_summarized_df %>%
  group_by(Group, Phylum_Grouped) %>% # 按组和处理后的门进行聚合
  summarise(MeanValue = mean(Abundance), .groups = 'drop') %>% # 计算平均丰度
  # 重命名列以匹配旧代码中的 'variable' 和 'value'
  rename(variable = Group, value = MeanValue, tax = Phylum_Grouped) %>%
  as.data.frame()
# 移除 tax 列中的 'p__' 前缀
f_mean_phylum_abundance$tax <- gsub("^p__", "", f_mean_phylum_abundance$tax)
# 首先获取所有唯一的门名称（包括 'Other'）
f_current_phylum <- unique(f_mean_phylum_abundance$tax)
# 将 'Other' 分离出来
f_other_category_phylum <- f_current_phylum[f_current_phylum == "Other" | f_current_phylum == "others"]
f_other_category_phylum <- f_other_category_phylum[1] # 确保只取一个 'Other' 或 'others'
# 排除 'Other'/'others' 后，对剩余的门进行字母排序
f_non_other_phylum <- f_current_phylum[f_current_phylum != "Other" & f_current_phylum != "others"]
f_sorted_phylum <- sort(f_non_other_phylum)
# 将排序后的门和 'Other' 组合起来，确保 'Other' 在最后
f_tax_order_phylum <- c(f_sorted_phylum)
if (length(f_other_category_phylum) > 0) {
  f_tax_order_phylum <- c(f_tax_order_phylum, f_other_category_phylum)
}
# 将 'tax' 列转换为因子，并应用排序
f_mean_phylum_abundance$tax <- factor(f_mean_phylum_abundance$tax, levels = f_tax_order_phylum)# 为绘图中的门和组设置顺序
# 确保 'Other' 类别在堆叠图和图例的最底部显示
if ("Other" %in% f_tax_order_phylum) {
  f_tax_order_phylum <- c(f_tax_order_phylum[f_tax_order_phylum != "Other"], "Other")
} else if ("others" %in% f_tax_order_phylum) { # 兼容你可能使用 'others' 的情况
  f_tax_order_phylum <- c(f_tax_order_phylum[f_tax_order_phylum != "others"], "others")
}
# 将 'tax' 和 'variable' 列转换为因子，并应用排序
f_mean_phylum_abundance$tax <- factor(f_mean_phylum_abundance$tax, levels = f_tax_order_phylum)
f_mean_phylum_abundance$variable <- factor(f_mean_phylum_abundance$variable, levels = unique(f_mean_phylum_abundance$variable))
# 定义自定义颜色板
# 获取绘图数据中实际存在的唯一 Phylum_Grouped 类别数量
f_unique_phyla_in_plot <- levels(f_mean_phylum_abundance$tax)
f_num_unique_phyla <- length(f_unique_phyla_in_plot)
# 你的自定义颜色列表
my_manual_colors_phylum <- c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
                             "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
                             "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                             "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                             "#cccccc") # 为 "Other" 添加一个灰色
# 检查颜色数量是否足够，如果不够则发出警告或进行调整
if (length(my_manual_colors_phylum) < f_num_unique_phyla) {
  warning("提供的颜色数量不足以覆盖所有门类别。某些颜色可能被重用或缺失。")
  # 可以在这里添加代码，自动生成更多颜色，例如使用 RColorBrewer
  # library(RColorBrewer)
  # my_manual_colors_phylum <- c(brewer.pal(n = f_num_unique_phyla - 1, name = "Set3"), "#cccccc") # 示例
} else if (length(my_manual_colors_phylum) > f_num_unique_phyla) {
  my_manual_colors_phylum <- my_manual_colors_phylum[1:f_num_unique_phyla] # 如果颜色太多则截取
}
# 确保 'Other' 类别被正确映射为灰色 (再次检查，以防裁剪后颜色被改变)
if ("Other" %in% f_unique_phyla_in_plot) {
  other_pos <- which(f_unique_phyla_in_plot == "Other")
  my_manual_colors_phylum[other_pos] <- "#cccccc"
}
# 将颜色名称映射到排序后的门类别，以确保正确的颜色分配
names(my_manual_colors_phylum) <- f_unique_phyla_in_plot
# 绘制堆叠柱状图 
f_mean_phylum_abundance$variable <- factor(f_mean_phylum_abundance$variable,
                                           levels = c("JRG", "JJG", "TZG", "PAG",
                                                      "JRN", "JJN", "TZN", "PAN"))
habitat_labeller <- c(
  "JRG" = "Jurong\nRhizosphere\nSoil",
  "JJG" = "Jingjiang\nRhizosphere\nSoil",
  "TZG" = "Tongzhou\nRhizosphere\nSoil",
  "PAG" = "Panan\nRhizosphere\nSoil",
  "JRN" = "Jurong\nBulb",
  "JJN" = "Jingjiang\nBulb",
  "TZN" = "Tongzhou\nBulb",
  "PAN" = "Panan\nBulb"
)
phylum_stack_fungi <- ggplot(f_mean_phylum_abundance,
                             aes(x = variable, # 直接使用 'variable' (即 'Group')
                                 y = value,
                                 fill = factor(tax, levels = f_tax_order_phylum))) + # rev(f_tax_order_phylum) 用于图例顺序
  geom_bar(stat = "identity", position = "fill", width = 0.7) + # 'position = "fill"' 将Y轴缩放到0-1
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0,0)) + # Y轴标签为百分比
  coord_cartesian(ylim = c(0,1)) + # 强制 Y 轴范围为 0 到 1
  xlab("") + # X轴标签
  ylab("Fungal relative abundance") + # Y轴标签
  theme_classic() + # 经典主题
  guides(fill = guide_legend(title = "Phylum")) + # 图例标题
  theme(legend.key.size = unit(0.4, "cm")) + # 图例键大小
  theme(text = element_text(family = "sans", size = 8)) + # 基础字体大小
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.background = element_blank(), # 移除面板背景
    panel.grid = element_blank(), # 移除网格线
    axis.text.y = element_text(size = 12, colour = "black", family = "sans", angle = 0),
    axis.text.x = element_text(size = 12, colour = "black", family = "sans", angle = 45, hjust = 1), # 旋转 X 轴标签
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # 标题居中加粗
    axis.line = element_line(size = 0.1, colour = "black") # 坐标轴线
  ) +
  scale_fill_manual(values = my_manual_colors_phylum) + # 应用自定义颜色
  scale_x_discrete(labels = habitat_labeller, # 应用自定义标签
                   guide = guide_axis(angle = 90))
phylum_stack_fungi
# ggsave("D:/study/master/Extended_Data_tables/Figure_3/1a_phylum_stack_fungi.png", plot = phylum_stack_fungi, width = 10, height = 8, dpi = 600, bg = "transparent")
# 将数据转换为相对丰度
f_merged_genus_ra <- transform_sample_counts(f_merged, function(x) x / sum(x))
#聚合到属 (Genus) 水平
f_merged_genus <- tax_glom(f_merged_genus_ra, taxrank = "Genus", NArm = FALSE)
# 提取并处理数据框，筛选平均相对丰度 > 1% 的门 
f_genus_ra_df <- psmelt(f_merged_genus)
# 确保 Genus 列是字符类型
f_genus_ra_df$Genus <- as.character(f_genus_ra_df$Genus)
# 定义一个函数，用于判断 Genus 是否为“特殊”类型（未培养、未分类、无等级或NA）
is_special_genus <- function(genus_name) {
  is.na(genus_name) | # 检查是否为NA
    stringr::str_detect(genus_name, "uncultured") |
    stringr::str_detect(genus_name, "unclassified") |
    stringr::str_detect(genus_name, "norank")
}
# 如果 Genus 属于“特殊”类型，则将其标记为“Other_Special_Genus”。
f_genus_pre_grouped_df <- f_genus_ra_df %>%
  mutate(Genus_Pre_Grouped = ifelse(is_special_genus(Genus), "Other_Special_Genus", Genus))
# 计算当前 Genus 级别（已包含 "Other_Special_Genus"）的平均相对丰度
f_genus_mean_ra <- f_genus_pre_grouped_df %>%
  group_by(Genus_Pre_Grouped) %>% # 注意这里是对 Genus_Pre_Grouped 列分组
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选出最终的“Top Phyla”列表
f_top_phyla_final <- f_genus_mean_ra %>%
  filter(MeanAbundance > 0.005) %>%
  pull(Genus_Pre_Grouped)
#最终的 Genus 分组：将所有不符合条件的门归入“Other”
f_genus_final_grouped_df <- f_genus_pre_grouped_df %>%
  mutate(Genus_Grouped = ifelse(
    Genus_Pre_Grouped %in% f_top_phyla_final & Genus_Pre_Grouped != "Other_Special_Genus",
    Genus_Pre_Grouped, # 如果是 Top Phyla 且不是特殊的“Other”标记，则保留其名称
    "Other"             # 否则，归类为最终的“Other”
  ))
# 新计算每个样本中（新）分组门的相对丰度
f_genus_summarized_df <- f_genus_final_grouped_df %>%
  group_by(Sample, Group, Genus_Grouped) %>% # 确保 'Group' 列存在于 psmelt 结果中
  summarise(Abundance = sum(Abundance), .groups = 'drop')

# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
f_genus_summarized_df <- f_genus_summarized_df %>%
  filter(!is.na(Group))
# 查看最终分组后的门类及其平均丰度
f_genus_summarized_df %>%
  group_by(Genus_Grouped) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(MeanAbundance)) %>%
  print()
# 计算每个门在所有样本中的平均相对丰度
f_genus_mean_ra <- f_genus_ra_df %>%
  group_by(Genus) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()
# 筛选平均相对丰度大于 0的门 
f_top_phyla <- f_genus_mean_ra %>%
  filter(MeanAbundance > 0) %>%
  pull(Genus)
# 过滤掉 'Group' 列为 NA 的行（如果存在），防止绘图警告
f_genus_summarized_df <- f_genus_summarized_df %>%
  filter(!is.na(Group))
# 准备绘图数据 (计算每个组的平均相对丰度)
f_mean_genus_abundance <- f_genus_summarized_df %>%
  group_by(Group, Genus_Grouped) %>% # 按组和处理后的门进行聚合
  summarise(MeanValue = mean(Abundance), .groups = 'drop') %>% # 计算平均丰度
  # 重命名列以匹配旧代码中的 'variable' 和 'value'
  rename(variable = Group, value = MeanValue, tax = Genus_Grouped) %>%
  as.data.frame()
# 移除 tax 列中的 'g__' 前缀
f_mean_genus_abundance$tax <- gsub("^g__", "", f_mean_genus_abundance$tax)
# 替换特定的长属名
f_mean_genus_abundance$tax <- recode(f_mean_genus_abundance$tax,
                                     "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" = "Rhizobium group",
                                     "Burkholderia-Caballeronia-Paraburkholderia" = "Burkholderia group",
                                     .default = as.character(f_mean_genus_abundance$tax))
# 首先获取所有唯一的门名称（包括 'Other'）
f_current_genus <- unique(f_mean_genus_abundance$tax)
# 将 'Other' 分离出来
f_other_category_genus <- f_current_genus[f_current_genus == "Other" | f_current_genus == "others"]
f_other_category_genus <- f_other_category_genus[1] # 确保只取一个 'Other' 或 'others'
# 排除 'Other'/'others' 后，对剩余的门进行字母排序
f_non_other_genus <- f_current_genus[f_current_genus != "Other" & f_current_genus != "others"]
f_sorted_genus <- sort(f_non_other_genus)
# 将排序后的门和 'Other' 组合起来，确保 'Other' 在最后
f_tax_order_genus <- c(f_sorted_genus)
if (length(f_other_category_genus) > 0) {
  f_tax_order_genus <- c(f_tax_order_genus, f_other_category_genus)
}
# 将 'tax' 列转换为因子，并应用排序
f_mean_genus_abundance$tax <- factor(f_mean_genus_abundance$tax, levels = f_tax_order_genus)# 为绘图中的门和组设置顺序
# 确保 'Other' 类别在堆叠图和图例的最底部显示
if ("Other" %in% f_tax_order_genus) {
  f_tax_order_genus <- c(f_tax_order_genus[f_tax_order_genus != "Other"], "Other")
} else if ("others" %in% tax_order) { # 兼容你可能使用 'others' 的情况
  f_tax_order_genus <- c(f_tax_order_genus[f_tax_order_genus != "others"], "others")
}
# 将 'tax' 和 'variable' 列转换为因子，并应用排序
f_mean_genus_abundance$tax <- factor(f_mean_genus_abundance$tax, levels = f_tax_order_genus)
f_mean_genus_abundance$variable <- factor(f_mean_genus_abundance$variable, levels = unique(f_mean_genus_abundance$variable))
# 定义自定义颜色板
# 获取绘图数据中实际存在的唯一 Genus_Grouped 类别数量
f_unique_genus_in_plot <- levels(f_mean_genus_abundance$tax)
f_num_unique_genus <- length(f_unique_genus_in_plot)
# 你的自定义颜色列表
my_manual_colors_genus <- c("#d2da93","#5196d5","#00ceff","#ff630d","#35978b",
                            "#e5acd7","#77aecd","#ec8181","#dfc6a5","#e50719",
                            "#d27e43","#8a4984","#fe5094","#8d342e","#f94e54",
                            "#ffad00","#36999d","#00fc8d","#b64aa0","#9b82e1",
                            "#F7CAD0","#F8E19B","#90E0C2","#00B4D8","#7A4D9B",
                            "#cccccc") # 为 "Other" 添加一个灰色
# 检查颜色数量是否足够，如果不够则发出警告或进行调整
if (length(my_manual_colors_genus) < f_num_unique_genus) {
  warning("提供的颜色数量不足以覆盖所有门类别。某些颜色可能被重用或缺失。")
} else if (length(my_manual_colors_genus) > f_num_unique_genus) {
  my_manual_colors_genus <- my_manual_colors_genus[1:f_num_unique_genus] # 如果颜色太多则截取
}
# 确保 'Other' 类别被正确映射为灰色 (再次检查，以防裁剪后颜色被改变)
if ("Other" %in% f_unique_genus_in_plot) {
  other_pos <- which(f_unique_genus_in_plot == "Other")
  my_manual_colors_genus[other_pos] <- "#cccccc"
}
# 将颜色名称映射到排序后的门类别，以确保正确的颜色分配
names(my_manual_colors_genus) <- f_unique_genus_in_plot
# 绘制堆叠柱状图 
f_mean_genus_abundance$variable <- factor(f_mean_genus_abundance$variable,
                                          levels = c("JRG", "JJG", "TZG", "PAG",
                                                     "JRN", "JJN", "TZN", "PAN"))
habitat_labeller <- c(
  "JRG" = "Jurong\nRhizosphere\nSoil",
  "JJG" = "Jingjiang\nRhizosphere\nSoil",
  "TZG" = "Tongzhou\nRhizosphere\nSoil",
  "PAG" = "Panan\nRhizosphere\nSoil",
  "JRN" = "Jurong\nBulb",
  "JJN" = "Jingjiang\nBulb",
  "TZN" = "Tongzhou\nBulb",
  "PAN" = "Panan\nBulb"
)
genus_stack_fungi <- ggplot(f_mean_genus_abundance,
                            aes(x = variable, # 直接使用 'variable' (即 'Group')
                                y = value,
                                fill = factor(tax, levels = f_tax_order_genus))) + # rev(f_tax_order_genus) 用于图例顺序
  geom_bar(stat = "identity", position = "fill", width = 0.7) + # 'position = "fill"' 将Y轴缩放到0-1
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01), expand = c(0,0)) + # Y轴标签为百分比
  coord_cartesian(ylim = c(0,1)) + # 强制 Y 轴范围为 0 到 1
  xlab("") + # X轴标签
  ylab("Fungal relative abundance") + # Y轴标签
  theme_classic() + # 经典主题
  guides(fill = guide_legend(title = "Genus")) + # 图例标题
  theme(legend.key.size = unit(0.4, "cm")) + # 图例键大小
  theme(text = element_text(family = "sans", size = 8)) + # 基础字体大小
  theme(
    legend.text = element_text(size = 12 , face = "italic"),
    legend.title = element_text(size = 12),
    panel.background = element_blank(), # 移除面板背景
    panel.grid = element_blank(), # 移除网格线
    axis.text.y = element_text(size = 12, colour = "black", family = "sans", angle = 0),
    axis.text.x = element_text(size = 12, colour = "black", family = "sans", angle = 45, hjust = 1), # 旋转 X 轴标签
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), # 标题居中加粗
    axis.line = element_line(size = 0.1, colour = "black") # 坐标轴线
  ) +
  scale_fill_manual(values = my_manual_colors_genus) + # 应用自定义颜色
  scale_x_discrete(labels = habitat_labeller, # 应用自定义标签
                   guide = guide_axis(angle = 90))
genus_stack_fungi
# ggsave("D:/study/master/Extended_Data_tables/Figure_3/1a_genus_stack_fungi.png", plot = genus_stack_fungi, width = 10, height = 8, dpi = 600, bg = "transparent")

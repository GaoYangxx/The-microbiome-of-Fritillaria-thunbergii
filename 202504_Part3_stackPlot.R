#图3b
install.packages("devtools")
install.packages(c("tidyverse", "RColorBrewer"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("phyloseq", "indicspecies"))
devtools::install_github("benjjneb/speedyseq")
devtools::install_github("joey711/phyloseqCompanion")
library(indicspecies)
library(speedyseq)
library(phyloseq)
library(tidyverse)
library(phyloseqCompanion)
library(RColorBrewer)
#细菌堆叠图
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
b_ASV <- as.matrix(otu_table(b_merged, taxa_are_rows=T))
b_ASV_df <- b_ASV[rowSums(b_ASV[])>0,]
#核心asv（在大于2个居群生态位中相对丰度大于0.1%）
# 计算 ASV 相对丰度
b_ASV_rel <- sweep(b_ASV_df, 2, colSums(b_ASV_df), "/")
# 先把样本、niche 和 ASV 丰度组合成一个长表
b_rel_long <- melt(as.matrix(b_ASV_rel), varnames = c("ASV", "Sample"))
sample_niches <- metadata$Niche#样本到生态位（niche）的映射向量
names(sample_niches) <- metadata$Sample.ID
b_rel_long$Niche <- sample_niches[b_rel_long$Sample]
# 按 ASV 和 Niche 聚合：计算平均相对丰度
b_asv_niche_mean <- b_rel_long %>%
  group_by(ASV, Niche) %>%
  summarise(mean_abundance = mean(value), .groups = "drop")
# 对每个 ASV，看有多少个 niche平均相对丰度超过阈值
b_asv_niche_count <- b_asv_niche_mean %>%
  filter(mean_abundance >= 0.001) %>%
  group_by(ASV) %>%
  summarise(niche_count = n())
# 筛选核心 ASV：在 =2 个 niche 中都有平均丰度 >0.1%
b_core_asvs_names <- b_asv_niche_count$ASV[b_asv_niche_count$niche_count == 2]
# 提取核心 ASV 的相对丰度数据
b_core_asv_data <- b_ASV_rel[b_core_asvs_names, ]
# 输出核心 ASV 数量
cat("核心细菌 ASV 数量:", length(b_core_asvs_names), "\n")
#生成核心ASV的phyloseq
b_merge_core_abundance <- b_ASV_df[rownames(b_ASV_df) %in% rownames(b_core_asv_data),]
b_mca_table <- otu_table(b_merge_core_abundance, taxa_are_rows=T)
b_merged_NOMIS_core_ab<- merge_phyloseq(b_mca_table, b_tax, metadata) 
#最特殊asv（只在一个分组中出现）
metadata_nomis <- sample.data.frame(b_merged)#提取元数据
b_asv_df <- as.data.frame(t(otu_table(b_ASV_df, taxa_are_rows=T))) #提取转置asv数量表
b_asvdfmelt <- melt(as.matrix(b_asv_df))#长宽表格转换
b_asvdfmelt <- b_asvdfmelt[b_asvdfmelt$value >0,]
b_control_asv <- b_asvdfmelt
b_control_asv$Group <- vapply((b_control_asv$Var1), function(x) metadata_nomis$Group[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))# 添加分组信息
b_controlasv_endemic <- b_control_asv %>%
  distinct(Group, Var2) %>%         # 去重：只保留一个Group-ASV组合
  group_by(Var2) %>%
  tally(name = "group_count") %>%   # 统计每个 ASV 出现在多少个 Group 中
  filter(group_count == 1)          # 只保留那些只在一个 Group 中出现的
b_merge_controlasv_abundance <- b_ASV_df[rownames(b_ASV_df) %in% b_controlasv_endemic$Var2,]
b_controlasv_table <- otu_table(b_merge_controlasv_abundance, taxa_are_rows=T)
b_merged_NOMIS_controlasv_ab<- merge_phyloseq(b_controlasv_table, b_tax, metadata) 
#只出现在一个 Group 中所有样本的asv
#从原始数据中提取只出现在一个 Group 中 ASV 的样本信息
b_asv_group_info <- b_control_asv %>%
  filter(Var2 %in% b_controlasv_endemic$Var2) %>%
  group_by(Var2, Group) %>%
  summarise(sample_count = n_distinct(Var1), .groups = "drop") %>%#统计每个 ASV 在这个 group 中 有丰度的样本数量（非零）
  filter(sample_count == 6)  # 只保留在某个组中6个样本都出现的ASV
#提取这些 ASV 的丰度表
b_final_asv_ids <- b_asv_group_info$Var2
b_final_asv_abundance <- b_ASV_df[rownames(b_ASV_df) %in% b_final_asv_ids, ]
length(b_final_asv_ids)
#构建新的 phyloseq 对象（可选）
b_final_asv_table <- otu_table(as.matrix(b_final_asv_abundance), taxa_are_rows = TRUE)
b_final_ps <- merge_phyloseq(b_final_asv_table, b_tax, metadata)
#查看结果（按需提取特定分组的样本等）
View(as.data.frame(otu_table(b_final_ps)))
#比较特殊asv（只在一个生态位中出现）
b_specific_asv <- b_asvdfmelt
b_specific_asv$Niche <- vapply((b_specific_asv$Var1), function(x) metadata_nomis$Niche[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))# 添加生态位信息
b_specific_endemic <- b_specific_asv %>%
  distinct(Niche, Var2) %>%         # 去重：只保留一个Niche-ASV组合
  group_by(Var2) %>%
  tally(name = "niche_count") %>%   # 统计每个 ASV 出现在多少个 Niche 中
  filter(niche_count == 1)          # 只保留那些只在一个 Niche 中出现的
b_merge_specific_abundance <- b_ASV_df[rownames(b_ASV_df) %in% b_specific_endemic$Var2,]
b_specific_table <- otu_table(b_merge_specific_abundance, taxa_are_rows=T)
b_merged_NOMIS_specific_ab<- merge_phyloseq(b_specific_table, b_tax, metadata) 
# 指示性asv，运行 multipatt 进行指示性分析，找出哪些 ASV 是某些 niche 的指示物种
b_indicator_multipatt <- multipatt(t(b_ASV_df), metadata$Niche, 
                                   func = "r.g", control = how(nperm = 999))
summary(b_indicator_multipatt)
# 获取显著的指示 ASV（p < 0.01）
b_indicator_names <- b_indicator_multipatt$sign[which(b_indicator_multipatt$sign$p.value < 0.01), ]
print(b_indicator_names)
#生成指示性ASV的phyloseq
b_merge_indicator_abundance <- b_ASV_df[rownames(b_ASV_df) %in% rownames(b_indicator_names),]
b_indicator_table <- otu_table(b_merge_indicator_abundance, taxa_are_rows=T)
b_merged_NOMIS_indicator_ab<- merge_phyloseq(b_indicator_table, b_tax, metadata) 
#统一读取
b_dat	<- as.data.frame(otu_table(b_merged, taxa_are_rows=T))
b_cores <- as.data.frame(otu_table(b_merged_NOMIS_core_ab, taxa_are_rows=T)) 
b_endemic <- as.data.frame(otu_table(b_merged_NOMIS_specific_ab, taxa_are_rows=T))#比较特殊
b_indicator <- as.data.frame(otu_table(b_merged_NOMIS_indicator_ab, taxa_are_rows=T))
b_uniq <- as.data.frame(otu_table(b_merged_NOMIS_controlasv_ab, taxa_are_rows=T))#最特殊
b_ASV <- rownames(b_dat)
b_ASV <- as.data.frame(b_ASV)
colnames(b_ASV) <- "ASV"
b_cores$ASV<-rownames(b_cores)
rownames(b_cores)<-NULL
b_cores <- b_cores %>%
  select(ASV)%>%# 只保留 ASV 这一列
  mutate(type = "core")#增加一列 type，所有行都标记为 "core"
b_endemic$ASV<-rownames(b_endemic)
rownames(b_endemic)<-NULL
b_endemic <- b_endemic %>%
  select(ASV)%>%
  mutate(type = "endemic")
colnames(b_endemic) <- c("ASV", "type")
b_indicator$ASV<-rownames(b_indicator)
rownames(b_indicator)<-NULL
b_indicator <- b_indicator %>%
  select(ASV)%>%
  mutate(type = "indicator")
colnames(b_indicator) <- c("ASV","type")
b_uniq$ASV<-rownames(b_uniq)
rownames(b_uniq)<-NULL
b_uniq <- b_uniq %>%
  select(ASV)%>%
  mutate(type = "uniq")
colnames(b_uniq) <- c("ASV","type")
b_dat_numb <- rbind(b_cores, b_endemic, b_uniq, b_indicator)#合并
b_ASV_other <- b_ASV %>%#筛选出那些没有出现在四类标签里的ASV
  filter(!(ASV %in% b_dat_numb$ASV))%>%
  mutate(type = "Other")
#最特殊asv
length(unique(row.names(b_uniq)))/length(unique(row.names(b_ASV_df)))# 占 总 ASV 数量 的比例
b_merge_uniq <- merge(b_uniq, b_asvdfmelt, by.x="ASV",by.y="Var2")#合并reads 数量
b_merge_uniq$Group <- vapply((b_merge_uniq$Var1), function(x) metadata_nomis$Group[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))#加列
colnames(b_merge_uniq) <- c("ASV","type","SampleID","nb_count","Group")#asv名字、出现在所有样本中的样本次数、样本名、asv出现在某一样本的reads数量、组名
b_uniq_Group_prev <- b_merge_uniq %>%
  group_by(Group) %>%             # 按照 "Group" 列进行分组
  summarize(prev = n())           # 计算每个分组的行数（即该组中只在一个样本出现的ASV数量）
b_uniq_one_plot <- b_merge_uniq %>% group_by(Group) 
b_uniq_one_plot<- b_uniq_one_plot[c("Group","ASV","type")]
b_uniq_one_plot$Color <- "Specific and unique"
# 计算每个 Group 的最特殊 ASV 数量 (prev) 占所有最特殊 ASV 总数量的比例
b_prop_unique_Group<- b_uniq_Group_prev%>%
  group_by(Group)%>% 
  summarize(prop=prev/length(unique(row.names(b_uniq)))) 
colnames(b_prop_unique_Group) <- c("group_range","prop_unique")
#其他asv
length(unique(row.names(b_ASV_other)))/length(unique(row.names(b_ASV_df)))# 占 总 ASV 数量 的比例
b_merge_other <- merge(b_ASV_other, b_asvdfmelt, by.x="ASV",by.y="Var2")#合并reads 数量
b_merge_other$Group <- vapply((b_merge_other$Var1), function(x) metadata_nomis$Group[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))#加列
colnames(b_merge_other) <- c("ASV","type","SampleID","nb_count","Group")#asv名字、出现在所有样本中的样本次数、样本名、asv出现在某一样本的reads数量、组名
b_other_Group_prev <- b_merge_other %>%
  group_by(Group) %>%             # 按照 "Group" 列进行分组
  summarize(prev = n_distinct(ASV))           # 计算每个分组的其他asv种类
b_other_one_plot <- b_merge_other %>% group_by(Group) 
b_other_one_plot<- b_other_one_plot[c("Group","ASV","type")]
b_other_one_plot$Color <- "Other"
# 计算每个 Group 的其他 ASV 数量 (prev) 占所有其他ASV 总数量的比例
b_prop_other_Group<- b_other_Group_prev%>%
  group_by(Group)%>% 
  summarize(prop=prev/length(unique(row.names(b_ASV_other)))) 
colnames(b_prop_other_Group) <- c("group_range","prop_other")
#合并最特殊和其他asv
b_df_full <- rbind(b_uniq_one_plot, b_other_one_plot)
# 使用 distinct() 保留 Group 和 ASV 组合唯一的行，.keep_all = TRUE 会保留该唯一行对应的原始数据框中的所有其他列
b_df_full_clean <- b_df_full %>%
  distinct(Group, ASV, .keep_all = TRUE)
niveaux <- c("JRG", "JJG", "TZG", "PAG","JRN", "JJN", "TZN", "PAN")
b_df_full_clean$Group<- factor(b_df_full_clean$Group, levels = niveaux) 
b_df_full_clean$Group <- fct_rev(b_df_full_clean$Group) # 将刚才设置的 Group 因子水平顺序反转
# 定义堆叠图中 Color 列的因子水平顺序
stacking_levels_order <- c("Other", "Specific and unique")
b_df_full_clean$Color <- factor(b_df_full_clean$Color, levels = stacking_levels_order) 
# 定义图例中 Color 列的因子水平顺序
legend_key_order_labels <- c("Specific and unique", "Other") 
plot_colors <- c("#606060", "#7FA3D6") # 假设 Other 用灰色, Specific 用蓝色
names(plot_colors) <- stacking_levels_order # 为颜色命名以匹配因子水平
# 生成柱状图，使用设定的堆叠顺序和图例顺序
stack_plot_bacteria <- ggplot(b_df_full_clean, aes(fill = Color, y = Group)) +
  geom_bar(position = "stack", stat = "count") + # stat="count" 根据行数堆叠
  labs(
    x = "Number of specific bacterial ASVs", # 设置横坐标标题
    y = NULL                       # 移除纵坐标标题 (或者 y = "")
  ) +
  theme_classic() +
  # 使用 scale_fill_manual 手动设定填充颜色和图例项目顺序
  scale_fill_manual(
    values = plot_colors, #指定堆叠图中显示的分类顺序
    breaks = legend_key_order_labels # 指定图例中显示的分类标签及其顺序
  ) +
  scale_x_continuous(
    limits = c(0, 23000), # 保持轴的范围到 25000
    breaks = seq(0, 20000, by = 5000), # 保持刻度和标签只到 20000
    expand = expansion(mult = c(0, 0)) # 新增或修改 expand 参数，左侧扩展为 0
  ) +
  theme(
    # 移除纵坐标轴刻度标签和刻度线
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(), # 确保 Y 轴线仍然存在 (theme_classic 默认包含)
    # 设置横坐标轴标题和刻度标签的大小
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14), 
    axis.ticks.x = element_line(), # 确保刻度线本身是可见的 (theme_classic 默认包含)
    axis.ticks.length.x = unit(-0.25, "cm"), # 设置 x 轴刻度线朝内，并指定长度
    # 移除图例标题
    legend.title = element_blank(),
    # 设置图例文字大小
    legend.text = element_text(size = 14),
    # 设置图例位置在绘图区域内部的右下角
    legend.position = c(1, 0.1),
    legend.justification = c(1, 0) # 设置图例框的右下角对齐到坐标(1, 0)
  )  
stack_plot_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3b_stack_bacteria.png", plot = stack_plot_bacteria, width = 8, height = 6, dpi = 600, bg = "transparent")
#每个分组中最特殊asv的占每个分组中所有asv的比例
#识别出那些只在一个特定分组中出现过的 ASV（不论它们在那个分组的多少个样本中出现），然后统计每个分组中有多少个这样的“组特有”ASV。先计算（去除）样本数，后计算（保留）分组数。sum(n())和n()一样。
b_controlasv_end<- b_control_asv %>% group_by(Group,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)%>%
  ungroup()%>% group_by(Group) %>%summarize(number=n())
#统计每个分组中总共有多少种 ASV（不论这些 ASV 是否在其他分组中也出现过）。
b_controlasv_total<- b_control_asv %>% group_by(Group,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(Group)%>%summarize(sumii=sum(n()))
b_prop_unique_all_Group <- merge(b_controlasv_end, b_controlasv_total, by="Group")
b_prop_unique_all_Group$prop_unique_all <- b_prop_unique_all_Group$number/b_prop_unique_all_Group$sumii
b_prop_unique_all_Group
b_prop_unique_all_Group$Group <- factor(b_prop_unique_all_Group$Group, levels = niveaux)#画图顺序
b_prop_unique_all_Group <- b_prop_unique_all_Group[order(b_prop_unique_all_Group$Group), ]# 数据框顺序
b_prop_unique_all_Group
sprintf("%.1f%%", b_prop_unique_all_Group$prop_unique_all * 100)
#真菌堆叠图
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
# 合并为 phyloseq 对象
f_ASV <- otu_table(f_ASV, taxa_are_rows = TRUE)
f_tax <- tax_table(as.matrix(f_tax))
f_merged <- merge_phyloseq(f_ASV, f_tax, metadata)
#去除非真菌
f_merged <- subset_taxa(f_merged, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
f_merged <- subset_taxa(f_merged, (Kingdom!="d__Archaea") | is.na(Kingdom))
f_merged <- subset_taxa(f_merged, (Order!="o__Chloroplast") )
f_merged <- subset_taxa(f_merged, (Family!="f__Mitochondria"))
f_merged <- subset_taxa(f_merged, (Family!="NA"))
# 提取 ASV 表
f_ASV <- as.matrix(otu_table(f_merged, taxa_are_rows=T))
f_ASV_df <- f_ASV[rowSums(f_ASV[])>0,]
#核心asv（在大于2个居群生态位中相对丰度大于0.1%）
# 计算 ASV 相对丰度
f_ASV_rel <- sweep(f_ASV_df, 2, colSums(f_ASV_df), "/")
# 先把样本、niche 和 ASV 丰度组合成一个长表
f_rel_long <- melt(as.matrix(f_ASV_rel), varnames = c("ASV", "Sample"))
sample_niches <- metadata$Niche#样本到生态位（niche）的映射向量
names(sample_niches) <- metadata$Sample.ID
f_rel_long$Niche <- sample_niches[f_rel_long$Sample]
# 按 ASV 和 Niche 聚合：计算平均相对丰度
f_asv_niche_mean <- f_rel_long %>%
  group_by(ASV, Niche) %>%
  summarise(mean_abundance = mean(value), .groups = "drop")
# 对每个 ASV，看有多少个 niche平均相对丰度超过阈值
f_asv_niche_count <- f_asv_niche_mean %>%
  filter(mean_abundance >= 0.001) %>%
  group_by(ASV) %>%
  summarise(niche_count = n())
# 筛选核心 ASV：在 =2 个 niche 中都有平均丰度 >0.1%
f_core_asvs_names <- f_asv_niche_count$ASV[f_asv_niche_count$niche_count == 2]
# 提取核心 ASV 的相对丰度数据
f_core_asv_data <- f_ASV_rel[f_core_asvs_names, ]
# 输出核心 ASV 数量
cat("核心真菌 ASV 数量:", length(f_core_asvs_names), "\n")
#生成核心ASV的phyloseq
f_merge_core_abundance <- f_ASV_df[rownames(f_ASV_df) %in% rownames(f_core_asv_data),]
f_mca_table <- otu_table(f_merge_core_abundance, taxa_are_rows=T)
f_merged_NOMIS_core_ab<- merge_phyloseq(f_mca_table, f_tax, metadata) 
#最特殊asv（只在一个分组中出现）
metadata_nomis <- sample.data.frame(f_merged)#提取元数据
f_asv_df <- as.data.frame(t(otu_table(f_ASV_df, taxa_are_rows=T))) #提取转置asv数量表
f_asvdfmelt <- melt(as.matrix(f_asv_df))#长宽表格转换
f_asvdfmelt <- f_asvdfmelt[f_asvdfmelt$value >0,]
f_control_asv <- f_asvdfmelt
f_control_asv$Group <- vapply((f_control_asv$Var1), function(x) metadata_nomis$Group[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))# 添加分组信息
f_controlasv_endemic <- f_control_asv %>%
  distinct(Group, Var2) %>%         # 去重：只保留一个Group-ASV组合
  group_by(Var2) %>%
  tally(name = "group_count") %>%   # 统计每个 ASV 出现在多少个 Group 中
  filter(group_count == 1)          # 只保留那些只在一个 Group 中出现的
f_merge_controlasv_abundance <- f_ASV_df[rownames(f_ASV_df) %in% f_controlasv_endemic$Var2,]
f_controlasv_table <- otu_table(f_merge_controlasv_abundance, taxa_are_rows=T)
f_merged_NOMIS_controlasv_ab<- merge_phyloseq(f_controlasv_table, f_tax, metadata) 
#只出现在一个 Group 中所有样本的asv
#从原始数据中提取只出现在一个 Group 中 ASV 的样本信息
f_asv_group_info <- f_control_asv %>%
  filter(Var2 %in% f_controlasv_endemic$Var2) %>%
  group_by(Var2, Group) %>%
  summarise(sample_count = n_distinct(Var1), .groups = "drop") %>%#统计每个 ASV 在这个 group 中 有丰度的样本数量（非零）
  filter(sample_count == 6)  # 只保留在某个组中6个样本都出现的ASV
#提取这些 ASV 的丰度表
f_final_asv_ids <- f_asv_group_info$Var2
f_final_asv_abundance <- f_ASV_df[rownames(f_ASV_df) %in% f_final_asv_ids, ]
length(f_final_asv_ids)
#构建新的 phyloseq 对象（可选）
f_final_asv_table <- otu_table(as.matrix(f_final_asv_abundance), taxa_are_rows = TRUE)
f_final_ps <- merge_phyloseq(f_final_asv_table, f_tax, metadata)
#查看结果（按需提取特定分组的样本等）
View(as.data.frame(otu_table(f_final_ps)))
#比较特殊asv（只在一个生态位中出现）
f_specific_asv <- f_asvdfmelt
f_specific_asv$Niche <- vapply((f_specific_asv$Var1), function(x) metadata_nomis$Niche[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))# 添加生态位信息
f_specific_endemic <- f_specific_asv %>%
  distinct(Niche, Var2) %>%         # 去重：只保留一个Niche-ASV组合
  group_by(Var2) %>%
  tally(name = "niche_count") %>%   # 统计每个 ASV 出现在多少个 Niche 中
  filter(niche_count == 1)          # 只保留那些只在一个 Niche 中出现的
f_merge_specific_abundance <- f_ASV_df[rownames(f_ASV_df) %in% f_specific_endemic$Var2,]
f_specific_table <- otu_table(f_merge_specific_abundance, taxa_are_rows=T)
f_merged_NOMIS_specific_ab<- merge_phyloseq(f_specific_table, f_tax, metadata) 
# 指示性asv，运行 multipatt 进行指示性分析，找出哪些 ASV 是某些 niche 的指示物种
f_indicator_multipatt <- multipatt(t(f_ASV_df), metadata$Niche, 
                                   func = "r.g", control = how(nperm = 999))
summary(f_indicator_multipatt)
# 获取显著的指示 ASV（p < 0.01）
f_indicator_names <- f_indicator_multipatt$sign[which(f_indicator_multipatt$sign$p.value < 0.01), ]
print(f_indicator_names)
#生成指示性ASV的phyloseq
f_merge_indicator_abundance <- f_ASV_df[rownames(f_ASV_df) %in% rownames(f_indicator_names),]
f_indicator_table <- otu_table(f_merge_indicator_abundance, taxa_are_rows=T)
f_merged_NOMIS_indicator_ab<- merge_phyloseq(f_indicator_table, f_tax, metadata) 
#统一读取
f_dat	<- as.data.frame(otu_table(f_merged, taxa_are_rows=T))
f_cores <- as.data.frame(otu_table(f_merged_NOMIS_core_ab, taxa_are_rows=T)) 
f_endemic <- as.data.frame(otu_table(f_merged_NOMIS_specific_ab, taxa_are_rows=T))#比较特殊
f_indicator <- as.data.frame(otu_table(f_merged_NOMIS_indicator_ab, taxa_are_rows=T))
f_uniq <- as.data.frame(otu_table(f_merged_NOMIS_controlasv_ab, taxa_are_rows=T))#最特殊
f_ASV <- rownames(f_dat)
f_ASV <- as.data.frame(f_ASV)
colnames(f_ASV) <- "ASV"
f_cores$ASV<-rownames(f_cores)
rownames(f_cores)<-NULL
f_cores <- f_cores %>%
  select(ASV)%>%# 只保留 ASV 这一列
  mutate(type = "core")#增加一列 type，所有行都标记为 "core"
f_endemic$ASV<-rownames(f_endemic)
rownames(f_endemic)<-NULL
f_endemic <- f_endemic %>%
  select(ASV)%>%
  mutate(type = "endemic")
colnames(f_endemic) <- c("ASV", "type")
f_indicator$ASV<-rownames(f_indicator)
rownames(f_indicator)<-NULL
f_indicator <- f_indicator %>%
  select(ASV)%>%
  mutate(type = "indicator")
colnames(f_indicator) <- c("ASV","type")
f_uniq$ASV<-rownames(f_uniq)
rownames(f_uniq)<-NULL
f_uniq <- f_uniq %>%
  select(ASV)%>%
  mutate(type = "uniq")
colnames(f_uniq) <- c("ASV","type")
f_dat_numb <- rbind(f_cores, f_endemic, f_uniq, f_indicator)#合并
f_ASV_other <- f_ASV %>%#筛选出那些没有出现在四类标签里的ASV
  filter(!(ASV %in% f_dat_numb$ASV))%>%
  mutate(type = "Other")
#最特殊asv
length(unique(row.names(f_uniq)))/length(unique(row.names(f_ASV_df)))# 占 总 ASV 数量 的比例
f_merge_uniq <- merge(f_uniq, f_asvdfmelt, by.x="ASV",by.y="Var2")#合并reads 数量
f_merge_uniq$Group <- vapply((f_merge_uniq$Var1), function(x) metadata_nomis$Group[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))#加列
colnames(f_merge_uniq) <- c("ASV","type","SampleID","nf_count","Group")#asv名字、出现在所有样本中的样本次数、样本名、asv出现在某一样本的reads数量、组名
f_uniq_Group_prev <- f_merge_uniq %>%
  group_by(Group) %>%             # 按照 "Group" 列进行分组
  summarize(prev = n())           # 计算每个分组的行数（即该组中只在一个样本出现的ASV数量）
f_uniq_one_plot <- f_merge_uniq %>% group_by(Group) 
f_uniq_one_plot<- f_uniq_one_plot[c("Group","ASV","type")]
f_uniq_one_plot$Color <- "Specific and unique"
# 计算每个 Group 的最特殊 ASV 数量 (prev) 占所有最特殊 ASV 总数量的比例
f_prop_unique_Group<- f_uniq_Group_prev%>%
  group_by(Group)%>% 
  summarize(prop=prev/length(unique(row.names(f_uniq)))) 
colnames(f_prop_unique_Group) <- c("group_range","prop_unique")
#其他asv
length(unique(row.names(f_ASV_other)))/length(unique(row.names(f_ASV_df)))# 占 总 ASV 数量 的比例
f_merge_other <- merge(f_ASV_other, f_asvdfmelt, by.x="ASV",by.y="Var2")#合并reads 数量
f_merge_other$Group <- vapply((f_merge_other$Var1), function(x) metadata_nomis$Group[metadata_nomis$Sample.ID == x], FUN.VALUE = character(1))#加列
colnames(f_merge_other) <- c("ASV","type","SampleID","nf_count","Group")#asv名字、出现在所有样本中的样本次数、样本名、asv出现在某一样本的reads数量、组名
f_other_Group_prev <- f_merge_other %>%
  group_by(Group) %>%             # 按照 "Group" 列进行分组
  summarize(prev = n_distinct(ASV))           # 计算每个分组的其他asv种类
f_other_one_plot <- f_merge_other %>% group_by(Group) 
f_other_one_plot<- f_other_one_plot[c("Group","ASV","type")]
f_other_one_plot$Color <- "Other"
# 计算每个 Group 的其他 ASV 数量 (prev) 占所有其他ASV 总数量的比例
f_prop_other_Group<- f_other_Group_prev%>%
  group_by(Group)%>% 
  summarize(prop=prev/length(unique(row.names(f_ASV_other)))) 
colnames(f_prop_other_Group) <- c("group_range","prop_other")
#合并最特殊和其他asv
f_df_full <- rbind(f_uniq_one_plot, f_other_one_plot)
# 使用 distinct() 保留 Group 和 ASV 组合唯一的行，.keep_all = TRUE 会保留该唯一行对应的原始数据框中的所有其他列
f_df_full_clean <- f_df_full %>%
  distinct(Group, ASV, .keep_all = TRUE)
niveaux <- c("JRG", "JJG", "TZG", "PAG","JRN", "JJN", "TZN", "PAN")
f_df_full_clean$Group<- factor(f_df_full_clean$Group, levels = niveaux) 
f_df_full_clean$Group <- fct_rev(f_df_full_clean$Group) # 将刚才设置的 Group 因子水平顺序反转
# 定义堆叠图中 Color 列的因子水平顺序
stacking_levels_order <- c("Other", "Specific and unique")
f_df_full_clean$Color <- factor(f_df_full_clean$Color, levels = stacking_levels_order) 
# 定义图例中 Color 列的因子水平顺序
legend_key_order_labels <- c("Specific and unique", "Other") 
plot_colors <- c("#606060", "#7FA3D6") # 假设 Other 用灰色, Specific 用蓝色
names(plot_colors) <- stacking_levels_order # 为颜色命名以匹配因子水平
# 生成柱状图，使用设定的堆叠顺序和图例顺序
stack_plot_fungi <- ggplot(f_df_full_clean, aes(fill = Color, y = Group)) +
  geom_bar(position = "stack", stat = "count") + # stat="count" 根据行数堆叠
  labs(
    x = "Number of specific fungal ASVs", # 设置横坐标标题
    y = NULL                       # 移除纵坐标标题 (或者 y = "")
  ) +
  theme_classic() +
  # 使用 scale_fill_manual 手动设定填充颜色和图例项目顺序
  scale_fill_manual(
    values = plot_colors, #指定堆叠图中显示的分类顺序
    breaks = legend_key_order_labels # 指定图例中显示的分类标签及其顺序
  ) +
  scale_x_continuous(
    limits = c(0, 3500), # 保持轴的范围到 4000
    breaks = seq(0, 3000, by = 1000), # 保持刻度和标签只到 3000
    expand = expansion(mult = c(0, 0)) # 新增或修改 expand 参数，左侧扩展为 0
  ) +
  theme(
    # 移除纵坐标轴刻度标签和刻度线
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(), # 确保 Y 轴线仍然存在 (theme_classic 默认包含)
    # 设置横坐标轴标题和刻度标签的大小
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14), 
    axis.ticks.x = element_line(), # 确保刻度线本身是可见的 (theme_classic 默认包含)
    axis.ticks.length.x = unit(-0.25, "cm"), # 设置 x 轴刻度线朝内，并指定长度
    # 移除图例标题
    legend.title = element_blank(),
    # 设置图例文字大小
    legend.text = element_text(size = 14),
    # 设置图例位置在绘图区域内部的右下角
    legend.position = c(1, 0.1),
    legend.justification = c(1, 0) # 设置图例框的右下角对齐到坐标(1, 0)
  )  
stack_plot_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3b_stack_fungi.png", plot = stack_plot_fungi, width = 8, height = 6, dpi = 600, bg = "transparent")
#每个分组中最特殊asv的占每个分组中所有asv的比例
#识别出那些只在一个特定分组中出现过的 ASV（不论它们在那个分组的多少个样本中出现），然后统计每个分组中有多少个这样的“组特有”ASV。先计算（去除）样本数，后计算（保留）分组数。sum(n())和n()一样。
f_controlasv_end<- f_control_asv %>% group_by(Group,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)%>%
  ungroup()%>% group_by(Group) %>%summarize(number=n())
#统计每个分组中总共有多少种 ASV（不论这些 ASV 是否在其他分组中也出现过）。
f_controlasv_total<- f_control_asv %>% group_by(Group,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(Group)%>%summarize(sumii=sum(n()))
f_prop_unique_all_Group <- merge(f_controlasv_end, f_controlasv_total, by="Group")
f_prop_unique_all_Group$prop_unique_all <- f_prop_unique_all_Group$number/f_prop_unique_all_Group$sumii
f_prop_unique_all_Group
f_prop_unique_all_Group$Group <- factor(f_prop_unique_all_Group$Group, levels = niveaux)#画图顺序
f_prop_unique_all_Group <- f_prop_unique_all_Group[order(f_prop_unique_all_Group$Group), ]# 数据框顺序
f_prop_unique_all_Group
sprintf("%.1f%%", f_prop_unique_all_Group$prop_unique_all * 100)

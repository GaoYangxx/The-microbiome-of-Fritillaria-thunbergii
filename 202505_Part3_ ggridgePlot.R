#图3d
install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "phyloseqCompanion"))
install.packages(c("tidyverse", "ggridges"))
install.packages("adiv")
install.packages("vegan")
install.packages("indicspecies")
install.packages("RColorBrewer") 
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(RColorBrewer) 
library(indicspecies)
library(speedyseq)
library(phyloseq)
library(phyloseqCompanion)
library(tidyverse)
library(ggridges)
library(adiv)
library(vegan)
#细菌山峦图
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
b_indicator_table_filtered <- subset(b_indicator_table, !(row.names(b_indicator_table) %in% c(row.names(b_endemic), row.names(b_cores))))# 过滤
b_merge_indicator_phylo <- merge_phyloseq(b_indicator_table_filtered, b_tax, metadata)
b_indicator_ASV <- row.names(b_indicator_table_filtered)
b_Specific_ASV <- row.names(b_endemic)
b_core_ASV <- row.names(b_cores)
b_df_core_Specific_indicator <- c(b_indicator_ASV, b_Specific_ASV, b_core_ASV)#合并
b_df_core_Specific_indicator_unique <- unique(b_df_core_Specific_indicator)# 去重
b_nomis_asv_count_unique <- b_ASV_df[rownames(b_ASV_df) %in% b_df_core_Specific_indicator_unique,]# ASV 计数
b_asv_nomis_unique <- otu_table(b_nomis_asv_count_unique, taxa_are_rows=T)
b_nomis_asv_unique_phylo <- merge_phyloseq(b_asv_nomis_unique, b_tax, metadata)
#将ASV计数聚合到属（Genus）级别
b_unique_Genus_taxglom <- tax_glom(b_nomis_asv_unique_phylo, taxrank=rank_names(b_nomis_asv_unique_phylo)[6], NArm=F)
#将属级别的计数转换成相对丰度
b_transf_unique = transform_sample_counts(b_unique_Genus_taxglom, function(x) x / sum(x))
#门相对丰度
b_unique_phylum_taxglom <- tax_glom(b_nomis_asv_unique_phylo, taxrank=rank_names(b_nomis_asv_unique_phylo)[2], NArm=F)
b_transf_unique_phylum = transform_sample_counts(b_unique_phylum_taxglom, function(x) x / sum(x))
#科相对丰度
b_unique_family_taxglom <- tax_glom(b_nomis_asv_unique_phylo, taxrank=rank_names(b_nomis_asv_unique_phylo)[5], NArm=F)
b_transf_family = transform_sample_counts(b_unique_family_taxglom, function(x) x / sum(x))
#获取总相对丰度最高的14个Family的ID
b_TopASV_f <- names(sort(taxa_sums(b_transf_family), TRUE)[1:14])
#从属级别相对丰度表中，保留属于这前14个科的属
b_top25_NOMIS_f <- prune_taxa(b_TopASV_f, b_transf_unique)
b_top25_NOMIS_f <- prune_taxa(taxa_sums(b_top25_NOMIS_f)>0, b_top25_NOMIS_f)
b_top_Genus<-as.data.frame(tax_table(b_top25_NOMIS_f))# 提取最终选定属的分类信息表
#提取属、科、门级别相对丰度数据（OTU表）和分类信息
b_asv_table_all_genus <- otu_table(b_transf_unique, taxa_are_rows=T)
b_tax_table_all_genus <- tax_table(b_transf_unique)
b_asv_table_all_family <- otu_table(b_transf_family, taxa_are_rows=T)
b_tax_table_all_family <- tax_table(b_transf_family)
b_asv_table_all_phylum <- otu_table(b_transf_unique_phylum, taxa_are_rows=T)
b_tax_table_all_phylum <- tax_table(b_transf_unique_phylum)
#核心asv
sample_data(b_merged_NOMIS_core_ab)$Group <- as.factor(sample_data(b_merged_NOMIS_core_ab)$Group)
sample_data(b_merged_NOMIS_core_ab)$Sample.ID <- as.factor(sample_data(b_merged_NOMIS_core_ab)$Sample.ID)
sample_data(b_merged_NOMIS_core_ab)$Origin <- as.factor(sample_data(b_merged_NOMIS_core_ab)$Origin)
sample_data(b_merged_NOMIS_core_ab)$Niche <- as.factor(sample_data(b_merged_NOMIS_core_ab)$Niche)
#计算核心物种在每个样本内的相对丰度
b_core_RA = transform_sample_counts(b_merged_NOMIS_core_ab, function(x) x / sum(x))
#按"Group"列合并样本
b_core_RA = merge_samples(b_core_RA, "Group")
#重新计算合并后各组的相对丰度
b_core_RA = transform_sample_counts(b_core_RA, function(x) x / sum(x))
#将phyloseq对象转换成长格式数据框
b_data_core <- psmelt(b_core_RA) 
b_data_core$Genus <-as.character(b_data_core$Genus)
#筛选核心asv中相对丰度最高的10个属，并计算其总丰度
b_sumtot_core <- b_data_core %>%
  group_by(Genus) %>% 
  summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% b_top_Genus$Genus) %>%
  filter(!(Genus %in% c("", "g__uncultured"))) %>%
  filter(!grepl("g__unclassified|uncultured|unidentified|metagenome", Genus))
#将不在筛选列表中的核心asv属名替换为"Other"
b_data_core$Genus[!(b_data_core$Genus %in% b_sumtot_core$Genus)] <- "Other"
#添加一个新列"core"，并赋值为"core"
b_data_core$core <- "core"
#计算b_data_core数据框中所有相对丰度的总和
b_data_core$totalAbundance <- sum(b_data_core$Abundance)
#按属和"core"列分组，计算每个组的总相对丰度（除以总Abundance），并去重
b_data_core_mod <- b_data_core%>%
  group_by(Genus, core)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()
#选择颜色
n <- 20
#获取RColorBrewer包中所有定性（qualitative）调色板的信息
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
#生成并合并所有定性调色板中的颜色到一个长向量（虽然这个向量在后续没有被使用）
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#画堆叠柱状图
b_barplot_biogeo_core <- ggplot(data=b_data_core_mod, aes(x=core, y=abundance, fill=Genus))
b_barplot_biogeo_core <- b_barplot_biogeo_core + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
b_barplot_biogeo_core<- b_barplot_biogeo_core+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
#比较特殊asv（框架同上）
sample_data(b_merged_NOMIS_specific_ab)$Group <- as.factor(sample_data(b_merged_NOMIS_specific_ab)$Group)
sample_data(b_merged_NOMIS_specific_ab)$Sample.ID <- as.factor(sample_data(b_merged_NOMIS_specific_ab)$Sample.ID)
sample_data(b_merged_NOMIS_specific_ab)$Origin <- as.factor(sample_data(b_merged_NOMIS_specific_ab)$Origin)
sample_data(b_merged_NOMIS_specific_ab)$Niche <- as.factor(sample_data(b_merged_NOMIS_specific_ab)$Niche)
b_Specific_RA = transform_sample_counts(b_merged_NOMIS_specific_ab, function(x) x / sum(x))
b_Specific_RA = merge_samples(b_Specific_RA, "Group")
b_Specific_RA = transform_sample_counts(b_Specific_RA, function(x) x / sum(x))
b_data_Specific <- psmelt(b_Specific_RA) 
b_data_Specific$Genus <-as.character(b_data_Specific$Genus) 
b_sumtot_Specific <-
  b_data_Specific %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% b_top_Genus$Genus) %>% filter(!(Genus %in% c(""," g__uncultured"))) %>%
  filter(!grepl("g__unclassified|uncultured|unidentified|metagenome", Genus)) 
b_data_Specific$Genus[!(b_data_Specific$Genus %in% b_sumtot_Specific$Genus)] <- "Other"
b_data_Specific$Specific <- "Specific"
b_data_Specific$totalAbundance <- sum(b_data_Specific$Abundance)
b_data_Specific_mod <- b_data_Specific%>%
  group_by(Genus, Specific)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
b_barplot_biogeo_Specific <- ggplot(data= b_data_Specific_mod, aes(x=Specific, y=abundance, fill=Genus))
b_barplot_biogeo_Specific <- b_barplot_biogeo_Specific + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
b_barplot_biogeo_Specific<- b_barplot_biogeo_Specific + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
#指示性asv（框架同上）
sample_data(b_merged_NOMIS_indicator_ab)$Group <- as.factor(sample_data(b_merged_NOMIS_indicator_ab)$Group)
sample_data(b_merged_NOMIS_indicator_ab)$Sample.ID <- as.factor(sample_data(b_merged_NOMIS_indicator_ab)$Sample.ID)
sample_data(b_merged_NOMIS_indicator_ab)$Origin <- as.factor(sample_data(b_merged_NOMIS_indicator_ab)$Origin)
sample_data(b_merged_NOMIS_indicator_ab)$Niche <- as.factor(sample_data(b_merged_NOMIS_indicator_ab)$Niche)
b_indicator_RA = transform_sample_counts(b_merged_NOMIS_indicator_ab, function(x) x / sum(x))
b_indicator_RA = merge_samples(b_indicator_RA, "Group")
b_indicator_RA = transform_sample_counts(b_indicator_RA, function(x) x / sum(x))
b_data_indicator <- psmelt(b_indicator_RA) 
b_data_indicator$Genus <-as.character(b_data_indicator$Genus)
b_sumtot_indicator <-
  b_data_indicator %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% b_top_Genus$Genus) %>%filter(!(Genus %in% c(""," g__uncultured")))%>%
  filter(!grepl("g__unclassified|uncultured|unidentified|metagenome", Genus)) 
b_data_indicator$Genus[!(b_data_indicator$Genus %in% b_sumtot_indicator$Genus)] <- "Other"
b_data_indicator$indicator <- "indicator"
b_data_indicator$totalAbundance <- sum(b_data_indicator$Abundance)
b_data_indicator_mod <- b_data_indicator%>%
  group_by(Genus, indicator)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
b_barplot_biogeo_indicator <- ggplot(data= b_data_indicator_mod, aes(x=indicator, y=abundance, fill=Genus))
b_barplot_biogeo_indicator <- b_barplot_biogeo_indicator + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
b_barplot_biogeo_indicator <- b_barplot_biogeo_indicator + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
#合并之前绘制的三个柱状图
ggarrange(b_barplot_biogeo_Specific, b_barplot_biogeo_core, b_barplot_biogeo_indicator, ncol = 3, nrow = 1, common.legend=T)
#过滤并准备比较特殊、核心、指示性物种的数据用于绘制山峦图
b_dataset_Specific_filter <- as.data.frame(b_data_Specific[b_data_Specific$Abundance >0,])
b_rename_Specific <- rename(b_dataset_Specific_filter, category = Specific)# 重命名
b_dataset_core_filter <- as.data.frame(b_data_core[b_data_core$Abundance >0,])
b_rename_core <- rename(b_dataset_core_filter, category = core)
b_dataset_indicator_filter <- as.data.frame(b_data_indicator[b_data_indicator$Abundance >0,])
b_rename_indicator <- rename(b_dataset_indicator_filter, category = indicator)
b_binddataset <- rbind(b_rename_Specific, b_rename_core, b_rename_indicator)
set.seed(3467)
#将 category 列转换为因子，并指定水平的顺序
b_binddataset$category <- factor(b_binddataset$category, levels = c("core", "Specific", "indicator"))
#获取数据中所有独特的原始属名
b_original_genus_levels <- unique(b_binddataset$Genus)
#创建一个命名向量，定义原始属名到显示标签的映射
b_genus_labels_map <- b_original_genus_levels # 复制原始属名作为起始
names(b_genus_labels_map) <- b_original_genus_levels # 将向量命名为其原始值
# 应用通用规则：移除 "g__" 前缀
b_genus_labels_map <- gsub("g__", "", b_genus_labels_map)
# 应用特定规则：将特定长名称替换为 "Rhizobium group"
b_long_original_name <- "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
if (b_long_original_name %in% names(b_genus_labels_map)) { # 检查该名称是否存在于数据中
  b_genus_labels_map[b_long_original_name] <- "Rhizobium group"
}
#在绘制图形前，将丰度转换为百分比
b_binddataset <- b_binddataset %>%
  mutate(Abundance_pct = Abundance * 100) 
# 定义 X 轴的刻度位置 (百分比值)
b_percentage_breaks <- c(0.001, 0.1, 10) # 定义常见的百分比刻度点
# 定义 X 轴刻度位置对应的标签
b_percentage_labels <- c(
  expression(10^-3), # 0.001 处显示 10^-3 的数学形式
  "0.1",           # 0.1 处显示文本 "0.1"
  "10"            # 10 处显示文本 "10"
)
#绘制山峦图，在 scale_y_discrete 中使用命名向量作为标签，X轴使用百分比丰度，自定义X轴刻度和标签
ggridge_bacteria <- ggplot(b_binddataset, aes(x = Abundance_pct, y = fct_reorder(Genus, Abundance, .desc = F), fill = category)) + # 将X轴映射改为 Abundance_pct
  geom_density_ridges(scale = 1, alpha = 0.8) + # 添加山峦图（密度曲线）图层
  # scale_fill_cyclical(values = c("blue", "green","red"))+
  theme_bw() + # 使用黑白主题作为基础
  scale_x_continuous(
    trans = "log10", # X轴对数转换
    breaks = b_percentage_breaks, # 设置刻度的位置
    labels = b_percentage_labels # 设置刻度位置对应的标签
  ) +
  scale_fill_brewer(
    palette = "Dark2", # 填充颜色使用 Dark2 调色板
    name = NULL, # 移除图例标题
    labels = c("Core", "Specific",
               "Indicator") # 设置图例项的显示文本
  ) +
  # 修改 Y 轴比例尺，使用命名向量作为标签 (保持不变)
  scale_y_discrete(labels = b_genus_labels_map) + # 使用创建好的命名向量
  labs(x = "Bacterial relative abundance (%)", # 设置横坐标标题 (保持不变，标题已经包含%)
       y = NULL) + # 将纵坐标标题设置为空 (保持不变)
  theme(
    panel.border = element_blank(), # 移除面板周围的边框 (保持不变)
    panel.grid.minor = element_blank(), # 移除次要网格线 (保持不变)
    #panel.grid.major = element_blank(), # 移除主要网格线 (如果您想保留主要网格线，请删除或注释掉这行)
    # 添加移除刻度线
    axis.ticks.x = element_blank(), # 移除横坐标刻度线
    axis.ticks.y = element_blank(), # 移除纵坐标刻度线
    axis.title.x = element_text(size = 16),# 修改横坐标标题大小
    #修改 Y 轴文本外观和对齐
    axis.text.y = element_text(face = "italic", # 设置字体为斜体
                               vjust = -0.5, # 设置垂直对齐方式为居中
                               hjust = 0, # 设置为左对齐
                               margin = margin(r = unit(-20, "pt")),# 设置右侧边距为0点，使其紧贴Y轴
                               size = 14), 
    # 修改 X 轴文本外观和对齐
    axis.text.x = element_text(margin = margin(t = unit(-20, "pt")),#设置顶部边距为0点，使其紧贴X轴
                               size = 14,
                               vjust = 0 ), # 设置为底端对齐 
    legend.text = element_text(size = 14), # 设置图例文本大小
    legend.spacing.y = unit(0.5, "cm") # 设置图例项之间的垂直间距，例如 0.5 厘米
  )
ggridge_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3d_ggridge_bacteria.png", plot = ggridge_bacteria, width = 6, height = 6, dpi = 600, bg = "transparent")
#找到相对丰度最高的门科属
bacteria_ASV <- read.csv(file = "D:/study/master/meiji/bacteria_ASV.csv",
                         sep=",",header=TRUE,check.names = FALSE)
# 定义文件的绝对路径（使用正斜杠）
metadata <- read_tsv("D:/study/master/metadata.tsv")
colnames(metadata)[1]<-"SampleID"
# 查看元数据
head(metadata)
# 读取 ASV 表（特征表）
b_ASV <- bacteria_ASV[, c("ASV",metadata$SampleID)]
b_ASV <- b_ASV %>% column_to_rownames(var = "ASV")
head(b_ASV)
# 读取 taxonomy（分类信息）
b_tax <- bacteria_ASV[, 2:9]
b_tax <- b_tax %>% column_to_rownames(var = "ASV")
head(b_tax)
# 提取 ASV 表
b_ASV_df <- b_ASV[rowSums(b_ASV[])>0,]
b_ASV_df$asv<-rownames(b_ASV_df)
b_tax$asv<-rownames(b_tax)
#在每个 Genus 分组，将同一个属下的所有 ASV 的计数在每个样本中加起来，得到每个属在每个样本的总计数。
b_dat_m <- b_ASV_df %>%
  left_join(b_tax)%>%#合并
  group_by(Genus)%>%#按属分组
  summarise_if(is.numeric, sum)%>%# 找到所有数值型的列，并对这些列计算它们的 sum
  na.omit()%>%
  pivot_longer(cols = !Genus)%>%# 从“宽”格式转换为“长”格式
  left_join(metadata, by = c("name" = "SampleID"))#合并
b_sum_all <- sum(b_dat_m$value)
b_num_all <- length(unique(b_dat_m$name))
#计算属相对丰度（单个属的菌/总菌）
b_dat_abun <- b_dat_m %>%
  group_by(Genus)%>%
  summarise(Abundance = sum(value)/b_sum_all)
colnames(b_dat_abun) <- c("Genus", "Abundance")
b_tax_sel <- b_tax%>%
  select(-c(asv,Species))%>%
  distinct()%>%#移除重复的行
  filter(Genus %in% b_dat_abun$Genus)
#结合属、相对丰度、分类信息
b_dat_final <- b_dat_abun%>%
  left_join(b_tax_sel)
#计算每个门的总丰度，并选出总丰度最高的门。
b_phylum_top <- b_dat_final %>%
  group_by(Phylum)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))
#找出相对丰度最高的科
b_family_top <- b_dat_final %>%
  group_by(Family)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))
#找出相对丰度最高的属
b_dat_final %>%
  dplyr::arrange(desc(Abundance)) -> b_genus_top
View(b_phylum_top)
View(b_family_top)
View(b_genus_top)
#真菌山峦图
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
f_indicator_table_filtered <- subset(f_indicator_table, !(row.names(f_indicator_table) %in% c(row.names(f_endemic), row.names(f_cores))))# 过滤
f_merge_indicator_phylo <- merge_phyloseq(f_indicator_table_filtered, f_tax, metadata)
f_indicator_ASV <- row.names(f_indicator_table_filtered)
f_Specific_ASV <- row.names(f_endemic)
f_core_ASV <- row.names(f_cores)
f_df_core_Specific_indicator <- c(f_indicator_ASV, f_Specific_ASV, f_core_ASV)#合并
f_df_core_Specific_indicator_unique <- unique(f_df_core_Specific_indicator)# 去重
f_nomis_asv_count_unique <- f_ASV_df[rownames(f_ASV_df) %in% f_df_core_Specific_indicator_unique,]# ASV 计数
f_asv_nomis_unique <- otu_table(f_nomis_asv_count_unique, taxa_are_rows=T)
f_nomis_asv_unique_phylo <- merge_phyloseq(f_asv_nomis_unique, f_tax, metadata)
#将ASV计数聚合到属（Genus）级别
f_unique_Genus_taxglom <- tax_glom(f_nomis_asv_unique_phylo, taxrank=rank_names(f_nomis_asv_unique_phylo)[6], NArm=F)
#将属级别的计数转换成相对丰度
f_transf_unique = transform_sample_counts(f_unique_Genus_taxglom, function(x) x / sum(x))
#门相对丰度
f_unique_phylum_taxglom <- tax_glom(f_nomis_asv_unique_phylo, taxrank=rank_names(f_nomis_asv_unique_phylo)[2], NArm=F)
f_transf_unique_phylum = transform_sample_counts(f_unique_phylum_taxglom, function(x) x / sum(x))
#科相对丰度
f_unique_family_taxglom <- tax_glom(f_nomis_asv_unique_phylo, taxrank=rank_names(f_nomis_asv_unique_phylo)[5], NArm=F)
f_transf_family = transform_sample_counts(f_unique_family_taxglom, function(x) x / sum(x))
#获取总相对丰度最高的14个Family的ID
f_TopASV_f <- names(sort(taxa_sums(f_transf_family), TRUE)[1:13])
#从属级别相对丰度表中，保留属于这前14个科的属
f_top25_NOMIS_f <- prune_taxa(f_TopASV_f, f_transf_unique)
f_top25_NOMIS_f <- prune_taxa(taxa_sums(f_top25_NOMIS_f)>0, f_top25_NOMIS_f)
f_top_Genus<-as.data.frame(tax_table(f_top25_NOMIS_f))# 提取最终选定属的分类信息表
#提取属、科、门级别相对丰度数据（OTU表）和分类信息
f_asv_table_all_genus <- otu_table(f_transf_unique, taxa_are_rows=T)
f_tax_table_all_genus <- tax_table(f_transf_unique)
f_asv_table_all_family <- otu_table(f_transf_family, taxa_are_rows=T)
f_tax_table_all_family <- tax_table(f_transf_family)
f_asv_table_all_phylum <- otu_table(f_transf_unique_phylum, taxa_are_rows=T)
f_tax_table_all_phylum <- tax_table(f_transf_unique_phylum)
#核心asv
sample_data(f_merged_NOMIS_core_ab)$Group <- as.factor(sample_data(f_merged_NOMIS_core_ab)$Group)
sample_data(f_merged_NOMIS_core_ab)$Sample.ID <- as.factor(sample_data(f_merged_NOMIS_core_ab)$Sample.ID)
sample_data(f_merged_NOMIS_core_ab)$Origin <- as.factor(sample_data(f_merged_NOMIS_core_ab)$Origin)
sample_data(f_merged_NOMIS_core_ab)$Niche <- as.factor(sample_data(f_merged_NOMIS_core_ab)$Niche)
#计算核心物种在每个分组内的相对丰度
f_core_RA = merge_samples(f_merged_NOMIS_core_ab, "Group")
#计算合并后各组的相对丰度
f_core_RA = transform_sample_counts(f_core_RA, function(x) x / sum(x))
#将phyloseq对象转换成长格式数据框
f_data_core <- psmelt(f_core_RA) 
f_data_core$Genus <-as.character(f_data_core$Genus)
#筛选核心asv中相对丰度最高的10个属，并计算其总丰度
f_sumtot_core <- f_data_core %>%
  group_by(Genus) %>% 
  summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% f_top_Genus$Genus) %>%
  filter(!(Genus %in% c("", "g__uncultured"))) %>%
  filter(!grepl("g__unclassified|uncultured|unidentified|metagenome", Genus))
#将不在筛选列表中的核心asv属名替换为"Other"
f_data_core$Genus[!(f_data_core$Genus %in% f_sumtot_core$Genus)] <- "Other"
#添加一个新列"core"，并赋值为"core"
f_data_core$core <- "core"
#计算f_data_core数据框中所有相对丰度的总和
f_data_core$totalAbundance <- sum(f_data_core$Abundance)
#按属和"core"列分组，计算每个组的总相对丰度（除以总Abundance），并去重
f_data_core_mod <- f_data_core%>%
  group_by(Genus, core)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()
#选择颜色
n <- 20
#获取RColorBrewer包中所有定性（qualitative）调色板的信息
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
#生成并合并所有定性调色板中的颜色到一个长向量（虽然这个向量在后续没有被使用）
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#画堆叠柱状图
f_barplot_biogeo_core <- ggplot(data=f_data_core_mod, aes(x=core, y=abundance, fill=Genus))
f_barplot_biogeo_core <- f_barplot_biogeo_core + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
f_barplot_biogeo_core<- f_barplot_biogeo_core+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
#比较特殊asv（框架同上）
sample_data(f_merged_NOMIS_specific_ab)$Group <- as.factor(sample_data(f_merged_NOMIS_specific_ab)$Group)
sample_data(f_merged_NOMIS_specific_ab)$Sample.ID <- as.factor(sample_data(f_merged_NOMIS_specific_ab)$Sample.ID)
sample_data(f_merged_NOMIS_specific_ab)$Origin <- as.factor(sample_data(f_merged_NOMIS_specific_ab)$Origin)
sample_data(f_merged_NOMIS_specific_ab)$Niche <- as.factor(sample_data(f_merged_NOMIS_specific_ab)$Niche)
f_Specific_RA = transform_sample_counts(f_merged_NOMIS_specific_ab, function(x) x / sum(x))
f_Specific_RA = merge_samples(f_Specific_RA, "Group")
f_Specific_RA = transform_sample_counts(f_Specific_RA, function(x) x / sum(x))
f_data_Specific <- psmelt(f_Specific_RA) 
f_data_Specific$Genus <-as.character(f_data_Specific$Genus) 
f_sumtot_Specific <-
  f_data_Specific %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% f_top_Genus$Genus) %>% filter(!(Genus %in% c(""," g__uncultured"))) %>%
  filter(!grepl("g__unclassified|uncultured|unidentified|metagenome", Genus)) 
f_data_Specific$Genus[!(f_data_Specific$Genus %in% f_sumtot_Specific$Genus)] <- "Other"
f_data_Specific$Specific <- "Specific"
f_data_Specific$totalAbundance <- sum(f_data_Specific$Abundance)
f_data_Specific_mod <- f_data_Specific%>%
  group_by(Genus, Specific)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
f_barplot_biogeo_Specific <- ggplot(data= f_data_Specific_mod, aes(x=Specific, y=abundance, fill=Genus))
f_barplot_biogeo_Specific <- f_barplot_biogeo_Specific + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
f_barplot_biogeo_Specific<- f_barplot_biogeo_Specific + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
#指示性asv（框架同上）
sample_data(f_merged_NOMIS_indicator_ab)$Group <- as.factor(sample_data(f_merged_NOMIS_indicator_ab)$Group)
sample_data(f_merged_NOMIS_indicator_ab)$Sample.ID <- as.factor(sample_data(f_merged_NOMIS_indicator_ab)$Sample.ID)
sample_data(f_merged_NOMIS_indicator_ab)$Origin <- as.factor(sample_data(f_merged_NOMIS_indicator_ab)$Origin)
sample_data(f_merged_NOMIS_indicator_ab)$Niche <- as.factor(sample_data(f_merged_NOMIS_indicator_ab)$Niche)
f_indicator_RA = transform_sample_counts(f_merged_NOMIS_indicator_ab, function(x) x / sum(x))
f_indicator_RA = merge_samples(f_indicator_RA, "Group")
f_indicator_RA = transform_sample_counts(f_indicator_RA, function(x) x / sum(x))
f_data_indicator <- psmelt(f_indicator_RA) 
f_data_indicator$Genus <-as.character(f_data_indicator$Genus)
f_sumtot_indicator <-
  f_data_indicator %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% f_top_Genus$Genus) %>%filter(!(Genus %in% c(""," g__uncultured")))%>%
  filter(!grepl("g__unclassified|uncultured|unidentified|metagenome", Genus)) 
f_data_indicator$Genus[!(f_data_indicator$Genus %in% f_sumtot_indicator$Genus)] <- "Other"
f_data_indicator$indicator <- "indicator"
f_data_indicator$totalAbundance <- sum(f_data_indicator$Abundance)
f_data_indicator_mod <- f_data_indicator%>%
  group_by(Genus, indicator)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
f_barplot_biogeo_indicator <- ggplot(data= f_data_indicator_mod, aes(x=indicator, y=abundance, fill=Genus))
f_barplot_biogeo_indicator <- f_barplot_biogeo_indicator + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
f_barplot_biogeo_indicator <- f_barplot_biogeo_indicator + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
#合并之前绘制的三个柱状图
ggarrange(f_barplot_biogeo_Specific, f_barplot_biogeo_core, f_barplot_biogeo_indicator, ncol = 3, nrow = 1, common.legend=T)
#过滤并准备比较特殊、核心、指示性物种的数据用于绘制山峦图
f_dataset_Specific_filter <- as.data.frame(f_data_Specific[f_data_Specific$Abundance >0,])
f_rename_Specific <- rename(f_dataset_Specific_filter, category = Specific)# 重命名
f_dataset_core_filter <- as.data.frame(f_data_core[f_data_core$Abundance >0,])
f_rename_core <- rename(f_dataset_core_filter, category = core)
f_dataset_indicator_filter <- as.data.frame(f_data_indicator[f_data_indicator$Abundance >0,])
f_rename_indicator <- rename(f_dataset_indicator_filter, category = indicator)
f_binddataset <- rbind(f_rename_Specific, f_rename_core, f_rename_indicator)
set.seed(3467)
#将 category 列转换为因子，并指定水平的顺序
f_binddataset$category <- factor(f_binddataset$category, levels = c("core", "Specific", "indicator"))
#获取数据中所有独特的原始属名
f_original_genus_levels <- unique(f_binddataset$Genus)
#创建一个命名向量，定义原始属名到显示标签的映射
f_genus_labels_map <- f_original_genus_levels # 复制原始属名作为起始
names(f_genus_labels_map) <- f_original_genus_levels # 将向量命名为其原始值
# 应用通用规则：移除 "g__" 前缀
f_genus_labels_map <- gsub("g__", "", f_genus_labels_map)
# 应用特定规则：将特定长名称替换为 "Rhizobium group"
f_long_original_name <- "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
if (f_long_original_name %in% names(f_genus_labels_map)) { # 检查该名称是否存在于数据中
  f_genus_labels_map[f_long_original_name] <- "Rhizobium group"
}
#在绘制图形前，将丰度转换为百分比
f_binddataset <- f_binddataset %>%
  mutate(Abundance_pct = Abundance * 100) 
# 定义 X 轴的刻度位置 (百分比值)
f_percentage_breaks <- c(0.001, 0.1, 10) # 定义常见的百分比刻度点
# 定义 X 轴刻度位置对应的标签
f_percentage_labels <- c(
  expression(10^-3), # 0.001 处显示 10^-3 的数学形式
  "0.1",           # 0.1 处显示文本 "0.1"
  "10"            # 10 处显示文本 "10"
)
#使用 Abundance 列计算每个属（包括 Other）的丰度中位数，用于后续排序
f_genus_medians_for_sorting <- f_binddataset %>%
  group_by(Genus) %>%
  summarise(median_abundance = median(Abundance)) %>% # 计算中位数
  ungroup()
# 确定最终的 Y 轴属显示顺序，将 "Other" 放在最后
ordered_genera_excluding_other <- f_genus_medians_for_sorting %>%
  filter(Genus != "Other") %>% # 排除 "Other" 类别
  arrange(median_abundance) %>% # 按中位数升序排序
  pull(Genus) # 提取排序后的属名向量
# 创建最终的顺序：排序好的属 + "Other"
final_genus_order <- c("Other", ordered_genera_excluding_other)
# 将 f_binddataset 中的 Genus 列转换为因子，并指定其水平顺序
f_binddataset$Genus <- factor(f_binddataset$Genus, levels = final_genus_order)
#绘制山峦图，在 scale_y_discrete 中使用命名向量作为标签，X轴使用百分比丰度，自定义X轴刻度和标签
ggridge_fungi <- ggplot(f_binddataset, aes(x = Abundance_pct, y = Genus, fill = category)) + # 将X轴映射改为 Abundance_pct
  geom_density_ridges(scale = 1, alpha = 0.8) + # 添加山峦图（密度曲线）图层
  # scale_fill_cyclical(values = c("blue", "green","red"))+
  theme_bw() + # 使用黑白主题作为基础
  scale_x_continuous(
    trans = "log10", # X轴对数转换
    breaks = f_percentage_breaks, # 设置刻度的位置
    labels = f_percentage_labels # 设置刻度位置对应的标签
  ) +
  scale_fill_brewer(
    palette = "Dark2", # 填充颜色使用 Dark2 调色板
    name = NULL, # 移除图例标题
    labels = c("Core", "Specific",
               "Indicator") # 设置图例项的显示文本
  ) +
  # 修改 Y 轴比例尺，使用命名向量作为标签 (保持不变)
  scale_y_discrete(labels = f_genus_labels_map) + # 使用创建好的命名向量
  labs(x = "Fungal relative abundance (%)", # 设置横坐标标题 (保持不变，标题已经包含%)
       y = NULL) + # 将纵坐标标题设置为空 (保持不变)
  theme(
    panel.border = element_blank(), # 移除面板周围的边框 (保持不变)
    panel.grid.minor = element_blank(), # 移除次要网格线 (保持不变)
    #panel.grid.major = element_blank(), # 移除主要网格线 (如果您想保留主要网格线，请删除或注释掉这行)
    # 添加移除刻度线
    axis.ticks.x = element_blank(), # 移除横坐标刻度线
    axis.ticks.y = element_blank(), # 移除纵坐标刻度线
    axis.title.x = element_text(size = 16),# 修改横坐标标题大小
    #修改 Y 轴文本外观和对齐
    axis.text.y = element_text(face = "italic", # 设置字体为斜体
                               vjust = -0.5, # 设置垂直对齐方式为居中
                               hjust = 0, # 设置为左对齐
                               margin = margin(r = unit(-40, "pt")),# 设置右侧边距为0点，使其紧贴Y轴
                               size = 14), 
    # 修改 X 轴文本外观和对齐
    axis.text.x = element_text(margin = margin(t = unit(-20, "pt")),#设置顶部边距为0点，使其紧贴X轴
                               size = 14,
                               vjust = 0 ), # 设置为底端对齐 
    legend.text = element_text(size = 14), # 设置图例文本大小
    legend.spacing.y = unit(0.5, "cm") # 设置图例项之间的垂直间距，例如 0.5 厘米
  )
ggridge_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3d_ggridge_fungi.png", plot = ggridge_fungi, width = 6, height = 6, dpi = 600, bg = "transparent")
#找到相对丰度最高的门科属
fungi_ASV <- read.csv(file = "D:/study/master/meiji/fungi_ASV.csv",
                      sep=",",header=TRUE,check.names = FALSE)
# 定义文件的绝对路径（使用正斜杠）
metadata <- read_tsv("D:/study/master/metadata.tsv")
colnames(metadata)[1]<-"SampleID"
# 查看元数据
head(metadata)
# 读取 ASV 表（特征表）
f_ASV <- fungi_ASV[, c("ASV",metadata$SampleID)]
f_ASV <- f_ASV %>% column_to_rownames(var = "ASV")
head(f_ASV)
# 读取 taxonomy（分类信息）
f_tax <- fungi_ASV[, 2:9]
f_tax <- f_tax %>% column_to_rownames(var = "ASV")
head(f_tax)
# 提取 ASV 表
f_ASV_df <- f_ASV[rowSums(f_ASV[])>0,]
f_ASV_df$asv<-rownames(f_ASV_df)
f_tax$asv<-rownames(f_tax)
#在每个 Genus 分组，将同一个属下的所有 ASV 的计数在每个样本中加起来，得到每个属在每个样本的总计数。
f_dat_m <- f_ASV_df %>%
  left_join(f_tax)%>%#合并
  group_by(Genus)%>%#按属分组
  summarise_if(is.numeric, sum)%>%# 找到所有数值型的列，并对这些列计算它们的 sum
  na.omit()%>%
  pivot_longer(cols = !Genus)%>%# 从“宽”格式转换为“长”格式
  left_join(metadata, by = c("name" = "SampleID"))#合并
f_sum_all <- sum(f_dat_m$value)
f_num_all <- length(unique(f_dat_m$name))
#计算属相对丰度（单个属的菌/总菌）
f_dat_abun <- f_dat_m %>%
  group_by(Genus)%>%
  summarise(Abundance = sum(value)/f_sum_all)
colnames(f_dat_abun) <- c("Genus", "Abundance")
f_tax_sel <- f_tax%>%
  select(-c(asv,Species))%>%
  distinct()%>%#移除重复的行
  filter(Genus %in% f_dat_abun$Genus)
#结合属、相对丰度、分类信息
f_dat_final <- f_dat_abun%>%
  left_join(f_tax_sel)
#计算每个门的总丰度，并选出总丰度最高的门。
f_phylum_top <- f_dat_final %>%
  group_by(Phylum)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))
#找出相对丰度最高的科
f_family_top <- f_dat_final %>%
  group_by(Family)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))
#找出相对丰度最高的属
f_dat_final %>%
  dplyr::arrange(desc(Abundance)) -> f_genus_top
View(f_phylum_top)
View(f_family_top)
View(f_genus_top)

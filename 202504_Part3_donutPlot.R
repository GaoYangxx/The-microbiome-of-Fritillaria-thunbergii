#图3a
install.packages(c("RColorBrewer", "ggplot2", "tidyverse", "reshape2", "ggpubr", " dplyr"))
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
#计算每类有多少种asv
b_datToPlot <- b_dat_numb %>%
  arrange(ASV, type) %>%#先对 ASV 和 type 排序，确保拼接顺序一致
  group_by(ASV)%>%# 每个ASV为一组
  summarise(type = str_c(type, collapse="_"))%>%#将同一个ASV的多个 type 合并成一个字符串，比如 "core_indicator"
  ungroup()%>%#取消分组状态
  bind_rows(b_ASV_other)%>%#合并"Other" 的asv
  group_by(type)%>%# 按 type 统计每一类ASV的种类数量
  summarise(sum = n())
b_datToPlot$type <- factor(b_datToPlot$type, levels = c("core","core_indicator", "indicator","endemic_indicator", "endemic_uniq", "endemic", "Other"))#设置分类顺序与颜色
colors <- c("#b10026","#006400", "#F2CA7F",
            "#4AC0FF","#2F609F","#30308B","#737373")
names(colors) <- c("core","core_indicator", "indicator",
                   "endemic_indicator", "endemic_uniq", "endemic", "Other")
b_datToPlot <- b_datToPlot[order(b_datToPlot$type),]#排序
b_datToPlot$fraction = b_datToPlot$sum / sum(b_datToPlot$sum)#计算每一类所占百分比
b_datToPlot$ymax = cumsum(b_datToPlot$fraction)#计算累积百分比
b_datToPlot$ymin = c(0, head(b_datToPlot$ymax, n=-1))# 计算累积最小值
b_datToPlot$labelPosition <- (b_datToPlot$ymax + b_datToPlot$ymin) / 2#计算标签的中心点
b_p1 <- ggplot(b_datToPlot, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect() + # 画出环形每一块矩形
  #geom_label( x=4, aes(y=labelPosition, label=type), size=6) + # 标签放在每一段中间
  scale_fill_manual(
    values = colors,
    labels = c("Core", "Core and indicator", "Indicator",
               "Specific and indicator", "Specific and unique", "Specific", "Other")
  )+
  coord_polar(theta="y") +#矩形再转换为极坐标系的环形图（甜甜圈图）
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) + # 控制环的厚度
  theme_void() +
  theme(
    legend.title = element_blank(),# 取消图例标题
    legend.text = element_text(size = 14) # , legend.position = "none"#取消图例
  )+ 
  guides(fill = guide_legend(nrow = 1)) 
b_p1
b_dat_sum <- rowSums(b_dat[1:48]) # 对每个 ASV 求所有样本丰度之和
b_dat_sum <- as.data.frame(b_dat_sum)
b_dat_sum$ASV <- rownames(b_dat_sum) # 加上 ASV 名称
#计算每类的相对丰度
b_datToPlotCov <- b_dat_numb %>%
  arrange(ASV, type) %>%
  group_by(ASV)%>%
  summarise(type = str_c(type, collapse="_"))%>% # 合并同一个ASV多个type
  ungroup()%>%
  bind_rows(b_ASV_other)%>% 
  left_join(b_dat_sum)%>%# 添加丰度信息（b_dat_sum 是每个 ASV 的总丰度）
  filter(!(b_dat_sum == 0))%>% # 移除丰度为0的ASV
  mutate(perc = b_dat_sum / sum(b_dat_sum))%>%#计算相对丰度
  group_by(type)%>%
  summarise(fraction = sum(perc)) # 统计每种类型的丰度总占比
sum(b_datToPlotCov$fraction)
b_datToPlotCov$type <- factor(b_datToPlotCov$type, levels = c("core","core_indicator", "indicator","endemic_indicator", "endemic_uniq", "endemic", "Other"))
b_datToPlotCov <- b_datToPlotCov[order(b_datToPlotCov$type),]
# b_datToPlotCov$fraction = b_datToPlotCov$sum / sum(b_datToPlotCov$sum)
b_datToPlotCov$ymax = cumsum(b_datToPlotCov$fraction)
b_datToPlotCov$ymin = c(0, head(b_datToPlotCov$ymax, n=-1))
b_datToPlotCov$labelPosition <- (b_datToPlotCov$ymax + b_datToPlotCov$ymin) / 2
b_p2 <- ggplot(b_datToPlotCov, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect() +
  #geom_label(x=4, aes(y=labelPosition, label=type), size=6) +
  scale_fill_manual(
    values = colors,
    labels = c("Core", "Core and indicator", "Indicator",
               "Specific and indicator", "Specific and unique", "Specific", "Other")
  )+
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(
    legend.title = element_blank(),# 取消图例标题
    legend.text = element_text(size = 14)# , legend.position = "none"#取消图例
  )+ 
  guides(fill = guide_legend(nrow = 1)) 
b_p2
p3_bacteria <- ggarrange(
  b_p1, b_p2,
  ncol = 2, 
  common.legend = TRUE, 
  legend = "bottom"
)
p3_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3a_donut_bacteria.png", plot = p3_bacteria, width = 8, height = 6, dpi = 600, bg = "transparent")
#比较特殊在各生态位相对丰度，各生态位相对丰度的平均值为总相对丰度
b_specific2 <- as.data.frame(otu_table(b_merged_NOMIS_specific_ab))
b_asv2 <- as.data.frame(otu_table(b_merged))
G_samples <- rownames(metadata_nomis[metadata_nomis$Niche == "G", ])
N_samples <- rownames(metadata_nomis[metadata_nomis$Niche == "N", ])
# G组特异ASV的相对丰度
b_G_rel_ab <- sum(b_specific2[, G_samples]) / sum(b_asv2[, G_samples])
b_G_rel_ab
# N组特异ASV的相对丰度
b_N_rel_ab <- sum(b_specific2[, N_samples]) / sum(b_asv2[, N_samples])
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
#计算每类有多少种asv
f_datToPlot <- f_dat_numb %>%
  arrange(ASV, type) %>%#先对 ASV 和 type 排序，确保拼接顺序一致
  group_by(ASV)%>%# 每个ASV为一组
  summarise(type = str_c(type, collapse="_"))%>%#将同一个ASV的多个 type 合并成一个字符串，比如 "core_indicator"
  ungroup()%>%#取消分组状态
  bind_rows(f_ASV_other)%>%#合并"Other" 的asv
  group_by(type)%>%# 按 type 统计每一类ASV的种类数量
  summarise(sum = n())
f_datToPlot$type <- factor(f_datToPlot$type, levels = c("core","core_indicator", "indicator","endemic_indicator", "endemic_uniq", "endemic", "Other"))#设置分类顺序与颜色
colors <- c("#b10026","#006400", "#F2CA7F",
            "#4AC0FF","#2F609F","#30308B","#737373")
names(colors) <- c("core","core_indicator", "indicator",
                   "endemic_indicator", "endemic_uniq", "endemic", "Other")
f_datToPlot <- f_datToPlot[order(f_datToPlot$type),]#排序
f_datToPlot$fraction = f_datToPlot$sum / sum(f_datToPlot$sum)#计算每一类所占百分比
f_datToPlot$ymax = cumsum(f_datToPlot$fraction)#计算累积百分比
f_datToPlot$ymin = c(0, head(f_datToPlot$ymax, n=-1))# 计算累积最小值
f_datToPlot$labelPosition <- (f_datToPlot$ymax + f_datToPlot$ymin) / 2#计算标签的中心点
f_p1 <- ggplot(f_datToPlot, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect() + # 画出环形每一块矩形
  #geom_label( x=4, aes(y=labelPosition, label=type), size=6) + # 标签放在每一段中间
  scale_fill_manual(
    values = colors,
    labels = c("Core", "Core and indicator", "Indicator",
               "Specific and indicator", "Specific and unique", "Specific", "Other")
  )+
  coord_polar(theta="y") +#矩形再转换为极坐标系的环形图（甜甜圈图）
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) + # 控制环的厚度
  theme_void() +
  theme(
    legend.title = element_blank(),# 取消图例标题
    legend.text = element_text(size = 14) # , legend.position = "none"#取消图例
  )+ 
  guides(fill = guide_legend(nrow = 1)) 
f_p1
f_dat_sum <- rowSums(f_dat[1:48]) # 对每个 ASV 求所有样本丰度之和
f_dat_sum <- as.data.frame(f_dat_sum)
f_dat_sum$ASV <- rownames(f_dat_sum) # 加上 ASV 名称
#计算每类的相对丰度
f_datToPlotCov <- f_dat_numb %>%
  arrange(ASV, type) %>%
  group_by(ASV)%>%
  summarise(type = str_c(type, collapse="_"))%>% # 合并同一个ASV多个type
  ungroup()%>%
  bind_rows(f_ASV_other)%>% 
  left_join(f_dat_sum)%>%# 添加丰度信息（f_dat_sum 是每个 ASV 的总丰度）
  filter(!(f_dat_sum == 0))%>% # 移除丰度为0的ASV
  mutate(perc = f_dat_sum / sum(f_dat_sum))%>%#计算相对丰度
  group_by(type)%>%
  summarise(fraction = sum(perc)) # 统计每种类型的丰度总占比
sum(f_datToPlotCov$fraction)
f_datToPlotCov$type <- factor(f_datToPlotCov$type, levels = c("core","core_indicator", "indicator","endemic_indicator", "endemic_uniq", "endemic", "Other"))
f_datToPlotCov <- f_datToPlotCov[order(f_datToPlotCov$type),]
# f_datToPlotCov$fraction = f_datToPlotCov$sum / sum(f_datToPlotCov$sum)
f_datToPlotCov$ymax = cumsum(f_datToPlotCov$fraction)
f_datToPlotCov$ymin = c(0, head(f_datToPlotCov$ymax, n=-1))
f_datToPlotCov$labelPosition <- (f_datToPlotCov$ymax + f_datToPlotCov$ymin) / 2
f_p2 <- ggplot(f_datToPlotCov, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
  geom_rect() +
  #geom_label(x=4, aes(y=labelPosition, label=type), size=6) +
  scale_fill_manual(
    values = colors,
    labels = c("Core", "Core and indicator", "Indicator",
               "Specific and indicator", "Specific and unique", "Specific", "Other")
  )+
  coord_polar(theta="y") +
  #xlim(c(2, 4)) +
  xlim(c(1, 4)) +
  theme_void() +
  theme(
    legend.title = element_blank(),# 取消图例标题
    legend.text = element_text(size = 14)# , legend.position = "none"#取消图例
  )+ 
  guides(fill = guide_legend(nrow = 1)) 
f_p2
p3_fungi <- ggarrange(
  f_p1, f_p2,
  ncol = 2, 
  common.legend = TRUE, 
  legend = "bottom"
)
p3_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3a_donut_fungi.png", plot = p3_fungi, width = 8, height = 6, dpi = 600, bg = "transparent")
#比较特殊在各生态位相对丰度，各生态位相对丰度的平均值为总相对丰度
f_specific2 <- as.data.frame(otu_table(f_merged_NOMIS_specific_ab))
f_asv2 <- as.data.frame(otu_table(f_merged))
G_samples <- rownames(metadata_nomis[metadata_nomis$Niche == "G", ])
N_samples <- rownames(metadata_nomis[metadata_nomis$Niche == "N", ])
# G组特异ASV的相对丰度
f_G_rel_ab <- sum(f_specific2[, G_samples]) / sum(f_asv2[, G_samples])
f_G_rel_ab
# N组特异ASV的相对丰度
f_N_rel_ab <- sum(f_specific2[, N_samples]) / sum(f_asv2[, N_samples])
f_N_rel_ab

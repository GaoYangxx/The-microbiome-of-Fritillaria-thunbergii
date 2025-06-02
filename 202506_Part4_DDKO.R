#图4b
#https://zenodo.org/records/10035668/files/Tax4Fun2_1.1.5.tar.gz?download=1
#https://zenodo.org/records/10035668/files/Tax4Fun2_ReferenceData_v2.tar.gz?download=1
install.packages("D:/study/master/KEGG/Tax4Fun2_1.1.5.tar.gz", repos = NULL, type = "source")
install.packages("devtools")
install.packages("remotes")
remotes::install_github("kasperskytte/ampvis2") 
list_of_packages <- c("tidyverse", "mgcv", "data.table", "gt", "vegan", "randomForest", "smacof", "ampvis2", "Maaslin2", "readxl", "ggsci", "ggalluvial", "pheatmap", "ggpubr", "geosphere", "here", "qiime2R", "phyloseq", "Tax4Fun2", "permuco", "effectsize", "permutes", "buildmer", "fishualize")
for (package_name in list_of_packages) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name)
  }
}
devtools::install_github("brendanf/FUNGuildR")
library(tidyverse)
library(mgcv)
library(data.table)
library(gt)
library(vegan)
library(randomForest)
library(smacof)
library(ampvis2)
library(Maaslin2)
library(readxl)
library(ggsci)
library(ggalluvial)
library(pheatmap)
library(ggpubr)
library(geosphere)
library(here)
library(qiime2R)
library(phyloseq)
library(Tax4Fun2)
library(permuco)
library(effectsize)
library(permutes)
library(buildmer)
library(fishualize)
library(FUNGuildR)
#细菌功能多样性距离衰减
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
b_ASV_df_table <- as.data.frame(b_ASV_df) 
b_ASV_df_table <- rownames_to_column(b_ASV_df_table, var = "ID")
write.table(b_ASV_df_table, "D:/study/master/KEGG/b_ASV_df_table.txt", sep = "\t", quote = FALSE, row.names = TRUE)
#KEGG 参考基因组注释，fasta过大可以拆分再合并放在原处
runRefBlast(path_to_otus = 'D:/study/master/meiji/b_ASV_reps.fasta', 
            path_to_reference_data = 'D:/study/master/KEGG/Tax4Fun2_ReferenceData_v2',
            path_to_temp_folder = "D:/study/master/KEGG/Kelp_Ref99NR",
            database_mode = 'Ref99NR', 
            use_force = TRUE, 
            num_threads = 4)
ref_blast <- fread("D:/study/master/KEGG/ref_blast/ref_blast_md5.txt")
ref_blast$V1 <- b_rename_vector[ref_blast$V1] #替换md5值
write.table(ref_blast, "D:/study/master/KEGG/Kelp_Ref99NR/ref_blast.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#预测群落功能
makeFunctionalPrediction(path_to_otu_table = 'D:/study/master/KEGG/b_ASV_df_table.txt', 
                         path_to_reference_data = 'D:/study/master/KEGG/Tax4Fun2_ReferenceData_v2',
                         path_to_temp_folder = 'D:/study/master/KEGG/Kelp_Ref99NR', 
                         database_mode = 'Ref99NR',
                         normalize_by_copy_number = TRUE, #默认，用参考数据库中每个序列计算的16S rRNA拷贝数的平均值进行归一化
                         min_identity_to_reference = 0.97, 
                         normalize_pathways = FALSE)#默认，将把每个KO的相对丰度关联到它所属的每个路径上
# 加载数据
kegg <- read.table(file = "D:/study/master/KEGG/Kelp_Ref99NR/functional_prediction.txt", header = TRUE, sep = "\t", quote = "", comment.char = "")# 读取KEGG计数文件
b_otu_mat <- kegg[,-ncol(kegg)]
rownames(b_otu_mat) <- b_otu_mat[,1]
b_otu_mat <- b_otu_mat[,-1]
# 计算 Bray-Curtis 相异性矩阵
b_bray_curtis_matrix <- vegdist(t(b_otu_mat), method = "bray")
b_bray_curtis_matrix_mat <- as.matrix(b_bray_curtis_matrix)
# 计算 Sørensen 相异性矩阵
b_sorensen_matrix <- vegdist(t(b_otu_mat), method = "bray", binary=T)
b_sorensen_matrix_mat <- as.matrix(b_sorensen_matrix)
# 使用 geosphere 包计算地理距离
otu_mat_dist <- data.frame(metadata$longitude, metadata$latitude)#经纬度
rownames(otu_mat_dist) <- rownames(metadata)
otu_mat_dist$metadata.longitude <- as.numeric(otu_mat_dist$metadata.longitude)
otu_mat_dist$metadata.latitude <- as.numeric(otu_mat_dist$metadata.latitude)
geographic_distances <- distm(otu_mat_dist, fun=distGeo)
rownames(geographic_distances) <- rownames(otu_mat_dist)
colnames(geographic_distances) <- rownames(otu_mat_dist)
# 创建 Bray-Curtis 相异性与地理距离的数据框
b_bray_curtis_df <- data.frame(Dissimilarity = as.vector(b_bray_curtis_matrix_mat), GeographicDistance = as.vector(geographic_distances), Type="Bray-Curtis")
# 过滤 Bray-Curtis 相异性不为 0 的数据
b_braycurtis_df_filtered <- b_bray_curtis_df %>% filter(Dissimilarity != 0)
b_allbray_m<- as.matrix(b_bray_curtis_matrix)
diag(b_allbray_m)=NA # 将对角线设为 NA
b_allbray_diag<-t(matrix(t(b_allbray_m)[which(!is.na(b_allbray_m))],nrow=47,ncol=48)) # 处理 NA 值
min(b_allbray_diag)
geographic_distances <- distm(otu_mat_dist, fun=distGeo)
dist_geo_all <- as.matrix(geographic_distances)
diag(dist_geo_all)= NA # 将对角线设为 NA
dist_geo_all_diag<-t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=47,ncol=48)) # 处理 NA 值
# 进行 Mantel 检验
mantel(b_allbray_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
# # 创建 Bray-Curtis 相异性与地理距离的散点图 (已注释)
# b_bray_curtis_plot <- ggplot(b_bray_curtis_df, aes(x = log10(GeographicDistance), y = 1-Dissimilarity )) +
#    geom_point() +
#    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # 添加线性趋势线
#    stat_cor(method = "pearson", label.x = 0.8, label.y = 0.6, label.sep = 0.1, size = 4) +  # 添加相关系数
#    labs(x = "Geographic Distance", y = "Bray-Curtis Dissimilarity", title = "Bray-Curtis Dissimilarity vs. Geographic Distance")
# # 打印图表
# b_bray_curtis_plot
# 为 Sørensen 相异性与地理距离创建数据框
b_sorensen_df <- data.frame(Dissimilarity = as.vector(b_sorensen_matrix_mat), GeographicDistance = as.vector(geographic_distances), Type = "Sorensen")
# 过滤 Sørensen 相异性不为 0 的数据
b_sorensen_df_filtered <- b_sorensen_df %>% filter(Dissimilarity != 0)
# Mantel 检验
diag(b_sorensen_matrix_mat)=NA # 对角线设为 NA
b_sorensen_diag<-t(matrix(t(b_sorensen_matrix_mat)[which(!is.na(b_sorensen_matrix_mat))],nrow=47,ncol=48)) # 处理 NA 值
mantel(b_sorensen_diag, dist_geo_all_diag, method = "pearson", permutations = 9999, na.rm = TRUE)
# 合并 Bray-Curtis 和 Sørensen 数据框
b_combined_df <- rbind(b_braycurtis_df_filtered, b_sorensen_df_filtered)
#write.csv(b_combined_df, "202505_b_combined_dd_function.csv") # 可选: 写入 CSV 文件
get_fishcol <- fish(4, option="Scarus_quoyi") # 获取鱼类调色板颜色
custom_colors <- c("Bray-Curtis" = "#3F459BFF", "Sorensen" = "#009E9EFF") # 自定义颜色
ggplot(b_combined_df) + # 创建 ggplot 对象
  geom_point(aes(y=1-Dissimilarity, x = GeographicDistance/1000, color= Type), alpha=0.1) + # 散点图
  scale_color_manual(values = custom_colors)+ # 使用自定义颜色
  #scale_color_fish(option="Scarus_quoyi",discrete=T,direction=1)+ # 可选: 使用鱼类颜色
  ylab("Community Similarity") + xlab("Geographic distance (km)") + ggtitle("") + theme_bw() + # 设置标签和主题
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) # 自定义主题
# 确保地理距离是数值型
b_combined_df$GeographicDistance <- as.numeric(b_combined_df$GeographicDistance)
# 将地理距离转换为千米
b_combined_df$GeographicDistance_km <- b_combined_df$GeographicDistance / 1000
# 创建计算 R 平方值的函数
b_calculate_r_squared <- function(model) { return(format(summary(model)$r.squared, digits = 3)) }
# 对每个 Type 拟合线性回归模型
b_lm_models <- by(b_combined_df, b_combined_df$Type, function(sub_df) { lm(1 - sub_df$Dissimilarity ~ sub_df$GeographicDistance_km, data = sub_df) })
# 从线性回归模型中提取系数和 R 平方值
b_intercepts <- sapply(b_lm_models, function(model) coef(model)[1]) # 提取截距
b_slopes <- sapply(b_lm_models, function(model) coef(model)[2])   # 提取斜率
b_r_squared <- sapply(b_lm_models, b_calculate_r_squared)    # 计算 R 平方值
get_fishcol <- fish(4, option="Scarus_quoyi") # 获取鱼类调色板颜色 
custom_colors <- c("Bray-Curtis" = "#3F459BFF", "Sorensen" = "#009E9EFF") # 自定义颜色
b_combined_df_clean <- b_combined_df %>% filter(GeographicDistance !=0)
#stat_regline_equation可以计算回归方程
ddko_plot_bacteria <- ggplot(b_combined_df, aes(x = GeographicDistance/1000, y = 1 - Dissimilarity, color = Type)) + # 创建 ggplot 对象
  geom_point(alpha = 0.1) + # 散点图
  scale_color_manual(name = "",
                     values = custom_colors,
                     breaks = c("Sorensen", "Bray-Curtis"), # 明确指定 breaks 的顺序
                     labels = c("Sørensen", "Bray–Curtis")) + # 自定义颜色和图例标签
  stat_smooth(aes(y = 1 - Dissimilarity, x = GeographicDistance/1000, color = Type),
              method = "lm", formula = y ~ x, size = 1.2, se = FALSE, linetype = "solid") + # 添加线性回归平滑曲线
  ylab("Bacterial functional similarity") +
  xlab("Geographic distance (km)") +
  theme_bw() +
  #scale_x_log10() +
  scale_y_continuous(breaks = seq(0.8, 1, by = 0.1), # 设置 y 轴刻度，包括 1.0
                     limits = c(min(1 - b_combined_df_clean$Dissimilarity, na.rm = TRUE), 1)) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))
ddko_plot_bacteria
#ggsave("D:/study/master/Main_Figure_tables/Figure_4/4b_distancedecaykegg_bacteria.png", plot = ddko_plot_bacteria, width = 6.5, height = 6, dpi = 600, bg = "transparent")
# 检验斜率差异
b_dist_bray <- data.frame(sim_bc=as.vector(1- b_allbray_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers") # Bray-Curtis 数据框
b_dist_sor <- data.frame(sim_sor=as.vector(1- b_sorensen_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers") # Sørensen 数据框
b_merge_bcsor <- merge(b_dist_bray, b_dist_sor, by="row.names") # 合并数据框
b_anco_bcsor <- b_merge_bcsor[c("sim_bc","sim_sor","geo.dist.x")] # 选择相关列
b_manco <- melt(b_anco_bcsor, id=c("geo.dist.x")) # 转换为长格式
ggscatter(b_manco, x = "geo.dist.x", y = "value", color = "variable", add = "reg.line") + # 散点图和回归线
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)) # 添加回归方程和 R 平方值
b_manco$geo.dist.log <- log((b_manco$geo.dist.x)/1000) # 计算对数地理距离
#只保留地理距离 (geo_dist.x) 大于或等于 1000 米（即 1 千米）的数据点
b_manco <- b_manco %>% filter(geo.dist.log >= 0)
b_glm_dd_bcsor <- (lm(b_manco, formula = value ~ geo.dist.log + geo.dist.log:variable + variable)) # 线性模型 (可能违反假设)
# 非参数 ANCOVA
b_formula <- value ~ geo.dist.log + geo.dist.log:variable + variable # 模型公式
b_model <- lm(b_formula, data = b_manco) # 拟合线性模型
b_perm_test <- aovperm(b_model, nperm = 999, method="freedman_lane") # 执行置换检验
b_interaction_p_value <- b_perm_test$results$`geo.dist.log:variable`$p.value # 提取交互作用的 p 值
#真菌功能多样性距离衰减
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
f_ASV<-as.data.frame(f_ASV)
f_ASV$ASV <- rownames(f_ASV)
#转换 taxonomy 表为 metacoder 兼容格式
f_tax_df<-as.data.frame(f_tax)
f_tax_df$ASV <- rownames(f_tax_df)
f_tax_df$Taxonomy <- paste(      #参考parse_tax_data函数参数
  f_tax_df$Kingdom,          # 添加 Kingdom
  f_tax_df$Phylum,           # 添加 Phylum
  f_tax_df$Class,            # 添加 Class
  f_tax_df$Order,            # 添加 Order
  f_tax_df$Family,           # 添加 Family
  f_tax_df$Genus,            # 添加 Genus
  f_tax_df$Species,          # 添加 Species
  sep = ";"                 # 使用分号作为分隔符
)
f_tax_df <-f_tax_df[, c("ASV", "Taxonomy")]
f_tax_df$Taxonomy <- gsub("k__|p__|c__|o__|f__|g__|s__", "", f_tax_df$Taxonomy)
f_tax_df$Taxonomy <- gsub("_", " ", f_tax_df$Taxonomy)
#新数据库http://www.stbates.org/funguild_db_2.php不稳定
funguild_db <- get_funguild_db(db = "http://www.stbates.org/funguild_db.php")
funguild_result <- funguild_assign(f_tax_df, db = funguild_db)
View(funguild_result)
funguild <- left_join(f_ASV, funguild_result, by = "ASV")
#去除guild中的NA
#funguild <- funguild %>%
#  filter(!is.na(guild))
# 获取所有样本列的名称
sample_cols <- colnames(funguild)[1:48] # 假设前 48 列是样本列
# 选择样本列和 guild 列
guild_abundance <- funguild %>%
  select(all_of(sample_cols), guild)
# 将数据从宽格式转换为长格式，方便按样本计算总丰度
guild_long <- guild_abundance %>%
  pivot_longer(cols = -guild, names_to = "Sample", values_to = "Abundance")
# 按样本计算总丰度
total_abundance_per_sample <- guild_long %>%
  group_by(Sample) %>%
  summarise(TotalAbundance = sum(Abundance))
# 将总丰度信息合并回长格式的数据
guild_long <- guild_long %>%
  left_join(total_abundance_per_sample, by = "Sample")
# 按 Guild 和样本分组，计算相对丰度
relative_abundance <- guild_long %>%
  group_by(guild, Sample) %>%
  summarise(RelativeAbundance = sum(Abundance) / first(TotalAbundance), .groups = 'drop')
# 将结果转换回宽格式，Guild 作为行，样本作为列
relative_abundance_wide <- relative_abundance %>%
  pivot_wider(names_from = Sample, values_from = RelativeAbundance, values_fill = 0)
relative_abundance_wide <- relative_abundance_wide %>%
  mutate(go = paste0("guild", 1:n()))
# 加载数据
f_otu_mat <- relative_abundance_wide[,-1]
f_otu_mat <- as.data.frame(f_otu_mat)
rownames(f_otu_mat) <- f_otu_mat[, ncol(f_otu_mat)]
f_otu_mat <- f_otu_mat[,-ncol(f_otu_mat)]
# 计算 Bray-Curtis 相异性矩阵
f_bray_curtis_matrix <- vegdist(t(f_otu_mat), method = "bray")
f_bray_curtis_matrix_mat <- as.matrix(f_bray_curtis_matrix)
# 计算 Sørensen 相异性矩阵
f_sorensen_matrix <- vegdist(t(f_otu_mat), method = "bray", binary=T)
f_sorensen_matrix_mat <- as.matrix(f_sorensen_matrix)
# 使用 geosphere 包计算地理距离
otu_mat_dist <- data.frame(metadata$longitude, metadata$latitude)#经纬度
rownames(otu_mat_dist) <- rownames(metadata)
otu_mat_dist$metadata.longitude <- as.numeric(otu_mat_dist$metadata.longitude)
otu_mat_dist$metadata.latitude <- as.numeric(otu_mat_dist$metadata.latitude)
geographic_distances <- distm(otu_mat_dist, fun=distGeo)
rownames(geographic_distances) <- rownames(otu_mat_dist)
colnames(geographic_distances) <- rownames(otu_mat_dist)
# 创建 Bray-Curtis 相异性与地理距离的数据框
f_bray_curtis_df <- data.frame(Dissimilarity = as.vector(f_bray_curtis_matrix_mat), GeographicDistance = as.vector(geographic_distances), Type="Bray-Curtis")
# 过滤 Bray-Curtis 相异性不为 0 的数据
f_braycurtis_df_filtered <- f_bray_curtis_df %>% filter(Dissimilarity != 0)
f_allbray_m<- as.matrix(f_bray_curtis_matrix)
diag(f_allbray_m)=NA # 将对角线设为 NA
f_allbray_diag<-t(matrix(t(f_allbray_m)[which(!is.na(f_allbray_m))],nrow=47,ncol=48)) # 处理 NA 值
min(f_allbray_diag)
geographic_distances <- distm(otu_mat_dist, fun=distGeo)
dist_geo_all <- as.matrix(geographic_distances)
diag(dist_geo_all)= NA # 将对角线设为 NA
dist_geo_all_diag<-t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=47,ncol=48)) # 处理 NA 值
# 进行 Mantel 检验
mantel(f_allbray_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
# # 创建 Bray-Curtis 相异性与地理距离的散点图 (已注释)
# f_bray_curtis_plot <- ggplot(f_bray_curtis_df, aes(x = log10(GeographicDistance), y = 1-Dissimilarity )) +
#    geom_point() +
#    geom_smooth(method = "lm", se = FALSE, color = "blue") +  # 添加线性趋势线
#    stat_cor(method = "pearson", label.x = 0.8, label.y = 0.6, label.sep = 0.1, size = 4) +  # 添加相关系数
#    labs(x = "Geographic Distance", y = "Bray-Curtis Dissimilarity", title = "Bray-Curtis Dissimilarity vs. Geographic Distance")
# # 打印图表
# f_bray_curtis_plot
# 为 Sørensen 相异性与地理距离创建数据框
f_sorensen_df <- data.frame(Dissimilarity = as.vector(f_sorensen_matrix_mat), GeographicDistance = as.vector(geographic_distances), Type = "Sorensen")
# 过滤 Sørensen 相异性不为 0 的数据
f_sorensen_df_filtered <- f_sorensen_df %>% filter(Dissimilarity != 0)
# Mantel 检验
diag(f_sorensen_matrix_mat)=NA # 对角线设为 NA
f_sorensen_diag<-t(matrix(t(f_sorensen_matrix_mat)[which(!is.na(f_sorensen_matrix_mat))],nrow=47,ncol=48)) # 处理 NA 值
mantel(f_sorensen_diag, dist_geo_all_diag, method = "pearson", permutations = 9999, na.rm = TRUE)
# 合并 Bray-Curtis 和 Sørensen 数据框
f_combined_df <- rbind(f_braycurtis_df_filtered, f_sorensen_df_filtered)
#write.csv(f_combined_df, "202505_f_combined_dd_function.csv") # 可选: 写入 CSV 文件
get_fishcol <- fish(4, option="Scarus_quoyi") # 获取鱼类调色板颜色
custom_colors <- c("Bray-Curtis" = "#3F459BFF", "Sorensen" = "#009E9EFF") # 自定义颜色
ggplot(f_combined_df) + # 创建 ggplot 对象
  geom_point(aes(y=1-Dissimilarity, x = GeographicDistance/1000, color= Type), alpha=0.1) + # 散点图
  scale_color_manual(values = custom_colors)+ # 使用自定义颜色
  #scale_color_fish(option="Scarus_quoyi",discrete=T,direction=1)+ # 可选: 使用鱼类颜色
  ylab("Community Similarity") + xlab("Geographic distance (km)") + ggtitle("") + theme_bw() + # 设置标签和主题
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) # 自定义主题
# 确保地理距离是数值型
f_combined_df$GeographicDistance <- as.numeric(f_combined_df$GeographicDistance)
# 将地理距离转换为千米
f_combined_df$GeographicDistance_km <- f_combined_df$GeographicDistance / 1000
# 创建计算 R 平方值的函数
f_calculate_r_squared <- function(model) { return(format(summary(model)$r.squared, digits = 3)) }
# 对每个 Type 拟合线性回归模型
f_lm_models <- by(f_combined_df, f_combined_df$Type, function(suf_df) { lm(1 - suf_df$Dissimilarity ~ suf_df$GeographicDistance_km, data = suf_df) })
# 从线性回归模型中提取系数和 R 平方值
f_intercepts <- sapply(f_lm_models, function(model) coef(model)[1]) # 提取截距
f_slopes <- sapply(f_lm_models, function(model) coef(model)[2])   # 提取斜率
f_r_squared <- sapply(f_lm_models, f_calculate_r_squared)    # 计算 R 平方值
get_fishcol <- fish(4, option="Scarus_quoyi") # 获取鱼类调色板颜色 
custom_colors <- c("Bray-Curtis" = "#3F459BFF", "Sorensen" = "#009E9EFF") # 自定义颜色
f_combined_df_clean <- f_combined_df %>% filter(GeographicDistance !=0)
#stat_regline_equation可以计算回归方程
ddko_plot_fungi <- ggplot(f_combined_df, aes(x = GeographicDistance/1000, y = 1 - Dissimilarity, color = Type)) + # 创建 ggplot 对象
  geom_point(alpha = 0.1) + # 散点图
  scale_color_manual(name = "",
                     values = custom_colors,
                     breaks = c("Sorensen", "Bray-Curtis"), # 明确指定 breaks 的顺序
                     labels = c("Sørensen", "Bray–Curtis")) + # 自定义颜色和图例标签
  stat_smooth(aes(y = 1 - Dissimilarity, x = GeographicDistance/1000, color = Type),
              method = "lm", formula = y ~ x, size = 1.2, se = FALSE, linetype = "solid") + # 添加线性回归平滑曲线
  ylab("Fungal functional similarity") +
  xlab("Geographic distance (km)") +
  theme_bw() +
  #scale_x_log10() +
  scale_y_continuous(limits = c(min(1 - f_combined_df$Dissimilarity, na.rm = TRUE), 1)) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))
ddko_plot_fungi
#ggsave("D:/study/master/Main_Figure_tables/Figure_4/4b_distancedecaykegg_fungi.png", plot = ddko_plot_fungi, width = 6.5, height = 6, dpi = 600, bg = "transparent")
# 检验斜率差异
f_dist_bray <- data.frame(sim_bc=as.vector(1- f_allbray_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers") # Bray-Curtis 数据框
f_dist_sor <- data.frame(sim_sor=as.vector(1- f_sorensen_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers") # Sørensen 数据框
f_merge_bcsor <- merge(f_dist_bray, f_dist_sor, by="row.names") # 合并数据框
f_anco_bcsor <- f_merge_bcsor[c("sim_bc","sim_sor","geo.dist.x")] # 选择相关列
f_manco <- melt(f_anco_bcsor, id=c("geo.dist.x")) # 转换为长格式
ggscatter(f_manco, x = "geo.dist.x", y = "value", color = "variable", add = "reg.line") + # 散点图和回归线
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)) # 添加回归方程和 R 平方值
f_manco$geo.dist.log <- log((f_manco$geo.dist.x)/1000) # 计算对数地理距离
#只保留地理距离 (geo_dist.x) 大于或等于 1000 米（即 1 千米）的数据点
f_manco <- f_manco %>% filter(geo.dist.log >= 0)
f_glm_dd_bcsor <- (lm(f_manco, formula = value ~ geo.dist.log + geo.dist.log:variable + variable)) # 线性模型 (可能违反假设)
# 非参数 ANCOVA
f_formula <- value ~ geo.dist.log + geo.dist.log:variable + variable # 模型公式
f_model <- lm(f_formula, data = f_manco) # 拟合线性模型
f_perm_test <- aovperm(f_model, nperm = 999, method="freedman_lane") # 执行置换检验
f_interaction_p_value <- f_perm_test$results$`geo.dist.log:variable`$p.value # 提取交互作用的 p 值

#图4a
install.packages(c("ggplot2", "dplyr", "scales", "reshape2", "tibble", "tidyverse", 
                   "rstatix", "broom", "nortest", "vegan", "geosphere", "fANCOVA", 
                   "phyloseqCompanion", "performance", "ape", "readxl"))
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "qiime2R", "rbiom"))
devtools::install_github("fishualize/fishualize")
devtools::install_github("speedyseq/speedyseq")
devtools::install_github("statnet/fANCOVA")
library(speedyseq)
library(phyloseq)
library(phyloseqCompanion)
library(geosphere)
library(vegan)
library(rbiom)
library(scales)
library(ggplot2)
library(dplyr)
library(fishualize)
library(ggpubr)
library(reshape2)
library(performance)
library(fANCOVA)
library(tibble)
library(qiime2R)
library(ape)
library(readxl)
library(tidyverse)
library(rstatix)
library(broom)
library(nortest)
#细菌群落多样性距离衰减
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
#转置
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
NOMIS_df <- data.frame(metadata$longitude, metadata$latitude)#经纬度
NOMIS_df$metadata.longitude <- as.numeric(NOMIS_df$metadata.longitude)
NOMIS_df$metadata.latitude <- as.numeric(NOMIS_df$metadata.latitude)
#计算经纬度的地理距离矩阵
dist_geo_all <- distm(NOMIS_df, NOMIS_df, fun=distGeo)
dist_geo_all <- as.matrix(dist_geo_all)
diag(dist_geo_all)=NA#对角线元素设置为NA
#提取所有非对角线（即样本间）的距离，并将它们放入一个矩形矩阵中
dist_geo_all_diag <- t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=47,ncol=48))
min(dist_geo_all_diag)# 最小值
#计算物种丰度的Bray-Curtis相似性指数（框架同上）
b_vegan_matrix_all <-vegan_otu(b_merged)
b_allregion_bray <-vegdist(log1p(b_vegan_matrix_all), method="bray")# Bray-Curtis 距离矩阵
b_allregion_m <- as.matrix(b_allregion_bray)
diag(b_allregion_m)=NA
b_allregion_diag <-t(matrix(t(b_allregion_m)[which(!is.na(b_allregion_m))],nrow=47,ncol=48))
min(b_allregion_diag)
#Mantel检验两个距离矩阵相关性
mantel(b_allregion_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
#组合Bray-Curtis 相异性值、Bray-Curtis 相似性 (1 - 相异性)、地理距离、字符串
b_dist_all_bray <- data.frame(b_BC_dist_bc=as.vector(b_allregion_diag), b_BC_sim_bc=as.vector(1- b_allregion_diag), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(b_dist_all_bray, aes(x=geo_dist/1000, y= b_BC_sim_bc)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Bray-Curtis") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
#提取 Bray-Curtis回归模型的截距和斜率
b_dist_all_bray %>% 
  do({
    # 用 geo_dist (地理距离) 作为解释变量来预测 BC_sim_bc (Bray-Curtis 相似性)
    mod = lm(b_BC_sim_bc ~ geo_dist, data = b_dist_all_bray)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
# 加权 UniFrac（考虑物种丰度信息）距离矩阵（框架同上）
b_asv_table_unif<-otu_table(b_merged, taxa_are_rows=T)
b_unifrac_essai<-phyloseq::UniFrac(b_merged, weighted=T, normalized=TRUE, parallel=FALSE, fast=TRUE) 
b_allregion_unif<-as.matrix(b_unifrac_essai)
diag(b_allregion_unif)=NA
b_allregion_diag_unif<-t(matrix(t(b_allregion_unif)[which(!is.na(b_allregion_unif))],nrow=47,ncol=48))
mantel(b_allregion_diag_unif, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
b_dist_all_wunif <- data.frame(dis_wunif=as.vector(b_allregion_diag_unif), b_sim_wunif=as.vector(1- b_allregion_diag_unif), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(b_dist_all_wunif, aes(x=geo_dist/1000, y= b_sim_wunif)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Weighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
b_dist_all_wunif %>% 
  do({
    mod = lm(b_sim_wunif ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],#截距
               Slope = coef(mod)[2])#斜率
  })
#计算无加权（不考虑物种丰度） UniFrac 距离（框架同上）
b_asv_table_unif <-otu_table(b_merged, taxa_are_rows=T)
b_uwunifrac_essai<-phyloseq::UniFrac(b_merged, weighted=F, normalized=TRUE, parallel=FALSE, fast=TRUE) 
b_allregion_uwnif<-as.matrix(b_uwunifrac_essai)
diag(b_allregion_uwnif)=NA
b_allregion_diag_uwnif<-t(matrix(t(b_allregion_uwnif)[which(!is.na(b_allregion_uwnif))],nrow=47,ncol=48))
mantel(b_allregion_diag_uwnif, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
b_dist_all_uw <-data.frame(b_dis_uw=as.vector(b_allregion_diag_uwnif), b_sim_uw=as.vector(1- b_allregion_diag_uwnif), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(b_dist_all_uw, aes(x=geo_dist/1000, y= b_sim_uw)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - UnWeighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
b_dist_all_uw %>% 
  do({
    mod = lm(b_sim_uw ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
#计算 Sorensen 相异性（基于物种存在/缺失）矩阵，进行全距离衰减（框架同上）
b_vegan_matrix_all <-vegan_otu(b_merged)
b_allregion_sor <-vegdist(b_vegan_matrix_all, binary=T) 
b_allregion_sor_m <-as.matrix(b_allregion_sor)
diag(b_allregion_sor_m)=NA
b_allregion_sor_diag <-t(matrix(t(b_allregion_sor_m)[which(!is.na(b_allregion_sor_m))],nrow=47,ncol=48))
min(b_allregion_sor_diag)
mantel(b_allregion_sor_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
b_dist_all_sor <-data.frame(b_dis_sor=as.vector(b_allregion_sor_diag), b_sim_sor=as.vector(1- b_allregion_sor_diag), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(b_dist_all_sor, aes(x=geo_dist/1000, y= b_sim_sor)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Sorensen") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
b_dist_all_sor %>% 
  do({
    mod = lm(b_sim_sor ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
#组合不同衰减模式
#组合Bray-Curtis和Sorensen
b_melt_bc <- melt(b_allregion_diag) #转换为长格式，Bray-Curtis 相异性矩阵
b_melt_bc$dis <- "BC"
b_melt_sor <- melt(b_allregion_sor_diag)# Sorensen 相异性矩阵
b_melt_sor$dis <- "SOR"
melt_geo <- melt(dist_geo_all_diag)# 地理距离
colnames(melt_geo) <- c("Var1","Var2","dist_geo")
b_bind_bcsor <- rbind(b_melt_bc, b_melt_sor)# 按行合并（堆叠）
b_merge_dis_geo <- merge(b_bind_bcsor, melt_geo)# 根据共同的列来合并
b_bcsor_plot<- ggplot(b_merge_dis_geo) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Scarus_quoyi",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
b_bcsor_plot
#组合有加权和无加权UniFrac 数据
b_melt_unifrac <- melt(b_allregion_diag_unif) 
b_melt_unifrac$dis <- "wunifrac"
b_melt_uwunifrac <- melt(b_allregion_diag_uwnif)
b_melt_uwunifrac$dis <- "unwunifrac"
b_bind_unifrac<- rbind(b_melt_unifrac, b_melt_uwunifrac)
b_merge_dis_unif <- merge(b_bind_unifrac, melt_geo)
b_unifrac_plot<- ggplot(b_merge_dis_unif) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Trimma_lantana",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
b_unifrac_plot
#所有群落相异性（Bray-Curtis, Sorensen, 加权 UniFrac, 无加权 UniFrac）与地理距离组合
b_bind_all_dist <- rbind(b_melt_bc, b_melt_sor, b_melt_uwunifrac, b_melt_unifrac) 
b_merge_alldist_geo <- merge(b_bind_all_dist, melt_geo)
get_fishcol <- fish(4, option="Scarus_quoyi") # 获取鱼类调色板颜色
custom_colors <- c("wunifrac" = "#D53288FF", "BC" = "#3F459BFF", "SOR" = "#009E9EFF", "unwunifrac" = "#AD9C35FF") # 自定义颜色
# 计算相似性和 Log10 转换的地理距离 (千米) 并过滤
b_data_for_log_lm <- b_merge_alldist_geo %>%
  mutate(
    similarity = 1 - value, # 计算相似性
    geo_dist_km = dist_geo / 1000 # 将地理距离转换为千米
  ) %>%
  # 计算 Log10 地理距离，并过滤掉非有限结果 (-Inf, NaN)
  filter(geo_dist_km > 0) %>% # 确保距离大于 0
  mutate(log10_geo_dist_km = log10(geo_dist_km))
# 按距离度量类型 (dis) 分组，并拟合线性模型 (similarity ~ log10_geo_dist_km)
#    并提取每个模型的系数和 R-squared
b_model_summaries <- b_data_for_log_lm %>%
  group_by(dis) %>%
  do({
    # 将 data = cur_data() 修改为 data = .
    model = lm(similarity ~ log10_geo_dist_km, data = .)
    # 使用 tidy 和 glance 函数提取模型信息 (需要 broom 包)
    tidy_output <- broom::tidy(model)
    glance_output <- broom::glance(model)
    # 将系数和 R-squared 组合到一行数据框中
    data.frame(
      Metric = unique(.$dis), # 使用 . 引用当前组的数据框
      Intercept = tidy_output$estimate[1], # 截距
      # 确保 log10_geo_dist_km 是第二个系数
      Slope_log10_km = tidy_output$estimate[grep("log10_geo_dist_km", tidy_output$term)],
      R_squared = glance_output$r.squared # 模型的 R-squared
    )
  })
# 打印结果
print(b_model_summaries)
b_merge_alldist_geo_clean <- b_merge_alldist_geo %>% filter(dist_geo !=0)
alldist_plot_bacteria<- ggplot(b_merge_alldist_geo_clean) + 
  geom_point(aes(y=1-value, x = dist_geo/1000, color= dis), alpha=0.1) + 
  stat_smooth(aes(y=1-value, x = dist_geo/1000, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_manual(name = "",
                     values = custom_colors,
                     breaks = c("wunifrac", "BC", "SOR", "unwunifrac"),
                     labels = c("Weighted UniFrac", "Bray–Curtis", "Sørensen", "Unweighted UniFrac"))+
  ylab("Bacterial community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 16),       # 坐标轴标题大小
        axis.text = element_text(size = 14),        # 坐标轴刻度大小
        legend.text = element_text(size = 14))   # 图例内容大小
alldist_plot_bacteria 
# ggsave("D:/study/master/Main_Figure_tables/Figure_4/4a_distancedecay_bacteria.png", plot = alldist_plot_bacteria, width = 7, height = 6, dpi = 600, bg = "transparent")
b_merge_alldist_geo$dist_geo <- as.numeric(b_merge_alldist_geo$dist_geo)
b_merge_alldist_geo$dist_geo_km <- b_merge_alldist_geo$dist_geo / 1000#地理距离转换为公里
calculate_r_squared <- function(model) {#计算r平方
  return(format(summary(model)$r.squared, digits = 3))
}
#线性回归，群落相似性随地理距离的变化趋势
b_lm_models_tax <- by(b_merge_alldist_geo, b_merge_alldist_geo$dis, function(sub_df_tax) {
  lm(1 - sub_df_tax$value ~ sub_df_tax$dist_geo_km, data = sub_df_tax)
}) #每个距离拟合一个线性模型
b_intercepts <- sapply(b_lm_models_tax, function(model) coef(model)[1])#计算截距
b_slopes <- sapply(b_lm_models_tax, function(model) coef(model)[2])# 斜率
b_r_squared <- sapply(b_lm_models_tax, calculate_r_squared)# R2 值
#subset_size <- 35000#随机抽样
#set.seed(42) 
#subset_merge_alldist <- merge_alldist_geo[sample(nrow(merge_alldist_geo), subset_size), , drop = FALSE]
b_subset_merge_alldist <- b_merge_alldist_geo
# Bray-Curtis 和 Sorensen 差异性模式比较
b_merge_bcsor <- merge(b_dist_all_bray, b_dist_all_sor, by="row.names") 
b_anco_bcsor <- b_merge_bcsor[c("b_BC_sim_bc","b_sim_sor","geo_dist.x")]
b_manco <- melt(b_anco_bcsor, id=c("geo_dist.x"))
# Bray-Curtis 和 Sorensen 相似性随线性地理距离变化的散点图，并在图上直接显示了各自的线性回归方程和 R2 值。
ggscatter(
  b_manco, x = "geo_dist.x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )
b_manco$geo_dist_log <- log((b_manco$geo_dist.x)/1000)
#组合加权和无加权 UniFrac 数据
b_merge_unifrac <- merge(b_dist_all_wunif, b_dist_all_uw, by="row.names")
b_anco_unifrac <- b_merge_unifrac[c("b_sim_wunif","b_sim_uw","geo_dist.x")]
b_manco_unifrac <- melt(b_anco_unifrac, id=c("geo_dist.x")) 
b_manco_unifrac$geo_dist_log <- log((b_manco_unifrac$geo_dist.x)/1000)
#绘制加权与无加权 UniFrac 比较图
ggscatter(
  b_manco_unifrac, x = "geo_dist.x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )
b_manco_unifrac$geo_dist_log <- log((b_manco_unifrac$geo_dist.x)/1000)
#BC/Sorensen方差分析 (ANOVA)
#只保留地理距离 (geo_dist.x) 大于或等于 1000 米（即 1 千米）的数据点
b_manco_clean <- b_manco %>% filter(geo_dist_log >= 0)
b_manco_clean %>% anova_test(sqrt(value) ~ variable * geo_dist_log)
#因变量是原始的 value (相似性)，自变量是对数地理距离 (geo_dist_log) 和距离度量类型 (variable)
b_model_bcsor <- lm((value) ~ geo_dist_log + variable, data = b_manco_clean) #不带交互作用
b_model_metrics_bcsor <- augment(b_model_bcsor) %>%# 检查模型诊断指标
  select(-.hat, -.sigma, -.fitted) #移除 .hat (帽子值), .sigma (残差的标准误差), .fitted (拟合值)
head(b_model_metrics_bcsor, 3)
ad.test(b_model_bcsor$residuals) #Anderson-Darling 正态性检验
#分析 UniFrac 数据，检验加权和无加权 UniFrac 相似性随线性地理距离变化的斜率是否显著不同，这里的 UniFrac 分析使用了原始线性地理距离 geo_dist.x，而 BC/Sorensen 分析使用了对数转换的地理距离 geo_dist_log
b_manco_unifrac %>% anova_test(log10(value) ~variable*geo_dist.x)# 方差分析
#中心化变量并拟合包含交互作用的线性模型
# BC/Sorensen中心化分析
b_manco_clean$geo_dist_log_centered <- datawizard::standardize(b_manco_clean$geo_dist_log, center = TRUE, scale = FALSE) 
b_glm_dd_bcsor <- lm(sqrt(value) ~ geo_dist_log_centered * variable, data = b_manco_clean) # 拟合带交互作用的平方根相似性对中心化对数距离模型
check_model(b_glm_dd_bcsor)
summary(b_glm_dd_bcsor)
#处理 UniFrac 数据（过滤、对数转换、中心化）并拟合模型
b_manco_unifrac_clean <- b_manco_unifrac %>% filter(geo_dist_log >= 0)# 过滤掉近距离点
b_manco_unifrac_clean$geo_dist_log <- log((b_manco_unifrac_clean$geo_dist.x)/1000)
b_manco_unifrac_clean$geo_dist_log_centered <- datawizard::standardize(b_manco_unifrac_clean$geo_dist_log, center = TRUE, scale = FALSE)
b_glm_dd_wufrac <- lm(b_manco_unifrac_clean, formula = value ~ geo_dist_log_centered*variable)
check_model(b_glm_dd_wufrac)
summary(b_glm_dd_wufrac)
#真菌群落多样性距离衰减
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
#转置
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
NOMIS_df <- data.frame(metadata$longitude, metadata$latitude)#经纬度
NOMIS_df$metadata.longitude <- as.numeric(NOMIS_df$metadata.longitude)
NOMIS_df$metadata.latitude <- as.numeric(NOMIS_df$metadata.latitude)
#计算经纬度的地理距离矩阵
dist_geo_all <- distm(NOMIS_df, NOMIS_df, fun=distGeo)
dist_geo_all <- as.matrix(dist_geo_all)
diag(dist_geo_all)=NA#对角线元素设置为NA
#提取所有非对角线（即样本间）的距离，并将它们放入一个矩形矩阵中
dist_geo_all_diag <- t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=47,ncol=48))
min(dist_geo_all_diag)# 最小值
#计算物种丰度的Bray-Curtis相似性指数（框架同上）
f_vegan_matrix_all <-vegan_otu(f_merged)
f_allregion_bray <-vegdist(log1p(f_vegan_matrix_all), method="bray")# Bray-Curtis 距离矩阵
f_allregion_m <- as.matrix(f_allregion_bray)
diag(f_allregion_m)=NA
f_allregion_diag <-t(matrix(t(f_allregion_m)[which(!is.na(f_allregion_m))],nrow=47,ncol=48))
min(f_allregion_diag)
#Mantel检验两个距离矩阵相关性
mantel(f_allregion_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
#组合Bray-Curtis 相异性值、Bray-Curtis 相似性 (1 - 相异性)、地理距离、字符串
f_dist_all_bray <- data.frame(f_BC_dist_bc=as.vector(f_allregion_diag), f_BC_sim_bc=as.vector(1- f_allregion_diag), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(f_dist_all_bray, aes(x=geo_dist/1000, y= f_BC_sim_bc)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Bray-Curtis") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
#提取 Bray-Curtis回归模型的截距和斜率
f_dist_all_bray %>% 
  do({
    # 用 geo_dist (地理距离) 作为解释变量来预测 BC_sim_bc (Bray-Curtis 相似性)
    mod = lm(f_BC_sim_bc ~ geo_dist, data = f_dist_all_bray)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
# 加权 UniFrac（考虑物种丰度信息）距离矩阵（框架同上）
f_asv_table_unif<-otu_table(f_merged, taxa_are_rows=T)
f_unifrac_essai<-phyloseq::UniFrac(f_merged, weighted=T, normalized=TRUE, parallel=FALSE, fast=TRUE) 
f_allregion_unif<-as.matrix(f_unifrac_essai)
diag(f_allregion_unif)=NA
f_allregion_diag_unif<-t(matrix(t(f_allregion_unif)[which(!is.na(f_allregion_unif))],nrow=47,ncol=48))
mantel(f_allregion_diag_unif, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
f_dist_all_wunif <- data.frame(dis_wunif=as.vector(f_allregion_diag_unif), f_sim_wunif=as.vector(1- f_allregion_diag_unif), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(f_dist_all_wunif, aes(x=geo_dist/1000, y= f_sim_wunif)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Weighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
f_dist_all_wunif %>% 
  do({
    mod = lm(f_sim_wunif ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],#截距
               Slope = coef(mod)[2])#斜率
  })
#计算无加权（不考虑物种丰度） UniFrac 距离（框架同上）
f_asv_table_unif <-otu_table(f_merged, taxa_are_rows=T)
f_uwunifrac_essai<-phyloseq::UniFrac(f_merged, weighted=F, normalized=TRUE, parallel=FALSE, fast=TRUE) 
f_allregion_uwnif<-as.matrix(f_uwunifrac_essai)
diag(f_allregion_uwnif)=NA
f_allregion_diag_uwnif<-t(matrix(t(f_allregion_uwnif)[which(!is.na(f_allregion_uwnif))],nrow=47,ncol=48))
mantel(f_allregion_diag_uwnif, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
f_dist_all_uw <-data.frame(f_dis_uw=as.vector(f_allregion_diag_uwnif), f_sim_uw=as.vector(1- f_allregion_diag_uwnif), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(f_dist_all_uw, aes(x=geo_dist/1000, y= f_sim_uw)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - UnWeighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
f_dist_all_uw %>% 
  do({
    mod = lm(f_sim_uw ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
#计算 Sorensen 相异性（基于物种存在/缺失）矩阵，进行全距离衰减（框架同上）
f_vegan_matrix_all <-vegan_otu(f_merged)
f_allregion_sor <-vegdist(f_vegan_matrix_all, binary=T) 
f_allregion_sor_m <-as.matrix(f_allregion_sor)
diag(f_allregion_sor_m)=NA
f_allregion_sor_diag <-t(matrix(t(f_allregion_sor_m)[which(!is.na(f_allregion_sor_m))],nrow=47,ncol=48))
min(f_allregion_sor_diag)
mantel(f_allregion_sor_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
f_dist_all_sor <-data.frame(f_dis_sor=as.vector(f_allregion_sor_diag), f_sim_sor=as.vector(1- f_allregion_sor_diag), geo_dist=as.vector(dist_geo_all_diag), Method="All populations")
ggplot(f_dist_all_sor, aes(x=geo_dist/1000, y= f_sim_sor)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Sorensen") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()
f_dist_all_sor %>% 
  do({
    mod = lm(f_sim_sor ~ geo_dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
#组合不同衰减模式
#组合Bray-Curtis和Sorensen
f_melt_bc <- melt(f_allregion_diag) #转换为长格式，Bray-Curtis 相异性矩阵
f_melt_bc$dis <- "BC"
f_melt_sor <- melt(f_allregion_sor_diag)# Sorensen 相异性矩阵
f_melt_sor$dis <- "SOR"
melt_geo <- melt(dist_geo_all_diag)# 地理距离
colnames(melt_geo) <- c("Var1","Var2","dist_geo")
f_bind_bcsor <- rbind(f_melt_bc, f_melt_sor)# 按行合并（堆叠）
f_merge_dis_geo <- merge(f_bind_bcsor, melt_geo)# 根据共同的列来合并
f_bcsor_plot<- ggplot(f_merge_dis_geo) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Scarus_quoyi",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
f_bcsor_plot
#组合有加权和无加权UniFrac 数据
f_melt_unifrac <- melt(f_allregion_diag_unif) 
f_melt_unifrac$dis <- "wunifrac"
f_melt_uwunifrac <- melt(f_allregion_diag_uwnif)
f_melt_uwunifrac$dis <- "unwunifrac"
f_bind_unifrac<- rbind(f_melt_unifrac, f_melt_uwunifrac)
f_merge_dis_unif <- merge(f_bind_unifrac, melt_geo)
f_unifrac_plot<- ggplot(f_merge_dis_unif) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Trimma_lantana",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
f_unifrac_plot
#所有群落相异性（Bray-Curtis, Sorensen, 加权 UniFrac, 无加权 UniFrac）与地理距离组合
f_bind_all_dist <- rbind(f_melt_bc, f_melt_sor, f_melt_uwunifrac, f_melt_unifrac) 
f_merge_alldist_geo <- merge(f_bind_all_dist, melt_geo)
get_fishcol <- fish(4, option="Scarus_quoyi") # 获取鱼类调色板颜色
custom_colors <- c("wunifrac" = "#D53288FF", "BC" = "#3F459BFF", "SOR" = "#009E9EFF", "unwunifrac" = "#AD9C35FF") # 自定义颜色
# 计算相似性和 Log10 转换的地理距离 (千米) 并过滤
f_data_for_log_lm <- f_merge_alldist_geo %>%
  mutate(
    similarity = 1 - value, # 计算相似性
    geo_dist_km = dist_geo / 1000 # 将地理距离转换为千米
  ) %>%
  # 计算 Log10 地理距离，并过滤掉非有限结果 (-Inf, NaN)
  filter(geo_dist_km > 0) %>% # 确保距离大于 0
  mutate(log10_geo_dist_km = log10(geo_dist_km))
# 按距离度量类型 (dis) 分组，并拟合线性模型 (similarity ~ log10_geo_dist_km)
#    并提取每个模型的系数和 R-squared
f_model_summaries <- f_data_for_log_lm %>%
  group_by(dis) %>%
  do({
    # 将 data = cur_data() 修改为 data = .
    model = lm(similarity ~ log10_geo_dist_km, data = .)
    # 使用 tidy 和 glance 函数提取模型信息 (需要 broom 包)
    tidy_output <- broom::tidy(model)
    glance_output <- broom::glance(model)
    # 将系数和 R-squared 组合到一行数据框中
    data.frame(
      Metric = unique(.$dis), # 使用 . 引用当前组的数据框
      Intercept = tidy_output$estimate[1], # 截距
      # 确保 log10_geo_dist_km 是第二个系数
      Slope_log10_km = tidy_output$estimate[grep("log10_geo_dist_km", tidy_output$term)],
      R_squared = glance_output$r.squared # 模型的 R-squared
    )
  })
# 打印结果
print(f_model_summaries)
f_merge_alldist_geo_clean <- f_merge_alldist_geo %>% filter(dist_geo !=0)
alldist_plot_fungi<- ggplot(f_merge_alldist_geo_clean) + 
  geom_point(aes(y=1-value, x = dist_geo/1000, color= dis), alpha=0.1) + 
  stat_smooth(aes(y=1-value, x = dist_geo/1000, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_manual(name = "",
                     values = custom_colors,
                     breaks = c("wunifrac", "BC", "SOR", "unwunifrac"),
                     labels = c("Weighted UniFrac", "Bray–Curtis", "Sørensen", "Unweighted UniFrac"))+
  ylab("Fungal community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 16),       # 坐标轴标题大小
        axis.text = element_text(size = 14),        # 坐标轴刻度大小
        legend.text = element_text(size = 14))   # 图例内容大小
alldist_plot_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_4/4a_distancedecay_fungi.png", plot = alldist_plot_fungi, width = 7, height = 6, dpi = 600, bg = "transparent")
f_merge_alldist_geo$dist_geo <- as.numeric(f_merge_alldist_geo$dist_geo)
f_merge_alldist_geo$dist_geo_km <- f_merge_alldist_geo$dist_geo / 1000#地理距离转换为公里
calculate_r_squared <- function(model) {#计算r平方
  return(format(summary(model)$r.squared, digits = 3))
}
#线性回归，群落相似性随地理距离的变化趋势
f_lm_models_tax <- by(f_merge_alldist_geo, f_merge_alldist_geo$dis, function(suf_df_tax) {
  lm(1 - suf_df_tax$value ~ suf_df_tax$dist_geo_km, data = suf_df_tax)
}) #每个距离拟合一个线性模型
f_intercepts <- sapply(f_lm_models_tax, function(model) coef(model)[1])#计算截距
f_slopes <- sapply(f_lm_models_tax, function(model) coef(model)[2])# 斜率
f_r_squared <- sapply(f_lm_models_tax, calculate_r_squared)# R2 值
#subset_size <- 35000#随机抽样
#set.seed(42) 
#subset_merge_alldist <- merge_alldist_geo[sample(nrow(merge_alldist_geo), subset_size), , drop = FALSE]
f_subset_merge_alldist <- f_merge_alldist_geo
# Bray-Curtis 和 Sorensen 差异性模式比较
f_merge_bcsor <- merge(f_dist_all_bray, f_dist_all_sor, by="row.names") 
f_anco_bcsor <- f_merge_bcsor[c("f_BC_sim_bc","f_sim_sor","geo_dist.x")]
f_manco <- melt(f_anco_bcsor, id=c("geo_dist.x"))
# Bray-Curtis 和 Sorensen 相似性随线性地理距离变化的散点图，并在图上直接显示了各自的线性回归方程和 R2 值。
ggscatter(
  f_manco, x = "geo_dist.x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )
f_manco$geo_dist_log <- log((f_manco$geo_dist.x)/1000)
#组合加权和无加权 UniFrac 数据
f_merge_unifrac <- merge(f_dist_all_wunif, f_dist_all_uw, by="row.names")
f_anco_unifrac <- f_merge_unifrac[c("f_sim_wunif","f_sim_uw","geo_dist.x")]
f_manco_unifrac <- melt(f_anco_unifrac, id=c("geo_dist.x")) 
f_manco_unifrac$geo_dist_log <- log((f_manco_unifrac$geo_dist.x)/1000)
#绘制加权与无加权 UniFrac 比较图
ggscatter(
  f_manco_unifrac, x = "geo_dist.x", y = "value",
  color = "variable", add = "reg.line"
)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
  )
f_manco_unifrac$geo_dist_log <- log((f_manco_unifrac$geo_dist.x)/1000)
#BC/Sorensen方差分析 (ANOVA)
#只保留地理距离 (geo_dist.x) 大于或等于 1000 米（即 1 千米）的数据点
f_manco_clean <- f_manco %>% filter(geo_dist_log >= 0)
f_manco_clean %>% anova_test(sqrt(value) ~ variable * geo_dist_log)
#因变量是原始的 value (相似性)，自变量是对数地理距离 (geo_dist_log) 和距离度量类型 (variable)
f_model_bcsor <- lm((value) ~ geo_dist_log + variable, data = f_manco_clean) #不带交互作用
f_model_metrics_bcsor <- augment(f_model_bcsor) %>%# 检查模型诊断指标
  select(-.hat, -.sigma, -.fitted) #移除 .hat (帽子值), .sigma (残差的标准误差), .fitted (拟合值)
head(f_model_metrics_bcsor, 3)
ad.test(f_model_bcsor$residuals) #Anderson-Darling 正态性检验
#分析 UniFrac 数据，检验加权和无加权 UniFrac 相似性随线性地理距离变化的斜率是否显著不同，这里的 UniFrac 分析使用了原始线性地理距离 geo_dist.x，而 BC/Sorensen 分析使用了对数转换的地理距离 geo_dist_log
f_manco_unifrac_corrected <- f_manco_unifrac %>%
  mutate(value_corrected = ifelse(value <= 0, min(value[value > 0], na.rm = TRUE) / 10, value))
f_manco_unifrac_corrected %>% anova_test(log10(value_corrected) ~ variable * geo_dist.x)# 方差分析
#中心化变量并拟合包含交互作用的线性模型
# BC/Sorensen中心化分析
f_manco_clean$geo_dist_log_centered <- datawizard::standardize(f_manco_clean$geo_dist_log, center = TRUE, scale = FALSE) 
f_glm_dd_bcsor <- lm(sqrt(value) ~ geo_dist_log_centered * variable, data = f_manco_clean) # 拟合带交互作用的平方根相似性对中心化对数距离模型
check_model(f_glm_dd_bcsor)
summary(f_glm_dd_bcsor)
#处理 UniFrac 数据（过滤、对数转换、中心化）并拟合模型
f_manco_unifrac_clean <- f_manco_unifrac %>% filter(geo_dist_log >= 0)# 过滤掉近距离点
f_manco_unifrac_clean$geo_dist_log <- log((f_manco_unifrac_clean$geo_dist.x)/1000)
f_manco_unifrac_clean$geo_dist_log_centered <- datawizard::standardize(f_manco_unifrac_clean$geo_dist_log, center = TRUE, scale = FALSE)
f_glm_dd_wufrac <- lm(f_manco_unifrac_clean, formula = value ~ geo_dist_log_centered*variable)
check_model(f_glm_dd_wufrac)
summary(f_glm_dd_wufrac)

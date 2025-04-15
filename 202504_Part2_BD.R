#图2d β多样性
install.packages(c("phyloseq", "vegan", "ggplot2", "ecodist"))
install.packages("mvabund")
install.packages("tibble")
install.packages('devtools')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("remotes")
remotes::install_github("jbisanz/phyloseqCompanion") 
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(phyloseq)
library(phyloseqCompanion)
library(vegan)
library(ggplot2)
library(ecodist)
library(tibble)
library(devtools)
library(pairwiseAdonis)
library(mvabund)
#细菌贝塔
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
# 提取分类表
b_tax <- tax_table(b_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(b_merged))
#把 phyloseq 中的 OTU 表标准化为一个行为样本、列为 ASV （互换）的 matrix 格式
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
# 标准化单列数据（例如 latitude 列）
metadata_nmds <- metadata
metadata_nmds$Niche <- as.factor(metadata_nmds$Niche)
metadata_nmds$ele_sp <- scale(metadata_nmds$elevation)
metadata_nmds$ele_attribute <- "A"
metadata_nmds$ele_attribute <- ifelse(metadata_nmds$ele_sp < 1, "B", "A")
# metadata_nmds$ele_attribute <- ifelse((metadata_nmds$ele_sp) < 1, metadata_nmds$ele_attribute == "B", metadata_nmds$ele_attribute == "A")
metadata_nmds$ele_attribute <- as.factor(metadata_nmds$ele_attribute)
## 画NMDS图
b_asv_table_nmds <- as.matrix((otu_table(b_merged, taxa_are_rows=T)))
b_asv_table_nmds_f <- b_asv_table_nmds[rowSums(b_asv_table_nmds[])>0,]
b_nmds_bc_nomis <- metaMDS(t(log1p(b_asv_table_nmds_f)), distance = "bray", k = 2, trymax=999)# 运行 NMDS 分析
stressplot(b_nmds_bc_nomis) # stressplot 图，评估拟合质量：点落在对角线， NMDS 表现好
b_nmds_bc_nomis$stress# stress 越小越好（通常 < 0.2 可接受，< 0.1 很好）
b_data.scores = as.data.frame(scores(b_nmds_bc_nomis)$sites)# 提取每个样本的 NMDS 坐标
b_data.scores$Sample = metadata_nmds$Sample.ID#给每个点加上对应的样本 ID
b_data.scores$Group = metadata_nmds$Group
b_data.scores$Niche = metadata_nmds$Niche
head(b_data.scores)
colors_group <- c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                  "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
colors_niche <- c("#2E2A2BFF", "#CF4E9CFF")
group_levels <- c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN")
niche_levels <- c("G", "N")
# 确保 Group 是因子，且顺序正确
b_data.scores$Group <- factor(b_data.scores$Group, levels = group_levels)
b_data.scores$Niche <- factor(b_data.scores$Niche, levels = niche_levels)
nmds_bc_plot_bacteria = ggplot(b_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(fill = Group), shape = 21, color = "black" , alpha = 0.8 , show.legend = FALSE)+ 
  stat_ellipse(aes(x=NMDS1, y=NMDS2,color= Niche),type = "norm" , size = 1 , alpha = 0.8)+#椭圆轮廓
  scale_fill_manual(name = "Region", values = colors_group) +  #点
  scale_color_manual(# 椭圆轮廓线
    values = colors_niche,
    labels = c("Rhizosphere Soil", "Bulb")  # 替换图例文字
  ) +
  guides(
    color = guide_legend(
      title = NULL,                   # 去掉标题
      override.aes = list(size = 2.5) # 延长图例线条
    )
  ) +  
  theme(
    axis.text.y = element_text(colour = "black", size = 14), # 坐标轴刻度标签（数字）
    axis.text.x = element_text(colour = "black", size = 14), 
    legend.text = element_text(size = 14, colour ="black"), # 图例文字
    legend.position = c(0.01, 1),  # 图例放在左上角（相对坐标）
    legend.justification = c(0, 1),   # 对齐方式：左上
    legend.key.width = unit(1.8, "cm"), # 控制横线长度
    axis.title.y = element_text(size = 16), 
    axis.title.x = element_text(size = 16, colour = "black"), 
    legend.title = element_text(size = 14, colour = "black", face = "bold"), 
    panel.background = element_blank(), 
    panel.border = element_blank(),  # 去掉默认边框
    axis.line.x.bottom = element_line(color = "black"),  # 仅下边框
    axis.line.y.left   = element_line(color = "black"),  # 仅左边框
    axis.line.x.top    = element_blank(),  # 去掉上边框
    axis.line.y.right  = element_blank(),  # 去掉右边框
    legend.key = element_blank()
  ) + 
  labs(x = "Bacterial NMDS1", y = "Bacterial NMDS2", shape = "Type")
nmds_bc_plot_bacteria
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2d_beta_bacteria.png", plot = nmds_bc_plot_bacteria, width = 8, height = 6, dpi = 600, bg = "transparent")
#包括平滑线来显示纬度
# 生成一个空图框，type = "n" 不画出任何点
plot(x=b_data.scores$NMDS1, y= b_data.scores$NMDS2, type="n", xlim = c(-2.5, 2.9), ylim = c(-1.6, 1.8)) 
# 添加样本点
points(b_nmds_bc_nomis, display = "sites", cex = 2.3, pch=19, col=alpha(colors_group[factor(b_data.scores$Group)], 0.8)) 
# 添加变量的等高线曲面
ordisurf(b_nmds_bc_nomis, metadata_nmds$ele_sp, add = TRUE, col="blue", labcex=1) 
#群落组间变异度分析
b_vegan_matrix<- vegan_otu(b_merged)
#计算 Bray-Curtis 距离矩阵
b_bray <- vegdist(log1p(b_vegan_matrix), method="bray") 
#分析流程：组内、组间、组间成对比较。具体为，betadisper（单因素、置换检验）检查组内变异是否一致；若不显著，继续做 PERMANOVA，adonis2 检查 Group 或海拔等变量是否显著影响群落组成，后可进一步细分pairwise.adonis 做组间成对比较，找出差异具体在哪些组；若显著，继续做 多响应变量广义线性模型。
#检查不同 Group 之间 beta 多样性是否差异显著
b_bdisp_nomis<- betadisper(b_bray, metadata_nmds$Group, type=c("centroid"))
b_bdisp_nomis
#单因素方差分析（ANOVA）来检验不同 Group 之间 beta dispersion 是否存在显著差异
b_aov_bdisp <-anova(b_bdisp_nomis)
#置换检验 + 成对比较
permutest(b_bdisp_nomis, pairwise=T)
#用纬度（latitude）做分组变量分析
b_bdisp_nomis_lat<- betadisper(b_bray, metadata_nmds$latitude, type=c("centroid"))
b_bdisp_nomis_lat
b_aov_bdisp_lat <-anova(b_bdisp_nomis_lat)
permutest(b_bdisp_nomis_lat, pairwise=T)
#如果不显著，则使用 adonis 和成对 adonis，并在稿件中报告这些不同的值
b_ado <- adonis2(b_bray ~ Group, permutations = 999, method = "bray", data=metadata_nmds) #分组
b_ado_latitude <- adonis2(b_bray ~ ele_attribute, permutations = 999, method = "bray", data=metadata_nmds)
pairwise_rb <- pairwise.adonis(b_bray, metadata_nmds$Group, p.adjust.m="holm")#分组
pairwise_b_lat <- pairwise.adonis(b_bray, metadata_nmds$ele_attribute, p.adjust.m="holm")
#用 mvabund 包 来进行 群落组成差异分析
b_asv_pu <- t(otu_table(b_merged, taxa_are_rows=T))
b_ab <- mvabund(b_asv_pu)#转换为 mvabund 对象
b_asv_nb_group <- manyglm(b_ab ~ Group,#多响应变量广义线性模型，分组
                          data = metadata_nmds, family = 'negative binomial') 
b_nomis_avo_group <- anova(b_asv_nb_group, p.uni= "adjusted", nBoot = 99, pairwise.comp=metadata_nmds$Group, show.time=T) #进行置换检验
b_nomis_avo_group$pairwise.comp.table#组间两两比较
#b_nomis_manyglm_res_group <- b_nomis_avo_group$uni.p#提取每个 ASV 的显著性结果
#write.csv(b_nomis_avo_group$pairwise.comp.table, "D:/study/master/Main_Figure_tables/Figure_2/2d_b_results_manyglm_nomis_group.csv")
b_asv_nb_niche <- manyglm(b_ab ~ Niche,#生态位
                          data = metadata_nmds, family = 'negative binomial') 
b_nomis_avo_niche <- anova(b_asv_nb_niche, p.uni= "adjusted", nBoot = 99, show.time=T) #进行置换检验
b_nomis_avo_niche$table#组间两两比较
#b_nomis_manyglm_res_niche <- b_nomis_avo_niche$uni.p#提取每个 ASV 的显著性结果
#write.csv(b_nomis_avo_niche$table, "D:/study/master/Main_Figure_tables/Figure_2/2d_b_results_manyglm_nomis_niche.csv")
#真菌贝塔
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
#去除非细菌
f_merged <- subset_taxa(f_merged, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
f_merged <- subset_taxa(f_merged, (Kingdom!="d__Archaea") | is.na(Kingdom))
f_merged <- subset_taxa(f_merged, (Order!="o__Chloroplast") )
f_merged <- subset_taxa(f_merged, (Family!="f__Mitochondria"))
f_merged <- subset_taxa(f_merged, (Family!="NA"))
# 提取 ASV 表
f_ASV <- as.matrix(otu_table(f_merged, taxa_are_rows=T))
f_ASV_df <- f_ASV[rowSums(f_ASV[])>0,]
# 提取分类表
f_tax <- tax_table(f_merged)
# 提取样本元数据
metadata <- as_tibble(sample_data(f_merged))
#把 phyloseq 中的 OTU 表标准化为一个行为样本、列为 ASV （互换）的 matrix 格式
vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
# 标准化单列数据（例如 latitude 列）
metadata_nmds <- metadata
metadata_nmds$Niche <- as.factor(metadata_nmds$Niche)
metadata_nmds$ele_sp <- scale(metadata_nmds$elevation)
metadata_nmds$ele_attribute <- "A"
metadata_nmds$ele_attribute <- ifelse(metadata_nmds$ele_sp < 1, "B", "A")
# metadata_nmds$ele_attribute <- ifelse((metadata_nmds$ele_sp) < 1, metadata_nmds$ele_attribute == "B", metadata_nmds$ele_attribute == "A")
metadata_nmds$ele_attribute <- as.factor(metadata_nmds$ele_attribute)
## 画NMDS图
f_asv_table_nmds <- as.matrix((otu_table(f_merged, taxa_are_rows=T)))
f_asv_table_nmds_f <- f_asv_table_nmds[rowSums(f_asv_table_nmds[])>0,]
f_nmds_bc_nomis <- metaMDS(t(log1p(f_asv_table_nmds_f)), distance = "bray", k = 2, trymax=999)# 运行 NMDS 分析
stressplot(f_nmds_bc_nomis) # stressplot 图，评估拟合质量：点落在对角线， NMDS 表现好
f_nmds_bc_nomis$stress# stress 越小越好（通常 < 0.2 可接受，< 0.1 很好）
f_data.scores = as.data.frame(scores(f_nmds_bc_nomis)$sites)# 提取每个样本的 NMDS 坐标
f_data.scores$Sample = metadata_nmds$Sample.ID#给每个点加上对应的样本 ID
f_data.scores$Group = metadata_nmds$Group
f_data.scores$Niche = metadata_nmds$Niche
head(f_data.scores)
colors_group <- c("#8C57A2FF","#3EBCB6","#82581FFF","#2F509EFF",
                  "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
colors_niche <- c("#2E2A2BFF", "#CF4E9CFF")
group_levels <- c("JRG", "JJG", "TZG", "PAG", "JRN", "JJN", "TZN", "PAN")
niche_levels <- c("G", "N")
# 确保 Group 是因子，且顺序正确
f_data.scores$Group <- factor(f_data.scores$Group, levels = group_levels)
f_data.scores$Niche <- factor(f_data.scores$Niche, levels = niche_levels)
nmds_bc_plot_fungi = ggplot(f_data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(fill = Group), shape = 21, color = "black" , alpha = 0.8 , show.legend = FALSE)+ 
  stat_ellipse(aes(x=NMDS1, y=NMDS2,color= Niche),type = "norm" , size = 1 , alpha = 0.8)+#椭圆轮廓
  scale_fill_manual(name = "Region", values = colors_group) +  #点
  scale_color_manual(# 椭圆轮廓线
    values = colors_niche,
    labels = c("Rhizosphere Soil", "Bulb")  # 替换图例文字
  ) +
  guides(
    color = guide_legend(
      title = NULL,                   # 去掉标题
      override.aes = list(size = 2.5) # 延长图例线条
    )
  ) +  
  theme(
    axis.text.y = element_text(colour = "black", size = 14), # 坐标轴刻度标签（数字）
    axis.text.x = element_text(colour = "black", size = 14), 
    legend.text = element_text(size = 14, colour ="black"), # 图例文字
    legend.position = c(0.01, 1),  # 图例放在左上角（相对坐标）
    legend.justification = c(0, 1),   # 对齐方式：左上
    legend.key.width = unit(1.8, "cm"), # 控制横线长度
    axis.title.y = element_text(size = 16), 
    axis.title.x = element_text(size = 16, colour = "black"), 
    legend.title = element_text(size = 14, colour = "black", face = "bold"), 
    panel.background = element_blank(), 
    panel.border = element_blank(),  # 去掉默认边框
    axis.line.x.bottom = element_line(color = "black"),  # 仅下边框
    axis.line.y.left   = element_line(color = "black"),  # 仅左边框
    axis.line.x.top    = element_blank(),  # 去掉上边框
    axis.line.y.right  = element_blank(),  # 去掉右边框
    legend.key = element_blank()
  ) + 
  labs(x = "Fungal NMDS1", y = "Fungal NMDS2", shape = "Type")
nmds_bc_plot_fungi
# ggsave("D:/study/master/Main_Figure_tables/Figure_2/2d_beta_fungi.png", plot = nmds_bc_plot_fungi, width = 8, height = 6, dpi = 600, bg = "transparent")
#包括平滑线来显示纬度
# 生成一个空图框，type = "n" 不画出任何点
plot(x=f_data.scores$NMDS1, y= f_data.scores$NMDS2, type="n", xlim = c(-2.5, 2.9), ylim = c(-1.6, 1.8)) 
# 添加样本点
points(f_nmds_bc_nomis, display = "sites", cex = 2.3, pch=19, col=alpha(colors_group[factor(f_data.scores$Group)], 0.8)) 
# 添加变量的等高线曲面
ordisurf(f_nmds_bc_nomis, metadata_nmds$ele_sp, add = TRUE, col="blue", labcex=1) 
#群落组间变异度分析
f_vegan_matrix<- vegan_otu(f_merged)
#计算 Bray-Curtis 距离矩阵
f_bray <- vegdist(log1p(f_vegan_matrix), method="bray") 
#分析流程：组内、组间、组间成对比较。具体为，betadisper（单因素、置换检验）检查组内变异是否一致；若不显著，继续做 PERMANOVA，adonis2 检查 Group 或海拔等变量是否显著影响群落组成，后可进一步细分pairwise.adonis 做组间成对比较，找出差异具体在哪些组；若显著，继续做 多响应变量广义线性模型。
#检查不同 Group 之间 beta 多样性是否差异显著
f_bdisp_nomis<- betadisper(f_bray, metadata_nmds$Group, type=c("centroid"))
f_bdisp_nomis
#单因素方差分析（ANOVA）来检验不同 Group 之间 beta dispersion 是否存在显著差异
f_aov_bdisp <-anova(f_bdisp_nomis)
#置换检验 + 成对比较
permutest(f_bdisp_nomis, pairwise=T)
#用纬度（latitude）做分组变量分析
f_bdisp_nomis_lat<- betadisper(f_bray, metadata_nmds$latitude, type=c("centroid"))
f_bdisp_nomis_lat
f_aov_bdisp_lat <-anova(f_bdisp_nomis_lat)
permutest(f_bdisp_nomis_lat, pairwise=T)
#如果不显著，则使用 adonis 和成对 adonis，并在稿件中报告这些不同的值
f_ado <- adonis2(f_bray ~ Group, permutations = 999, method = "bray", data=metadata_nmds) #分组
f_ado_latitude <- adonis2(f_bray ~ ele_attribute, permutations = 999, method = "bray", data=metadata_nmds)
pairwise_rb <- pairwise.adonis(f_bray, metadata_nmds$Group, p.adjust.m="holm")#分组
pairwise_f_lat <- pairwise.adonis(f_bray, metadata_nmds$ele_attribute, p.adjust.m="holm")
#用 mvabund 包 来进行 群落组成差异分析
f_asv_pu <- t(otu_table(f_merged, taxa_are_rows=T))
f_ab <- mvabund(f_asv_pu)#转换为 mvabund 对象
f_asv_nb_group <- manyglm(f_ab ~ Group,#多响应变量广义线性模型，分组
                          data = metadata_nmds, family = 'negative binomial') 
f_nomis_avo_group <- anova(f_asv_nb_group, p.uni= "adjusted", nBoot = 99, pairwise.comp=metadata_nmds$Group, show.time=T) #进行置换检验
f_nomis_avo_group$pairwise.comp.table#组间两两比较
#f_nomis_manyglm_res_group <- f_nomis_avo_group$uni.p#提取每个 ASV 的显著性结果
#write.csv(f_nomis_avo_group$pairwise.comp.table, "D:/study/master/Main_Figure_tables/Figure_2/2d_f_results_manyglm_nomis_group.csv")
f_asv_nb_niche <- manyglm(f_ab ~ Niche,#生态位
                          data = metadata_nmds, family = 'negative binomial') 
f_nomis_avo_niche <- anova(f_asv_nb_niche, p.uni= "adjusted", nBoot = 99, show.time=T) #进行置换检验
f_nomis_avo_niche$table#组间两两比较
#f_nomis_manyglm_res_niche <- f_nomis_avo_niche$uni.p#提取每个 ASV 的显著性结果
#write.csv(f_nomis_avo_niche$table, "D:/study/master/Main_Figure_tables/Figure_2/2d_f_results_manyglm_nomis_niche.csv")

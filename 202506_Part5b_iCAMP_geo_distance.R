#图5b 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")
install.packages(c("tidyverse", "readr", "readxl", "data.table", "phyloseq",
                   "qiime2R", "ape", "magrittr", "dplyr", "tidyr", "iCAMP",
                   "phytools", "ggplot2", "ggnewscale", "castor", "viridis",
                   "phangorn", "geosphere"))
library(tidyverse) 
library(readr) 
library(readxl) 
library(data.table) 
library(phyloseq)  
library(qiime2R)  
library(ape)  
library(magrittr)  
library(dplyr)
library(tidyr) 
library(iCAMP)
library(ggtree)
library(ggtreeExtra)
library(phytools)
library(ggplot2)
library(ggnewscale)
library(castor)
library(viridis)
library(phangorn)
library(geosphere)
#细菌icamp距离衰减
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
# 计算每个 ASV 的总丰度
b_asv_sums <- taxa_sums(b_merged)
# 计算每个 ASV 在多少个样本中存在
# (这里通过计算每个ASV有多少个非零计数来得到)
b_asv_prevalence <- apply(otu_table(b_merged), 1, function(x) sum(x > 0))
# 根据设定的阈值筛选 ASV
# 创建一个逻辑向量，TRUE 表示 ASV 满足筛选条件
b_keep_taxa <- (b_asv_sums > 3) & (b_asv_prevalence > 3)
# 使用 subset_taxa 函数对 phyloseq 对象进行筛选
b_merged <- phyloseq::subset_taxa(b_merged, b_keep_taxa)
# 从 phyloseq 对象中提取 ASV 表和样本数据
b_ASV_all <- as.data.frame(otu_table(b_merged)) # 确保 ASV 在行，样本在列
metadata <- data.frame(sample_data(b_merged)) # 样本 x 元数据
# 筛选分类信息和系统发育树
# 确保你原始的 `b_tax ` (ASV x Taxonomy) 和 `b_phylo_tree` 是可用的
b_filtered_ASV_ids_grouped <- rownames(b_ASV_all[rowSums(b_ASV_all) > 0, ])
b_filtered_tax_grouped <- b_tax[rownames(b_tax) %in% b_filtered_ASV_ids_grouped, ]
b_filtered_tree_grouped <- drop.tip(b_phylo_tree,
                                    setdiff(b_phylo_tree$tip.label, b_filtered_ASV_ids_grouped))
# 设置一个新的输出目录用于保存分组后的数据
output_wd_all <- "D:/study/master/iCAMP/bacteria_all_data" # 新的输出路径
if (!dir.exists(output_wd_all)) {
  dir.create(output_wd_all, recursive = TRUE)
}
# 保存分组后的 ASV 表格
write.table(b_ASV_all, file = paste0(output_wd_all, "/Asv_table_all.txt"), sep = "\t", quote = FALSE)
# 保存筛选后的分类信息
write.table(b_filtered_tax_grouped, file = paste0(output_wd_all, "/classification_all.txt"), sep = "\t", quote = FALSE)
# 保存筛选后的系统发育树
write.tree(b_filtered_tree_grouped, file = paste0(output_wd_all, "/tree_all.nwk"))
# 设置统一的输入和输出路径
input_output_wd <- output_wd_all # 现在输入就是分组后的数据目录
# 设置 iCAMP 分析参数
prefix <- "all_bacteria" # 使用新的前缀
# 设置输入文件路径
com.file <- paste0(input_output_wd, "/Asv_table_all.txt")
tree.file <- paste0(input_output_wd, "/tree_all.nwk")
clas.file <- paste0(input_output_wd, "/classification_all.txt")
# 设置输出文件夹路径 (与输入路径相同，或在其下创建子文件夹)
save.wd <- paste0(input_output_wd, "/iCAMP_results_all_data") # 为结果创建一个子文件夹
if (!dir.exists(save.wd)) {
  dir.create(save.wd, recursive = TRUE)
}
#  设置 iCAMP 分析参数 
prefix <- "all_bacteria" # 使用统一的前缀表示对所有数据进行分析
rand.time <- 1000
nworker <- 48
memory.G <- 200 # 确保内存设置足够大，以处理整个数据集
#  读取输入数据 
# comm 矩阵通常需要是 (samples x ASVs) 形式
comm <- t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
tree <- read.tree(file = tree.file)
clas <- read.table(clas.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#  匹配物种ID 
# 确保 OTU 表、分类信息和系统发育树的物种 ID 一致
spid.check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tree = tree))
comm <- spid.check$comm
clas <- spid.check$clas
tree <- spid.check$tree
#  计算系统发育距离矩阵 
# wd 参数指向 pdist.big 存储中间文件的目录，这里可以直接使用 save.wd
pd.big <- iCAMP::pdist.big(tree = tree, wd = save.wd, nworker = nworker, memory.G = memory.G)
#  执行 iCAMP 分析 
# iCAMP 分析参数设置
bin.size.limit <- 48
sig.index <- "Confidence"
# 执行 iCAMP 分析，基于丰度和系统发育距离分析群落构建
icres <- iCAMP::icamp.big(
  comm = comm,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  rand = rand.time,
  tree = tree,
  prefix = prefix,
  ds = 0.2,
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "both",
  sig.index = sig.index,
  bin.size.limit = bin.size.limit,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = FALSE,
  ignore.zero = TRUE,
  output.wd = save.wd,
  correct.special = TRUE,
  unit.sum = rowSums(comm), # 针对所有样本总丰度计算 unit.sum
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)
#  保存 iCAMP 分析结果 
# 保存物种的系统发育 bin 和基于 Bray-Curtis 的分析结果
bins <- icres$detail$taxabin$sp.bin
res1 <- icres$CbMPDiCBraya
write.csv(res1, paste0(save.wd, "/", prefix, ".res.csv"), row.names = FALSE)
write.csv(bins, paste0(save.wd, "/", prefix, ".bins.csv"), row.names = FALSE)
# 按系统发育 bin 汇总 iCAMP 结果
icbin <- iCAMP::icamp.bins(
  icamp.detail = icres$detail,
  clas = clas,
  silent = FALSE,
  boot = FALSE,
  rand.time = rand.time,
  between.group = FALSE # 通常在不分组时设为 FALSE
)
# 保存 iCAMP 汇总结果
save(icbin, file = paste0(save.wd, "/", prefix, ".iCAMP.Summary.rda"))
# 保存每个组、每个 bin 和每次 turnover 的生态过程重要性
# 注意：在不分组分析中，icbin$Pt 可能只包含一个整体的 Pt，或者 Ptuv 会包含所有样本对
write.csv(icbin$Pt, file = paste0(save.wd, "/", prefix, ".ProcessImportance_Overall.csv"), row.names = FALSE)
# icbin$Ptk 和 icbin$Ptuv 会包含所有样本对的信息，这里根据实际输出命名
write.csv(icbin$Ptk, file = paste0(save.wd, "/", prefix, ".ProcessImportance_EachBin_Overall.csv"), row.names = FALSE)
write.csv(icbin$Ptuv, file = paste0(save.wd, "/", prefix, ".ProcessImportance_EachTurnover_Overall.csv"), row.names = FALSE)
# 保存每个组中每个 bin 对每个过程的贡献
write.csv(icbin$BPtk, file = paste0(save.wd, "/", prefix, ".BinContributeToProcess_Overall.csv"), row.names = FALSE)
# 保存每个 bin 中的物种及其最高分类级别
write.csv(data.frame(ID = rownames(icbin$Class.Bin), icbin$Class.Bin, stringsAsFactors = FALSE),
          file = paste0(save.wd, "/", prefix, ".Taxon_Bin.csv"), row.names = FALSE)
# 保存每个 bin 的最高分类级别
write.csv(icbin$Bin.TopClass, file = paste0(save.wd, "/", prefix, ".Bin_TopTaxon.csv"), row.names = FALSE) 
# 读取包含iCAMP结果的CSV文件
b_iCAMP <- read.csv("D:/study/master/iCAMP/bacteria_all_data/iCAMP_results_all_data/all_bacteria.ProcessImportance_EachTurnover_Overall.csv", head=T)
b_iCAMP <- b_iCAMP %>% filter(Method == "CbMPDiCbraya")
b_iCAMP$cut <- paste(b_iCAMP$samp1, b_iCAMP$samp2, sep = "")
#经纬度
NOMIS_df <- data.frame(
  longitude = as.numeric(metadata$longitude),
  latitude = as.numeric(metadata$latitude),
  Sample = rownames(metadata)
)
#计算经纬度的地理距离矩阵
# 从坐标数据框中提取组名和坐标
sample_names <- NOMIS_df$Sample
sample_longitudes <- NOMIS_df$longitude
sample_latitudes <- NOMIS_df$latitude
# 创建一个空的数据框，用于存储每对组的距离
pairwise_geo_distances <- data.frame(
  Sample1 = character(),
  Sample2 = character(),
  Geo_Distance_Meters = numeric(), # 距离单位是米
  stringsAsFactors = FALSE
)
# 循环遍历所有唯一的组对并计算距离
num_samples <- length(sample_names)
for (i in 1:(num_samples - 1)) {
  for (j in (i + 1):num_samples) {
    # 获取第一个组的坐标
    sample1_name <- sample_names[i]
    coords1 <- c(sample_longitudes[i], sample_latitudes[i])
    # 获取第二个组的坐标
    sample2_name <- sample_names[j]
    coords2 <- c(sample_longitudes[j], sample_latitudes[j])
    # 计算地理距离
    distance <- geosphere::distGeo(p1 = coords1, p2 = coords2)
    # 将结果添加到数据框中
    pairwise_geo_distances <- rbind(pairwise_geo_distances, data.frame(
      Sample1 = sample1_name,
      Sample2 = sample2_name,
      Geo_Distance_Meters = distance
    ))
  }
}
pairwise_geo_distances$cut <- paste(pairwise_geo_distances$Sample1, pairwise_geo_distances$Sample2, sep = "")
b_iCAMP_geo <- b_iCAMP %>%
  # 使用 left_join 合并，以 'cut' 列为共同键
  left_join(pairwise_geo_distances, by = "cut") %>%
  # 选择你最终希望保留的列
  dplyr::select(
    samp1,
    samp2,
    HeS,
    HoS,
    DL,
    HD,
    DR,
    cut,
    Geo_Distance_Meters
  )
b_iCAMP_geo <- b_iCAMP_geo %>%
  filter(Geo_Distance_Meters != 0)
# 创建ggplot对象
p_bacteria <- ggplot() +
  # 添加地理距离与HoS散点图
  geom_point(data= b_iCAMP_geo, aes(Geo_Distance_Meters/1000, HoS*100, color = "HoS"), col=rgb(0,0,1,0.25), cex=2) +
  # 添加地理距离与HoS线性拟合平滑曲线
  geom_smooth(data=b_iCAMP_geo, aes(Geo_Distance_Meters/1000, HoS*100, color = "HoS"), method="lm") +
  # 添加地理距离与DR散点图
  geom_point(data=b_iCAMP_geo, aes(Geo_Distance_Meters/1000, DR*100, color = "DR"), col=rgb(1,0,0,0.25), cex=2) +
  # 添加地理距离与DR线性拟合平滑曲线
  geom_smooth(data=b_iCAMP_geo, aes(Geo_Distance_Meters/1000, DR*100, color = "DR"), method="lm") +
  # 手动设置颜色和图例标签
  scale_color_manual(
    values = c("HoS" = "blue", "DR" = "red"), # 定义 HoS 为蓝色，DR 为红色
    labels = c("HoS", "DR") # 定义图例显示为 "HoS" 和 "DR"
  ) +
  # 使用白色背景主题
  theme_bw() +
  # 对x轴进行log10转换
  scale_x_continuous(trans='log10') +
  # 设置坐标轴标签
  labs(x="Geographic distance (km)", y="Bacterial contribution (%)")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),            #移除图例标题
        legend.position = c(0.05, 0.95),           #设置图例位置 (x, y 坐标)
        legend.justification = c("left", "top"),   #图例对齐方式
        axis.title = element_text(size = 16),       # 坐标轴标题大小
        axis.text = element_text(size = 14),        # 坐标轴刻度大小
        legend.text = element_text(size = 14))   # 图例内容大小
# 显示图形
p_bacteria
#ggsave("D:/study/master/Main_Figure_tables/Figure_5/5b_iCAMP_geo_distance_bacteria.png", plot = p_bacteria, width = 8, height = 8, dpi = 600, bg = "transparent")
# 对 HoS 进行线性回归分析
# 注意：自变量是 log10 转换后的地理距离
b_model_hos <- lm(I(HoS * 100) ~ log10(Geo_Distance_Meters / 1000), data = b_iCAMP_geo)
# 查看 HoS 模型的统计摘要
summary(b_model_hos)
# 对 DR 进行线性回归分析
# 注意：自变量是 log10 转换后的地理距离
b_model_dr <- lm(I(DR * 100) ~ log10(Geo_Distance_Meters / 1000), data = b_iCAMP_geo)
# 查看 DR 模型的统计摘要
summary(b_model_dr)
#真菌icamp距离衰减
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
# 计算每个 ASV 的总丰度
f_asv_sums <- taxa_sums(f_merged)
# 计算每个 ASV 在多少个样本中存在
# (这里通过计算每个ASV有多少个非零计数来得到)
f_asv_prevalence <- apply(otu_table(f_merged), 1, function(x) sum(x > 0))
# 根据设定的阈值筛选 ASV
# 创建一个逻辑向量，TRUE 表示 ASV 满足筛选条件
f_keep_taxa <- (f_asv_sums > 3) & (f_asv_prevalence > 3)
# 使用 subset_taxa 函数对 phyloseq 对象进行筛选
f_merged <- phyloseq::subset_taxa(f_merged, f_keep_taxa)
# 从 phyloseq 对象中提取 ASV 表和样本数据
f_ASV_all <- as.data.frame(otu_table(f_merged)) # 确保 ASV 在行，样本在列
metadata <- data.frame(sample_data(f_merged)) # 样本 x 元数据
# 筛选分类信息和系统发育树
# 确保你原始的 `f_tax ` (ASV x Taxonomy) 和 `f_phylo_tree` 是可用的
f_filtered_ASV_ids_grouped <- rownames(f_ASV_all[rowSums(f_ASV_all) > 0, ])
f_filtered_tax_grouped <- f_tax[rownames(f_tax) %in% f_filtered_ASV_ids_grouped, ]
f_filtered_tree_grouped <- drop.tip(f_phylo_tree,
                                    setdiff(f_phylo_tree$tip.label, f_filtered_ASV_ids_grouped))
# 设置一个新的输出目录用于保存分组后的数据
output_wd_all <- "D:/study/master/iCAMP/fungi_all_data" # 新的输出路径
if (!dir.exists(output_wd_all)) {
  dir.create(output_wd_all, recursive = TRUE)
}
# 保存分组后的 ASV 表格
write.table(f_ASV_all, file = paste0(output_wd_all, "/Asv_table_all.txt"), sep = "\t", quote = FALSE)
# 保存筛选后的分类信息
write.table(f_filtered_tax_grouped, file = paste0(output_wd_all, "/classification_all.txt"), sep = "\t", quote = FALSE)
# 保存筛选后的系统发育树
write.tree(f_filtered_tree_grouped, file = paste0(output_wd_all, "/tree_all.nwk"))
# 设置统一的输入和输出路径
input_output_wd <- output_wd_all # 现在输入就是分组后的数据目录
# 设置 iCAMP 分析参数
prefix <- "all_fungi" # 使用新的前缀
# 设置输入文件路径
com.file <- paste0(input_output_wd, "/Asv_table_all.txt")
tree.file <- paste0(input_output_wd, "/tree_all.nwk")
clas.file <- paste0(input_output_wd, "/classification_all.txt")
# 设置输出文件夹路径 (与输入路径相同，或在其下创建子文件夹)
save.wd <- paste0(input_output_wd, "/iCAMP_results_all_data") # 为结果创建一个子文件夹
if (!dir.exists(save.wd)) {
  dir.create(save.wd, recursive = TRUE)
}
#  设置 iCAMP 分析参数 
prefix <- "all_fungi" # 使用统一的前缀表示对所有数据进行分析
rand.time <- 1000
nworker <- 48
memory.G <- 200 # 确保内存设置足够大，以处理整个数据集
#  读取输入数据 
# comm 矩阵通常需要是 (samples x ASVs) 形式
comm <- t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
tree <- read.tree(file = tree.file)
clas <- read.table(clas.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#  匹配物种ID 
# 确保 OTU 表、分类信息和系统发育树的物种 ID 一致
spid.check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tree = tree))
comm <- spid.check$comm
clas <- spid.check$clas
tree <- spid.check$tree
#  计算系统发育距离矩阵 
# wd 参数指向 pdist.big 存储中间文件的目录，这里可以直接使用 save.wd
pd.big <- iCAMP::pdist.big(tree = tree, wd = save.wd, nworker = nworker, memory.G = memory.G)
#  执行 iCAMP 分析 
# iCAMP 分析参数设置
bin.size.limit <- 48
sig.index <- "Confidence"
# 执行 iCAMP 分析，基于丰度和系统发育距离分析群落构建
icres <- iCAMP::icamp.big(
  comm = comm,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  rand = rand.time,
  tree = tree,
  prefix = prefix,
  ds = 0.2,
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "both",
  sig.index = sig.index,
  bin.size.limit = bin.size.limit,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = FALSE,
  ignore.zero = TRUE,
  output.wd = save.wd,
  correct.special = TRUE,
  unit.sum = rowSums(comm), # 针对所有样本总丰度计算 unit.sum
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)
#  保存 iCAMP 分析结果 
# 保存物种的系统发育 bin 和基于 Bray-Curtis 的分析结果
bins <- icres$detail$taxabin$sp.bin
res1 <- icres$CbMPDiCBraya
write.csv(res1, paste0(save.wd, "/", prefix, ".res.csv"), row.names = FALSE)
write.csv(bins, paste0(save.wd, "/", prefix, ".bins.csv"), row.names = FALSE)
# 按系统发育 bin 汇总 iCAMP 结果
icbin <- iCAMP::icamp.bins(
  icamp.detail = icres$detail,
  clas = clas,
  silent = FALSE,
  boot = FALSE,
  rand.time = rand.time,
  between.group = FALSE # 通常在不分组时设为 FALSE
)
# 保存 iCAMP 汇总结果
save(icbin, file = paste0(save.wd, "/", prefix, ".iCAMP.Summary.rda"))
# 保存每个组、每个 bin 和每次 turnover 的生态过程重要性
# 注意：在不分组分析中，icbin$Pt 可能只包含一个整体的 Pt，或者 Ptuv 会包含所有样本对
write.csv(icbin$Pt, file = paste0(save.wd, "/", prefix, ".ProcessImportance_Overall.csv"), row.names = FALSE)
# icbin$Ptk 和 icbin$Ptuv 会包含所有样本对的信息，这里根据实际输出命名
write.csv(icbin$Ptk, file = paste0(save.wd, "/", prefix, ".ProcessImportance_EachBin_Overall.csv"), row.names = FALSE)
write.csv(icbin$Ptuv, file = paste0(save.wd, "/", prefix, ".ProcessImportance_EachTurnover_Overall.csv"), row.names = FALSE)
# 保存每个组中每个 bin 对每个过程的贡献
write.csv(icbin$BPtk, file = paste0(save.wd, "/", prefix, ".BinContributeToProcess_Overall.csv"), row.names = FALSE)
# 保存每个 bin 中的物种及其最高分类级别
write.csv(data.frame(ID = rownames(icbin$Class.Bin), icbin$Class.Bin, stringsAsFactors = FALSE),
          file = paste0(save.wd, "/", prefix, ".Taxon_Bin.csv"), row.names = FALSE)
# 保存每个 bin 的最高分类级别
write.csv(icbin$Bin.TopClass, file = paste0(save.wd, "/", prefix, ".Bin_TopTaxon.csv"), row.names = FALSE) 
# 读取包含iCAMP结果的CSV文件
f_iCAMP <- read.csv("D:/study/master/iCAMP/fungi_all_data/iCAMP_results_all_data/all_fungi.ProcessImportance_EachTurnover_Overall.csv", head=T)
f_iCAMP <- f_iCAMP %>% filter(Method == "CbMPDiCbraya")
f_iCAMP$cut <- paste(f_iCAMP$samp1, f_iCAMP$samp2, sep = "")
#经纬度
NOMIS_df <- data.frame(
  longitude = as.numeric(metadata$longitude),
  latitude = as.numeric(metadata$latitude),
  Sample = rownames(metadata)
)
#计算经纬度的地理距离矩阵
# 从坐标数据框中提取组名和坐标
sample_names <- NOMIS_df$Sample
sample_longitudes <- NOMIS_df$longitude
sample_latitudes <- NOMIS_df$latitude
# 创建一个空的数据框，用于存储每对组的距离
pairwise_geo_distances <- data.frame(
  Sample1 = character(),
  Sample2 = character(),
  Geo_Distance_Meters = numeric(), # 距离单位是米
  stringsAsFactors = FALSE
)
# 循环遍历所有唯一的组对并计算距离
num_samples <- length(sample_names)
for (i in 1:(num_samples - 1)) {
  for (j in (i + 1):num_samples) {
    # 获取第一个组的坐标
    sample1_name <- sample_names[i]
    coords1 <- c(sample_longitudes[i], sample_latitudes[i])
    # 获取第二个组的坐标
    sample2_name <- sample_names[j]
    coords2 <- c(sample_longitudes[j], sample_latitudes[j])
    # 计算地理距离
    distance <- geosphere::distGeo(p1 = coords1, p2 = coords2)
    # 将结果添加到数据框中
    pairwise_geo_distances <- rbind(pairwise_geo_distances, data.frame(
      Sample1 = sample1_name,
      Sample2 = sample2_name,
      Geo_Distance_Meters = distance
    ))
  }
}
pairwise_geo_distances$cut <- paste(pairwise_geo_distances$Sample1, pairwise_geo_distances$Sample2, sep = "")
f_iCAMP_geo <- f_iCAMP %>%
  # 使用 left_join 合并，以 'cut' 列为共同键
  left_join(pairwise_geo_distances, by = "cut") %>%
  # 选择你最终希望保留的列
  dplyr::select(
    samp1,
    samp2,
    HeS,
    HoS,
    DL,
    HD,
    DR,
    cut,
    Geo_Distance_Meters
  )
f_iCAMP_geo <- f_iCAMP_geo %>%
  filter(Geo_Distance_Meters != 0)
# 创建ggplot对象
p_fungi <- ggplot() +
  # 添加地理距离与HoS散点图
  geom_point(data= f_iCAMP_geo, aes(Geo_Distance_Meters/1000, HoS*100, color = "HoS"), col=rgb(0,0,1,0.25), cex=2) +
  # 添加地理距离与HoS线性拟合平滑曲线
  geom_smooth(data=f_iCAMP_geo, aes(Geo_Distance_Meters/1000, HoS*100, color = "HoS"), method="lm") +
  # 添加地理距离与DR散点图
  geom_point(data=f_iCAMP_geo, aes(Geo_Distance_Meters/1000, DR*100, color = "DR"), col=rgb(1,0,0,0.25), cex=2) +
  # 添加地理距离与DR线性拟合平滑曲线
  geom_smooth(data=f_iCAMP_geo, aes(Geo_Distance_Meters/1000, DR*100, color = "DR"), method="lm") +
  # 手动设置颜色和图例标签
  scale_color_manual(
    values = c("HoS" = "blue", "DR" = "red"), # 定义 HoS 为蓝色，DR 为红色
    labels = c("HoS", "DR") # 定义图例显示为 "HoS" 和 "DR"
  ) +
  # 使用白色背景主题
  theme_bw() +
  # 对x轴进行log10转换
  scale_x_continuous(trans='log10') +
  # 设置坐标轴标签
  labs(x="Geographic distance (km)", y="Fungal contribution (%)")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),            #移除图例标题
        legend.position = c(0.05, 0.95),           #设置图例位置 (x, y 坐标)
        legend.justification = c("left", "top"),   #图例对齐方式
        axis.title = element_text(size = 16),       # 坐标轴标题大小
        axis.text = element_text(size = 14),        # 坐标轴刻度大小
        legend.text = element_text(size = 14))   # 图例内容大小
# 显示图形
p_fungi
#ggsave("D:/study/master/Main_Figure_tables/Figure_5/5b_iCAMP_geo_distance_fungi.png", plot = p_fungi, width = 8, height = 8, dpi = 600, bg = "transparent")
# 对 HoS 进行线性回归分析
# 注意：自变量是 log10 转换后的地理距离
f_model_hos <- lm(I(HoS * 100) ~ log10(Geo_Distance_Meters / 1000), data = f_iCAMP_geo)
# 查看 HoS 模型的统计摘要
summary(f_model_hos)
# 对 DR 进行线性回归分析
# 注意：自变量是 log10 转换后的地理距离
f_model_dr <- lm(I(DR * 100) ~ log10(Geo_Distance_Meters / 1000), data = f_iCAMP_geo)
# 查看 DR 模型的统计摘要
summary(f_model_dr)

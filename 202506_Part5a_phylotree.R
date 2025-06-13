#图5a 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")
install.packages(c("tidyverse", "readr", "readxl", "data.table", "phyloseq",
                   "qiime2R", "ape", "magrittr", "dplyr", "tidyr", "iCAMP",
                   "phytools", "ggplot2", "ggnewscale", "castor", "viridis",
                   "phangorn"))
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
#细菌系统发育树
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
b_tree <- read_qza("D:/study/master/meiji/b_rooted-tree_3_3_uncl_na/b_final_rooted_tree.qza")
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
# 提取组信息
groups <- unique(metadata$Group)
# 遍历每个组
for (group in groups) {
  # 构建当前组的工作目录
  group_wd <- paste0("D:/study/master/iCAMP/bacteria/", group)
  # 如果目录不存在则创建
  if (!dir.exists(group_wd)) {
    dir.create(group_wd, recursive = TRUE)
  }
  # 获取当前组的样本ID
  group_samples <- rownames(subset(metadata, Group == group))
  # 筛选当前组的ASV表格
  group_ASV <- b_ASV[, colnames(b_ASV) %in% group_samples]
  # 获取当前组中存在的ASV的ID（去除丰度为0的ASV）
  group_ASV_ids <- rownames(group_ASV[rowSums(group_ASV) > 0, ])
  # 筛选当前组的分类信息
  group_tax <- b_tax[rownames(b_tax) %in% group_ASV_ids, ]
  # 筛选当前组的系统发育树（去除不在当前组中的ASV）
  group_tree <- drop.tip(b_phylo_tree, setdiff(b_phylo_tree$tip.label, group_ASV_ids))
  # 保存当前组的ASV表格
  write.table(group_ASV, file = paste0(group_wd, "/Asv_table.txt"), sep = "\t", quote = FALSE)
  # 保存当前组的分类信息
  write.table(group_tax, file = paste0(group_wd, "/classification.txt"), sep = "\t", quote = FALSE)
  # 保存当前组的系统发育树
  write.tree(group_tree, file = paste0(group_wd, "/tree.nwk"))
}
# 获取所有组的文件夹
dirs <- list.dirs("D:/study/master/iCAMP/bacteria", full.names = TRUE, recursive = FALSE)
# 遍历每个组的文件夹进行iCAMP分析
for (group_wd in dirs) {
  # 设置输入文件路径
  com.file <- paste0(group_wd, "/Asv_table.txt")
  tree.file <- paste0(group_wd, "/tree.nwk")
  clas.file <- paste0(group_wd, "/classification.txt")
  # 设置输出文件夹路径
  save.wd <- paste0(group_wd, "/out2")
  if (!dir.exists(save.wd)) {
    dir.create(save.wd)
  }
  # 设置分析参数
  prefix <- basename(group_wd) # 使用组名作为输出前缀
  rand.time <- 1000
  nworker <- 48
  memory.G <- 200
  # 读取输入数据
  comm <- t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  tree <- read.tree(file = tree.file)
  clas <- read.table(clas.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # 匹配物种ID
  spid.check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tree = tree))
  comm <- spid.check$comm
  clas <- spid.check$clas
  tree <- spid.check$tree
  # 计算系统发育距离矩阵
  pd.big <- iCAMP::pdist.big(tree = tree, wd = save.wd, nworker = nworker, memory.G = memory.G)
  # iCAMP分析参数设置
  bin.size.limit <- 48
  sig.index <- "Confidence"
  # 执行iCAMP分析，基于丰度和系统发育距离分析群落构建
  icres <- iCAMP::icamp.big(comm = comm, pd.desc = pd.big$pd.file,
                            pd.spname = pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree = tree, prefix = prefix,
                            ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "both", sig.index = sig.index, bin.size.limit = bin.size.limit,
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE,
                            detail.save = TRUE, qp.save = FALSE, detail.null = FALSE,
                            ignore.zero = TRUE, output.wd = save.wd, correct.special = TRUE,
                            unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975, omit.option = "no", meta.ab = NULL)
  # 保存物种的系统发育bin和基于Bray-Curtis的分析结果
  bins <- icres$detail$taxabin$sp.bin
  res1 <- icres$CbMPDiCBraya
  write.csv(res1, paste0(group_wd, "/", prefix, ".res.csv"))
  write.csv(bins, paste0(group_wd, "/", prefix, ".bins.csv"))
  # 按系统发育bin汇总iCAMP结果
  icbin <- iCAMP::icamp.bins(icamp.detail = icres$detail, clas = clas,
                             silent = FALSE, boot = FALSE,
                             rand.time = rand.time, between.group = FALSE)
  # 保存iCAMP汇总结果
  save(icbin, file = paste0(group_wd, "/", prefix, ".iCAMP.Summary.rda"))
  # 保存每个组、每个bin和每次turnover的生态过程重要性
  write.csv(icbin$Pt, file = paste0(group_wd, "/", prefix, ".ProcessImportance_EachGroup.csv"), row.names = FALSE)
  write.csv(icbin$Ptk, file = paste0(group_wd, "/", prefix, ".ProcessImportance_EachBin_EachGroup.csv"), row.names = FALSE)
  write.csv(icbin$Ptuv, file = paste0(group_wd, "/", prefix, ".ProcessImportance_EachTurnover.csv"), row.names = FALSE)
  # 保存每个组中每个bin对每个过程的贡献
  write.csv(icbin$BPtk, file = paste0(group_wd, "/", prefix, ".BinContributeToProcess_EachGroup.csv"), row.names = FALSE)
  # 保存每个bin中的物种及其最高分类级别
  write.csv(data.frame(ID = rownames(icbin$Class.Bin), icbin$Class.Bin, stringsAsFactors = FALSE),
            file = paste0(group_wd, "/", prefix, ".Taxon_Bin.csv"), row.names = FALSE)
  # 保存每个bin的最高分类级别
  write.csv(icbin$Bin.TopClass, file = paste0(group_wd, "/", prefix, ".Bin_TopTaxon.csv"), row.names = FALSE)
}
# 初始化存储所有组合iCAMP结果的数据框
iCAMP_all_combined <- data.frame()
# 遍历每个组的文件夹
for (group_wd in dirs) {
  # 获取组名
  prefix <- basename(group_wd)
  # 读取每个bin的OTU信息
  taxon_bin_data <- read.csv(paste0(group_wd, "/", prefix, ".Taxon_Bin.csv"), header = TRUE, stringsAsFactors = FALSE)
  # 读取每个bin的生态过程信息
  mntd_data <- read.csv(paste0(group_wd, "/", prefix, ".ProcessImportance_EachBin_EachGroup.csv"), header = FALSE, stringsAsFactors = FALSE)
  # 选取Taxon_Bin的前三列(Bin, Top.Taxon, Top.Taxon.Frac)
  taxon_bin_data <- taxon_bin_data[, 1:3]
  # 转置生态过程数据并转换为数据框
  mntd_data <- as.data.frame(t(mntd_data[, -c(1:3)]))
  # 选取Bin列和CbMNTDiCbraya列(假设是第14列)
  mntd_data <- mntd_data[, c(1, 14)]
  # 设置列名
  colnames(mntd_data) <- as.character(mntd_data[1, ])
  # 移除第一行(已作为列名)
  mntd_data <- mntd_data[-1, ]
  # 设置最终列名
  colnames(mntd_data) <- c("Bin", "process")
  # 统一Bin的格式
  mntd_data$Bin <- gsub("bin", "Bin", mntd_data$Bin)
  # 合并taxon_bin_data和mntd_data
  iCAMP_data <- merge(taxon_bin_data, mntd_data, by = "Bin", all.x = TRUE)
  # 添加组名列
  iCAMP_data$group <- prefix
  # 合并到总数据框
  iCAMP_all_combined <- rbind(iCAMP_all_combined, iCAMP_data)
}
# 保存合并后的结果
write.table(iCAMP_all_combined, file = "D:/study/master/iCAMP/bacteria/iCAMP_all_combined.csv",
            sep = ",", row.names = FALSE, quote = FALSE)
# ASV丰度转换为数据框
b_ASVs <- as.data.frame(b_ASV_df)
# 获取系统发育树并中点根化
b_tree <- phy_tree(b_merged)
b_tree <- midpoint(b_tree)
# 计算平均相对丰度
b_mean.RA <- data.frame(rowMeans(b_ASVs))
b_mean.RA$ASVs <- rownames(b_ASVs)
b_mean.RA$log.RA <- log1p(b_mean.RA$rowMeans.b_ASVs.)
colnames(b_mean.RA) <- c("mean.RA", "ASVs", "log.RA")
# 读取iCAMP合并结果
b_all.ASVs <- read.csv("D:/study/master/iCAMP/bacteria/iCAMP_all_combined.csv")
# 筛选DL、HoS和DR过程的ASV
b_all.ASVs <- b_all.ASVs[b_all.ASVs$process == "DL" | b_all.ASVs$process == "HoS" | b_all.ASVs$process == "DR",]
# 分类信息转换为数据框
b_tax <- as.data.frame(b_tax)
b_tax$ASV <- rownames(b_tax)
b_tax <- b_tax %>% filter(ASV %in% b_tree$tip.label)
# 筛选特定门的ASV
metadata <- as.data.frame(metadata)
b_ASV_df <- as.data.frame(b_ASV_df)
b_ASV_df$ASV <- rownames(b_ASV_df)
b_tax$ASV <- rownames(b_tax)
#在每个 Genus 分组，将同一个属下的所有 ASV 的计数在每个样本中加起来，得到每个属在每个样本的总计数。
b_dat_m <- b_ASV_df %>%
  left_join(b_tax)%>%#合并
  filter(Class != "c__bacteriap25") %>%
  group_by(Genus)%>%#按属分组
  filter(!(Genus == "g__uncultured" | stringr::str_detect(Genus, "g__unclassified")))%>%
  summarise_if(is.numeric, sum)%>%# 找到所有数值型的列，并对这些列计算它们的 sum
  na.omit()%>%
  pivot_longer(cols = !Genus)%>%# 从“宽”格式转换为“长”格式
  left_join(metadata, by = c("name" = "Sample.ID"))#合并
b_sum_all <- sum(b_dat_m$value)
#计算属相对丰度（单个属的菌/总菌）
b_dat_abun <- b_dat_m %>%
  group_by(Genus)%>%
  summarise(Abundance = sum(value)/b_sum_all)
colnames(b_dat_abun) <- c("Genus", "Abundance")
b_tax_sel <- b_tax%>%
  select(-c(ASV,Species))%>%
  distinct()%>%#移除重复的行
  filter(Genus %in% b_dat_abun$Genus)
#结合属、相对丰度、分类信息
b_dat_final <- b_dat_abun%>%
  left_join(b_tax_sel)
#优势菌门
b_phylumToPlot <- b_dat_final %>%
  group_by(Phylum)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(10)
b_phylumToPlot
p_prot_asvs <- b_tax[b_tax$Phylum == "p__Proteobacteria", ]
p_acti_asvs <- b_tax[b_tax$Phylum == "p__Actinobacteriota", ]
p_bact_asvs <- b_tax[b_tax$Phylum == "p__Bacteroidota", ]
p_firm_asvs <- b_tax[b_tax$Phylum == "p__Firmicutes", ]
p_acid_asvs <- b_tax[b_tax$Phylum == "p__Acidobacteriota", ]
p_myxo_asvs <- b_tax[b_tax$Phylum == "p__Myxococcota", ]
p_nitr_asvs <- b_tax[b_tax$Phylum == "p__Nitrospirota", ]
p_chlo_asvs <- b_tax[b_tax$Phylum == "p__Chloroflexi", ]
p_verr_asvs <- b_tax[b_tax$Phylum == "p__Verrucomicrobiota", ]
p_ento_asvs <- b_tax[b_tax$Phylum == "p__Entotheonellaeota", ]
# 创建ASV、Family和Phylum信息数据框
b_taxa_df <- data.frame(ASV = b_tree$tip.label, Family = rep('Others', length(b_tree$tip.label)),
                        Phylum = rep('Others', length(b_tree$tip.label)))
# 根据门信息更新Phylum列
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_prot_asvs)] <- 'Proteobacteria'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_acti_asvs)] <- 'Actinobacteriota'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_bact_asvs)] <- 'Bacteroidota'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_firm_asvs)] <- 'Firmicutes'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_acid_asvs)] <- 'Acidobacteriota'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_myxo_asvs)] <- 'Myxococcota'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_nitr_asvs)] <- 'Nitrospirota'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_chlo_asvs)] <- 'Chloroflexi'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_verr_asvs)] <- 'Verrucomicrobiota'
b_taxa_df$Phylum[b_taxa_df$ASV %in% rownames(p_ento_asvs)] <- 'Entotheonellaeota'
# 绘制带有分类学、热图和平均相对丰度的树
b_p = ggtree(b_tree, layout="fan", size=0.25, open.angle=10)
# 定义你想要的顺序
desired_group_order <- c("JRN", "JJN", "TZN", "PAN", "JRG", "JJG", "TZG", "PAG")
# 重新排序 f_all.ASVs$group 的因子水平
b_all.ASVs$group <- factor(b_all.ASVs$group, levels = desired_group_order)
# 添加门水平分类注释
b_p2 <- b_p +
  geom_fruit(data= b_taxa_df, geom=geom_tile,
             mapping=aes(y=ASV, fill=Phylum),
             width = 0.05, offset = 0.02) +
  scale_fill_manual(values = c('#46CF6B', '#bbfa7b', '#FA907B','#84CCFA',
                               '#fafa7b','#cf46aa', '#3234D9','#529C7E','white',
                               '#C42D50','#7bbbfa') , guide = "none") + new_scale_fill() +
  # 添加iCAMP过程热图注释
  geom_fruit(
    data= b_all.ASVs,
    geom=geom_tile,
    mapping=aes(y=ID, x=group, fill=process),
    offset=0.08, pwidth=0.4
  ) +
  scale_fill_manual(
    values=c("#3303df", "#dfac03","#df3e03"),
    name = NULL,
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  )
# 添加平均相对丰度条形图注释
b_p3 <- b_p2 +
  new_scale_fill()+
  geom_fruit(data= b_mean.RA,
             geom=geom_bar,
             mapping=aes(y=ASVs, x=log.RA, fill=log.RA),
             offset = 0.07, orientation='y',
             stat="identity",
             pwidth = 0.4) +#相对总宽度
  scale_fill_viridis(option="D", guide = "none")
p4_bacteria <- b_p3 +
  theme(
    # 设置图例位置，c(x, y) 向量表示在主图区域内的比例位置 (0到1)。
    # 对于圆图的右上角，大概是 x 接近 1，y 接近 1 的位置。
    legend.position = c(0.8, 0.85), # 示例：放置在图的右上方
    legend.justification = c(1, 1), # 右上对齐
    legend.box.just = "right",      # 图例框右对齐
    legend.background = element_rect(fill = "white", color = NA), # 可选：设置图例背景，避免被图覆盖
    legend.margin = margin(0, 0, 0, 0), # 可选：调整图例边距
    legend.text = element_text(size = 16)
  )
p4_bacteria
#ggsave("D:/study/master/Main_Figure_tables/Figure_5/5a_phy_tree_bacteria.png", plot = p4_bacteria, width = 15, height = 15, dpi = 600, bg = "transparent")
# 自定义标签和颜色（与你主图保持一致）
legend_df <- data.frame(
  x = rep(seq(1, 4), each = 2),  # 4 列，每列 2 行
  y = rep(2:1, 4),               # 行方向
  label = c("Jurong Bulb", "Jurong Rhizosphere Soil",
            "Jingjiang Bulb", "Jingjiang Rhizosphere Soil",
            "Tongzhou Bulb", "Tongzhou Rhizosphere Soil",
            "Panan Bulb", "Panan Rhizosphere Soil"),
  color = c("#E5614CFF", "#8C57A2FF",   # Jurong Bulb, Soil
            "#97A1A7FF", "#3EBCB6",    # Jingjiang Bulb, Soil
            "#DC9445FF", "#82581FFF",  # Tongzhou Bulb, Soil
            "#bee183", "#2F509EFF")    # Panan Bulb, Soil
  #label = c("Jurong Rhizosphere Soil", "Jingjiang Rhizosphere Soil", "Tongzhou Rhizosphere Soil", "Panan Rhizosphere Soil",
  #"Jurong Bulb", "Jingjiang Bulb", "Tongzhou Bulb", "Panan Bulb"),
  #color = c("#8C57A2", "#3EBCB6", "#82581F", "#2F509E", "#E5614C", "#97A1A7", "#DC9445", "#BEE183")
)
# 图例图层（圆点 + 标签）
legend_plot <- ggplot(legend_df, aes(x, y)) +
  geom_point(aes(color = color), size = 4.5, shape = 16, show.legend = FALSE) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.05, color = "black", size = 4.5) +
  scale_color_identity() +
  theme_void() +
  xlim(1, 5) + ylim(0.5, 2.5)
legend_plot
#ggsave("D:/study/master/Main_Figure_tables/Figure_5/5a_point_legend.png", plot = legend_plot, width = 12, height = 0.8, dpi = 2000, bg = "transparent")
# 定义图例数据框，点之间更接近
legend_df_single_row <- data.frame(
  # x 坐标：进一步缩小总跨度，使其更紧凑
  x = seq(from = 1, to = 2.75, length.out = 8), # 总跨度从 4.5 缩小到 2.75 (大约一半)
  # y 坐标：全部设置为 1，以便在同一行
  y = rep(1, 8),
  color = c("#E5614CFF", "#97A1A7FF", "#DC9445FF", "#bee183", # Bulb 颜色
            "#8C57A2FF", "#3EBCB6", "#82581FFF", "#2F509EFF")  # Soil 颜色
)
# 创建单行图例图层
legend_plot_single_row <-
  ggplot(legend_df_single_row, aes(x, y)) +
  geom_point(aes(color = color), size = 6, shape = 16, show.legend = FALSE) +
  scale_color_identity() +
  theme_void() + # 移除所有轴和背景元素
  # 调整 xlim 以适应新的、更紧凑的点排列
  xlim(0.8, 2.95) + # 根据新的 x 坐标范围进行调整
  # 保持 ylim 用于单行，并居中
  ylim(0.8, 1.2) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm")) # 最小化绘图边距
legend_plot_single_row
# ggsave("D:/study/master/Main_Figure_tables/Figure_5/5a_single_row_legend.png", plot = legend_plot_single_row, width = 3, height = 0.5, dpi = 2000, bg = "transparent") 
#推导主要生态过程并计算百分比
# 假设 dirs 变量已经包含了所有组的文件夹路径
dirs <- list.dirs("D:/study/master/iCAMP/bacteria", full.names = TRUE, recursive = FALSE)
# 初始化一个空的数据框，用于存储所有组的 ProcessImportance_EachTurnover 结果
b_all_turnovers_combined <- data.frame()
# 遍历每个组的文件夹
for (group_wd in dirs) {
  prefix <- basename(group_wd)
  turnover_file <- paste0(group_wd, "/", prefix, ".ProcessImportance_EachTurnover.csv") 
  if (file.exists(turnover_file)) {
    current_turnover_data <- read.csv(turnover_file, header = TRUE, stringsAsFactors = FALSE)
    # 定义包含过程得分的列名
    process_score_cols_to_check <- c("HeS", "HoS", "DL", "HD", "DR") 
    # 检查所有必需的过程得分列是否存在
    if (all(process_score_cols_to_check %in% colnames(current_turnover_data))) {
      current_turnover_data$group <- prefix
      # 不再需要将 Method 列赋给 Process
      # 'Process' 列将在下一步从得分中推导出来
      b_all_turnovers_combined <- rbind(b_all_turnovers_combined, current_turnover_data)
    } else {
      warning(paste0("警告: 文件 ", turnover_file, " 中未找到预期的过程得分列 (HeS, HoS, DL, HD, DR)。请检查文件头。跳过此文件。"))
    }
  } else {
    warning(paste0("警告: 未找到文件: ", turnover_file, "。跳过此文件。"))
  }
}
#现在，为每个群落对推导主要过程
if (nrow(b_all_turnovers_combined) > 0) {
  # 定义包含过程得分的列名
  process_score_cols <- c("HeS", "HoS", "DL", "HD", "DR") 
  # 创建一个新的 'Process' 列，通过查找每行中得分最高的列名来赋值
  # apply 函数会逐行操作：对于每一行，找到最大值的列的索引，然后返回该列的名称
  b_all_turnovers_combined$Process <- apply(
    b_all_turnovers_combined[, process_score_cols], 1, # 针对这些列，按行操作
    function(x) {
      # 找到最大值的索引
      max_idx <- which.max(x)
      # 返回该索引对应的列名
      return(names(x)[max_idx])
    }
  ) 
  # 筛选出你感兴趣的特定过程 (DL, HoS, DR)
  # 我们这里是根据论文结论来筛选，即扩散限制 (DL)、同质选择 (HoS) 和生态漂移 (DR)
  relevant_processes_data <- b_all_turnovers_combined[
    b_all_turnovers_combined$Process %in% c("DL", "HoS", "DR"),
  ] 
  # 检查过滤后是否有相关数据
  if (nrow(relevant_processes_data) == 0) {
    cat("警告：过滤后，没有找到 'DL', 'HoS', 'DR' 这三种过程的群落对。这可能意味着其他过程是主要的驱动力。\n")
    # 如果你好奇其他过程的分类情况，可以打印出所有被分配到的过程类型
    cat("所有群落对分配到的独特过程类型有：\n")
    print(unique(b_all_turnovers_combined$Process))
    # 如果没有相关数据，就直接返回，避免后续统计报错
    return(NULL) 
  } 
  # 统计每个过程的出现次数
  process_counts <- table(relevant_processes_data$Process) 
  # 计算这些相关过程的总群落对数量
  total_relevant_pairs <- sum(process_counts) 
  # 计算百分比
  # 为防止某个过程没有出现而导致计算错误，我们使用 if 语句进行检查，如果不存在则设为 0
  percentage_DL <- if ("DL" %in% names(process_counts)) (process_counts["DL"] / total_relevant_pairs) * 100 else 0
  percentage_HoS <- if ("HoS" %in% names(process_counts)) (process_counts["HoS"] / total_relevant_pairs) * 100 else 0
  percentage_DR <- if ("DR" %in% names(process_counts)) (process_counts["DR"] / total_relevant_pairs) * 100 else 0 
  # 打印结果
  cat("基于群落对的生态组装过程百分比:\n")
  cat("扩散限制 (DL): ", round(percentage_DL, 2), "% 的群落对\n")
  cat("同质选择 (HoS): ", round(percentage_HoS, 2), "% 的群落对\n")
  cat("生态漂移 (DR): ", round(percentage_DR, 2), "% 的群落对\n")
  cat("\n原始计数:\n")
  print(process_counts) 
} else {
  cat("错误：未能成功读取或处理任何 ProcessImportance_EachTurnover.csv 文件。请检查文件路径和内容。\n")
}
# --- 确保 b_asv_process_abundance 包含 'group' 列 ---
# 假设 b_all.ASVs_processed 已包含 'group' 列，来源于 iCAMP_all_combined.csv。
# 如果 b_all.ASVs_processed 中没有 'group' 列，需在此合并 group 信息。
# 例如，若 b_all.ASVs 仅有 ID 和 process，而 group 信息在 metadata 中，则需先将 group 信息加入 b_all.ASVs_processed 中。
# 此处假设 b_all.ASVs_processed 已包含 'group' 列。
b_all.ASVs_processed <- b_all.ASVs
colnames(b_all.ASVs_processed)[colnames(b_all.ASVs_processed) == "ID"] <- "ASVs" # 根据代码，ID列是ASV标识。
# 合并 b_all.ASVs_processed（包含ASV主导过程和group）与 b_mean.RA（包含ASV平均相对丰度）。
# 确保 b_asv_process_abundance 现在有 'ASVs', 'process', 'mean.RA', 和 'group' 列。
b_asv_process_abundance <- merge(b_all.ASVs_processed, b_mean.RA, by = "ASVs", all.x = TRUE)
# 移除合并后可能产生的NA值（通常不会有，但安全起见）。
b_asv_process_abundance <- na.omit(b_asv_process_abundance)
# --- 核心分组分析：按组计算各过程主导ASV的总相对丰度及其百分比贡献 ---
# 针对论文范例中“在每个山脉范围内”的分析。
b_process_abundance_summary_by_group <- b_asv_process_abundance %>%
  dplyr::group_by(group, process) %>% # 关键修改：先按组分组，再按过程分组。
  dplyr::summarise(TotalRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>% # 计算该组该过程下所有ASV的总相对丰度。
  dplyr::ungroup() %>% # 第一次ungroup()，为下一步按group计算百分比。
  dplyr::group_by(group) %>% # 再次分组，确保百分比是组内的。
  dplyr::mutate(PercentageContribution = (TotalRelativeAbundance / sum(TotalRelativeAbundance)) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, desc(PercentageContribution)) # 按组和百分比降序排列。
print("\n--- 基于ASV平均相对丰度，各组各生态组装过程主导的ASV贡献百分比 ---")
print(b_process_abundance_summary_by_group)
# --- 按组识别特定过程下的关键属 ---
# 确保 b_tax 数据框可以用于合并（ASVs作为一列）。
b_tax_df_for_merge <- as.data.frame(b_tax)
b_tax_df_for_merge$ASVs <- rownames(b_tax_df_for_merge) # 修改为 'ASVs' 以匹配合并。
# 合并分类学信息到包含ASV过程和丰度的数据框（b_asv_process_abundance_taxa 已包含 group 列）。
b_asv_process_abundance_taxa <- merge(b_asv_process_abundance, b_tax_df_for_merge, by = "ASVs", all.x = TRUE)
print("\n--- 按组识别特定过程下的关键属 ---")
# 示例：查找在“同质选择（HoS）”主导下，各组丰度最高的属。
b_top_hos_genera_by_group <- b_asv_process_abundance_taxa %>%
  dplyr::filter(process == "HoS") %>%
  dplyr::group_by(group, Genus) %>% # 按组和属分组。
  dplyr::summarise(GenusRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>%
  dplyr::arrange(group, desc(GenusRelativeAbundance)) %>%
  dplyr::slice_head(n = 5) # 对每个组取前5个。
print("\n同质选择 (HoS) 主导下的各组主要属（按丰度降序）:")
print(b_top_hos_genera_by_group)
# 示例：查找在“扩散限制（DL）”主导下，各组丰度最高的属。
b_top_dl_genera_by_group <- b_asv_process_abundance_taxa %>%
  dplyr::filter(process == "DL") %>%
  dplyr::group_by(group, Genus) %>% # 按组和属分组。
  dplyr::summarise(GenusRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>%
  dplyr::arrange(group, desc(GenusRelativeAbundance)) %>%
  dplyr::slice_head(n = 5)
print("\n扩散限制 (DL) 主导下的各组主要属（按丰度降序）:")
print(b_top_dl_genera_by_group)
# 示例：查找在“生态漂移（DR）”主导下，各组丰度最高的属。
b_top_dr_genera_by_group <- b_asv_process_abundance_taxa %>%
  dplyr::filter(process == "DR") %>%
  dplyr::group_by(group, Genus) %>% # 按组和属分组。
  dplyr::summarise(GenusRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>%
  dplyr::arrange(group, desc(GenusRelativeAbundance)) %>%
  dplyr::slice_head(n = 5)
print("\n生态漂移 (DR) 主导下的各组主要属（按丰度降序）:")
print(b_top_dr_genera_by_group)
# --- 计算各过程贡献百分比的中位数和IQR ---
# 计算同质选择 (HoS) 的贡献百分比中位数和IQR。
b_hos_contributions <- b_process_abundance_summary_by_group %>%
  dplyr::filter(process == "HoS") %>%
  dplyr::pull(PercentageContribution) # 提取 PercentageContribution 列。
b_median_hos <- median(b_hos_contributions, na.rm = TRUE)
b_q1_hos <- quantile(b_hos_contributions, probs = 0.25, na.rm = TRUE)
b_q3_hos <- quantile(b_hos_contributions, probs = 0.75, na.rm = TRUE)
cat(paste0("同质选择 (HoS) 主导的ASV贡献了: ", round(b_median_hos, 2), "% ",
           "(IQR: ", round(b_q1_hos, 2), "–", round(b_q3_hos, 2), "%)",
           " 的相对丰度，在每个组中。\n"))
# 计算生态漂移 (DR) 的贡献百分比中位数和IQR。
b_dr_contributions <- b_process_abundance_summary_by_group %>%
  dplyr::filter(process == "DR") %>%
  dplyr::pull(PercentageContribution) # 提取 PercentageContribution 列。
b_median_dr <- median(b_dr_contributions, na.rm = TRUE)
b_q1_dr <- quantile(b_dr_contributions, probs = 0.25, na.rm = TRUE)
b_q3_dr <- quantile(b_dr_contributions, probs = 0.75, na.rm = TRUE)
cat(paste0("生态漂移 (DR) 主导的ASV贡献了: ", round(b_median_dr, 2), "% ",
           "(IQR: ", round(b_q1_dr, 2), "–", round(b_q3_dr, 2), "%)",
           " 的相对丰度，在每个组中。\n"))
#真菌系统发育树
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
f_tree <- read_qza("D:/study/master/meiji/f_rooted-tree_3_3_uncl_na/f_final_rooted_tree.qza")
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
# 提取组信息
groups <- unique(metadata$Group)
# 遍历每个组
for (group in groups) {
  # 构建当前组的工作目录
  group_wd <- paste0("D:/study/master/iCAMP/fungi/", group)
  # 如果目录不存在则创建
  if (!dir.exists(group_wd)) {
    dir.create(group_wd, recursive = TRUE)
  }
  # 获取当前组的样本ID
  group_samples <- rownames(subset(metadata, Group == group))
  # 筛选当前组的ASV表格
  group_ASV <- f_ASV[, colnames(f_ASV) %in% group_samples]
  # 获取当前组中存在的ASV的ID（去除丰度为0的ASV）
  group_ASV_ids <- rownames(group_ASV[rowSums(group_ASV) > 0, ])
  # 筛选当前组的分类信息
  group_tax <- f_tax[rownames(f_tax) %in% group_ASV_ids, ]
  # 筛选当前组的系统发育树（去除不在当前组中的ASV）
  group_tree <- drop.tip(f_phylo_tree, setdiff(f_phylo_tree$tip.label, group_ASV_ids))
  # 保存当前组的ASV表格
  write.table(group_ASV, file = paste0(group_wd, "/Asv_table.txt"), sep = "\t", quote = FALSE)
  # 保存当前组的分类信息
  write.table(group_tax, file = paste0(group_wd, "/classification.txt"), sep = "\t", quote = FALSE)
  # 保存当前组的系统发育树
  write.tree(group_tree, file = paste0(group_wd, "/tree.nwk"))
}
# 获取所有组的文件夹
dirs <- list.dirs("D:/study/master/iCAMP/fungi", full.names = TRUE, recursive = FALSE)
# 遍历每个组的文件夹进行iCAMP分析
for (group_wd in dirs) {
  # 设置输入文件路径
  com.file <- paste0(group_wd, "/Asv_table.txt")
  tree.file <- paste0(group_wd, "/tree.nwk")
  clas.file <- paste0(group_wd, "/classification.txt")
  # 设置输出文件夹路径
  save.wd <- paste0(group_wd, "/out2")
  if (!dir.exists(save.wd)) {
    dir.create(save.wd)
  }
  # 设置分析参数
  prefix <- basename(group_wd) # 使用组名作为输出前缀
  rand.time <- 1000
  nworker <- 48
  memory.G <- 200
  # 读取输入数据
  comm <- t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  tree <- read.tree(file = tree.file)
  clas <- read.table(clas.file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # 匹配物种ID
  spid.check <- match.name(cn.list = list(comm = comm), rn.list = list(clas = clas), tree.list = list(tree = tree))
  comm <- spid.check$comm
  clas <- spid.check$clas
  tree <- spid.check$tree
  # 计算系统发育距离矩阵
  pd.big <- iCAMP::pdist.big(tree = tree, wd = save.wd, nworker = nworker, memory.G = memory.G)
  # iCAMP分析参数设置
  bin.size.limit <- 48
  sig.index <- "Confidence"
  # 执行iCAMP分析，基于丰度和系统发育距离分析群落构建
  icres <- iCAMP::icamp.big(comm = comm, pd.desc = pd.big$pd.file,
                            pd.spname = pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree = tree, prefix = prefix,
                            ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "both", sig.index = sig.index, bin.size.limit = bin.size.limit,
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE,
                            detail.save = TRUE, qp.save = FALSE, detail.null = FALSE,
                            ignore.zero = TRUE, output.wd = save.wd, correct.special = TRUE,
                            unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut = 0.975, omit.option = "no", meta.ab = NULL)
  # 保存物种的系统发育bin和基于Bray-Curtis的分析结果
  bins <- icres$detail$taxabin$sp.bin
  res1 <- icres$CbMPDiCBraya
  write.csv(res1, paste0(group_wd, "/", prefix, ".res.csv"))
  write.csv(bins, paste0(group_wd, "/", prefix, ".bins.csv"))
  # 按系统发育bin汇总iCAMP结果
  icbin <- iCAMP::icamp.bins(icamp.detail = icres$detail, clas = clas,
                             silent = FALSE, boot = FALSE,
                             rand.time = rand.time, between.group = FALSE)
  # 保存iCAMP汇总结果
  save(icbin, file = paste0(group_wd, "/", prefix, ".iCAMP.Summary.rda"))
  # 保存每个组、每个bin和每次turnover的生态过程重要性
  write.csv(icbin$Pt, file = paste0(group_wd, "/", prefix, ".ProcessImportance_EachGroup.csv"), row.names = FALSE)
  write.csv(icbin$Ptk, file = paste0(group_wd, "/", prefix, ".ProcessImportance_EachBin_EachGroup.csv"), row.names = FALSE)
  write.csv(icbin$Ptuv, file = paste0(group_wd, "/", prefix, ".ProcessImportance_EachTurnover.csv"), row.names = FALSE)
  # 保存每个组中每个bin对每个过程的贡献
  write.csv(icbin$BPtk, file = paste0(group_wd, "/", prefix, ".BinContributeToProcess_EachGroup.csv"), row.names = FALSE)
  # 保存每个bin中的物种及其最高分类级别
  write.csv(data.frame(ID = rownames(icbin$Class.Bin), icbin$Class.Bin, stringsAsFactors = FALSE),
            file = paste0(group_wd, "/", prefix, ".Taxon_Bin.csv"), row.names = FALSE)
  # 保存每个bin的最高分类级别
  write.csv(icbin$Bin.TopClass, file = paste0(group_wd, "/", prefix, ".Bin_TopTaxon.csv"), row.names = FALSE)
}
# 初始化存储所有组合iCAMP结果的数据框
iCAMP_all_combined <- data.frame()
# 遍历每个组的文件夹
for (group_wd in dirs) {
  # 获取组名
  prefix <- basename(group_wd)
  # 读取每个bin的OTU信息
  taxon_bin_data <- read.csv(paste0(group_wd, "/", prefix, ".Taxon_Bin.csv"), header = TRUE, stringsAsFactors = FALSE)
  # 读取每个bin的生态过程信息
  mntd_data <- read.csv(paste0(group_wd, "/", prefix, ".ProcessImportance_EachBin_EachGroup.csv"), header = FALSE, stringsAsFactors = FALSE)
  # 选取Taxon_Bin的前三列(Bin, Top.Taxon, Top.Taxon.Frac)
  taxon_bin_data <- taxon_bin_data[, 1:3]
  # 转置生态过程数据并转换为数据框
  mntd_data <- as.data.frame(t(mntd_data[, -c(1:3)]))
  # 选取Bin列和CbMNTDiCbraya列(假设是第14列)
  mntd_data <- mntd_data[, c(1, 14)]
  # 设置列名
  colnames(mntd_data) <- as.character(mntd_data[1, ])
  # 移除第一行(已作为列名)
  mntd_data <- mntd_data[-1, ]
  # 设置最终列名
  colnames(mntd_data) <- c("Bin", "process")
  # 统一Bin的格式
  mntd_data$Bin <- gsub("bin", "Bin", mntd_data$Bin)
  # 合并taxon_bin_data和mntd_data
  iCAMP_data <- merge(taxon_bin_data, mntd_data, by = "Bin", all.x = TRUE)
  # 添加组名列
  iCAMP_data$group <- prefix
  # 合并到总数据框
  iCAMP_all_combined <- rbind(iCAMP_all_combined, iCAMP_data)
}
# 保存合并后的结果
write.table(iCAMP_all_combined, file = "D:/study/master/iCAMP/fungi/iCAMP_all_combined.csv",
            sep = ",", row.names = FALSE, quote = FALSE)
# ASV丰度转换为数据框
f_ASVs <- as.data.frame(f_ASV_df)
# 获取系统发育树并中点根化
f_tree <- phy_tree(f_merged)
f_tree <- midpoint(f_tree)
# 计算平均相对丰度
f_mean.RA <- data.frame(rowMeans(f_ASVs))
f_mean.RA$ASVs <- rownames(f_ASVs)
f_mean.RA$log.RA <- log1p(f_mean.RA$rowMeans.f_ASVs.)
colnames(f_mean.RA) <- c("mean.RA", "ASVs", "log.RA")
# 读取iCAMP合并结果
f_all.ASVs <- read.csv("D:/study/master/iCAMP/fungi/iCAMP_all_combined.csv")
# 筛选DL、HoS和DR过程的ASV
f_all.ASVs <- f_all.ASVs[f_all.ASVs$process == "HeS" | f_all.ASVs$process == "HoS" | f_all.ASVs$process == "DR",]
# 分类信息转换为数据框
f_tax <- as.data.frame(f_tax)
f_tax$ASV <- rownames(f_tax)
f_tax <- f_tax %>% filter(ASV %in% f_tree$tip.label)
# 筛选特定门的ASV
metadata <- as.data.frame(metadata)
f_ASV_df <- as.data.frame(f_ASV_df)
f_ASV_df$ASV <- rownames(f_ASV_df)
f_tax$ASV <- rownames(f_tax)
#在每个 Genus 分组，将同一个属下的所有 ASV 的计数在每个样本中加起来，得到每个属在每个样本的总计数。
f_dat_m <- f_ASV_df %>%
  left_join(f_tax)%>%#合并
  filter(Class != "c__fungip25") %>%
  group_by(Genus)%>%#按属分组
  filter(!(Genus == "g__uncultured" | stringr::str_detect(Genus, "g__unclassified")))%>%
  summarise_if(is.numeric, sum)%>%# 找到所有数值型的列，并对这些列计算它们的 sum
  na.omit()%>%
  pivot_longer(cols = !Genus)%>%# 从“宽”格式转换为“长”格式
  left_join(metadata, by = c("name" = "Sample.ID"))#合并
f_sum_all <- sum(f_dat_m$value)
#计算属相对丰度（单个属的菌/总菌）
f_dat_abun <- f_dat_m %>%
  group_by(Genus)%>%
  summarise(Abundance = sum(value)/f_sum_all)
colnames(f_dat_abun) <- c("Genus", "Abundance")
f_tax_sel <- f_tax%>%
  select(-c(ASV,Species))%>%
  distinct()%>%#移除重复的行
  filter(Genus %in% f_dat_abun$Genus)
#结合属、相对丰度、分类信息
f_dat_final <- f_dat_abun%>%
  left_join(f_tax_sel)
#优势菌门
f_phylumToPlot <- f_dat_final %>%
  group_by(Phylum)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(10)
f_phylumToPlot
p_asco_asvs <- f_tax[f_tax$Phylum == "p__Ascomycota", ]
p_mort_asvs <- f_tax[f_tax$Phylum == "p__Mortierellomycota", ]
p_basi_asvs <- f_tax[f_tax$Phylum == "p__Basidiomycota", ]
p_olpi_asvs <- f_tax[f_tax$Phylum == "p__Olpidiomycota", ]
p_muco_asvs <- f_tax[f_tax$Phylum == "p__Mucoromycota", ]
p_chyt_asvs <- f_tax[f_tax$Phylum == "p__Chytridiomycota", ]
p_zoop_asvs <- f_tax[f_tax$Phylum == "p__Zoopagomycota", ]
p_kick_asvs <- f_tax[f_tax$Phylum == "p__Kickxellomycota", ]
p_basid_asvs <- f_tax[f_tax$Phylum == "p__Basidiobolomycota", ] 
p_calc_asvs <- f_tax[f_tax$Phylum == "p__Calcarisporiellomycota", ]
# 创建ASV、Family和Phylum信息数据框
f_taxa_df <- data.frame(ASV = f_tree$tip.label, Family = rep('Others', length(f_tree$tip.label)),
                        Phylum = rep('Others', length(f_tree$tip.label)))
# 根据门信息更新Phylum列
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_asco_asvs)] <- 'Ascomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_mort_asvs)] <- 'Mortierellomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_basi_asvs)] <- 'Basidiomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_olpi_asvs)] <- 'Olpidiomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_muco_asvs)] <- 'Mucoromycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_chyt_asvs)] <- 'Chytridiomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_zoop_asvs)] <- 'Zoopagomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_kick_asvs)] <- 'Kickxellomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_basid_asvs)] <- 'Basidiobolomycota'
f_taxa_df$Phylum[f_taxa_df$ASV %in% rownames(p_calc_asvs)] <- 'Calcarisporiellomycota'
# 绘制带有分类学、热图和平均相对丰度的树
f_p = ggtree(f_tree, layout="fan", size=0.25, open.angle=10)
# 定义你想要的顺序
desired_group_order <- c("JRN", "JJN", "TZN", "PAN", "JRG", "JJG", "TZG", "PAG")
# 重新排序 f_all.ASVs$group 的因子水平
f_all.ASVs$group <- factor(f_all.ASVs$group, levels = desired_group_order)
# 添加门水平分类注释
f_p2 <- f_p +
  geom_fruit(data= f_taxa_df, geom=geom_tile,
             mapping=aes(y=ASV, fill=Phylum),
             width = 0.1, offset = 0.02) +
  scale_fill_manual(values = c('#C42D50', '#bbfa7b', '#FA907B','#84CCFA',
                               '#cf46aa','#fafa7b', '#3234D9','white','#529C7E',
                               '#46CF6B', '#7bbbfa'), guide = "none") + new_scale_fill() +
  # 添加iCAMP过程热图注释
  geom_fruit(
    data= f_all.ASVs,
    geom=geom_tile,
    mapping=aes(y=ID, x=group, fill=process),
    offset=0.08, pwidth=0.4
  ) +
  scale_fill_manual(
    values=c("#dfac03", "#006400", "#df3e03"),
    name = NULL,
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  )
# 添加平均相对丰度条形图注释
f_p3 <- f_p2 +
  new_scale_fill()+
  geom_fruit(data= f_mean.RA,
             geom=geom_bar,
             mapping=aes(y=ASVs, x=log.RA, fill=log.RA),
             offset = 0.07, orientation='y',
             stat="identity",
             pwidth = 0.4) +#相对总宽度
  scale_fill_viridis(option="D", guide = "none")
p4_fungi <- f_p3 +
  theme(
    # 设置图例位置，c(x, y) 向量表示在主图区域内的比例位置 (0到1)。
    # 对于圆图的右上角，大概是 x 接近 1，y 接近 1 的位置。
    legend.position = c(0.8, 0.85), # 示例：放置在图的右上方
    legend.justification = c(1, 1), # 右上对齐
    legend.box.just = "right",      # 图例框右对齐
    legend.background = element_rect(fill = "white", color = NA), # 可选：设置图例背景，避免被图覆盖
    legend.margin = margin(0, 0, 0, 0), # 可选：调整图例边距
    legend.text = element_text(size = 16)
  )
p4_fungi
#ggsave("D:/study/master/Main_Figure_tables/Figure_5/5a_phy_tree_fungi.png", plot = p4_fungi, width = 15, height = 15, dpi = 600, bg = "transparent")
# 假设 dirs 变量已经包含了所有组的文件夹路径
dirs <- list.dirs("D:/study/master/iCAMP/fungi", full.names = TRUE, recursive = FALSE)
# 初始化一个空的数据框，用于存储所有组的 ProcessImportance_EachTurnover 结果
f_all_turnovers_combined <- data.frame()
# 遍历每个组的文件夹
for (group_wd in dirs) {
  prefix <- basename(group_wd)
  turnover_file <- paste0(group_wd, "/", prefix, ".ProcessImportance_EachTurnover.csv") 
  if (file.exists(turnover_file)) {
    current_turnover_data <- read.csv(turnover_file, header = TRUE, stringsAsFactors = FALSE)
    # 定义包含过程得分的列名
    process_score_cols_to_check <- c("HeS", "HoS", "DL", "HD", "DR") 
    # 检查所有必需的过程得分列是否存在
    if (all(process_score_cols_to_check %in% colnames(current_turnover_data))) {
      current_turnover_data$group <- prefix
      # 'Process' 列将在下一步从得分中推导出来
      f_all_turnovers_combined <- rbind(f_all_turnovers_combined, current_turnover_data)
    } else {
      warning(paste0("警告: 文件 ", turnover_file, " 中未找到预期的过程得分列 (HeS, HoS, DL, HD, DR)。请检查文件头。跳过此文件。"))
    }
  } else {
    warning(paste0("警告: 未找到文件: ", turnover_file, "。跳过此文件。"))
  }
}
# 现在，为每个群落对推导主要过程
if (nrow(f_all_turnovers_combined) > 0) {
  # 定义包含过程得分的列名
  process_score_cols <- c("HeS", "HoS", "DL", "HD", "DR") 
  # 创建一个新的 'Process' 列，通过查找每行中得分最高的列名来赋值
  # apply 函数会逐行操作：对于每一行，找到最大值的列的索引，然后返回该列的名称
  f_all_turnovers_combined$Process <- apply(
    f_all_turnovers_combined[, process_score_cols], 1, # 针对这些列，按行操作
    function(x) {
      # 找到最大值的索引
      max_idx <- which.max(x)
      # 返回该索引对应的列名
      return(names(x)[max_idx])
    }
  )
  # 筛选出你感兴趣的特定过程 (DR, HeS, HoS)
  relevant_processes_data <- f_all_turnovers_combined[
    f_all_turnovers_combined$Process %in% c("DR", "HeS", "HoS"),
  ]
  # 检查过滤后是否有相关数据
  if (nrow(relevant_processes_data) == 0) {
    cat("警告：过滤后，没有找到 'DR', 'HeS', 'HoS' 这三种过程的群落对。这可能意味着其他过程是主要的驱动力。\n")
    # 如果你好奇其他过程的分类情况，可以打印出所有被分配到的过程类型
    cat("所有群落对分配到的独特过程类型有：\n")
    print(unique(f_all_turnovers_combined$Process))
    # 如果没有相关数据，就直接返回，避免后续统计报错
    return(NULL) 
  }
  # 统计每个过程的出现次数
  process_counts <- table(relevant_processes_data$Process) 
  # 计算这些相关过程的总群落对数量
  total_relevant_pairs <- sum(process_counts) 
  # 计算百分比
  # 为防止某个过程没有出现而导致计算错误，我们使用 if 语句进行检查，如果不存在则设为 0
  percentage_DR <- if ("DR" %in% names(process_counts)) (process_counts["DR"] / total_relevant_pairs) * 100 else 0
  percentage_HeS <- if ("HeS" %in% names(process_counts)) (process_counts["HeS"] / total_relevant_pairs) * 100 else 0
  percentage_HoS <- if ("HoS" %in% names(process_counts)) (process_counts["HoS"] / total_relevant_pairs) * 100 else 0
  # 打印结果
  cat("基于群落对的生态组装过程百分比:\n")
  cat("生态漂移 (DR): ", round(percentage_DR, 2), "% 的群落对\n")
  cat("异质选择 (HeS): ", round(percentage_HeS, 2), "% 的群落对\n")
  cat("同质选择 (HoS): ", round(percentage_HoS, 2), "% 的群落对\n")
  cat("\n原始计数:\n")
  print(process_counts) 
} else {
  cat("错误：未能成功读取或处理任何 ProcessImportance_EachTurnover.csv 文件。请检查文件路径和内容。\n")
}
# --- 确保 f_asv_process_abundance 包含 'group' 列 ---
# 假设 f_all.ASVs_processed 已包含 'group' 列，来源于 iCAMP_all_combined.csv。
# 如果 f_all.ASVs_processed 中没有 'group' 列，需在此合并 group 信息。
# 例如，若 f_all.ASVs 仅有 ID 和 process，而 group 信息在 metadata 中，则需先将 group 信息加入 f_all.ASVs_processed 中。
# 此处假设 f_all.ASVs_processed 已包含 'group' 列。
f_all.ASVs_processed <- f_all.ASVs
colnames(f_all.ASVs_processed)[colnames(f_all.ASVs_processed) == "ID"] <- "ASVs" # ASV 标识列。
# 合并 f_all.ASVs_processed（包含 ASV 主导过程和 group）与 f_mean.RA（包含 ASV 平均相对丰度）。
f_asv_process_abundance <- merge(f_all.ASVs_processed, f_mean.RA, by = "ASVs", all.x = TRUE)
# 移除合并后可能产生的 NA 值（通常不会有，但安全起见）。
f_asv_process_abundance <- na.omit(f_asv_process_abundance)
# --- 核心分组分析：按组计算各过程主导 ASV 的总相对丰度及其百分比贡献 ---
# 针对论文范例中“在每个山脉范围内”的分析。
f_process_abundance_summary_by_group <- f_asv_process_abundance %>%
  dplyr::group_by(group, process) %>% # 按组分组，再按过程分组。
  dplyr::summarise(TotalRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>% # 计算该组该过程下所有 ASV 的总相对丰度。
  dplyr::ungroup() %>% # 第一次 ungroup()，为下一步按 group 计算百分比。
  dplyr::group_by(group) %>% # 再次分组，确保百分比是组内的。
  dplyr::mutate(PercentageContribution = (TotalRelativeAbundance / sum(TotalRelativeAbundance)) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, desc(PercentageContribution)) # 按组和百分比降序排列。
print("\n--- 基于 ASV 平均相对丰度，各组各生态组装过程主导的 ASV 贡献百分比 ---")
print(f_process_abundance_summary_by_group)
# --- 按组识别特定过程下的关键属 ---
# 确保 f_tax 数据框可以用于合并（ASVs 作为一列）。
f_tax_df_for_merge <- as.data.frame(f_tax)
f_tax_df_for_merge$ASVs <- rownames(f_tax_df_for_merge) # 匹配 'ASVs' 列名。
# 合并分类学信息到包含 ASV 过程和丰度的数据框（f_asv_process_abundance_taxa 已包含 group 列）。
f_asv_process_abundance_taxa <- merge(f_asv_process_abundance, f_tax_df_for_merge, by = "ASVs", all.x = TRUE)
print("\n--- 按组识别特定过程下的关键属 ---")
# 示例：查找在“同质选择（HoS）”主导下，各组丰度最高的属。
f_top_hos_genera_by_group <- f_asv_process_abundance_taxa %>%
  dplyr::filter(process == "HoS") %>%
  dplyr::group_by(group, Genus) %>% # 按组和属分组。
  dplyr::summarise(GenusRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>%
  dplyr::arrange(group, desc(GenusRelativeAbundance)) %>%
  dplyr::slice_head(n = 5) # 对每个组取前 5 个。
print("\n同质选择 (HoS) 主导下的各组主要属（按丰度降序）:")
print(f_top_hos_genera_by_group)
# 示例：查找在“生态漂移（DR）”主导下，各组丰度最高的属。
f_top_dr_genera_by_group <- f_asv_process_abundance_taxa %>%
  dplyr::filter(process == "DR") %>%
  dplyr::group_by(group, Genus) %>% # 按组和属分组。
  dplyr::summarise(GenusRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>%
  dplyr::arrange(group, desc(GenusRelativeAbundance)) %>%
  dplyr::slice_head(n = 5)
print("\n生态漂移 (DR) 主导下的各组主要属（按丰度降序）:")
print(f_top_dr_genera_by_group)
# 示例：查找在“异质选择（HeS）”主导下，各组丰度最高的属。
f_top_hes_genera_by_group <- f_asv_process_abundance_taxa %>%
  dplyr::filter(process == "HeS") %>%
  dplyr::group_by(group, Genus) %>% # 按组和属分组。
  dplyr::summarise(GenusRelativeAbundance = sum(mean.RA), .groups = 'drop_last') %>%
  dplyr::arrange(group, desc(GenusRelativeAbundance)) %>%
  dplyr::slice_head(n = 5)
print("\n异质选择 (HeS) 主导下的各组主要属（按丰度降序）:")
print(f_top_hes_genera_by_group)
# --- 计算各过程贡献百分比的中位数和 IQR (HoS, DR, HeS) ---
# 计算同质选择 (HoS) 的贡献百分比中位数和 IQR。
f_hos_contributions <- f_process_abundance_summary_by_group %>%
  dplyr::filter(process == "HoS") %>%
  dplyr::pull(PercentageContribution) # 提取 PercentageContribution 列。
f_median_hos <- median(f_hos_contributions, na.rm = TRUE)
f_q1_hos <- quantile(f_hos_contributions, probs = 0.25, na.rm = TRUE)
f_q3_hos <- quantile(f_hos_contributions, probs = 0.75, na.rm = TRUE)
cat(paste0("同质选择 (HoS) 主导的 ASV 贡献了: ", round(f_median_hos, 2), "% ",
           "(IQR: ", round(f_q1_hos, 2), "–", round(f_q3_hos, 2), "%)",
           " 的相对丰度，在每个组中。\n"))
# 计算生态漂移 (DR) 的贡献百分比中位数和 IQR。
f_dr_contributions <- f_process_abundance_summary_by_group %>%
  dplyr::filter(process == "DR") %>%
  dplyr::pull(PercentageContribution) # 提取 PercentageContribution 列。
f_median_dr <- median(f_dr_contributions, na.rm = TRUE)
f_q1_dr <- quantile(f_dr_contributions, probs = 0.25, na.rm = TRUE)
f_q3_dr <- quantile(f_dr_contributions, probs = 0.75, na.rm = TRUE)
cat(paste0("生态漂移 (DR) 主导的 ASV 贡献了: ", round(f_median_dr, 2), "% ",
           "(IQR: ", round(f_q1_dr, 2), "–", round(f_q3_dr, 2), "%)",
           " 的相对丰度，在每个组中。\n"))
# 计算异质选择 (HeS) 的贡献百分比中位数和 IQR。
f_hes_contributions <- f_process_abundance_summary_by_group %>%
  dplyr::filter(process == "HeS") %>%
  dplyr::pull(PercentageContribution) # 提取 PercentageContribution 列。
f_median_hes <- median(f_hes_contributions, na.rm = TRUE)
f_q1_hes <- quantile(f_hes_contributions, probs = 0.25, na.rm = TRUE)
f_q3_hes <- quantile(f_hes_contributions, probs = 0.75, na.rm = TRUE)
cat(paste0("异质选择 (HeS) 主导的 ASV 贡献了: ", round(f_median_hes, 2), "% ",
           "(IQR: ", round(f_q1_hes, 2), "–", round(f_q3_hes, 2), "%)",
           " 的相对丰度，在每个组中。\n"))

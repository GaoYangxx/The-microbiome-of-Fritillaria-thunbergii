#图5c
install.packages(c("picante", "adespatial", "speedyseq", "phyloseq", "ape",
                   "progress", "tibble", "qiime2R", "readxl", "dplyr", "tidyr"))
library(picante)
library(adespatial)
library(speedyseq)
library(phyloseq)
library(ape)
library(progress)
library(tibble)
library(qiime2R)
library(readxl)
library(dplyr) 
library(tidyr)
#细菌系统发育β多样性
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
# 读取ASV丰度表
b_ASVs <- data.frame(b_ASV)
# 读取样本元数据表
samp <- data.frame(metadata)
# 读取系统发育树
tree <- phy_tree(b_merged)
# beta多样性的系统发育结构
# 计算每个区域的平均ASV丰度
regs <- c()
for(j in levels(factor(samp$Group))){ # 遍历每个区域
  reg.j<-samp[samp$Group==j,] # 筛选当前区域的样本元数据
  n.j<-rowMeans(b_ASVs[reg.j$Sample.ID]) # 计算当前区域每个ASV的平均丰度
  regs<-cbind(regs,n.j) # 将当前区域的平均丰度添加到regs
}
regs <- data.frame(regs) # 将regs转换为数据框
colnames(regs) <- levels(factor(samp$Group)) # 设置regs的列名为区域名称
# 计算每个区域的平均最近类群距离(MNTD)
MNTDs <- c()
for(k in levels(factor(samp$Group))){ # 遍历每个区域
  reg.k<-samp[samp$Group==k,] # 筛选当前区域的样本元数据
  n.k<-b_ASVs[reg.k$Sample.ID] # 获取当前区域的ASV丰度
  n.k<-n.k[rowSums(n.k)>0,] # 移除在当前区域中丰度为0的ASV
  tree.k<-prune.sample(t(n.k), tree) # 修剪树，只保留当前区域的ASV
  mntd.k<-mntd(t(n.k), cophenetic(tree.k)) # 计算当前区域的MNTD
  mntd.k<-data.frame(mntd.k) # 将MNTD转换为数据框
  mntd.k$reg<-k # 添加区域标签
  MNTDs<-rbind(MNTDs, mntd.k) # 合并MNTD结果
}
# 绘制MNTD箱线图
boxplot(MNTDs$mntd.k~MNTDs$reg, las=2, xlab="", ylab="MNTD")
# 添加所有区域MNTD平均值的水平线
abline(h=mean(na.omit(as.numeric(MNTDs$mntd.k))), col="red", lwd=2)
# 跨区域的beta多样性谱
b_across.part <- c()
b.part.across.all <- beta.div.comp(t(regs), coef="S", quant=T) # 计算跨区域的beta多样性分解
b_across.part <- rbind(b_across.part, c(b.part.across.all$part, nrow(regs),0)) # 存储初始结果
pb <- progress_bar$new(total = length(seq(0.005,0.2,0.005))) # 初始化进度条
regs.phy <- phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree)) # 创建phyloseq对象
for(h in seq(0.005,0.2,0.005)){ # 遍历不同的系统发育聚合深度
  glom.h<-speedyseq::tree_glom(regs.phy, h) # 在当前深度聚合系统发育树
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
  b.part.across<-beta.div.comp(t(OTU.h), coef="S", quant=T) # 计算当前聚合水平的beta多样性分解
  b_across.part<-rbind(b_across.part, c(b.part.across$part, nrow(OTU.h), h)) # 存储结果
  pb$tick() # 更新进度条
}
b_across.part <- data.frame(b_across.part) # 转换为数据框
colnames(b_across.part)[6:7] <- c("tips","phy.agglom") # 设置列名
# write.csv(b_across.part, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_across.part.csv", row.names = FALSE)
b_across.part <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_across.part.csv", sep=",", header=TRUE, check.names = FALSE)
# 绘制beta多样性随系统发育聚合深度的变化
plot(b_across.part$phy.agglom, b_across.part$BDtotal,type="b",lwd=1.75, cex=1, xlim=c(0,0.2), ylim=c(0,0.5), xlab="phylogenetic agglomeration depth", ylab="beta diversity")
# 添加MNTD平均值的垂直线
abline(v=mean(na.omit(MNTDs$mntd.k)), lty=1, lwd=2)
# 区域内的beta多样性谱
samp2 <- samp 
b_regs.parts <- c()
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本元数据
  reg.a<-b_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除在当前区域中丰度为0的ASV
  b.part.a<-beta.div.comp(t(reg.a), coef="S", quant=T) # 计算区域内的beta多样性分解
  b_regs.parts<-rbind(b_regs.parts, c(b.part.a$part, nrow(reg.a), nrow(reg.a)/nrow(reg.a), h=0, a))} # 存储结果
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本元数据
  reg.a<-b_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除在当前区域中丰度为0的ASV
  reg.a.tree<-prune.sample(t(reg.a),tree) # 修剪树，只保留当前区域的ASV
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree)) # 创建phyloseq对象
  for(h in seq(0.005,0.2,0.005)){ # 遍历不同的系统发育聚合深度
    glom.h<-speedyseq::tree_glom(reg.a.phy, h) # 在当前深度聚合系统发育树
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
    b.part<-beta.div.comp(t(OTU.h), coef="S", quant=T) # 计算当前聚合水平的beta多样性分解
    b_regs.parts<-rbind(b_regs.parts, c(b.part$part, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a)) # 存储结果
    print(c(a,h))} # 打印当前区域和聚合深度
}
# 将区域内beta多样性分解结果转换为数据框
b_regs.parts <- data.frame(b_regs.parts)
# 设置b_regs.parts的列名
colnames(b_regs.parts)[6:9] <- c("n_tips","perc_tips","phy_glom","region")
# write.csv(b_regs.parts, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_regs.parts.csv", row.names = FALSE)
b_regs.parts <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_regs.parts.csv", sep=",", header=TRUE, check.names = FALSE)
# 设置文件路径和名称
file_path <- "D:/study/master/Main_Figure_tables/Figure_5/5c_phylogenetic_beta_bacteria.png"
# 打开PNG图形设备
png(filename = file_path,
    width = 8,        # 图片宽度，单位是英寸 (默认)
    height = 8,       # 图片高度，单位是英寸 (默认)
    units = "in",     # 明确指定单位为英寸
    res = 600,        # 分辨率，即 DPI (dots per inch)
    bg = "transparent" # 设置背景为透明
)
#调整字体大小
par(mar = c(5, 5, 2, 2) + 0.1, cex.lab = 1.4, cex.axis = 1.2) 
# 绘制区域内beta多样性随系统发育聚合深度的变化
plot(b_regs.parts$phy_glom, b_regs.parts$BDtotal, type = "n",
     ylim = c(0, 6),
     xlab = "Bacterial phylogenetic agglomeration depth (a.u.)",
     ylab = "Bacterial beta-diversity change (%)",
     axes = FALSE, # 禁用默认的坐标轴
     bty = "l"    #设置图的边框类型为 "l" (只显示左和下)
)
# 重新绘制坐标轴和刻度
# x轴
axis(side = 1) # 绘制x轴 (底部)
# y轴，并设置 las = 2 使刻度文本垂直
axis(side = 2, las = 2) # 绘制y轴 (左侧)，las = 2 让刻度标签垂直
# 添加图的盒子，确保只有左和下
box(bty = "l") # 重新绘制边框，只保留L型
# 遍历每个区域，绘制beta多样性变化
for(i in levels(factor(b_regs.parts$region))){
  reg.i <- b_regs.parts[b_regs.parts$region == i,]
  points(reg.i$phy_glom[1:40], abs(diff(as.numeric(reg.i$BDtotal)))/as.numeric(reg.i$BDtotal)[1]*100,
         pch = 16, cex = 1.1, col = "grey")
}
# 绘制跨区域beta多样性变化
points(b_across.part$phy.agglom[1:40], abs(diff(b_across.part$BDtotal)/b_across.part$BDtotal[1])*100,
       pch = 16, cex = 1.1, col = "red", type = "b")
# 添加垂直线
mean(MNTDs$mntd.k)
abline(v = 0.1732129, lty = 5, col = "darkgrey")
# 拟合跨区域beta多样性变化的指数模型
x <- as.numeric(b_across.part$phy.agglom)[1:40]
y <- abs(diff(b_across.part$BDtotal)/b_across.part$BDtotal[1])*100
m2 <- nls(y ~ a * exp(b * x), start = list(a = 8, b = -50))
# 绘制拟合的指数曲线
lines(x, predict(m2), col = "red", lwd = 3)
# 关闭图形设备，将图保存到文件
dev.off()
## 系统发育独特性分析
#跨区域
b_across.unique <-c(sum(rowSums(regs>0)==1), nrow(regs),1,0) # 计算跨区域的独特ASV数量
regs.phy<-phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree)) # 创建phyloseq对象
for(h in seq(0.005,0.2,0.005)){ # 遍历系统发育聚合深度
  glom.h<-speedyseq::tree_glom(regs.phy, h) # 聚合系统发育树
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
  b_across.unique<-rbind(b_across.unique, c(sum(rowSums(OTU.h>0)==1), nrow(OTU.h),nrow(OTU.h)/nrow(regs),h)) # 计算并存储独特ASV数量
  print(c(h))}
b_across.unique <-data.frame(b_across.unique) # 转换为数据框
colnames(b_across.unique) <-c("uniques","n_tips","n_total","h") # 设置列名
# write.csv(b_across.unique, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_across.unique.csv", row.names = TRUE)
b_across.unique <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_across.unique.csv", sep=",", header=TRUE, row.names = 1, check.names = FALSE)
#区域内
b_regs.uniques <-c()
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本
  reg.a<-b_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除丰度为0的ASV
  reg.a.uni<-sum(rowSums(reg.a>0)==1) # 计算当前区域的独特ASV数量
  b_regs.uniques<-rbind(b_regs.uniques, c(reg.a.uni, nrow(reg.a), 1, 0, a))} # 存储结果
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本
  reg.a<-b_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除丰度为0的ASV
  reg.a.tree<-prune.sample(t(reg.a),tree) # 修剪系统发育树
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree)) # 创建phyloseq对象
  for(h in seq(0.005,0.2,0.005)){ # 遍历系统发育聚合深度
    glom.h<-speedyseq::tree_glom(reg.a.phy, h) # 聚合系统发育树
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
    reg.uni<-sum(rowSums(OTU.h>0)==1) # 计算当前聚合水平的独特ASV数量
    b_regs.uniques<-rbind(b_regs.uniques, c(reg.uni, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a)) # 存储结果
    print(c(a,h))}
}
b_regs.uniques <- data.frame(b_regs.uniques) # 转换为数据框
colnames(b_regs.uniques) <-c("uniq_tips","all_tips", "perc_tips","h","Group") # 设置列名
# write.csv(b_regs.uniques, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_regs.uniques.csv", row.names = FALSE)
b_regs.uniques <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_b_regs.uniques.csv", sep=",", header=TRUE, check.names = FALSE)
# 绘制独特类群比例随聚合水平的变化
plot(b_regs.uniques$h, b_regs.uniques$perc_tips, type="n", xlab="agglomeration level",ylab="unique taxa (%)", ylim=c(0,100))
for(z in levels(factor(b_regs.uniques$Group))){ # 遍历每个区域
  b_regs.uniques.z<-b_regs.uniques[b_regs.uniques$Group==z,]
  points(b_regs.uniques.z$h, as.numeric(b_regs.uniques.z$perc_tips)*100,pch=16, col="grey", cex=1)} #绘制每个区域的点
# 绘制跨区域的独特类群比例
points(b_across.unique$h, b_across.unique$uniques/max(b_across.unique$uniques)*100, type="b", col="red", pch=16)
# 计算所有区域的平均 MNTD
b_avg_mntd <- mean(na.omit(MNTDs$mntd.k))
print(paste("所有区域的平均 MNTD:", round(b_avg_mntd, 7))) # 确保这个值是 0.1732129
# 获取起始点 (聚合深度最小) 的 Beta 多样性
b_min_phy_glom <- min(b_regs.parts$phy_glom)
# 筛选出最小聚合深度的数据
b_start_bd_data <- b_regs.parts %>%
  filter(phy_glom == b_min_phy_glom)
# 计算起始点的中位数和 IQR
b_start_median <- median(b_start_bd_data$BDtotal, na.rm = TRUE)
b_start_iqr <- quantile(b_start_bd_data$BDtotal, probs = c(0.25, 0.75), na.rm = TRUE)
cat(paste0("起始点 (phy_glom=", round(b_min_phy_glom, 3), ") 的 Beta 多样性:\n",
           "  中位数: ", round(b_start_median, 3), "\n",
           "  IQR: ", round(b_start_iqr[1], 3), " - ", round(b_start_iqr[2], 3), "\n\n"))
# 获取结束点 (聚合深度最大) 的 Beta 多样性
b_max_phy_glom <- max(b_regs.parts$phy_glom)
# 筛选出最大聚合深度的数据
b_end_bd_data <- b_regs.parts %>%
  filter(phy_glom == b_max_phy_glom)
# 计算结束点的中位数和 IQR
b_end_median <- median(b_end_bd_data$BDtotal, na.rm = TRUE)
b_end_iqr <- quantile(b_end_bd_data$BDtotal, probs = c(0.25, 0.75), na.rm = TRUE)
cat(paste0("结束点 (phy_glom=", round(b_max_phy_glom, 3), ") 的 Beta 多样性:\n",
           "  中位数: ", round(b_end_median, 3), "\n",
           "  IQR: ", round(b_end_iqr[1], 3), " - ", round(b_end_iqr[2], 3), "\n\n"))
# 计算在平均 MNTD 之前的下降百分比
# 找到最接近平均 MNTD 的那个 phy_glom 值
b_phy_glom_at_avg_mntd <- b_regs.parts %>%
  filter(phy_glom <= b_avg_mntd) %>%
  pull(phy_glom) %>%
  max() # 取小于等于 b_avg_mntd 中的最大 glom 值，作为切割点
# 筛选出从起始点到 b_phy_glom_at_avg_mntd 这段范围内的数据
b_data_before_avg_mntd <- b_regs.parts %>%
  filter(phy_glom <= b_phy_glom_at_avg_mntd)
# 对每个区域计算其在该段范围内的 Beta 多样性下降百分比
# 从每个区域的起始 Beta 多样性（在 b_min_phy_glom 处）到在 b_phy_glom_at_avg_mntd 处的 Beta 多样性
b_descent_percentages <- b_data_before_avg_mntd %>%
  group_by(region) %>%
  summarise(
    b_bd_at_start = BDtotal[which.min(phy_glom)], # 每个区域在起始点的BDtotal
    b_bd_at_avg_mntd_depth = BDtotal[which.max(phy_glom)], # 每个区域在MNTD平均深度点的BDtotal
    # 计算下降百分比
    b_descent_percent = ((b_bd_at_start - b_bd_at_avg_mntd_depth) / b_bd_at_start) * 100
  ) %>%
  ungroup() %>%
  pull(b_descent_percent) # 提取所有区域的下降百分比
# 计算下降百分比的中位数和 IQR
b_descent_median <- median(b_descent_percentages, na.rm = TRUE)
b_descent_iqr <- quantile(b_descent_percentages, probs = c(0.25, 0.75), na.rm = TRUE)
cat(paste0("区域内 Beta 多样性在系统发育距离短于平均 MNTD (", round(b_avg_mntd, 7), ") 时的下降百分比:\n",
           "  中位数: ", round(b_descent_median, 1), "% \n", # 保留一位小数
           "  IQR: ", round(b_descent_iqr[1], 1), "% - ", round(b_descent_iqr[2], 1), "% \n"))
#真菌系统发育β多样性
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
# 读取ASV丰度表
f_ASVs <- data.frame(f_ASV)
# 读取样本元数据表
samp <- data.frame(metadata)
# 读取系统发育树
tree <- phy_tree(f_merged)
# beta多样性的系统发育结构
# 计算每个区域的平均ASV丰度
regs <- c()
for(j in levels(factor(samp$Group))){ # 遍历每个区域
  reg.j<-samp[samp$Group==j,] # 筛选当前区域的样本元数据
  n.j<-rowMeans(f_ASVs[reg.j$Sample.ID]) # 计算当前区域每个ASV的平均丰度
  regs<-cbind(regs,n.j) # 将当前区域的平均丰度添加到regs
}
regs <- data.frame(regs) # 将regs转换为数据框
colnames(regs) <- levels(factor(samp$Group)) # 设置regs的列名为区域名称
# 计算每个区域的平均最近类群距离(MNTD)
MNTDs <- c()
for(k in levels(factor(samp$Group))){ # 遍历每个区域
  reg.k<-samp[samp$Group==k,] # 筛选当前区域的样本元数据
  n.k<-f_ASVs[reg.k$Sample.ID] # 获取当前区域的ASV丰度
  n.k<-n.k[rowSums(n.k)>0,] # 移除在当前区域中丰度为0的ASV
  tree.k<-prune.sample(t(n.k), tree) # 修剪树，只保留当前区域的ASV
  mntd.k<-mntd(t(n.k), cophenetic(tree.k)) # 计算当前区域的MNTD
  mntd.k<-data.frame(mntd.k) # 将MNTD转换为数据框
  mntd.k$reg<-k # 添加区域标签
  MNTDs<-rbind(MNTDs, mntd.k) # 合并MNTD结果
}
# 绘制MNTD箱线图
boxplot(MNTDs$mntd.k~MNTDs$reg, las=2, xlab="", ylab="MNTD")
# 添加所有区域MNTD平均值的水平线
abline(h=mean(na.omit(as.numeric(MNTDs$mntd.k))), col="red", lwd=2)
# 跨区域的beta多样性谱
f_across.part <- c()
b.part.across.all <- beta.div.comp(t(regs), coef="S", quant=T) # 计算跨区域的beta多样性分解
f_across.part <- rbind(f_across.part, c(b.part.across.all$part, nrow(regs),0)) # 存储初始结果
pb <- progress_bar$new(total = length(seq(0.005,0.2,0.005))) # 初始化进度条
regs.phy <- phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree)) # 创建phyloseq对象
for(h in seq(0.005,0.2,0.005)){ # 遍历不同的系统发育聚合深度
  glom.h<-speedyseq::tree_glom(regs.phy, h) # 在当前深度聚合系统发育树
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
  b.part.across<-beta.div.comp(t(OTU.h), coef="S", quant=T) # 计算当前聚合水平的beta多样性分解
  f_across.part<-rbind(f_across.part, c(b.part.across$part, nrow(OTU.h), h)) # 存储结果
  pb$tick() # 更新进度条
}
f_across.part <- data.frame(f_across.part) # 转换为数据框
colnames(f_across.part)[6:7] <- c("tips","phy.agglom") # 设置列名
# write.csv(f_across.part, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_across.part.csv", row.names = FALSE)
f_across.part <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_across.part.csv", sep=",", header=TRUE, check.names = FALSE)
# 绘制beta多样性随系统发育聚合深度的变化
plot(f_across.part$phy.agglom, f_across.part$BDtotal,type="b",lwd=1.75, cex=1, xlim=c(0,0.2), ylim=c(0,0.5), xlab="phylogenetic agglomeration depth", ylab="beta diversity")
# 添加MNTD平均值的垂直线
abline(v=mean(na.omit(MNTDs$mntd.k)), lty=1, lwd=2)
# 区域内的beta多样性谱
samp2 <- samp 
f_regs.parts <- c()
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本元数据
  reg.a<-f_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除在当前区域中丰度为0的ASV
  b.part.a<-beta.div.comp(t(reg.a), coef="S", quant=T) # 计算区域内的beta多样性分解
  f_regs.parts<-rbind(f_regs.parts, c(b.part.a$part, nrow(reg.a), nrow(reg.a)/nrow(reg.a), h=0, a))} # 存储结果
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本元数据
  reg.a<-f_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除在当前区域中丰度为0的ASV
  reg.a.tree<-prune.sample(t(reg.a),tree) # 修剪树，只保留当前区域的ASV
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree)) # 创建phyloseq对象
  for(h in seq(0.005,0.2,0.005)){ # 遍历不同的系统发育聚合深度
    glom.h<-speedyseq::tree_glom(reg.a.phy, h) # 在当前深度聚合系统发育树
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
    b.part<-beta.div.comp(t(OTU.h), coef="S", quant=T) # 计算当前聚合水平的beta多样性分解
    f_regs.parts<-rbind(f_regs.parts, c(b.part$part, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a)) # 存储结果
    print(c(a,h))} # 打印当前区域和聚合深度
}
# 将区域内beta多样性分解结果转换为数据框
f_regs.parts <- data.frame(f_regs.parts)
# 设置f_regs.parts的列名
colnames(f_regs.parts)[6:9] <- c("n_tips","perc_tips","phy_glom","region")
# write.csv(f_regs.parts, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_regs.parts.csv", row.names = FALSE)
f_regs.parts <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_regs.parts.csv", sep=",", header=TRUE, check.names = FALSE)
# 设置文件路径和名称
file_path <- "D:/study/master/Main_Figure_tables/Figure_5/5c_phylogenetic_beta_fungi.png"
# 打开PNG图形设备
png(filename = file_path,
    width = 8,        # 图片宽度，单位是英寸 (默认)
    height = 8,       # 图片高度，单位是英寸 (默认)
    units = "in",     # 明确指定单位为英寸
    res = 600,        # 分辨率，即 DPI (dots per inch)
    bg = "transparent" # 设置背景为透明
)
#调整字体大小
par(mar = c(5, 5, 2, 2) + 0.1, cex.lab = 1.4, cex.axis = 1.2)
# 绘制区域内beta多样性随系统发育聚合深度的变化
plot(f_regs.parts$phy_glom, f_regs.parts$BDtotal, type = "n",
     ylim = c(0, 6.5),
     xlab = "Fungal phylogenetic agglomeration depth (a.u.)",
     ylab = "Fungal beta-diversity change (%)",
     axes = FALSE, # 禁用默认的坐标轴
     bty = "l"    #设置图的边框类型为 "l" (只显示左和下)
)
# 重新绘制坐标轴和刻度
# x轴
axis(side = 1) # 绘制x轴 (底部)
# y轴，并设置 las = 2 使刻度文本垂直
axis(side = 2, las = 2) # 绘制y轴 (左侧)，las = 2 让刻度标签垂直
# 添加图的盒子，确保只有左和下
box(bty = "l") # 重新绘制边框，只保留L型
# 遍历每个区域，绘制beta多样性变化
for(i in levels(factor(f_regs.parts$region))){
  reg.i <- f_regs.parts[f_regs.parts$region == i,]
  points(reg.i$phy_glom[1:40], abs(diff(as.numeric(reg.i$BDtotal)))/as.numeric(reg.i$BDtotal)[1]*100,
         pch = 16, cex = 1.1, col = "grey")
}
# 绘制跨区域beta多样性变化
points(f_across.part$phy.agglom[1:40], abs(diff(f_across.part$BDtotal)/f_across.part$BDtotal[1])*100,
       pch = 16, cex = 1.1, col = "red", type = "b")
# 添加垂直线
mean(MNTDs$mntd.k)
abline(v = 0.1708501, lty = 5, col = "darkgrey")
# 拟合跨区域beta多样性变化的指数模型
x <- as.numeric(f_across.part$phy.agglom)[1:40]
y <- abs(diff(f_across.part$BDtotal)/f_across.part$BDtotal[1])*100
m2 <- nls(y ~ a * exp(b * x), start = list(a = 8, b = -50))
# 绘制拟合的指数曲线
lines(x, predict(m2), col = "red", lwd = 3)
# 关闭图形设备，将图保存到文件
dev.off()
## 系统发育独特性分析
#跨区域
f_across.unique <-c(sum(rowSums(regs>0)==1), nrow(regs),1,0) # 计算跨区域的独特ASV数量
regs.phy<-phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree)) # 创建phyloseq对象
for(h in seq(0.005,0.2,0.005)){ # 遍历系统发育聚合深度
  glom.h<-speedyseq::tree_glom(regs.phy, h) # 聚合系统发育树
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
  f_across.unique<-rbind(f_across.unique, c(sum(rowSums(OTU.h>0)==1), nrow(OTU.h),nrow(OTU.h)/nrow(regs),h)) # 计算并存储独特ASV数量
  print(c(h))}
f_across.unique <-data.frame(f_across.unique) # 转换为数据框
colnames(f_across.unique) <-c("uniques","n_tips","n_total","h") # 设置列名
# write.csv(f_across.unique, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_across.unique.csv", row.names = TRUE)
f_across.unique <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_across.unique.csv", sep=",", header=TRUE, row.names = 1, check.names = FALSE)
#区域内
f_regs.uniques <-c()
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本
  reg.a<-f_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除丰度为0的ASV
  reg.a.uni<-sum(rowSums(reg.a>0)==1) # 计算当前区域的独特ASV数量
  f_regs.uniques<-rbind(f_regs.uniques, c(reg.a.uni, nrow(reg.a), 1, 0, a))} # 存储结果
for(a in levels(factor(samp2$Group))){ # 遍历每个区域
  samp.a<-samp2[samp2$Group==a,] # 筛选当前区域的样本
  reg.a<-f_ASVs[samp.a$Sample.ID] # 获取当前区域的ASV丰度
  reg.a<-reg.a[rowSums(reg.a)>0,] # 移除丰度为0的ASV
  reg.a.tree<-prune.sample(t(reg.a),tree) # 修剪系统发育树
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree)) # 创建phyloseq对象
  for(h in seq(0.005,0.2,0.005)){ # 遍历系统发育聚合深度
    glom.h<-speedyseq::tree_glom(reg.a.phy, h) # 聚合系统发育树
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix") # 获取聚合后的OTU表
    reg.uni<-sum(rowSums(OTU.h>0)==1) # 计算当前聚合水平的独特ASV数量
    f_regs.uniques<-rbind(f_regs.uniques, c(reg.uni, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a)) # 存储结果
    print(c(a,h))}
}
f_regs.uniques <- data.frame(f_regs.uniques) # 转换为数据框
colnames(f_regs.uniques) <-c("uniq_tips","all_tips", "perc_tips","h","Group") # 设置列名
# write.csv(f_regs.uniques, file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_regs.uniques.csv", row.names = FALSE)
f_regs.uniques <- read.csv(file = "D:/study/master/Main_Figure_tables/Figure_5/5c_f_regs.uniques.csv", sep=",", header=TRUE, check.names = FALSE)
# 绘制独特类群比例随聚合水平的变化
plot(f_regs.uniques$h, f_regs.uniques$perc_tips, type="n", xlab="agglomeration level",ylab="unique taxa (%)", ylim=c(0,100))
for(z in levels(factor(f_regs.uniques$Group))){ # 遍历每个区域
  f_regs.uniques.z<-f_regs.uniques[f_regs.uniques$Group==z,]
  points(f_regs.uniques.z$h, as.numeric(f_regs.uniques.z$perc_tips)*100,pch=16, col="grey", cex=1)} #绘制每个区域的点
# 绘制跨区域的独特类群比例
points(f_across.unique$h, f_across.unique$uniques/max(f_across.unique$uniques)*100, type="b", col="red", pch=16)
# 计算所有区域的平均 MNTD
f_avg_mntd <- mean(na.omit(MNTDs$mntd.k))
print(paste("所有区域的平均 MNTD:", round(f_avg_mntd, 7))) # 确保这个值是 0.1708501
# 获取起始点 (聚合深度最小) 的 Beta 多样性
f_min_phy_glom <- min(f_regs.parts$phy_glom)
# 筛选出最小聚合深度的数据
f_start_bd_data <- f_regs.parts %>%
  filter(phy_glom == f_min_phy_glom)
# 计算起始点的中位数和 IQR
f_start_median <- median(f_start_bd_data$BDtotal, na.rm = TRUE)
f_start_iqr <- quantile(f_start_bd_data$BDtotal, probs = c(0.25, 0.75), na.rm = TRUE)
cat(paste0("起始点 (phy_glom=", round(f_min_phy_glom, 3), ") 的 Beta 多样性:\n",
           "  中位数: ", round(f_start_median, 3), "\n",
           "  IQR: ", round(f_start_iqr[1], 3), " - ", round(f_start_iqr[2], 3), "\n\n"))
# 获取结束点 (聚合深度最大) 的 Beta 多样性
f_max_phy_glom <- max(f_regs.parts$phy_glom)
# 筛选出最大聚合深度的数据
f_end_bd_data <- f_regs.parts %>%
  filter(phy_glom == f_max_phy_glom)
# 计算结束点的中位数和 IQR
f_end_median <- median(f_end_bd_data$BDtotal, na.rm = TRUE)
f_end_iqr <- quantile(f_end_bd_data$BDtotal, probs = c(0.25, 0.75), na.rm = TRUE)
cat(paste0("结束点 (phy_glom=", round(f_max_phy_glom, 3), ") 的 Beta 多样性:\n",
           "  中位数: ", round(f_end_median, 3), "\n",
           "  IQR: ", round(f_end_iqr[1], 3), " - ", round(f_end_iqr[2], 3), "\n\n"))
# 计算在平均 MNTD 之前的下降百分比
# 找到最接近平均 MNTD 的那个 phy_glom 值
f_phy_glom_at_avg_mntd <- f_regs.parts %>%
  filter(phy_glom <= f_avg_mntd) %>%
  pull(phy_glom) %>%
  max() # 取小于等于 f_avg_mntd 中的最大 glom 值，作为切割点
# 筛选出从起始点到 f_phy_glom_at_avg_mntd 这段范围内的数据
f_data_before_avg_mntd <- f_regs.parts %>%
  filter(phy_glom <= f_phy_glom_at_avg_mntd)
# 对每个区域计算其在该段范围内的 Beta 多样性下降百分比
# 从每个区域的起始 Beta 多样性（在 f_min_phy_glom 处）到在 f_phy_glom_at_avg_mntd 处的 Beta 多样性
f_descent_percentages <- f_data_before_avg_mntd %>%
  group_by(region) %>%
  summarise(
    f_bd_at_start = BDtotal[which.min(phy_glom)], # 每个区域在起始点的BDtotal
    f_bd_at_avg_mntd_depth = BDtotal[which.max(phy_glom)], # 每个区域在MNTD平均深度点的BDtotal
    # 计算下降百分比
    f_descent_percent = ((f_bd_at_start - f_bd_at_avg_mntd_depth) / f_bd_at_start) * 100
  ) %>%
  ungroup() %>%
  pull(f_descent_percent) # 提取所有区域的下降百分比
# 计算下降百分比的中位数和 IQR
f_descent_median <- median(f_descent_percentages, na.rm = TRUE)
f_descent_iqr <- quantile(f_descent_percentages, probs = c(0.25, 0.75), na.rm = TRUE)
cat(paste0("区域内 Beta 多样性在系统发育距离短于平均 MNTD (", round(f_avg_mntd, 7), ") 时的下降百分比:\n",
           "  中位数: ", round(f_descent_median, 1), "% \n", # 保留一位小数
           "  IQR: ", round(f_descent_iqr[1], 1), "% - ", round(f_descent_iqr[2], 1), "% \n"))

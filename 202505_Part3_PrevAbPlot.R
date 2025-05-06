#图3c 
install.packages(c("tidyverse", "ggrepel", "RColorBrewer", "readr", "dplyr", "ggplot2", "scales", "stringr"))
install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("remotes")
remotes::install_github("kevinwolz/hisafer")
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(phyloseq)
library(readr)
library(hisafer)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
#细菌常见的分类学水平散点图，展示属的丰度和普遍性，并按纲 Class 着色，用文本标注了丰度最高的几个属，线展示了累积丰度随普遍性增加而增加的趋势
#读取数据
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
  filter(Class != "c__bacteriap25") %>%
  group_by(Genus)%>%#按属分组
  filter(!(Genus == "g__uncultured" | stringr::str_detect(Genus, "g__unclassified")))%>%
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
#计算每个属的普遍性（单个属出现的样本/总样本）
b_dat_prev <- b_dat_m %>%
  group_by(Genus, name)%>%
  filter(value > 0)%>%
  ungroup()%>%
  group_by(Genus)%>%
  summarise(Prevalence = n()/b_num_all)
b_tax_sel <- b_tax%>%
  select(-c(asv,Species))%>%
  distinct()%>%#移除重复的行
  filter(Genus %in% b_dat_abun$Genus)
#结合属、相对丰度、普遍性、分类信息
b_dat_final <- b_dat_abun%>%
  left_join(b_dat_prev)%>%
  left_join(b_tax_sel)
#计算每个纲的总丰度，并选出总丰度最高的 9 个纲。
b_classToPlot <- b_dat_final %>%
  group_by(Class)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(9)
#找出相对丰度最高的 8 个属，用于后续在图中标注
b_dat_final %>%
  dplyr::arrange(desc(Abundance)) %>%
  .[1:8, ] -> b_text_size
# 修改 b_text_size 数据框，创建新的标签列
b_text_size <- b_text_size %>%
  mutate(Label_Genus = str_replace(Genus, "g__", "")) %>% # 移除 "g__" 前缀
  mutate(Label_Genus = ifelse(
    Label_Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", # 如果是这个长名字
    "Rhizobium group", # 替换为 "Rhizobium group"
    Label_Genus # 否则保持移除 "g__" 后的名字
  ))
#将不在前 9 个纲的属的 Class 设置为 "Other"
b_dat_final <- b_dat_final%>%
  mutate(Class = if_else(Class %in% b_classToPlot$Class, Class, "Other"))
b_colors <- c(brewer.pal(9, "Set1"), "black")#调色板颜色
names(b_colors) <- c(b_classToPlot$Class, "Other")
#属的相对丰度和普遍性散点图，用颜色区分了主要分类（纲），并突出了丰度最高的几个属的名称
b_p1 <- ggplot(b_dat_final, aes(x = Prevalence, y = log10(Abundance), color = Class))+
  geom_jitter(width = 0.05, height = 0.05) +  # 增加水平方向和垂直方向的抖动幅度
  geom_point()+
  scale_y_continuous(limits = c(-7,0))+#, breaks = c(1e-5, 1e-4, 1e-3, 1e-2,1e-1, 1))+
  scale_color_manual(values = b_colors)+
  geom_text_repel(data = b_text_size, aes(label = Label_Genus), size=3)+
  theme_minimal()
b_p1
#ggsave_fitmax("PrevalenceAbundanceNOMIS_Genus.pdf", maxwidth = 10, p1)
#按精确的普遍性值分组，然后计算每个普遍性值下的所有属的 Abundance 总和
b_dat_sum <- b_dat_final %>%
  select(Abundance, Prevalence)%>%
  group_by(Prevalence)%>%
  reframe(sum = sum(Abundance))
#属的每个普遍性和每个精确普遍值的所有属总丰度线图
b_p2 <- ggplot(b_dat_sum, aes(x = Prevalence, y = log10(sum)))+
  geom_line()+
  scale_y_continuous(limits = c(-7,0))+
  theme_minimal()+
  geom_smooth(method = "gam")
b_p2
breaks <- (0:10)/10
#将 Prevalence 列的数值范围分割成 40 个等宽的区间（或箱），并将每个 Prevalence 值归入对应的区间
b_dat_cut <- b_dat_final%>%
  mutate( ints = cut(Prevalence ,breaks = 40)) %>% 
  group_by(ints) %>% 
  summarise(sum = sum(Abundance))
b_dat_cut$Prevalence <- as.numeric(b_dat_cut$ints)
#属的每个普遍性箱和每个箱内所有属总丰度线图
b_p3 <- ggplot(b_dat_cut, aes(x = as.numeric(ints), y = log10(sum)))+
  geom_line()+
  scale_y_continuous(limits = c(-7,0))+
  theme_bw()
b_p3
#在每个 Prevalence 精确值的分组内，计算所有具有该普遍性值的属的 Abundance 总和
b_temp <- b_dat_final %>%
  group_by(Prevalence)%>%
  summarise(sum = sum(Abundance))
b_temp$cumsum <- cumsum(b_temp$sum) * 100
# 定义抖动幅度 (使用你原来 geom_jitter 的值)
jitter_width <- 0.05
jitter_height <- 0.05
# 计算条件抖动后的坐标，并添加到数据框中
b_dat_final_jittered <- b_dat_final %>%
  mutate(
    # 计算 log10(Abundance) 先，方便后续计算
    Abundance_log10 = log10(Abundance),
    # 计算条件抖动后的 x 坐标: 如果 Prevalence < 1 则抖动
    Prevalence_jittered = ifelse(
      Prevalence < 1, # 条件：如果 Prevalence 小于 1
      Prevalence + runif(n(), -jitter_width / 2, jitter_width / 2), # 如果是，则加上随机抖动
      Prevalence # 如果不是 (即 Prevalence = 1)，则不加抖动
    ),
    # 计算条件抖动后的 y 坐标: 如果 Prevalence < 1 则抖动
    Abundance_jittered = ifelse(
      Prevalence < 1, # 条件：如果 Prevalence 小于 1
      Abundance_log10 + runif(n(), -jitter_height / 2, jitter_height / 2)+2, # 如果是，则加上随机抖动
      Abundance_log10+2 # 如果不是，则不加抖动
    )
  )
#去除负值
b_dat_final_jittered$Prevalence_jittered[b_dat_final_jittered$Prevalence_jittered < 0] <- 0.002
#抖动标签
b_dat_final_jittered %>%
  dplyr::arrange(desc(Abundance)) %>%
  .[1:8, ] -> b_text_size_jittered 
b_text_size_jittered <- b_text_size_jittered %>%
  mutate(Label_Genus = str_replace(Genus, "g__", "")) %>% # 移除 "g__" 前缀
  mutate(Label_Genus = ifelse(
    Label_Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", # 如果是这个长名字
    "Rhizobium group", # 替换为 "Rhizobium group"
    Label_Genus # 否则保持移除 "g__" 后的名字
  ))
# 创建一个数据框用于 Cumulative abundance 标签的位置和文本
b_cumulative_label_data <- data.frame(
  x_pos = 0.2, # 选择一个 Prevalence 值靠左的位置 (0到1之间)
  y_pos = 1,  # 选择一个 Abundance_log10 scale 上靠上的位置 (-5到2之间)
  text_label = "Cumulative abundance"
)
#属散点图和累积丰度线图，cumsum 列的每个值代表所有普遍性小于或等于当前行普遍性的属的总丰度
rp_p4_bacteria <- ggplot(b_dat_final_jittered, aes(x = Prevalence_jittered, y = Abundance_jittered))+
  #geom_jitter(aes(color = Class), width = 0.05, height = 0.05) +  # 增加水平方向和垂直方向的抖动幅度
  geom_point(aes(color = Class))+
  scale_x_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, by = 0.25), # 设置刻度位置，例如 0, 0.25, 0.5, 0.75, 1
    labels = scales::percent_format(suffix = "") ,# 将刻度标签格式化为百分比 (0% 到 100%)
    expand = expansion(mult = c(0, 0)) # 添加：X轴紧贴数据范围
  ) +
  scale_y_continuous(
    limits = c(-5, 2.4), # 保持 Y 轴范围
    breaks = seq(-5, 2, by = 1), # 保持刻度位置为 -7, -6, ..., 0
    labels = function(y) { #  自定义标签格式化函数开始
      # y 是当前刻度的值向量 (-5, -4, ..., 2)
      formatted_labels <- character(length(y)) # 创建一个空的字符向量来存储结果
      # 找到需要显示为 10 的次方的刻度（例如 <= -2）
      # 注意：这里我们根据您的描述是 -5 到 -2 次方保持原样
      power_of_10_indices <- which(y >= -5 & y <= -2)
      if (length(power_of_10_indices) > 0) {
        formatted_labels[power_of_10_indices] <- scales::math_format(10^.x)(y[power_of_10_indices])
      }
      # 找到需要显示为十进制的刻度（例如 >= -1）
      # 根据您的描述是 -1 到 2 次方用十进制
      decimal_indices <- which(y >= -1 & y <= 2)
      if (length(decimal_indices) > 0) {
        # 计算 10 的这些次方值
        decimal_values <- 10^(y[decimal_indices])
        # 将这些值格式化为字符串
        # scales::number_format() 可以用来控制小数位数，但对于 0.1, 1, 10, 100 直接转字符即可
        formatted_labels[decimal_indices] <- as.character(decimal_values)
      }
      return(formatted_labels) # 返回格式化后的标签向量
    } ,#  自定义标签格式化函数结束
    expand = expansion(mult = c(0, 0)) #  添加：Y轴紧贴数据范围
  ) +
  scale_color_manual(
    values = b_colors, # 使用之前定义的颜色向量
    name = NULL,      #  设置 name = NULL 来去掉图例标题
    labels = function(x) { stringr::str_replace(x, "c__", "") }, #  使用 labels 函数来修改标签文本
    guide = guide_legend(
      nrow = 4, 
      ncol = 3,
      override.aes = list(size = 4.5)
    )
  )+
  geom_text_repel(
    data = b_text_size_jittered,
    aes(x = Prevalence_jittered, y = Abundance_jittered, label = Label_Genus, color = Class),
    size = 4,       # 设置标签大小为 14
    fontface = "italic", # 设置标签为斜体
    max.overlaps = Inf,
    show.legend = FALSE # 确保文本标签不生成图例
  ) + # 尽最大可能找到一个无重叠的布局
  geom_line(data = b_temp, aes(x = Prevalence, y = log10(cumsum)))+
  # 添加 Cumulative abundance 标签
  geom_text(
    data = b_cumulative_label_data, # 使用专门的数据框
    aes(x = x_pos, y = y_pos, label = text_label), # 映射位置和文本
    size = 4.5, # 设置标签大小
    color = "black", # 设置标签颜色，可以根据需要更改
    hjust = 0, # 将文本的左边缘对齐到 x_pos
    vjust = 0.3 # 将文本的垂直中心对齐到 y_pos
  ) +
  # geom_smooth(data = b_dat_sum, aes(x = Prevalence, y = log10(sum)), method = "gam", se = F)+
  theme_classic()+
  #  添加主题元素来控制图例位置
  theme(
    legend.position = c(1, 0.01),      # 将图例锚点放在图面板的右下角
    legend.justification = c(1, 0) ,# 将图例框的右下角对齐到锚点
    legend.byrow = TRUE,
    axis.title.x = element_text(size = 16), # 坐标轴标题（x、y）
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14), # 坐标轴刻度标签（数字）
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14), # 图例文字
    legend.title = element_text(size = 14)
  )+
  labs(
    x = "Bacterial prevalence (%)", # X 轴标题
    y = "Bacterial relative abundance (%)" # Y 轴标题
  )
rp_p4_bacteria 
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3c_rp_bacteria.png", plot = rp_p4_bacteria, width = 8, height = 6, dpi = 600, bg = "transparent")
#优势菌门
b_phylumToPlot <- b_dat_final %>%
  group_by(Phylum)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(9)
#真菌常见的分类学水平散点图，展示属的丰度和普遍性，并按纲 Class 着色，用文本标注了丰度最高的几个属，线展示了累积丰度随普遍性增加而增加的趋势
#读取数据
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
  filter(Class != "c__fungip25") %>%
  group_by(Genus)%>%#按属分组
  filter(!(Genus == "g__uncultured" | stringr::str_detect(Genus, "g__unclassified")))%>%
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
#计算每个属的普遍性（单个属出现的样本/总样本）
f_dat_prev <- f_dat_m %>%
  group_by(Genus, name)%>%
  filter(value > 0)%>%
  ungroup()%>%
  group_by(Genus)%>%
  summarise(Prevalence = n()/f_num_all)
f_tax_sel <- f_tax%>%
  select(-c(asv,Species))%>%
  distinct()%>%#移除重复的行
  filter(Genus %in% f_dat_abun$Genus)
#结合属、相对丰度、普遍性、分类信息
f_dat_final <- f_dat_abun%>%
  left_join(f_dat_prev)%>%
  left_join(f_tax_sel)
#计算每个纲的总丰度，并选出总丰度最高的 9 个纲。
f_classToPlot <- f_dat_final %>%
  group_by(Class)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(9)
#找出相对丰度最高的 8 个属，用于后续在图中标注
f_dat_final %>%
  dplyr::arrange(desc(Abundance)) %>%
  .[1:8, ] -> f_text_size
# 修改 f_text_size 数据框，创建新的标签列
f_text_size <- f_text_size %>%
  mutate(Label_Genus = str_replace(Genus, "g__", "")) %>% # 移除 "g__" 前缀
  mutate(Label_Genus = ifelse(
    Label_Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", # 如果是这个长名字
    "Rhizobium group", # 替换为 "Rhizobium group"
    Label_Genus # 否则保持移除 "g__" 后的名字
  ))
#将不在前 9 个纲的属的 Class 设置为 "Other"
f_dat_final <- f_dat_final%>%
  mutate(Class = if_else(Class %in% f_classToPlot$Class, Class, "Other"))
f_colors <- c(brewer.pal(9, "Set1"), "black")#调色板颜色
names(f_colors) <- c(f_classToPlot$Class, "Other")
#属的相对丰度和普遍性散点图，用颜色区分了主要分类（纲），并突出了丰度最高的几个属的名称
f_p1 <- ggplot(f_dat_final, aes(x = Prevalence, y = log10(Abundance), color = Class))+
  geom_jitter(width = 0.05, height = 0.05) +  # 增加水平方向和垂直方向的抖动幅度
  geom_point()+
  scale_y_continuous(limits = c(-7,0))+#, breaks = c(1e-5, 1e-4, 1e-3, 1e-2,1e-1, 1))+
  scale_color_manual(values = f_colors)+
  geom_text_repel(data = f_text_size, aes(label = Label_Genus), size=3)+
  theme_minimal()
f_p1
#ggsave_fitmax("PrevalenceAbundanceNOMIS_Genus.pdf", maxwidth = 10, p1)
#按精确的普遍性值分组，然后计算每个普遍性值下的所有属的 Abundance 总和
f_dat_sum <- f_dat_final %>%
  select(Abundance, Prevalence)%>%
  group_by(Prevalence)%>%
  reframe(sum = sum(Abundance))
#属的每个普遍性和每个精确普遍值的所有属总丰度线图
f_p2 <- ggplot(f_dat_sum, aes(x = Prevalence, y = log10(sum)))+
  geom_line()+
  scale_y_continuous(limits = c(-7,0))+
  theme_minimal()+
  geom_smooth(method = "gam")
f_p2
breaks <- (0:10)/10
#将 Prevalence 列的数值范围分割成 40 个等宽的区间（或箱），并将每个 Prevalence 值归入对应的区间
f_dat_cut <- f_dat_final%>%
  mutate( ints = cut(Prevalence ,breaks = 40)) %>% 
  group_by(ints) %>% 
  summarise(sum = sum(Abundance))
f_dat_cut$Prevalence <- as.numeric(f_dat_cut$ints)
#属的每个普遍性箱和每个箱内所有属总丰度线图
f_p3 <- ggplot(f_dat_cut, aes(x = as.numeric(ints), y = log10(sum)))+
  geom_line()+
  scale_y_continuous(limits = c(-7,0))+
  theme_bw()
f_p3
#在每个 Prevalence 精确值的分组内，计算所有具有该普遍性值的属的 Abundance 总和
f_temp <- f_dat_final %>%
  group_by(Prevalence)%>%
  summarise(sum = sum(Abundance))
f_temp$cumsum <- cumsum(f_temp$sum) * 100
# 定义抖动幅度 (使用你原来 geom_jitter 的值)
jitter_width <- 0.05
jitter_height <- 0.05
# 计算条件抖动后的坐标，并添加到数据框中
f_dat_final_jittered <- f_dat_final %>%
  mutate(
    # 计算 log10(Abundance) 先，方便后续计算
    Abundance_log10 = log10(Abundance),
    # 计算条件抖动后的 x 坐标: 如果 Prevalence < 1 则抖动
    Prevalence_jittered = ifelse(
      Prevalence < 1, # 条件：如果 Prevalence 小于 1
      Prevalence + runif(n(), -jitter_width / 2, jitter_width / 2), # 如果是，则加上随机抖动
      Prevalence # 如果不是 (即 Prevalence = 1)，则不加抖动
    ),
    # 计算条件抖动后的 y 坐标: 如果 Prevalence < 1 则抖动
    Abundance_jittered = ifelse(
      Prevalence < 1, # 条件：如果 Prevalence 小于 1
      Abundance_log10 + runif(n(), -jitter_height / 2, jitter_height / 2)+2, # 如果是，则加上随机抖动
      Abundance_log10+2 # 如果不是，则不加抖动
    )
  )
#去除负值
f_dat_final_jittered$Prevalence_jittered[f_dat_final_jittered$Prevalence_jittered < 0] <- 0.002
#抖动标签
f_dat_final_jittered %>%
  dplyr::arrange(desc(Abundance)) %>%
  .[1:8, ] -> f_text_size_jittered 
f_text_size_jittered <- f_text_size_jittered %>%
  mutate(Label_Genus = str_replace(Genus, "g__", "")) %>% # 移除 "g__" 前缀
  mutate(Label_Genus = ifelse(
    Label_Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", # 如果是这个长名字
    "Rhizobium group", # 替换为 "Rhizobium group"
    Label_Genus # 否则保持移除 "g__" 后的名字
  ))
# 创建一个数据框用于 Cumulative abundance 标签的位置和文本
f_cumulative_label_data <- data.frame(
  x_pos = 0.2, # 选择一个 Prevalence 值靠左的位置 (0到1之间)
  y_pos = 1.5,  # 选择一个 Abundance_log10 scale 上靠上的位置 (-5到2之间)
  text_label = "Cumulative abundance"
)
#属散点图和累积丰度线图，cumsum 列的每个值代表所有普遍性小于或等于当前行普遍性的属的总丰度
rp_p4_fungi <- ggplot(f_dat_final_jittered, aes(x = Prevalence_jittered, y = Abundance_jittered))+
  #geom_jitter(aes(color = Class), width = 0.05, height = 0.05) +  # 增加水平方向和垂直方向的抖动幅度
  geom_point(aes(color = Class))+
  scale_x_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, by = 0.25), # 设置刻度位置，例如 0, 0.25, 0.5, 0.75, 1
    labels = scales::percent_format(suffix = "") ,# 将刻度标签格式化为百分比 (0% 到 100%)
    expand = expansion(mult = c(0, 0)) # 添加：X轴紧贴数据范围
  ) +
  scale_y_continuous(
    limits = c(-5, 2.4), # 保持 Y 轴范围
    breaks = seq(-5, 2, by = 1), # 保持刻度位置为 -7, -6, ..., 0
    labels = function(y) { #  自定义标签格式化函数开始
      # y 是当前刻度的值向量 (-5, -4, ..., 2)
      formatted_labels <- character(length(y)) # 创建一个空的字符向量来存储结果
      # 找到需要显示为 10 的次方的刻度（例如 <= -2）
      # 注意：这里我们根据您的描述是 -5 到 -2 次方保持原样
      power_of_10_indices <- which(y >= -5 & y <= -2)
      if (length(power_of_10_indices) > 0) {
        formatted_labels[power_of_10_indices] <- scales::math_format(10^.x)(y[power_of_10_indices])
      }
      # 找到需要显示为十进制的刻度（例如 >= -1）
      # 根据您的描述是 -1 到 2 次方用十进制
      decimal_indices <- which(y >= -1 & y <= 2)
      if (length(decimal_indices) > 0) {
        # 计算 10 的这些次方值
        decimal_values <- 10^(y[decimal_indices])
        # 将这些值格式化为字符串
        # scales::number_format() 可以用来控制小数位数，但对于 0.1, 1, 10, 100 直接转字符即可
        formatted_labels[decimal_indices] <- as.character(decimal_values)
      }
      return(formatted_labels) # 返回格式化后的标签向量
    } ,#  自定义标签格式化函数结束
    expand = expansion(mult = c(0, 0)) #  添加：Y轴紧贴数据范围
  ) +
  scale_color_manual(
    values = f_colors, # 使用之前定义的颜色向量
    name = NULL,      #  设置 name = NULL 来去掉图例标题
    labels = function(x) { stringr::str_replace(x, "c__", "") }, #  使用 labels 函数来修改标签文本
    guide = guide_legend(
      nrow = 4, 
      ncol = 3,
      override.aes = list(size = 4.5)
    )
  )+
  geom_text_repel(
    data = f_text_size_jittered,
    aes(x = Prevalence_jittered, y = Abundance_jittered, label = Label_Genus, color = Class),
    size = 4,       # 设置标签大小为 14
    fontface = "italic", # 设置标签为斜体
    max.overlaps = Inf,
    show.legend = FALSE # 确保文本标签不生成图例
  ) + # 尽最大可能找到一个无重叠的布局
  geom_line(data = f_temp, aes(x = Prevalence, y = log10(cumsum)))+
  # 添加 Cumulative abundance 标签
  geom_text(
    data = f_cumulative_label_data, # 使用专门的数据框
    aes(x = x_pos, y = y_pos, label = text_label), # 映射位置和文本
    size = 4.5, # 设置标签大小
    color = "black", # 设置标签颜色，可以根据需要更改
    hjust = 0, # 将文本的左边缘对齐到 x_pos
    vjust = 0.3 # 将文本的垂直中心对齐到 y_pos
  ) +
  # geom_smooth(data = f_dat_sum, aes(x = Prevalence, y = log10(sum)), method = "gam", se = F)+
  theme_classic()+
  #  添加主题元素来控制图例位置
  theme(
    legend.position = c(1, 0.01),      # 将图例锚点放在图面板的右下角
    legend.justification = c(1, 0) ,# 将图例框的右下角对齐到锚点
    legend.byrow = TRUE,
    axis.title.x = element_text(size = 16), # 坐标轴标题（x、y）
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14), # 坐标轴刻度标签（数字）
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14), # 图例文字
    legend.title = element_text(size = 14)
  )+
  labs(
    x = "Fungal prevalence (%)", # X 轴标题
    y = "Fungal relative abundance (%)" # Y 轴标题
  )
rp_p4_fungi 
# ggsave("D:/study/master/Main_Figure_tables/Figure_3/3c_rp_fungi.png", plot = rp_p4_fungi, width = 8, height = 6, dpi = 600, bg = "transparent")
#优势菌门
f_phylumToPlot <- f_dat_final %>%
  group_by(Phylum)%>%
  summarise(sum= sum(Abundance))%>%
  arrange(desc(sum))%>%
  top_n(9)

#图1a
options(timeout = 600)  # 设置超时时间为600秒
install.packages("raster")
install.packages("httr")
install.packages("geodata")
install.packages("sf")
install.packages("fmsb")
install.packages("scales")
install.packages("ggplot2")
install.packages("terra")
install.packages("png")
install.packages("grid")
install.packages("dplyr")
library(raster) 
library(httr)
library(geodata)
library(sf)
library(fmsb)
library(scales)
library(ggplot2)
library(terra)
library(png)
library(grid)
library(dplyr)
# ---------------- 1. 获取坐标和植物数据 ----------------
#输入坐标
pops <- data.frame(
  code = c("JR", "JJ", "TZ", "PA"),
  longitude = c(119.08889, 120.19389, 121.0225, 120.35167),  # 经度
  latitude = c(32.13889, 32.00111, 31.99611, 28.92445)  # 纬度
)
# 确保 pops 数据框中的经纬度是数值型
pops$longitude <- as.numeric(pops$longitude)
pops$latitude <- as.numeric(pops$latitude)
# 使用 sf 包转换为空间数据框
pops_sf <- st_as_sf(pops, coords = c("longitude", "latitude"), crs = 4326) 
pops_sf
# 检查 pops_sf 数据的 CRS
crs(pops_sf)
# 查看 pops_sf 的范围
ext(pops_sf)
#获取植物数据
rb_metadata <- read.delim(file = "D:/study/master/rb_metadata.tsv", sep = "\t", header = TRUE, check.names = FALSE) 
# 计算每个 Group 的数值型变量的均值
plant_data_mean_id <- rb_metadata %>%
  group_by(Origin) %>%  # 按 Origin 分组
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "keep") %>%
  mutate(Origin = factor(Origin, levels = c("JR", "TZ", "JJ", "PA"))) %>%  # 自定义排序
  arrange(Origin)  # 按自定义顺序排列
print(plant_data_mean_id)
# ---------------- 2. 绘制雷达图 ----------------
# 去除 ID 列，保留植物数据
plant_data_mean <- plant_data_mean_id %>%
  ungroup() %>%  # 取消分组
  select(-Origin)  # 删除 Origin 列
# 最小-最大归一化函数
normalize_data <- function(data) {
  return((data - min(data)) / (max(data) - min(data)) * 100)
}
# 对植物数据进行归一化
plant_data_mean_normalized <- as.data.frame(lapply(plant_data_mean, normalize_data))
# 定义雷达图绘制函数
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = NULL,vlcex = 0.7,
                                        caxislabels = NULL, title = NULL,...) {
  radarchart(
    data, axistype = 0,  # 设置为0来完全去掉轴
    # 自定义多边形
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # 自定义网格
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # 自定义轴标签
    axislabcol = "grey", 
    # 变量标签设为NULL，去掉变量名
    vlcex = 0.0007, vlabels = NULL,
    caxislabels = caxislabels, title = title, ...
  )
}
#外观性状雷达图
plant_appearance <- plant_data_mean[4:9]
plant_appearance_normalized <- plant_data_mean_normalized[4:9]
#保存图片
#png("D:/study/master/Main_Figure_tables/Figure_1/1a_radar_chart_appearance.png", width = 2000, height = 2000, res = 300, bg = "transparent")
# 设置图形布局
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(2, 2))  # 生成多个雷达图（根据采样点数量调整）
# 颜色定义
colors <- c("JR" = "#DC9445FF", 
            "JJ" = "#E5614CFF", 
            "TZ" = "#bee183", 
            "PA" = "#ADD8E6")
# 循环绘制每个采样点的雷达图
for(i in 1:nrow(plant_appearance_normalized)) {
  # 创建雷达图数据，仍然使用归一化的数据进行绘制
  radar_data <- rbind(rep(100, ncol(plant_appearance_normalized)),  # 最大值
                      rep(0, ncol(plant_appearance_normalized)),    # 最小值
                      plant_appearance_normalized[i, ])  # 当前采样点的归一化数据
  # 获取当前采样点的ID
  sample_id <- plant_data_mean_id$Origin[i]
  color <- colors[sample_id]
  # 绘制雷达图
  colnames(radar_data)<-rep("",ncol(radar_data))# 去除变量名标签
  create_beautiful_radarchart(
    data = radar_data,
    vlabels = NULL,  
    caxislabels = NULL,  # 去除轴标签
    title = NULL,        # 去掉标题
    color = color
  )
  # 添加每个变量的标签（显示原始数据值）
  for(j in 1:ncol(plant_appearance_normalized)) {
    # 获取每个点的位置
    angle <- (2 * pi * (j + 0.5)) / ncol(plant_appearance_normalized)#旋转角度
    # 计算标签位置
    x_pos <- cos(angle) * (as.numeric(plant_appearance_normalized[i, j])+50) /130   # 设置标签距离中心的距离
    y_pos <- sin(angle) * (as.numeric(plant_appearance_normalized[i, j])+50) /130   # 设置标签距离中心的距离
    # 在雷达图上添加原始数据值标签
    text(x_pos, y_pos, labels = sprintf("%.1f", plant_appearance[i, j]), col = "black", cex = 1)
  }
}
# 关闭图形设备，保存文件
#dev.off()
#植物化学成分雷达图
plant_phytochemicals <- plant_data_mean[10:16]
plant_phytochemicals_normalized <- plant_data_mean_normalized[10:16]
#保存图片
#png("D:/study/master/Main_Figure_tables/Figure_1/1a_radar_chart_phytochemicals.png", width = 2000, height = 2000, res = 300, bg = "transparent")
# 设置图形布局
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(2, 2))  # 生成多个雷达图（根据采样点数量调整）
# 颜色定义
colors <- c("JR" = "#DC9445FF", 
            "JJ" = "#E5614CFF", 
            "TZ" = "#bee183", 
            "PA" = "#ADD8E6")
# 循环绘制每个采样点的雷达图
for(i in 1:nrow(plant_phytochemicals_normalized)) {
  # 创建雷达图数据，仍然使用归一化的数据进行绘制
  radar_data <- rbind(rep(100, ncol(plant_phytochemicals_normalized)),  # 最大值
                      rep(0, ncol(plant_phytochemicals_normalized)),    # 最小值
                      plant_phytochemicals_normalized[i, ])  # 当前采样点的归一化数据
  # 获取当前采样点的ID
  sample_id <- plant_data_mean_id$Origin[i]
  color <- colors[sample_id]
  # 绘制雷达图
  colnames(radar_data)<-rep("",ncol(radar_data))# 去除变量名标签
  create_beautiful_radarchart(
    data = radar_data,
    vlabels = NULL,  
    caxislabels = NULL,  # 去除轴标签
    title = NULL,        # 去掉标题
    color = color
  )
  # 添加每个变量的标签（显示原始数据值）
  for(j in 1:ncol(plant_phytochemicals_normalized)) {
    # 获取每个点的位置
    angle <- (2 * pi * (j + 0.75)) / ncol(plant_phytochemicals_normalized) #旋转角度
    # 计算标签位置
    x_pos <- cos(angle) * (as.numeric(plant_phytochemicals_normalized[i, j])+50) /130   # 设置标签距离中心的距离
    y_pos <- sin(angle) * (as.numeric(plant_phytochemicals_normalized[i, j])+50) /130   # 设置标签距离中心的距离
    # 在雷达图上添加原始数据值标签
    text(x_pos, y_pos, labels = sprintf("%.1f", plant_phytochemicals[i, j]), col = "black", cex = 1)
  }
}
# 关闭图形设备，保存文件
#dev.off()
par(mfrow = c(1, 1))  # 恢复为单个图形
#定义只有变量名标签的空白雷达图函数
create_empty_radarchart <- function(data, color = "#00AFBB", #vlabels = NULL,
                                    vlcex = 0.7,
                                    caxislabels = NULL, title = NULL,...) {
  radarchart(
    data, axistype = 0,  # 设置为0来完全去掉轴
    # 自定义多边形
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 0, plty = 1,
    # 自定义网格
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # 自定义轴标签
    axislabcol = "grey", 
    # 变量标签设为NULL，默认变量名
    vlcex = 1.4, #vlabels = NULL,
    caxislabels = caxislabels, title = title, ...
  )
}
#绘制空白雷达图，忽略错误警告
#忽视错误
#try(create_empty_radarchart(data = climate_values_data_normalized), silent = TRUE) 
#png("D:/study/master/Main_Figure_tables/Figure_1/1a_empty_radar_chart_appearance.png", width = 2500, height = 2500, res = 300, bg = "transparent")
#try(create_empty_radarchart(data = plant_appearance_normalized, vlabels = c("Plant height (cm)","Leaf length\n(cm)","Leaf\nwidth\n(cm)",expression("Stem diameter (cm"^2*")"),"Bulb\nraw\nweight (g)","Bulb dry\nweight (g)")) , silent = TRUE)
#dev.off()
#png("D:/study/master/Main_Figure_tables/Figure_1/1a_empty_radar_chart_phytochemicals.png", width = 2500, height = 1500, res = 300, bg = "transparent")
#try(create_empty_radarchart(data = plant_phytochemicals_normalized, vlabels = c("Soluble protein content (mg/g)","Chlorophyll a\ncontent (mg/g)","Chlorophyll\nb content\n(mg/g)","Total chlorophyll\ncontent (mg/g)","Malondialdehyde\n(nmol/g)","Peimine\ncontent\n(mg/g)","Peiminine\ncontent (mg/g)")) , silent = TRUE)
#dev.off()
# ---------------- 3. 绘制中国（江苏省和浙江省）地图 ----------------
# 下载并加载 GeoJSON 数据
geojson_url <- "https://geo.datav.aliyun.com/areas_v3/bound/geojson?code=100000_full"
china_map <- st_read(geojson_url)
# 可视化地图并添加坐标点（统一黑色）
chinampa_plot <-ggplot(data = china_map) +
  geom_sf(fill = "gray85", color = "gray85") +  # 地图填充为灰色，边界为深灰色
  geom_sf(data = pops_sf, size = 2, color = "black") +  # 添加黑色点
  theme_void() +  # 移除经纬度和背景网格
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # 设置画布背景为白色
    plot.background = element_rect(fill = "white", color = NA)   # 设置整个画布为白色
  )
chinampa_plot
library(ggplot2)
library(sf)
# 下载并加载江苏省（320000）和浙江省（330000）的 GeoJSON 数据
jiangsu_geojson <- "https://geo.datav.aliyun.com/areas_v3/bound/geojson?code=320000_full"
zhejiang_geojson <- "https://geo.datav.aliyun.com/areas_v3/bound/geojson?code=330000_full"
# 读取地图数据
jiangsu_map <- st_read(jiangsu_geojson)
zhejiang_map <- st_read(zhejiang_geojson)
# 合并两省地图数据
china_jiangsu_zhejiang_map <- rbind(jiangsu_map, zhejiang_map)
# 可视化地图并添加采样点（黑色点）
china_jiangsu_zhejiang_map_plot <- ggplot() +
  geom_sf(data = china_jiangsu_zhejiang_map, fill = "gray85", color = "gray50") +  # 江苏浙江地图
  geom_sf(data = pops_sf, size = 2, color = "black") +  # 采样点
  theme_void() +  # 移除经纬度和网格
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)   
  )
# 显示地图
china_jiangsu_zhejiang_map_plot
#ggsave("D:/study/master/Main_Figure_tables/Figure_1/1a_china_jiangsu_zhejiang_map.png", plot = china_jiangsu_zhejiang_map_plot, width = 8, height = 8, dpi = 300, bg = "transparent")

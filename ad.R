# 加载所需的R包，如果没有安装则自动安装
if(!require(multtest)) install.packages("multtest")           # 用于多重假设检验
if(!require(Seurat)) install.packages("Seurat")               # 用于单细胞RNA测序数据分析
if(!require(dplyr)) install.packages("dplyr")                 # 用于数据处理和操作
if(!require(mindr)) install.packages("mindr")                 # 用于思维导图生成
if(!require(tidyverse)) install.packages("tidyverse")         # 数据科学工具集合
if(!require(data.table)) install.packages("data.table")       # 高效数据处理
if(!require(pheatmap)) install.packages("pheatmap")           # 热图绘制
if(!require(SingleR)) install.packages("SingleR")             # 单细胞注释工具
# 安装生物信息学相关包
if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")  # 用于功能富集分析
if(!require(org.Mm.eg.db)) BiocManager::install("org.Mm.eg.db")        # 小鼠基因注释数据库
if(!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")        # 人类基因注释数据库

# 加载其他必要的库
library(tidyr)      # 数据整理
library(harmony)    # 批次效应校正
library(ggpubr)     # 基于ggplot2的出版级图形
library(ggrepel)    # 避免标签重叠
library(Seurat)
library(celldex)
library(SingleR)
library(dplyr)


# 设置工作目录及创建结果文件夹
setwd('E:/AD/GSE147424_RAW')  # 修改为你的项目路径

# 创建result文件夹（如果不存在）
if (!dir.exists("../result")) {
  dir.create("../result")
  print("已创建result文件夹用于保存结果")
}

# 清空工作环境并执行垃圾回收
rm(list=ls())
gc()

# 读取样本文件
samples = list.files(".")  # 从当前目录读取文件
dir <- file.path('.', samples)
names(dir) <- samples

# 创建一个空列表存储Seurat对象
sc <- list()

# 循环读取每个样本文件并创建Seurat对象
for(i in 1:length(dir)){
  print(i)  # 打印当前处理的样本编号
  counts <- read.table(dir[i], sep=',', header=T, row.names=1)  # 读取计数矩阵
  sc[[i]] <- CreateSeuratObject(counts, min.cells=3)  # 创建Seurat对象，仅保留在至少3个细胞中表达的基因
}

# 为每个样本添加组别信息
sc[[1]]$Group <- "Lesional"       # 病变组
sc[[2]]$Group <- "Lesional"
sc[[3]]$Group <- "Non-lesional"   # 非病变组
sc[[4]]$Group <- "Healthy"        # 健康组
sc[[5]]$Group <- "Lesional"
sc[[6]]$Group <- "Healthy"
sc[[7]]$Group <- "Lesional"
sc[[8]]$Group <- "Healthy"
sc[[9]]$Group <- "Healthy"
sc[[10]]$Group <- "Healthy"
sc[[11]]$Group <- "Non-lesional"
sc[[12]]$Group <- "Healthy"
sc[[13]]$Group <- "Healthy"
sc[[14]]$Group <- "Non-lesional"
sc[[15]]$Group <- "Non-lesional"
sc[[16]]$Group <- "Non-lesional"
sc[[17]]$Group <- "Healthy"

# 下面是注释掉的代码，用于查找特定类型的基因
# grep('^MT',x=rownames(sc[[1]]@assays$RNA@data),value = T)  # 查找线粒体基因
# grep('^RP[SL]',x=rownames(sc[[1]]@assays$RNA@data),value = T)  # 查找核糖体基因
# grep('^HB[^(P)]',x=rownames(sc[[1]]@assays$RNA@data),value = T)  # 查找血红蛋白基因

# 合并所有样本的Seurat对象
sc2 <- merge(sc[[1]], y=c(sc[[2]], sc[[3]], sc[[4]], sc[[5]], sc[[6]], sc[[7]], sc[[8]], sc[[9]], sc[[10]], 
                          sc[[11]], sc[[12]], sc[[13]], sc[[14]], sc[[15]], sc[[16]], sc[[17]]))



# 检查并打印NA值的数量
na_count <- sum(is.na(sc2@meta.data$Group))
print(paste("Group中NA值的数量:", na_count))

sc2  # 查看合并后的对象
table(sc2@meta.data$Group)  # 查看各组样本数量

# 保存合并后的数据到result文件夹
saveRDS(sc2, '../result/step1_0822merged.rds')


s.integrated <- readRDS('../result/step1_0822merged.rds')
print("读取RDS文件后:")
print(table(is.na(s.integrated@meta.data$Group)))


# 修正Group标签大小写问题
s.integrated@meta.data$Group <- gsub("Non-lesional", "Non-Lesional", s.integrated@meta.data$Group)
print("修正标签后:")
print(table(s.integrated@meta.data$Group))

s.integrated <- PercentageFeatureSet(sc2,'^MT',col.name = 'percent_MT')
s.integrated <- subset(s.integrated, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent_MT < 25)
print("质控筛选后:")
print(table(is.na(s.integrated@meta.data$Group)))
table(s.integrated@meta.data$Group)

s.integrated <- NormalizeData(s.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
s.integrated <- FindVariableFeatures(s.integrated, selection.method = "vst", nfeatures = 2000)
s.integrated <- ScaleData(object = s.integrated,features = rownames(s.integrated))
s.integrated <- RunPCA(s.integrated,  features = VariableFeatures(object = s.integrated),reduction.name = "pca")
ElbowPlot(s.integrated) 

#####harmony去批次
s.integrated <- RunHarmony(s.integrated,group.by.vars='orig.ident',reduction.save = "harmony")
s.integrated <- RunTSNE(s.integrated, reduction = "harmony", dims = 1:15)
s.integrated <- FindNeighbors(s.integrated, reduction = "pca", dims = 1:15)
s.integrated <- FindClusters(s.integrated, resolution = 0.1)

#annoate
s.se <- HumanPrimaryCellAtlasData()
datas = GetAssayData(s.integrated,slot='data')
clusters.combined = s.integrated@meta.data$seurat_clusters
pred.s <- SingleR(test = datas,ref = s.se,labels= s.se$label.main,clusters = clusters.combined,
                  assay.type.test='logcounts',assay.type.ref='logcounts')
celltype = data.frame(cluster=rownames(pred.s),celltype=pred.s$labels)
for(i in 1:45311){
  index=s.integrated@meta.data$seurat_clusters[i]
  s.integrated@meta.data$celltype[i]=celltype[index,2]
}
print("注释后:")
print(table(is.na(s.integrated@meta.data$celltype)))



# 保存注释后的数据到result文件夹
saveRDS(s.integrated, '../result/step3_0822annoated.rds')

# 读取保存的数据
s.integrated <- readRDS('../result/step3_0822annoated.rds')
print("读取RDS文件后:")
print(table(is.na(s.integrated@meta.data$Group)))

###########第一个图：细胞类型tSNE可视化###########

# 使用tSNE降维展示细胞类型分布
###########First figure
p <- DimPlot(s.integrated,group.by = 'celltype',reduction='tsne',label=TRUE,pt.size=1)+theme(text=element_text(size=10))
p

# 保存图形到result文件夹
ggsave("../result/08SingleR.png", p, width = 7, height = 5, dpi = 300)
ggsave("../result/08TSNE.svg", p, width = 6, height = 6, dpi = 300)



###########Second figure
Idents(s.integrated) <- s.integrated@meta.data$celltype
markers <- FindAllMarkers(s.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p1 <- DoHeatmap(s.integrated, features = top5$gene,label=FALSE) + NoLegend() +theme(legend.position = "right",text=element_text(size=15),legend.text=element_text(size=15),axis.text.y=element_text(size=16))+ 
  scale_fill_gradient2(high = "#ff8c69", low = '#6495ed')
#scale_fill_discrete_qualitative(palette = "cold")

p1
#p+theme(text=element_text(size=16,  family="serif"))
# 保存热图到result文件夹
ggsave("../result/08cellmarker.png", p1, width = 9, height = 9, dpi = 300)
ggsave("../result/08cellmarker.bmp", p1, width = 9, height = 9, dpi = 300)
ggsave("../result/08cellmarker.svg", p1, width = 9, height = 9)




############Third figure cell type distribution
metadata <- FetchData(s.integrated, vars = c("orig.ident", 'Group', "celltype"))
print("FetchData后:")
print(table(is.na(metadata$Group)))
print(table(is.na(metadata$orig.ident)))
print(table(is.na(metadata$celltype)))

# 检查orig.ident的唯一值
print("orig.ident的唯一值:")
print(unique(metadata$orig.ident))

# 确保分组名称一致性
unique_groups <- unique(metadata$Group)
print("唯一的Group值:")
print(unique_groups)

# 正确设置因子水平
metadata$Group <- factor(metadata$Group, levels = unique_groups)

# 设置样本顺序，与图1保持一致
sample_order <- c("S4", "S6", "S9", "S10", "S12", "S13", "S17",
                  "S3", "S11", "S14", "S15", "S16",
                  "S1", "S2", "S5", "S7")

# 检查是否有样本不在sample_order中
all_samples <- unique(metadata$orig.ident)
missing_samples <- setdiff(all_samples, sample_order)
print("不在sample_order中的样本:")
print(missing_samples)

# 过滤掉不在sample_order中的样本
metadata_filtered <- metadata %>% filter(orig.ident %in% sample_order)
print("过滤后的样本数量:")
print(nrow(metadata_filtered))
print("过滤前的样本数量:")
print(nrow(metadata))

# 使用过滤后的数据计算各细胞类型的比例
counts <- metadata_filtered %>%
  group_by(orig.ident, Group, celltype) %>%
  summarize(Count = n(), .groups = 'drop') %>%
  ungroup()

# 设置样本顺序
counts$orig.ident <- factor(counts$orig.ident, levels = sample_order)

# 设置细胞类型顺序
cell_order <- c("DC", "Endothelial_cells", "Epithelial_cells", "Fibroblasts", 
                "Neurons", "Tissue_stem_cells", "T_cells", "Keratinocytes")
counts$celltype <- factor(counts$celltype, levels = cell_order)

# 计算每个样本内的细胞类型比例
counts <- counts %>%
  group_by(orig.ident) %>%
  mutate(prop = Count / sum(Count))

# 为图1中使用的颜色方案
colors <- c(
  "DC" = "#b4eeb4",            # 浅绿色
  "Endothelial_cells" = "#FFEBCD",  # 浅橙黄色
  "Epithelial_cells" = "#FFEBCD",   # 浅橙黄色
  "Fibroblasts" = "#ff8c69",       # 橙红色
  "Neurons" = "#C1CDCD",           # 灰色
  "Tissue_stem_cells" = "#FFFFDB", # 浅黄色
  "T_cells" = "#F0F8FF",          # 浅蓝色
  "Keratinocytes" = "#6495ed"      # 蓝色
)

# 修正分组标签以匹配数据中的实际值
correct_group_names <- unique_groups  # 使用实际数据中的Group名称
group_positions <- c(3.5, 10.5, 14)   # 保持位置不变
group_colors <- c("green3", "turquoise3", "gold")  # 保持颜色不变

# 创建图形
p <- ggplot(counts, aes(x = orig.ident, y = prop, fill = celltype)) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values = colors) +
  labs(x = "", y = "Celltype(Proportion)", fill = "Celltype") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = c(
      rep("green3", 7),  # Healthy 组样本颜色
      rep("turquoise3", 5), # Non-lesional 组样本颜色 (注意小写l)
      rep("gold", 4)     # Lesional 组样本颜色
    )),
    panel.grid.major.x = element_blank()
  )

# 添加分组标签，使用数据中的实际分组名称
group_labels <- data.frame(
  x = group_positions,  # 标签位置
  y = -0.1,             # y轴位置
  Group = correct_group_names,  # 使用实际数据中的分组名称
  color = group_colors
)

p <- p + geom_text(
  data = group_labels,
  aes(x = x, y = y, label = Group, color = Group),
  inherit.aes = FALSE,
  size = 5
) +
  scale_color_manual(values = setNames(group_colors, correct_group_names)) +
  guides(color = "none") +  # 隐藏颜色图例
  ylim(-0.15, 1)  # 扩展y轴范围以容纳标签

print(p)

# 保存图表
ggsave("../result/celltypedistribution.png", p2, width = 10, height = 6, dpi = 300)
ggsave("../result/08celltypedistribution去s8.svg", p2, width = 12, height = 6, dpi = 300)


# 检查是否存在NA值
table(is.na(metadata$Group))
table(is.na(metadata$orig.ident))

###########第四个图：角质形成细胞（KC）比例分析###########
dak <- filter(counts, celltype == "Keratinocytes")
dak$orig.ident <- factor(dak$orig.ident, levels = c("S4", "S6", "S9", "S10", "S12", "S13", "S17",
                                                    "S3", "S11", "S14", "S15","S16",
                                                    "S1", "S2", "S5", "S7"))
#factor(dak$Group)
dak$Group <- factor(dak$Group,levels = c("Healthy","Non-lesional","Lesional"))

# 绘制箱体图和样本点 画kc fraction

p3 <- ggplot(dak, aes(x = Group, y = prop)) +
  geom_boxplot(width = 0.3, fill = "white", color = "black") +  # 调节箱体宽度和颜色
  geom_jitter(width = 0.2,size=3,height = 0,aes(color = Group), alpha = 0.3) +  # 添加样本点并设置颜色设orig.ident就样本一个颜色
  labs(x = "Group", y="Percentage",title  = "Keratinocytes fraction") +
  geom_hline(yintercept = 0.1, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 0.15, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 0.05, lty = 4, col = "black", lwd = 0.8)+
  #geom_hline(yintercept = 4, lty = 4, col = "black", lwd = 0.8)+
  #geom_hline(yintercept = 5, lty = 4, col = "black", lwd = 0.8)+
  theme(        
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.45),
    text=element_text(size=16),legend.position="none")+
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA)

#my_comparisons <- list(c("healthy","Non-lesional"), c("healthy", "Lesional"),
#                       c("Non-lesional", "Lesional"))
p3 #+ stat_compare_means(comparisons = my_comparisons,method='t.test')

# 保存图形到result文件夹
ggsave("../result/08KCfraction.png", p3, width = 5, height = 5, dpi = 300)
ggsave("../result/08KCfraction.bmp", p3, width = 5, height = 5, dpi = 300)
ggsave("../result/08KCfraction去s8.svg", p3, width = 5, height = 5)


###########角质形成细胞亚群分析###########

KC_cluster <- subset(s.integrated, subset = celltype == "Keratinocytes")


gene_expression <- FetchData(KC_cluster, vars = c("CDC42"), slot = "data")
Group <- KC_cluster@meta.data$Group
data <- data.frame(gene_expression, Group)
data <- data[data$CDC42 != 0, ]
factor(data$Group)
data$Group <- factor(data$Group,levels = c("Healthy","Non-lesional","Lesional"))

#加p值
my_comparisons <- list(c("Healthy","Non-lesional"), c("Healthy", "Lesional"),
                       c("Non-lesional", "Lesional"))

p <- ggplot(data, aes(x = Group, y = CDC42,fill = Group)) +
  geom_boxplot(fill = "white", color = "black",width=0.2) +
  labs(x = "Group", y = "Gene Expression (CDC42) in Keratinocytes") +
  #geom_hline(yintercept = 1, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2.5, lty = 4, col = "black", lwd = 0.8)+
  #geom_hline(yintercept = 4, lty = 4, col = "black", lwd = 0.8)+
  #geom_hline(yintercept = 5, lty = 4, col = "black", lwd = 0.8)+
  scale_y_continuous(breaks = seq(1, 4, by = 0.5)) +
  stat_compare_means(comparisons = my_comparisons )+
  theme(        
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
    text=element_text(size=16))+
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA) +
  # 添加标题
  labs(title = "CDC42")
p

# 保存CDC42表达图到result文件夹
ggsave("../result/08CDC42_expression.png", p, width = 5, height = 5, dpi = 300)

##########角质形成细胞亚群聚类分析##########

# 数据预处理
KC_cluster <- ScaleData(KC_cluster)  # 数据标准化
KC_cluster <- FindVariableFeatures(KC_cluster, selection.method = "vst", nfeatures = 2000)  # 寻找高变异基因
KC_cluster <- RunPCA(object = KC_cluster, features = VariableFeatures(KC_cluster))  # 主成分分析
p_elbow <- ElbowPlot(KC_cluster)  # 查看主成分贡献率
ggsave("../result/08KC_ElbowPlot.png", p_elbow, width = 7, height = 5, dpi = 300)

# JackStraw分析确定有意义的主成分
KC_cluster <- JackStraw(KC_cluster, num.replicate = 100) 
KC_cluster <- ScoreJackStraw(KC_cluster, dims = 1:20) 
p_jackstraw <- JackStrawPlot(KC_cluster, dims = 1:20)
ggsave("../result/08KC_JackStrawPlot.png", p_jackstraw, width = 10, height = 7, dpi = 300)

# 构建KNN网络
KC_cluster <- FindNeighbors(KC_cluster, dims = 1:11)

# 使用合适的分辨率进行聚类（3个亚类）
res.used <- 0.045  # 分辨率参数
KC_cluster <- FindClusters(object = KC_cluster, verbose = T, resolution = res.used)

# tSNE降维可视化
KC_cluster <- RunTSNE(object = KC_cluster, dims = 1:11, do.fast = TRUE)

# 绘制tSNE图
kcd <- DimPlot(KC_cluster, reduction='tsne', label=TRUE) + 
  theme(text=element_text(size=10), legend.position = "none") + 
  scale_color_manual(values = c("lightblue", "grey", "orange"))
kcd  # 显示图形
ggsave("../result/08kcd.svg", kcd, width = 5, height = 5, dpi = 300)
ggsave("../result/08kcd.png", kcd, width = 5, height = 5, dpi = 300)

# 寻找亚群标志基因
KC_cluster.markers <- FindAllMarkers(KC_cluster, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

# 选择每个亚群的前20个差异基因
top20 <- KC_cluster.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# 绘制热图
p_heatmap1 <- DoHeatmap(KC_cluster, features = top20$gene) + 
  NoLegend() + 
  theme(legend.position = "right", legend.text=element_text(size=15), axis.text.y=element_text(size=6)) + 
  scale_fill_gradient2(high = "#ff8c69", low = '#6495ed')
ggsave("../result/08KC_heatmap_before_rename.png", p_heatmap1, width = 12, height = 9, dpi = 300)

# 给亚群命名
new.cluster.ids <- c("Basal", "Suprabasal", "Proliferating")
names(new.cluster.ids) <- levels(KC_cluster)
KC_cluster <- RenameIdents(KC_cluster, new.cluster.ids)
KC_cluster@meta.data$Subtype <- Idents(KC_cluster)

# 再次可视化
p_dimplot_renamed <- DimPlot(KC_cluster, reduction = "tsne", label=T)
ggsave("../result/08KC_clusters_renamed.png", p_dimplot_renamed, width = 6, height = 5, dpi = 300)

# 绘制新的热图
p_heatmap2 <- DoHeatmap(KC_cluster, features = top20$gene, label=FALSE) + 
  NoLegend() + 
  theme(legend.position = "right", text=element_text(size=15), 
        legend.text=element_text(size=15), axis.text.y=element_text(size=16)) + 
  scale_fill_gradient2(high = "#ff8c69", low = '#6495ed')
ggsave("../result/08KC_heatmap_after_rename.png", p_heatmap2, width = 12, height = 9, dpi = 300)

# 提取基底细胞亚群
Basal <- subset(KC_cluster, subset = seurat_clusters == "0")
Idents(Basal) <- Basal@meta.data$Group

# 移除特定样本（S8）
sample_to_remove <- "S8"
Basal <- subset(Basal, subset = orig.ident != sample_to_remove)

# 保存数据到result文件夹
saveRDS(Basal, '../result/step5_0822basalannoated.rds')
saveRDS(s.integrated, '../result/step4_0822kcsubannoated.rds')

######将Seurat对象转换为HDF5格式（用于Python分析）######
library(SeuratDisk)
Basal <- DietSeurat(Basal, scale.data = FALSE)  # 减小对象大小
SaveH5Seurat(Basal, filename = "../result/py0825basal", assays = "RNA")  # 保存为H5Seurat格式
Convert("../result/py0825basal.h5seurat", dest = "h5ad", overwrite = TRUE)  # 转换为H5AD格式（兼容Python的scanpy）

###################################
# 可视化特定基因表达
p_feature <- FeaturePlot(KC_cluster, features = c("S100A9"), reduction = 'tsne')
ggsave("../result/08KC_S100A9_expression.png", p_feature, width = 6, height = 5, dpi = 300)

# 保存角质形成细胞亚群数据到result文件夹
saveRDS(KC_cluster, '../result/step6_0826kcsub.rds')

# 导出为loom格式（兼容Velocity分析）
library(loomR)
library(SeuratDisk)
kc.loom <- as.loom(x=KC_cluster, filename='../result/08kc.loom', verbose=FALSE)
kc.loom$close_all()

# 加载保存的数据
KC_cluster <- readRDS('../result/step6_0826kcsub.rds')

# UMAP降维可视化
KC_cluster <- RunUMAP(object = KC_cluster, dims = 1:11, do.fast = TRUE)
p_umap <- DimPlot(KC_cluster, reduction = "umap", label=T)
ggsave("../result/08KC_umap.png", p_umap, width = 6, height = 5, dpi = 300)

# 可视化AQP3基因在UMAP上的表达
p_aqp3 <- FeaturePlot(KC_cluster, features = c("AQP3"), reduction='umap')
ggsave("../result/08KC_AQP3_umap.png", p_aqp3, width = 6, height = 5, dpi = 300)

# 移除特定样本（S8）
sample_to_remove <- "S8"
KC_cluster <- subset(KC_cluster, subset = orig.ident != sample_to_remove)
saveRDS(KC_cluster, '../result/step7_0826kcsub_no_S8.rds')

# 打印完成信息
print("分析完成！所有结果已保存至'result'文件夹")


##堆叠小提琴图
Idents(KC_cluster)
gene_list <- c('KRT1','KRT5','KRT14','KRT15','KRT17',
               'KRT19','AQP3','MGST1','S100A8','S100A9')
gene_exp <- FetchData(KC_cluster, vars = gene_list, slot = "data")
#Group <- Basal@meta.data$Group
data <- data.frame(gene_exp, KC_cluster@active.ident)
#data <- data[data$KRT1 != 0, ]

vp <- VlnPlot(KC_cluster, features = gene_list,
              stack=T,pt.size=0
              
)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(colour = 'black',size = 10,angle = 0),
        legend.position = 'none')
vp
ggsave("../result/08vp.svg", vp , width = 7.5, height = 3.5, dpi = 300)

###################################
# 从step3的数据开始读取
s.integrated <- readRDS('../result/step3_0822annoated.rds')

# 先移除S8样本，然后再提取子集
sample_to_remove <- "S8"
s.integrated <- subset(s.integrated, subset = orig.ident != sample_to_remove)

# 提取角质形成细胞子集
KC_cluster <- subset(s.integrated, subset = celltype == "Keratinocytes")

# 查看当前可用的群集
print(table(KC_cluster$seurat_clusters))

# 基于表格结果，选择群集2作为基底角质形成细胞
Basal <- subset(KC_cluster, subset = seurat_clusters == "2")

# 检查提取是否成功
dim(Basal)

#####################################
# CDC42表达分析部分
#####################################

# 获取CDC42基因表达数据
gene_expression <- FetchData(Basal, vars = c("CDC42"), slot = "data")
Group <- Basal@meta.data$Group
data <- data.frame(gene_expression, Group)
data <- data[data$CDC42 != 0, ]

# 设置分组顺序
data$Group <- factor(data$Group, levels = c("Healthy", "Non-lesional", "Lesional"))

# 设置比较组
my_comparisons <- list(c("Healthy","Non-lesional"), c("Non-lesional", "Lesional"), c("Healthy", "Lesional"))

# 创建绘图
pg <- ggplot(data, aes(x = Group, y = CDC42, fill = Group)) +
  geom_boxplot(fill = "white", color = "black", width = 0.4) +
  labs(x = "Group", y = "Gene Expression (CDC42) in Keratinocytes") +
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2.5, lty = 4, col = "black", lwd = 0.8) +
  scale_y_continuous(breaks = seq(1, 4, by = 0.5)) +
  stat_compare_means(comparisons = my_comparisons) +
  theme(        
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
    text = element_text(size = 16)
  ) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA) +
  labs(title = "CDC42")

# 显示绘图
print(pg)

# 保存绘图
ggsave("../result/08CDC42去s8.svg", pg, width = 5, height = 5.5, dpi = 300)

#####################################
# 差异表达基因分析部分
#####################################

# 检查Group变量
print(table(Basal$Group))

# 设置当前的ident为Group
Idents(Basal) <- "Group"

# 验证ident设置
print(table(Idents(Basal)))

# L&H差异表达分析
le.markers <- FindMarkers(Basal, assay = 'RNA', slot = 'data', 
                          ident.1 = "Lesional", ident.2 = "Healthy",
                          logfc.threshold = 0, min.pct = 0.25)
row_names <- rownames(le.markers)
le.markers$Gene <- row_names
write.csv(le.markers, '../result/0822去除s8kc类L&H_deg_data_fc.csv')

# L&NL差异表达分析
lnl.markers <- FindMarkers(Basal, assay = 'RNA', slot = 'data', 
                           ident.1 = "Lesional", ident.2 = "Non-lesional",
                           logfc.threshold = 0, min.pct = 0.25)
row_names <- rownames(lnl.markers)
lnl.markers$Gene <- row_names
write.csv(lnl.markers, '../result/0822去除s8kc类L&NL_deg_data.csv')

# NL&H差异表达分析
nlh.markers <- FindMarkers(Basal, assay = 'RNA', slot = 'data',
                           ident.1 = "Non-lesional", ident.2 = "Healthy",
                           logfc.threshold = 0, min.pct = 0.25)
row_names <- rownames(nlh.markers)
nlh.markers$Gene <- row_names
write.csv(nlh.markers, '../result/0822去除s8kc类NL&H_deg_data.csv')

# LNL&H差异表达分析
lnlh.markers <- FindMarkers(Basal, assay = 'RNA', slot = 'data', 
                            ident.1 = c("Lesional", "Non-lesional"), ident.2 = "Healthy",
                            logfc.threshold = 0, min.pct = 0.25)
row_names <- rownames(lnlh.markers)
lnlh.markers$Gene <- row_names
write.csv(lnlh.markers, '../result/0822去除s8kc类LNL&H_deg_data.csv')

###################################
# 火山图数据准备部分
###################################
# 提取各组的CDC42表达数据
h_exp <- data$CDC42[data$Group == "Healthy"]
nl_exp <- data$CDC42[data$Group == "Non-lesional"]
l_exp <- data$CDC42[data$Group == "Lesional"]
lnl_exp <- data$CDC42[data$Group == "Lesional"|data$Group == "Non-lesional"]

# 计算不同组间的p值和log2FC
# NL & H
p = t.test(nl_exp, h_exp)$p.value
log2fc = log2(mean(nl_exp) / mean(h_exp))

# L & H
p1 = t.test(l_exp, h_exp)$p.value
log2fc1 = log2(mean(l_exp) / mean(h_exp))

# LNL & H
p2 = t.test(lnl_exp, h_exp)$p.value
log2fc2 = log2(mean(lnl_exp) / mean(h_exp))

# L & NL
p3 = t.test(l_exp, nl_exp)$p.value
log2fc3 = log2(mean(l_exp) / mean(nl_exp))

###################################
# 批量计算所有基因的表达差异
###################################
# 初始化结果数据框
gene_list <- rownames(Basal@assays$RNA)
result_nl_h <- data.frame(Gene = character(0), PValue = numeric(0), Log2FC = numeric(0))
result_l_h <- data.frame(Gene = character(0), PValue = numeric(0), Log2FC = numeric(0))
result_lnl_h <- data.frame(Gene = character(0), PValue = numeric(0), Log2FC = numeric(0))
result_l_nl <- data.frame(Gene = character(0), PValue = numeric(0), Log2FC = numeric(0))

# 循环计算每个基因的差异表达
for (gene in gene_list) {
  gene_expression <- FetchData(Basal, vars = c(gene), slot = "data")
  Group <- Basal@meta.data$Group
  data <- data.frame(gene_expression, Group)
  data <- data[data[,1] != 0, ]
  
  h_exp <- data[,1][data$Group == "Healthy"]
  nl_exp <- data[,1][data$Group == "Non-lesional"]
  l_exp <- data[,1][data$Group == "Lesional"]
  lnl_exp <- data[,1][data$Group == "Lesional" | data$Group == "Non-lesional"]
  
  # 检查组间细胞数量是否足够进行t检验
  if (length(nl_exp) >= 2 && length(h_exp) >= 2) {
    p_nl_h <- t.test(nl_exp, h_exp)$p.value
    log2fc_nl_h <- log2(mean(nl_exp) / mean(h_exp))
  } else {
    p_nl_h <- NA
    log2fc_nl_h <- NA
  }
  
  if (length(l_exp) >= 2 && length(h_exp) >= 2) {
    p_l_h <- t.test(l_exp, h_exp)$p.value
    log2fc_l_h <- log2(mean(l_exp) / mean(h_exp))
  } else {
    p_l_h <- NA
    log2fc_l_h <- NA
  }
  
  if (length(lnl_exp) >= 2 && length(h_exp) >= 2) {
    p_lnl_h <- t.test(lnl_exp, h_exp)$p.value
    log2fc_lnl_h <- log2(mean(lnl_exp) / mean(h_exp))
  } else {
    p_lnl_h <- NA
    log2fc_lnl_h <- NA
  }
  
  if (length(l_exp) >= 2 && length(nl_exp) >= 2) {
    p_l_nl <- t.test(l_exp, nl_exp)$p.value
    log2fc_l_nl <- log2(mean(l_exp) / mean(nl_exp))
  } else {
    p_l_nl <- NA
    log2fc_l_nl <- NA
  }
  
  # 添加结果到数据框
  result_nl_h <- rbind(result_nl_h, data.frame(Gene = gene, PValue = p_nl_h, Log2FC = log2fc_nl_h))
  result_l_h <- rbind(result_l_h, data.frame(Gene = gene, PValue = p_l_h, Log2FC = log2fc_l_h))
  result_lnl_h <- rbind(result_lnl_h, data.frame(Gene = gene, PValue = p_lnl_h, Log2FC = log2fc_lnl_h))
  result_l_nl <- rbind(result_l_nl, data.frame(Gene = gene, PValue = p_l_nl, Log2FC = log2fc_l_nl))
}

# 保存t检验结果
write.csv(result_nl_h, file = '../result/08去0NL&H.csv')
write.csv(result_l_h, '../result/08去0L&H.csv')
write.csv(result_lnl_h, '../result/08去0LNL&H.csv')
write.csv(result_l_nl, '../result/08去0L&NL.csv')

###################################
# 火山图数据处理部分
###################################
# 读取保存的结果文件
result_nl_h <- read.csv('../result/08去0NL&H.csv', row.names = 1)
result_l_h <- read.csv('../result/08去0L&H.csv', row.names = 1)
result_l_nl <- read.csv('../result/08去0L&NL.csv', row.names = 1)
result_lnl_h <- read.csv('../result/08去0LNL&H.csv', row.names = 1)

# 数据清理和过滤
result_nl_h <- unique(na.omit(result_nl_h), by = "Gene")
result_l_h <- unique(na.omit(result_l_h), by = "Gene")
result_lnl_h <- unique(na.omit(result_lnl_h), by = "Gene")
result_l_nl <- unique(na.omit(result_l_nl), by = "Gene")

# 筛选p值显著的基因
re1 <- result_nl_h[result_nl_h$PValue < 0.01,]
re2 <- result_l_h[result_l_h$PValue < 0.01,]
re3 <- result_lnl_h[result_lnl_h$PValue < 0.01,]
re4 <- result_l_nl[result_l_nl$PValue < 0.01,]

# 查找不同组合中共同的差异表达基因
gene_intersect <- intersect(intersect(re2$Gene, re4$Gene), intersect(re3$Gene, re4$Gene))
gene_intersect1 <- intersect(intersect(re1$Gene, re2$Gene), re4$Gene)

# Rho家族基因列表
Rhogene <- c("RAC3", "RND2", "RHOT2", "RAC2", "RHOV", "RHOQ", "RHOF", "RHOD", "RHOC", "RHOA", 
             "RAC1", "RHOBTB2", "RHOU", "RHOJ", "RHOG", "RHOBTB1", "RHOH", "RND1", "RND3", 
             "RHOT1", "CDC42", "RHOB")

# 查找显著差异表达的Rho家族基因
fg <- intersect(gene_intersect1, Rhogene)
print(fg)

#############################20250212插入43个gene画富集分析
###################################
# 富集分析部分
###################################
# 加载必要的R包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# 基因ID转换
gene.df <- bitr(gene_intersect1, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)

# 定义富集分析结果可视化函数
erich2plot <- function(data4plot) {
  data4plot <- data4plot[order(data4plot$qvalue, decreasing = F)[1:20], ]
  data4plot$BgRatio <-
    apply(data4plot, 1, function(x) {
      as.numeric(strsplit(x[3], '/')[[1]][1])
    }) / apply(data4plot, 1, function(x) {
      as.numeric(strsplit(x[4], '/')[[1]][1])
    })
  
  p <- ggplot(data4plot, aes(BgRatio, Description))
  p <- p + geom_point()
  
  pbubble <- p + geom_point(aes(size = Count, color = -1 * log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low = "#90EE90", high = "red") + 
    labs(color = expression(-log[10](qvalue)), size = "observed.gene.count") +
    theme_bw() +
    theme(text = element_text(size = 16))
  
  return(pr)
}

# 进行KEGG通路富集分析
ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism = 'hsa',
                    pvalueCutoff = 0.05, pAdjustMethod = 'BH', qvalueCutoff = 0.2,
                    minGSSize = 10, maxGSSize = 500, use_internal_data = F)
ekegg <- setReadable(ekegg, 'org.Hs.eg.db', 'ENTREZID')

# 绘制并保存KEGG富集分析结果
pp <- erich2plot(ekegg@result) + 
  theme(axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 16))
print(pp)
ggsave("../result/KEGG.png", erich2plot(ekegg@result), width = 7, height = 5, dpi = 300)
ggsave("../result/KEGG0612.svg", pp, width = 12, height = 7, dpi = 300)
ggsave("../result/KEGG.bmp", erich2plot(ekegg@result), width = 7, height = 5, dpi = 300)

# 进行GO富集分析 (分子功能-MF)
ggoMF <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, 
                  minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
MF_plot <- erich2plot(ggoMF@result)
print(MF_plot)
ggsave("../result/GO_MF.svg", MF_plot, width = 12, height = 7, dpi = 300)

# 进行GO富集分析 (细胞组分-CC)
ggoCC <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, 
                  minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
CC_plot <- erich2plot(ggoCC@result)
print(CC_plot)
ggsave("../result/GO_CC.svg", CC_plot, width = 12, height = 7, dpi = 300)

# 进行GO富集分析 (生物过程-BP)
ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, 
                  minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
BP_plot <- erich2plot(ggoBP@result)
print(BP_plot)
ggsave("../result/GO_BP.svg", BP_plot, width = 18, height = 10, dpi = 300)


##########0613用AD比healthy marker  gene 去富集
df <- read.csv('../result/0822去除s8kc类LNL&H_deg_data.csv',row.names = 1)

gene.df <- bitr(df$Gene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count")+
    #ggtitle("KEGG") +
    theme_bw() +
    theme(text=element_text(size=16))
  pr
}


ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism='hsa',
                    pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)

ekegg <- setReadable(ekegg,'org.Hs.eg.db','ENTREZID')
pp <- erich2plot(ekegg@result)+ theme(axis.text.y=element_text(size=20)
                                      ,axis.title.x = element_blank()
                                      ,axis.title.y = element_blank()
                                      ,text=element_text(size=16))
pp

#########################################################

##cca intersected gene
cgene <- c('S100A8','S100A9','IFI27','IFITM3','CCL27','C19orf33',
           'NFIB','CRNDE','KRT16','KRT6A','SERPINB4','S100A7','IFITM1',
           'CDC42','TPPP3','CCL2','C12orf75','CTSC','TIMP1','TYMP')
intersect(cgene,l)
l
############################Fifth figure 火山图
labellist <-  c("MTCL1", "SRGN", "CDC42", "LGALS7B"
                ,"S100A7","ETV4","SRGN","CCL27")

re1$Significant <- ifelse(re1$PValue < 0.01 & abs(re1$Log2FC) >= 0.1, ifelse(re1$Log2FC > 0.1, "Up", "Down"), "Stable")

p5 <- ggplot(
  # 数据、映射、颜色
  re1, aes(x = Log2FC, y = -log10(PValue))) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("steelblue","grey", "red")) +
  # 注释
  geom_text_repel(
    data = subset(re1,re1$Gene %in% labellist),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # 辅助线
  #geom_vline(xintercept = c(-0.1, 0.1), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  #geom_hline(yintercept = 10, lty = 4, col = "black", lwd = 0.8) +
  #geom_hline(yintercept = 20, lty = 4, col = "black", lwd = 0.8) +
  # 坐标轴
  labs(x = "Fold Change",
       y = "-log10 (PValue)") +
  # 图例
  theme(
    legend.position = "None",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    text=element_text(size=16)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA) +
  # 添加标题
  labs(title = "Non-Lesional vs Healthy")

p5
ggsave("../result/08nlh0.5_1.svg", p5, width = 5, height = 5,  dpi = 300)



#########l&H
re2$Significant <- ifelse(re2$PValue < 0.01 & abs(re2$Log2FC) >= 0.1, 
                          ifelse(re2$Log2FC > 0.1, "Up", "Down"), "Stable")

labellist2 <-  c("LGALS7B", "KRT14", "CDC42", "GAS5"
                 ,"APOE","S100A8","KRT6A","CCL27")

p6 <- ggplot(
  # 数据、映射、颜色
  re2, aes(x = Log2FC, y = -log10(PValue))) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("steelblue","grey", "red")) +
  # 注释
  geom_text_repel(
    data = subset(re2,re2$Gene %in% labellist2),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # 辅助线
  #geom_vline(xintercept = c(-0.1, 0.1), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  #geom_hline(yintercept = 10, lty = 4, col = "black", lwd = 0.8) +
  #geom_hline(yintercept = 20, lty = 4, col = "black", lwd = 0.8) +
  # 坐标轴
  labs(x = "Fold Change",
       y = "-log10 (P value)") +
  # 图例
  theme(
    legend.position = "None",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    text=element_text(size=16)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA) +
  # 添加标题
  labs(title = "Lesional vs Healthy")

p6
ggsave("../result/08L&H0.5_1.svg", p6 , width = 5, height = 5, dpi = 300)


###l&NL
labellist3 <-  c("KRT17", "SLC9A3R2", "CDC42", "LRFN3"
                 ,"S100A8","KRT6A","KRT16")

re4$Significant <- ifelse(re4$PValue < 0.01 & abs(re4$Log2FC) >= 0.1, 
                          ifelse(re4$Log2FC > 0.1, "Up", "Down"), "Stable")
p7 <- ggplot(
  # 数据、映射、颜色
  re4, aes(x = Log2FC, y = -log10(PValue))) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("steelblue","grey", "red")) +
  # 注释
  geom_text_repel(
    data = subset(re4,re4$Gene %in% labellist3),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # 辅助线
  #geom_vline(xintercept = c(-0.1, 0.1), lty = 4, col = "black", lwd = 0.8) +
  #geom_hline(yintercept = 1, lty = 4, col = "black", lwd = 0.8) +
  #geom_hline(yintercept = 5, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  
  # 坐标轴
  labs(x = "Fold Change",
       y = "-log10 (P value)") +
  # 图例
  theme(
    legend.position = "None",
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    text=element_text(size=16) ) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 20, color = "black", fill = NA) +
  scale_y_continuous(limits = c(0, 20))+
  labs(title = "Lesional vs Non-Lesional")
p7

ggsave("../result/08L&NL0.5_1.svg", p7 , width = 5, height = 5, dpi = 300)


############Nine figure vn图
library(pheatmap)
re_nlh <- re1[abs(re1$Log2FC)>0.1,]

re_lh <- re2[abs(re2$Log2FC)>0.1,]

re_lnl <- re4[abs(re4$Log2FC)>0.1,]

l <- intersect(intersect(re_nlh$Gene,re_lh$Gene),re_lnl$Gene)
#23

############ten figure pheatmap

lgene_expression <- FetchData(Basal, vars = gene_intersect1, slot = "data")
lorig.ident <- Basal@meta.data$orig.ident
ldata <- data.frame(lgene_expression, lorig.ident)

##不考虑表达量为0的细胞
mean_without_zeros <- function(x) {
  mean_value <- mean(x[x != 0])
  ifelse(is.nan(mean_value), 0, mean_value)
}

lmedata <- ldata %>%
  group_by(lorig.ident) %>%
  summarize(across(all_of(gene_intersect1), mean_without_zeros, .names = "{.col}"))


lmedata <- as.data.frame(lmedata)
rownames(lmedata) <- lmedata[,1]
lmedata <- lmedata[,-1]
m <- t(lmedata)

sample_order <- c("S4", "S6", "S9", "S10", "S12", "S13", "S17",
                  "S3", "S11", "S14", "S15","S16",
                  "S1", "S2", "S5", "S7")
m <- m[, sample_order]
p10 <- pheatmap(m, 
                cluster_rows = TRUE, 
                cluster_cols = FALSE, 
                fontsize = 10,
                fontsize_row=10,
                fontsize_col=10,
                angle_col = c("90"),
                color = colorRampPalette(c("white" ,"#ffa500" ))(20))
p10
ggsave("../result/08pheat.svg", p10 , width = 7, height = 9, dpi = 300)

###241101重画heatmap
lgene_expression <- FetchData(Basal, vars = gene_intersect1, slot = "data")
lorig.ident <- Basal@meta.data$orig.ident
ldata <- data.frame(lgene_expression, lorig.ident)

mean_with_zeros <- function(x) {
  mean_value <- mean(x)
  ifelse(is.nan(mean_value), 0, mean_value)
}

lmedata <- ldata %>%
  group_by(lorig.ident) %>%
  summarize(across(all_of(gene_intersect1), mean_with_zeros, .names = "{.col}"))

lmedata <- as.data.frame(lmedata)
rownames(lmedata) <- lmedata[,1]
lmedata <- lmedata[,-1]
m <- t(lmedata)

sample_order <- c("S4", "S6", "S9", "S10", "S12", "S13", "S17",
                  "S3", "S11", "S14", "S15","S16",
                  "S1", "S2", "S5", "S7")
m <- m[, sample_order]
p10 <- pheatmap(m, 
                cluster_rows = TRUE, 
                cluster_cols = FALSE, 
                fontsize = 10,
                fontsize_row=10,
                fontsize_col=10,
                angle_col = c("90"),
                color = colorRampPalette(c("white" ,"#ffa500" ))(20))
p10
ggsave("../result/241101pheat.svg", p10 , width = 7, height = 9, dpi = 300)


target_gene <- "CDC42"  # 替换为实际基因名

# 提取目标基因在每个样本中的表达值
gene_expression_values <- lgene_expression[, target_gene]

# 创建一个数据框，包含样本标识和目标基因的表达值
gene_expression_df <- data.frame(orig.ident = lorig.ident, expression = gene_expression_values)

# 查看数据框
print(gene_expression_df)


############eleven figure BP
gene.df <- bitr(gene_intersect1,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count")+
    #ggtitle("KEGG") +
    theme_bw() +
    theme(text=element_text(size=16))
  pr
}


ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism='hsa',
                    pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)

ekegg <- setReadable(ekegg,'org.Hs.eg.db','ENTREZID')
pp <- erich2plot(ekegg@result)+ theme(axis.text.y=element_text(size=20)
                                      ,axis.title.x = element_blank()
                                      ,axis.title.y = element_blank()
                                      ,text=element_text(size=16))
pp
ggsave("KEGG.png", erich2plot(ekegg@result) , width = 7, height = 5, dpi = 300)
ggsave("08KEGG.svg", pp , width = 12, height = 7, dpi = 300)
ggsave("KEGG.bmp", erich2plot(ekegg@result) , width = 7, height = 5, dpi = 300)


head(ekegg)



#####go
##go

ggoMF <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)


p21 <- erich2plot(ggoMF@result)+ theme(axis.text.y=element_text(size=20)
                                       #,axis.title.x = element_blank()
                                       ,axis.title.y = element_blank()
                                       ,text=element_text(size=16))+labs(x='MF')
p21

ggsave("../redsult/08gomf.svg", p21 , width = 13.5, height = 7, dpi = 300)


ggoCC <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)

p22 <- erich2plot(ggoCC@result)+ theme(axis.text.y=element_text(size=20)
                                       #,axis.title.x = element_blank()
                                       ,axis.title.y = element_blank()
                                       ,text=element_text(size=16))+labs(x='CC')
p22
ggsave("..result/08gocc.svg", p22 , width = 15, height = 7, dpi = 300)


ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)

p12 <- erich2plot(ggoBP@result)+ theme(axis.text.y=element_text(size=20)
                                       ,axis.title.x = element_blank()
                                       ,axis.title.y = element_blank()
                                       ,text=element_text(size=16))

p12
ggsave("08gobp41gene.svg", p12 , width = 16, height = 7, dpi = 300)


goplot(ggoBP)
goplot(ggoCC)
goplot(ggoMF)

table(ggoBP)
write.csv(ggoBP@result,'../resultggobp.csv')
write.csv(ekegg@result,'../result/kegg.csv')


###################################20250228富集


result_nl_h <- read.csv('../result/08去0NL&H.csv',row.names = 1)
result_l_h <- read.csv('../result/08去0L&H.csv',row.names = 1)
result_l_nl <- read.csv('../result/08去0L&NL.csv',row.names = 1)
result_lnl_h <- read.csv('../result/08去0LNL&H.csv',row.names = 1)
#######pre valco plot
result_nl_h <- unique(na.omit(result_nl_h), by = "Gene")
result_l_h <- unique(na.omit(result_l_h), by = "Gene")
result_lnl_h <- unique(na.omit(result_lnl_h), by = "Gene")
result_l_nl <- unique(na.omit(result_l_nl), by = "Gene")

###filter p<0.01
re1 <- result_nl_h[result_nl_h$PValue<0.01,]
re2 <- result_l_h[result_l_h$PValue<0.01,]
re3 <- result_lnl_h[result_lnl_h$PValue<0.01,]
re4 <- result_l_nl[result_l_nl$PValue<0.01,]

filtered_df <- re4[re4$PValue < 0.05, ]  # 先筛选 p 值显著的
top50_combined <- filtered_df[order(abs(filtered_df$Log2FC), decreasing = TRUE), ][1:50, ]


#############################
gene.df <- bitr(top50_combined$Gene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)

erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count")+
    #ggtitle("KEGG") +
    theme_bw() +
    theme(text=element_text(size=16))
  pr
}


ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism='hsa',
                    pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)

ekegg <- setReadable(ekegg,'org.Hs.eg.db','ENTREZID')
pp <- erich2plot(ekegg@result)+ theme(axis.text.y=element_text(size=20)
                                      ,axis.title.x = element_blank()
                                      ,axis.title.y = element_blank()
                                      ,text=element_text(size=16))
pp

##go

ggoMF <- enrichGO(top50_combined$Gene, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
erich2plot(ggoMF@result)


ggoCC <- enrichGO(top50_combined$Gene, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)

erich2plot(ggoCC@result)

ggoBP <- enrichGO(top50_combined$Gene, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)

erich2plot(ggoBP@result)
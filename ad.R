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
sc2  # 查看合并后的对象
table(sc2@meta.data$Group)  # 查看各组样本数量

# 保存合并后的数据到result文件夹
saveRDS(sc2, '../result/step1_0822merged.rds')

s.integrated <- PercentageFeatureSet(sc2,'^MT',col.name = 'percent_MT')
s.integrated <- subset(s.integrated, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent_MT < 25)
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


# 计算线粒体基因的表达比例
s.integrated <- PercentageFeatureSet(sc2, '^MT', col.name = 'percent_MT')
# 根据基因数和线粒体基因表达比例过滤细胞
s.integrated <- subset(s.integrated, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent_MT < 25)
table(s.integrated@meta.data$Group)  # 查看过滤后各组细胞数量
# 数据标准化
s.integrated <- NormalizeData(s.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# 鉴定高变异基因
s.integrated <- FindVariableFeatures(s.integrated, selection.method = "vst", nfeatures = 2000)
# 数据缩放
s.integrated <- ScaleData(object = s.integrated, features = rownames(s.integrated))
# 主成分分析(PCA)
s.integrated <- RunPCA(s.integrated, features = VariableFeatures(object = s.integrated), reduction.name = "pca")
p_elbow <- ElbowPlot(s.integrated)  # 查看主成分贡献率，帮助确定使用多少个主成分
ggsave("../result/elbow_plot.png", p_elbow, width = 7, height = 5, dpi = 300)

# 使用Harmony进行批次效应校正
s.integrated <- RunHarmony(s.integrated, group.by.vars='orig.ident', reduction.save = "harmony")

# 使用校正后的数据进行t-SNE降维
s.integrated <- RunTSNE(s.integrated, reduction = "harmony", dims = 1:15)

# 寻找邻近点并聚类
s.integrated <- FindNeighbors(s.integrated, reduction = "pca", dims = 1:15)
s.integrated <- FindClusters(s.integrated, resolution = 0.1)  # resolution参数控制聚类粒度

# 保存校正后的数据到result文件夹
saveRDS(s.integrated, '../result/step2_0822harmonyed.rds')

# 加载参考数据集
s.se <- HumanPrimaryCellAtlasData()  # 来自celldex包 [[1]]

# 确保默认Assay是RNA
DefaultAssay(s.integrated) <- "RNA"

# 合并多层数据（关键步骤！）
s.integrated <- JoinLayers(s.integrated, layers = "data", assay = "RNA")  # [[2]][[3]]

# 提取标准化数据
datas <- GetAssayData(s.integrated, layer = "data")  # 此时不再报错

# 后续步骤（基因对齐、SingleR注释等）
common_genes <- intersect(rownames(datas), rownames(s.se))
datas <- datas[common_genes, ]
s.se <- s.se[common_genes, ]

pred.s <- SingleR(
  test = datas,
  ref = s.se,
  labels = s.se$label.main,
  clusters = s.integrated$seurat_clusters
)

# 将结果添加到Seurat对象
s.integrated <- AddMetaData(s.integrated, metadata = pred.s$labels, col.name = "celltype")

# 保存注释后的数据到result文件夹
saveRDS(s.integrated, '../result/step3_0822annoated.rds')

# 读取保存的数据
s.integrated <- readRDS('../result/step3_0822annoated.rds')


###########第一个图：细胞类型tSNE可视化###########

# 使用tSNE降维展示细胞类型分布


# 保存图形到result文件夹
ggsave("../result/08SingleR.png", p, width = 7, height = 5, dpi = 300)
ggsave("../result/08TSNE.svg", p, width = 6, height = 6, dpi = 300)


###########第二个图：细胞类型标志基因热图###########

# 设置分群依据为细胞类型
Idents(s.integrated) <- s.integrated@meta.data$celltype

# 寻找每个细胞类型的标志基因（仅考虑上调的基因）
markers <- FindAllMarkers(s.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 选择每个细胞类型表达量变化最大的前5个基因
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# 绘制热图
p1 <- DoHeatmap(s.integrated, features = top5$gene, label=FALSE) + 
  NoLegend() + 
  theme(legend.position = "right", text=element_text(size=15), 
        legend.text=element_text(size=15), axis.text.y=element_text(size=16)) + 
  scale_fill_gradient2(high = "#ff8c69", low = '#6495ed')  # 自定义颜色：从蓝色到橙色渐变

p1  # 显示热图

# 保存热图到result文件夹
ggsave("../result/08cellmarker.png", p1, width = 9, height = 9, dpi = 300)
ggsave("../result/08cellmarker.bmp", p1, width = 9, height = 9, dpi = 300)
ggsave("../result/08cellmarker.svg", p1, width = 9, height = 9)


############第三个图：细胞类型分布堆叠柱状图############

# 提取元数据，包括样本ID、分组和细胞类型
metadata <- FetchData(s.integrated, vars = c("orig.ident", 'Group', "celltype"))

# 统计每个样本中各细胞类型的数量
counts <- metadata %>%
  group_by(orig.ident, Group, celltype) %>%
  summarize(Count = n()) %>%
  ungroup()

# 自定义样本顺序
counts$orig.ident <- factor(counts$orig.ident, levels = c("S4", "S6", "S9", "S10", "S12", "S13", "S17",
                                                          "S3", "S11", "S14", "S15", "S16",
                                                          "S1", "S2", "S5", "S7"))

# 自定义细胞类型顺序
counts$celltype <- factor(counts$celltype, levels=c("DC", "Endothelial_cells", "Epithelial_cells", "Fibroblasts",
                                                    "Neurons", "Tissue_stem_cells", "T_cells", "Keratinocytes"))

# 在每个样本内计算各细胞类型的比例
counts$prop <- ave(counts$Count, counts$orig.ident, FUN = function(x) x/sum(x))

# 自定义颜色方案
col <- rev(c('#6495ed', '#F0F8FF', '#FFFFDB', '#C1CDCD', 
             '#ff8c69', '#EEC591', '#FFEBCD', '#b4eeb4'))

# 绘制堆叠柱状图
p2 <- ggplot(counts, aes(x = orig.ident, y = prop, fill = celltype)) +
  geom_col(position = 'stack', width = 0.6) +
  labs(x = "orig.ident", y = "Celltype(Proportion)", fill = "Celltype") +
  scale_fill_manual(values=col) + 
  theme_minimal(base_size = 20) +
  theme(text=element_text(size=16))

p2  # 显示图形

# 保存图形到result文件夹
ggsave("../result/celltypedistribution.png", p2, width = 7, height = 5, dpi = 300)
ggsave("../result/celltypedistribution.bmp", p2, width = 7, height = 5, dpi = 300)
ggsave("../result/08celltypedistribution去s8.svg", p2, width = 10, height = 5, dpi = 300)


###########第四个图：角质形成细胞（KC）比例分析###########

# 筛选出角质形成细胞的数据
dak <- filter(counts, celltype == "Keratinocytes")

# 自定义样本顺序
dak$orig.ident <- factor(dak$orig.ident, levels = c("S4", "S6", "S9", "S10", "S12", "S13", "S17",
                                                    "S3", "S11", "S14", "S15", "S16",
                                                    "S1", "S2", "S5", "S7"))

# 自定义分组顺序
dak$Group <- factor(dak$Group, levels = c("Healthy", "Non-lesional", "Lesional"))

# 绘制角质形成细胞比例的箱线图
p3 <- ggplot(dak, aes(x = Group, y = prop)) +
  geom_boxplot(width = 0.3, fill = "white", color = "black") +  # 设置箱体样式
  geom_jitter(width = 0.2, size=3, height = 0, aes(color = Group), alpha = 0.3) +  # 添加散点
  labs(x = "Group", y="Percentage", title = "Keratinocytes fraction") +
  # 添加参考线
  geom_hline(yintercept = 0.1, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 0.15, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 0.05, lty = 4, col = "black", lwd = 0.8) +
  theme(        
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.45),
    text=element_text(size=16), legend.position="none") +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA)

# 定义统计比较组（注释掉的代码）
# my_comparisons <- list(c("healthy","Non-lesional"), c("healthy", "Lesional"),
#                        c("Non-lesional", "Lesional"))
p3  # 显示图形
# 添加统计比较（注释掉的代码）
# p3 + stat_compare_means(comparisons = my_comparisons, method='t.test')

# 保存图形到result文件夹
ggsave("../result/08KCfraction.png", p3, width = 5, height = 5, dpi = 300)
ggsave("../result/08KCfraction.bmp", p3, width = 5, height = 5, dpi = 300)
ggsave("../result/08KCfraction去s8.svg", p3, width = 5, height = 5)


###########角质形成细胞亚群分析###########

# 提取角质形成细胞
KC_cluster <- subset(s.integrated, subset = celltype == "Keratinocytes")

# CDC42基因表达分析
gene_expression <- FetchData(KC_cluster, vars = c("CDC42"), slot = "data")
Group <- KC_cluster@meta.data$Group
data <- data.frame(gene_expression, Group)
data <- data[data$CDC42 != 0, ]  # 仅保留表达CDC42的细胞
data$Group <- factor(data$Group, levels = c("Healthy", "Non-lesional", "Lesional"))

# 定义统计比较组
my_comparisons <- list(c("Healthy", "Non-lesional"), c("Healthy", "Lesional"),
                       c("Non-lesional", "Lesional"))

# 绘制CDC42基因表达箱线图
p <- ggplot(data, aes(x = Group, y = CDC42, fill = Group)) +
  geom_boxplot(fill = "white", color = "black", width=0.2) +
  labs(x = "Group", y = "Gene Expression (CDC42) in Keratinocytes") +
  # 添加参考线
  geom_hline(yintercept = 2, lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = 2.5, lty = 4, col = "black", lwd = 0.8) +
  scale_y_continuous(breaks = seq(1, 4, by = 0.5)) +
  stat_compare_means(comparisons = my_comparisons) +  # 添加统计比较
  theme(        
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.45),
    text=element_text(size=16)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, color = "black", fill = NA) +
  labs(title = "CDC42")
p  # 显示图形

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

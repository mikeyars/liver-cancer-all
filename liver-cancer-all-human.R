# 加载必要的 R 包
library(Seurat)

# 加载数据
sample_paths <- c("g:/siglecell/humen/X1CA.matrix","g:/siglecell/humen/X2CA.matrix","g:/siglecell/humen/X3CA.matrix","g:/siglecell/humen/X4CA.matrix","g:/siglecell/humen/X1PCa.matrix","g:/siglecell/humen/X2PCa.matrix","g:/siglecell/humen/X3PCa.matrix")
seurat_list <- lapply(sample_paths, Read10X)

# 创建 Seurat 对象列表
seurat_list <- lapply(seurat_list, function(data) {
  # 创建 Seurat 对象
  seurat_object <- CreateSeuratObject(counts = data)
  
  # 过滤细胞、标准化等预处理步骤
  preprocess_sample <- function(seurat_object) {
    # 过滤细胞
    seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
    
    # 标准化数据
    seurat_object <- NormalizeData(seurat_object)
    
    # 寻找可变特征
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
    
    # 缩放数据
    seurat_object <- ScaleData(seurat_object)
    
    return(seurat_object)
  }
  
  # 对每个样本执行预处理步骤
  return(preprocess_sample(seurat_object))
})

# 输出预处理后的 Seurat 对象列表
seurat_list

# 获取整合数据的锚点（anchors）
anchors <- FindIntegrationAnchors(object.list = seurat_list)

# 整合数据，指定dims参数
integrated_seurat <- IntegrateData(anchor.features = anchors, dims = 1:20)

# 进行下游分析，例如主成分分析、聚类等
integrated_seurat <- RunPCA(integrated_seurat, npcs = 20)
integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:20)
integrated_seurat <- FindClusters(integrated_seurat)

# 运行UMAP降维
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:20)

# 可视化整合后的数据
TSNEPlot(integrated_seurat, group.by = "seurat_clusters")
UMAPPlot(integrated_seurat, group.by = "seurat_clusters")
# ...
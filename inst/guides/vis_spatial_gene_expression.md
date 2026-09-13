# `vis_compare_spatial_gene_expression_groups()` 使用说明

这个函数用于**已经完成整合/注释的空间转录组 Seurat 对象**，输入一个 gene vector，对多个 Visium 样本进行：

1. 全样本 gene expression 汇总；
2. 按临床变量（如 response、treatment、WHO subtype）进行 sample-level 统计比较；
3. 可选按 spatial domain / pathological region 分层；
4. 每个样本的 SpatialFeaturePlot；
5. 每个样本每个基因的 Moran's I 空间自相关；
6. 基因表达与 cell2location cell abundance 的 spot-level 相关，并在 sample level 汇总；
7. Excel + PDF + PNG 自动输出。

## 统计原则

默认把**sample/patient sample**作为独立统计单位。

对于每个 `sample × region × gene`，先在 spots 内计算：

- n_spots
- n_positive_spots
- pct_positive_spots
- mean_expression
- median_expression

然后使用 `sample_aggregation = "median"` 或 `"mean"` 得到每个样本的一个 expression 值，再进行临床组间 Wilcoxon / Kruskal-Wallis 检验。

这样可避免把同一张切片内的大量 spots 当作独立生物学重复。

---

## 1. 只看所有样本

```r
library(biomed)

gene2show <- c(
  "TACSTD2", "NECTIN4", "CD276",
  "MET", "DLL3", "CEACAM5"
)

res <- vis_compare_spatial_gene_expression_groups(
  vis,
  genes = gene2show,
  statistic_group = NULL,
  sample_col = "orig.ident",
  assay = "RNA",
  layer = "data",
  sample_aggregation = "median",
  output_dir = "spatial_gene_all_samples"
)
```

---

## 2. 按多个临床变量分别比较

```r
res <- vis_compare_spatial_gene_expression_groups(
  vis,
  genes = gene2show,
  statistic_group = c(
    "response",
    "treatment",
    "WHO_group"
  ),
  sample_col = "orig.ident",
  assay = "RNA",
  layer = "data",
  sample_aggregation = "median",
  output_dir = "spatial_gene_clinical"
)
```

每个 `statistic_group` 会独立分析，不会把多个临床变量拼成一个联合分组。

---

## 3. 按 spatial domain / pathological region 分层

例如 metadata 中已经有：

```r
table(vis$spatial_domain)
```

可以运行：

```r
res <- vis_compare_spatial_gene_expression_groups(
  vis,
  genes = gene2show,
  statistic_group = "response",
  sample_col = "orig.ident",
  split_by_region = "spatial_domain",
  assay = "RNA",
  layer = "data",
  output_dir = "spatial_gene_by_domain"
)
```

统计单位会变成：

```text
sample × spatial_domain × gene
```

例如分别比较 Tumor、Stroma、Immune-rich region 内 response vs non-response 的基因表达。

---

## 4. cell2location abundance 在 metadata 中

如果 cell2location 的结果已经写入：

```r
colnames(vis@meta.data)
```

例如：

```text
T_cell
B_cell
Macrophage
Fibroblast
Tumor
```

运行：

```r
res <- vis_compare_spatial_gene_expression_groups(
  vis,
  genes = gene2show,
  statistic_group = "response",
  sample_col = "orig.ident",

  cell2location_metadata_cols = c(
    "T_cell",
    "B_cell",
    "Macrophage",
    "Fibroblast",
    "Tumor"
  ),

  correlation_method = "spearman",
  min_spots_for_correlation = 20,

  output_dir = "spatial_gene_cell2location"
)
```

函数会在**每一个样本内部**计算：

```text
gene expression across spots
        vs
cell2location abundance across spots
```

得到：

```text
sample × gene × cell type × correlation
```

之后还可以按 `response` 等临床变量比较这些 per-sample correlation。

---

## 5. cell2location abundance 在单独 assay 中

假设：

```r
vis[["cell2location"]]
```

的 rows/features 是 cell types，columns 是 spots：

```r
res <- vis_compare_spatial_gene_expression_groups(
  vis,
  genes = gene2show,
  statistic_group = "response",
  sample_col = "orig.ident",

  cell2location_assay = "cell2location",
  cell2location_layer = "data",

  output_dir = "spatial_gene_cell2location"
)
```

可以通过：

```r
cell2location_features = c(
  "T_cell",
  "Macrophage",
  "Fibroblast"
)
```

只分析指定 cell types。

---

## 6. Moran's I

默认：

```r
run_moran = TRUE
```

用于判断每个基因在每张空间切片中是否具有空间聚集/空间自相关。

对于非常大的切片，默认：

```r
moran_max_spots = 8000
moran_large_sample = "skip"
```

避免 Moran's I 占用过多内存。

如果希望随机抽样后计算：

```r
moran_large_sample = "subsample"
```

---

## 7. 非标准 Seurat layer

如果使用：

```r
layer = "data"
layer = "counts"
layer = "scale.data"
```

函数优先调用 Seurat `SpatialFeaturePlot()`，可叠加组织图像。

如果使用自定义 layer，例如：

```r
layer = "scvi_normalized"
```

表达统计仍然可以运行；空间图在 `SpatialFeaturePlot()` 无法读取该 layer 时会自动回退为基于 tissue coordinates 的空间散点图。

---

## 主要输出

```text
spatial_gene_comparison/
├── SpatialGeneExpression_results.xlsx
├── 00_all_samples/
│   └── SampleHeatmap/
├── clinical_response/
│   ├── DotPlot/
│   └── SampleBoxplot/
├── clinical_treatment/
├── SpatialFeaturePlot/
├── SpatialRegionPlot/
├── MoranI/
└── Cell2locationCorrelation/
```

Excel 中主要包括：

- `Parameters`
- `Image_Sample_Map`
- `AllSample_Summary`
- `Sample_Gene_Summary`
- `MoranI`
- `C2L_Correlations`
- 每个临床变量的 Group Summary / Overall Test / Pairwise Test
- cell2location correlation 的临床组比较

## 建议

对于正式的组间统计，优先使用：

```r
sample_aggregation = "median"
```

并保持 patient/sample 为独立统计单位。

spot-level 数据更适合用于：

- SpatialFeaturePlot
- positive spot fraction
- Moran's I
- gene–cell type abundance spatial correlation
- spatial domain / pathological region 内部分析

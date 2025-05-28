library(pheatmap)

# 1. 读取 CSV 文件并转换为矩阵
cts <- read.csv("D:\\Bioana\\dbparsemergedv1\\grouped_trimmed.csv",
                sep = '\t',
                row.names = 1,
                check.names = FALSE)

mat <- as.matrix(cts)

# 2. 计算 Lactobacillus crispatus 和 Lactobacillus iners 的数量
col_names <- colnames(cts)

# 计算 cris_count，从 "CAZyme" 之后开始数
cris_count <- sum(grepl("^Lactobacillus crispatus", col_names[-1]))

# 计算 total_col（不包括 "CAZyme"）
total_col <- length(col_names)

# 计算 iners_count
iners_count <- total_col - cris_count

# 确定插入位置（最后一个 L. crispatus 之后）
insert_after <- cris_count + 1


cat("Lactobacillus crispatus 数量: ", cris_count, "\n")
cat("Lactobacillus iners 数量: ", iners_count, "\n")
cat("插入新列的位置: ", insert_after, "\n")

# 3. 插入额外的列
n_extra <- 3  # 需要插入的列数

# 检查原始矩阵的列数是否足够
if(ncol(mat) < insert_after){
  stop(paste("原始矩阵的列数少于", insert_after, "，无法插入额外的列。"))
}

# 定义新列的值（对于二元数据 0/1，新值设为 2）
new_val <- max(mat, na.rm = TRUE) + 1

# 创建新列矩阵，并命名为 "-"
new_cols <- matrix(new_val, nrow = nrow(mat), ncol = n_extra)
colnames(new_cols) <- rep("-", n_extra)

# 将新列插入到指定位置
mat_new <- cbind(mat[, 1:insert_after], new_cols, mat[, (insert_after + 1):ncol(mat)])

# 取消列聚类：直接组合左、中、右列
final_col_order <- colnames(mat_new)

# 按照定义的顺序重新排列矩阵的列
mat_final <- mat_new[, final_col_order]

# 4. 添加列分组注释
annotation_col <- data.frame(
  Species = c(rep("L. crispatus", cris_count), rep("divide", n_extra), rep("L. iners", iners_count))
)

# **去除列名重复**
make.unique.colnames <- function(colnames_vec) {
  return(make.unique(colnames_vec))
}

colnames(mat_final) <- make.unique.colnames(colnames(mat_final))

cat("mat_final 列数: ", ncol(mat_final), "\n")
cat("annotation_col 行数: ", nrow(annotation_col), "\n")
cat("annotation_col 前几行:\n")
print(head(annotation_col))
cat("mat_final 列名前几项:\n")
print(head(colnames(mat_final)))

# 确保列注释行数和矩阵列数一致
if (ncol(mat_final) != nrow(annotation_col)) {
  stop("列注释(annotation_col)的行数与矩阵(mat_final)的列数不匹配！")
}

# 现在可以安全地设置 `rownames(annotation_col)`
rownames(annotation_col) <- colnames(mat_final)


# 定义注释颜色
annotation_colors <- list(
  Species = c("L. crispatus" = "orange", "divide" = "red", "L. iners" = "grey")
)

# 5. **定义颜色调色板**
# 0: blue, 1: green, 2: red
color_palette <- c("blue", "green", "red")

# 定义断点，确保每个值有单独的区间
breaks <- c(-0.5, 0.5, 1.5, 2.5)

# 6. **生成热图并保存为 PDF 文件**
pdf("D:\\Bioana\\dbparsemergedv1\\grouped_trimmed.pdf", width = 40, height = 15)
pheatmap(mat_final,
         color = color_palette,       # 使用定义的颜色调色板
         breaks = breaks,             # 定义断点
         cluster_rows = FALSE,        # 不聚类行
         cluster_cols = FALSE,        # 不聚类列，保持原始顺序
         annotation_col = annotation_col,  # 添加列分组注释
         annotation_colors = annotation_colors,  # 定义注释颜色
         show_rownames = TRUE,
         show_colnames = FALSE,       # 不显示列名
         annotation_names_col = TRUE,  # 显示顶部标签行
         annotation_legend = TRUE,    # 显示图例
         fontsize = 25,               # 增加整体字体大小
         fontsize_row = 20,
         fontsize_col = 15,
         legend_breaks = c(0, 1, 2),  # 保证图例中显示所有的类别
         legend_labels = c("Absent", "Present", "Inserted"), # 自定义图例标签
         angle_col = 90)              # 设置列名角度为垂直
dev.off()

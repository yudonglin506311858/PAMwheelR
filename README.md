# PAMwheelR

PAMwheelR 是一个用于可视化 PAM (Protospacer Adjacent Motif) 序列的 R 包。它可以将不同长度和频率的 PAM 序列以轮状图(wheel diagram)的形式展示，支持任意层数的嵌套可视化。

## 安装

```r
# 安装依赖
if (!require("remotes")) install.packages("devtools")

# 从GitHub安装
devtools::install_github("https://github.com/yudonglin506311858/PAMwheelR.git")

```

## 示例数据
![image](https://github.com/user-attachments/assets/91a91edc-84e3-46b5-8a7b-bddb673fa10a)

## 示例使用

```r
library(PAMwheelR)
library(dplyr)
library(grDevices)
library(ggplot2)

# 生成随机PAM矩阵函数
generate_pam_matrix <- function(codon_length = 3, seed = NULL) {  
  if (!is.null(seed)) set.seed(seed)  
  bases <- c("A", "T", "C", "G")
  codons <- replicate(4^codon_length, {
    paste(sample(bases, codon_length, replace = TRUE), collapse = "")
  })
  ratios <- rexp(4^codon_length, rate = 1)
    data.frame(
    codon = codons,
    ratio = ratios
  )
}


pam_matrix <- generate_pam_matrix(codon_length = 2, seed = 123)
plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 2)

pam_matrix <- generate_pam_matrix(codon_length = 3, seed = 123)
plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 3)

pam_matrix <- generate_pam_matrix(codon_length = 4, seed = 123)
plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 4)

pam_matrix <- generate_pam_matrix(codon_length = 5, seed = 123)
plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 5)

```


![image](https://github.com/user-attachments/assets/aef7cea5-4c62-4cb5-9d9a-43df3e1d8674)

## 参考文章

https://pmc.ncbi.nlm.nih.gov/articles/PMC4826307/

https://www.science.org/doi/10.1126/sciadv.abb3350

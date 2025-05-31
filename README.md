# PAMwheelR

## 安装

```r
# 安装依赖
if (!require("remotes")) install.packages("devtools")

# 从GitHub安装
devtools::install_github("https://github.com/yudonglin506311858/PAMwheelR.git")

```

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

head(pam_matrix)

plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 2)

pam_matrix <- generate_pam_matrix(codon_length = 3, seed = 123)

head(pam_matrix)

plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 3)

pam_matrix <- generate_pam_matrix(codon_length = 4, seed = 123)

head(pam_matrix)

plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 4)

pam_matrix <- generate_pam_matrix(codon_length = 5, seed = 123)

```


![image](https://github.com/user-attachments/assets/aef7cea5-4c62-4cb5-9d9a-43df3e1d8674)


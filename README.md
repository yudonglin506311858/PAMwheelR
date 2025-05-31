# PAMwheelR

source("R/plot_pam_wheel.R")

library(dplyr)
library(grDevices)
library(ggplot2)

# 生成随机PAM矩阵函数
generate_pam_matrix <- function(codon_length = 3, seed = NULL) {
  # 设置随机种子（可选）
  if (!is.null(seed)) set.seed(seed)
  
  # 定义碱基集合
  bases <- c("A", "T", "C", "G")
  
  # 生成随机密码子
  codons <- replicate(4^codon_length, {
    paste(sample(bases, codon_length, replace = TRUE), collapse = "")
  })
  
  # 生成随机比例（指数分布，产生长尾效果）
  ratios <- rexp(4^codon_length, rate = 1)
  
  # 返回数据框
  data.frame(
    codon = codons,
    ratio = ratios
  )
}

# 示例使用
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


![image](https://github.com/user-attachments/assets/aef7cea5-4c62-4cb5-9d9a-43df3e1d8674)


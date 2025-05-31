#' Plot PAM Wheel Diagram (支持任意层数)
#'
#' @param pam_matrix A data frame with columns: codon (k-letter strings) and ratio (numeric)
#' @param inner_radius Relative radius of innermost circle (0-1)
#' @param n_layers Number of layers to display (default: 3)
#' @return ggplot object
#'
#' @examples
#' # Generate random PAM matrix
#' generate_pam_matrix <- function(codon_length = 3, seed = NULL) {
#'   if (!is.null(seed)) set.seed(seed)
#'   bases <- c("A", "T", "C", "G")
#'   codons <- replicate(4^codon_length, {
#'     paste(sample(bases, codon_length, replace = TRUE), collapse = "")
#'   })
#'   ratios <- rexp(4^codon_length, rate = 1)
#'   data.frame(codon = codons, ratio = ratios)
#' }
#'
#' # 2-layer example
#' pam_matrix <- generate_pam_matrix(codon_length = 2, seed = 123)
#' plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 2)
#'
#' # 3-layer example
#' pam_matrix <- generate_pam_matrix(codon_length = 3, seed = 123)
#' plot_pam_wheel(pam_matrix, inner_radius = 0.3, n_layers = 3)
#'
#' @export
plot_pam_wheel <- function(pam_matrix, inner_radius = 0.2, n_layers = 3) {
  # 检查输入格式
  if (!all(c("codon", "ratio") %in% colnames(pam_matrix))) {
    stop("Input must contain 'codon' and 'ratio' columns")
  }

  # 获取密码子长度
  codon_length <- unique(nchar(as.character(pam_matrix$codon)))
  if (length(codon_length) > 1) {
    stop("All codons must have the same length")
  }

  # 验证层数
  if (n_layers > codon_length) {
    stop("Number of layers cannot exceed codon length")
  }

  # 定义基础颜色（最内层）
  base_colors <- c(
    A = "#79c9cb",  # 浅蓝绿色
    C = "#a0d071",  # 浅绿色
    G = "#d17575",  # 浅红色
    T = "#9875b1"   # 浅紫色
  )

  # 定义最外层颜色
  outer_colors <- c(
    A = "#a6d7db",  # 浅蓝绿色（更浅）
    C = "#bbdd9f",  # 浅绿色（更浅）
    G = "#dea4a1",  # 浅红色（更浅）
    T = "#bb9fcb"   # 浅紫色（更浅）
  )

  # 为每个碱基创建渐变色生成器
  color_generators <- list(
    A = grDevices::colorRampPalette(c(base_colors["A"], outer_colors["A"])),
    C = grDevices::colorRampPalette(c(base_colors["C"], outer_colors["C"])),
    G = grDevices::colorRampPalette(c(base_colors["G"], outer_colors["G"])),
    T = grDevices::colorRampPalette(c(base_colors["T"], outer_colors["T"]))
  )

  # 提取碱基信息 - 更安全的方式
  pam_data <- pam_matrix %>%
    dplyr::mutate(codon = as.character(codon))

  # 添加所有需要的碱基列
  for (i in 1:n_layers) {
    base_col <- paste0("base", i)
    pam_data <- pam_data %>%
      dplyr::mutate(!!base_col := substr(codon, i, i))
  }

  # 动态生成层数据
  layer_list <- list()

  for (i in 1:n_layers) {
    layer_name <- paste0("layer", i)
    base_cols <- paste0("base", 1:i)

    # 计算分组键
    group_keys <- base_cols

    # 计算比例
    layer_data <- pam_data

    # 如果有父层，先按父层分组
    if (i > 1) {
      parent_cols <- paste0("base", 1:(i-1))
      layer_data <- layer_data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(parent_cols)))
    }

    # 计算当前层比例
    layer_data <- layer_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_keys)), .add = TRUE) %>%
      dplyr::summarise(total = sum(ratio), .groups = "drop")

    # 重新按父层分组计算比例
    if (i > 1) {
      layer_data <- layer_data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(parent_cols))) %>%
        dplyr::mutate(
          group_total = sum(total),
          prop = total / group_total,
          end = cumsum(prop),
          start = dplyr::lag(end, default = 0)
        ) %>%
        dplyr::ungroup()
    } else {
      layer_data <- layer_data %>%
        dplyr::mutate(
          group_total = sum(total),
          prop = total / group_total,
          end = cumsum(prop),
          start = dplyr::lag(end, default = 0)
        )
    }

    # 添加其他信息
    layer_data <- layer_data %>%
      dplyr::mutate(
        layer = layer_name,
        label = !!rlang::sym(paste0("base", i)),
        parent = if (i > 1) {
          apply(dplyr::select(., dplyr::all_of(parent_cols)), 1, paste0, collapse = "")
        } else {
          NA_character_
        },
        full_label = apply(dplyr::select(., dplyr::all_of(base_cols)), 1, paste0, collapse = ""),
        layer_num = i
      )

    # 添加颜色值
    if (i == 1) {
      layer_data <- layer_data %>%
        dplyr::mutate(
          color_value = base_colors[.data$base1]
        )
    } else {
      layer_data <- layer_data %>%
        dplyr::mutate(
          color_value = purrr::map2_chr(
            .data$base1, .data$full_label,
            ~ color_generators[[.x]](n_layers)[i]
          )
        )
    }

    layer_list[[i]] <- layer_data
  }

  # 合并所有层数据
  plot_data <- dplyr::bind_rows(layer_list) %>%
    dplyr::mutate(
      layer = factor(layer, levels = paste0("layer", 1:n_layers)),
      # 设置半径
      radius_min = inner_radius * (layer_num - 1),
      radius_max = inner_radius * layer_num
    )

  # 计算全局角度（确保嵌套对齐）
  plot_data_angles <- NULL

  for (i in 1:n_layers) {
    layer_name <- paste0("layer", i)
    layer_data <- plot_data %>% dplyr::filter(layer == layer_name)

    if (i == 1) {
      # 第一层角度计算
      layer_angles <- layer_data %>%
        dplyr::mutate(
          global_start = 2 * pi * start,
          global_end = 2 * pi * end
        )
    } else {
      # 基于父层计算角度
      parent_layer <- paste0("layer", i-1)
      parent_angles <- plot_data_angles %>%
        dplyr::filter(layer == parent_layer) %>%
        dplyr::select(full_label, parent_start = global_start, parent_end = global_end)

      layer_angles <- layer_data %>%
        dplyr::inner_join(
          parent_angles,
          by = c("parent" = "full_label")
        ) %>%
        dplyr::mutate(
          global_start = parent_start + start * (parent_end - parent_start),
          global_end = parent_start + end * (parent_end - parent_start)
        )
    }

    # 添加到角度数据
    plot_data_angles <- dplyr::bind_rows(plot_data_angles, layer_angles)
  }

  # 生成扇区多边形坐标的函数
  generate_sector <- function(start, end, r_min, r_max, group) {
    n_points <- 50
    seq_angles <- seq(start, end, length.out = n_points)

    # 外圆弧
    outer_x <- r_max * cos(seq_angles)
    outer_y <- r_max * sin(seq_angles)

    # 内圆弧（反向）
    inner_x <- r_min * cos(rev(seq_angles))
    inner_y <- r_min * sin(rev(seq_angles))

    data.frame(
      x = c(outer_x, inner_x),
      y = c(outer_y, inner_y),
      group = group
    )
  }

  # 为每个扇区创建多边形
  plot_data_angles <- plot_data_angles %>%
    dplyr::mutate(sector_id = dplyr::row_number())

  polygons <- purrr::pmap_dfr(
    list(
      plot_data_angles$global_start,
      plot_data_angles$global_end,
      plot_data_angles$radius_min,
      plot_data_angles$radius_max,
      plot_data_angles$sector_id
    ),
    generate_sector
  )

  # 添加扇区信息到多边形数据
  polygons <- polygons %>%
    dplyr::inner_join(
      plot_data_angles %>% dplyr::select(sector_id, layer, color_value, label, start, end, layer_num),
      by = c("group" = "sector_id")
    )

  # 创建标签数据
  labels <- plot_data_angles %>%
    dplyr::mutate(
      mid_angle = (global_start + global_end) / 2,
      # 调整标签位置：在半径范围内居中
      #label_r = (radius_min + radius_max) / 2,
      # 调整标签位置：第一层的标签向外偏移，其他层保持原来的计算方式
      label_r = dplyr::case_when(
        layer == "layer1" ~ (radius_min + radius_max) / 2 + 0.07,  # 第一层向外偏移0.1
        TRUE ~ (radius_min + radius_max) / 2  # 其他层保持原样
      ),
      label_x = label_r * cos(mid_angle),
      label_y = label_r * sin(mid_angle),
      # 统一标签大小
      label_size = 7,
      # 计算标签角度 (使文字垂直于圆心)
      label_angle = (mid_angle * 180 / pi) + 270,
      # 只显示足够大的扇区标签
      show_label = (end - start) > (0.1 / layer_num)  # 根据层数调整阈值
    )

  # 创建白色同心圆数据（带黑色边框）
  white_circle_radius <- inner_radius * 0.5
  circle_data <- data.frame(
    x = white_circle_radius * cos(seq(0, 2*pi, length.out = 100)),
    y = white_circle_radius * sin(seq(0, 2*pi, length.out = 100))
  )

  # 绘制图形
  p <- ggplot2::ggplot() +
    # 添加PAM序列扇区
    ggplot2::geom_polygon(
      data = polygons,
      ggplot2::aes(x = x, y = y, group = group, fill = color_value),
      color = "white", linewidth = 0.5
    ) +
    # 添加白色同心圆
    ggplot2::geom_polygon(
      data = circle_data,
      ggplot2::aes(x = x, y = y),
      fill = "white", color = "black", linewidth = 1
    ) +
    # 添加标签
    ggplot2::geom_text(
      data = labels %>% dplyr::filter(show_label),
      ggplot2::aes(
        x = label_x, y = label_y,
        label = label,
        angle = label_angle,
        size = label_size
      ),
      color = "black", show.legend = FALSE
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::labs(fill = "Base") +
    ggplot2::theme(
      legend.position = "bottom",
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
    )

  return(p)
}

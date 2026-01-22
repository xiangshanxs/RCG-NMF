library(Seurat)
library(rliger)
library(dplyr)

#### 第一步：导入单细胞数据 ####
sc_input = readRDS("./scRef_input_mainExample.RDS")
sc_count = sc_input$sc_count
sc_meta = sc_input$sc_meta
ct.varname = sc_input$ct.varname
sample.varname = sc_input$sample.varname

# 获取单细胞数据的所有基因
sc_genes <- rownames(sc_count)
cat("Number of genes in scRNA-seq:", length(sc_genes), "\n")

#### 第二步：导入和处理空间转录组数据 ####
base_dir <- "E:/PyWorkSpace/pycharm/DLPFC"
sample_ids <- c("151507", "151508", "151509", "151510")

expr_list <- list()
coord_list <- list()
label_list <- list()

for (sample in sample_ids) {
  cat("Processing sample:", sample, "\n")
  sample_dir <- file.path(base_dir, sample)
  h5_file <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
  metadata_file <- file.path(sample_dir, "metadata.csv")
  spatial_file <- file.path(sample_dir, "spatial/tissue_positions_list.csv")
  
  seurat_data <- Read10X_h5(file = h5_file)
  expr_matrix <- as.matrix(seurat_data)
  
  TC <- read.csv(metadata_file, header = TRUE)
  true_labels <- TC$ground_truth
  
  coord_matrix <- as.matrix(read.csv(spatial_file, header = FALSE))
  coord_matrix <- coord_matrix[, c(1, 5:6)]
  rownames(coord_matrix) <- coord_matrix[, 1]
  coord_matrix <- coord_matrix[, -1]
  colnames(coord_matrix) <- c("x", "y")
  
  valid_idx <- !is.na(true_labels)
  expr_matrix <- expr_matrix[, valid_idx]
  true_labels <- true_labels[valid_idx]
  TC <- TC[valid_idx, ]
  rownames(TC) <- TC[, 1]
  TC <- TC[, -1]
  coord_matrix <- coord_matrix[rownames(TC), ]
  
  coord_df <- as.data.frame(coord_matrix)
  coord_df$x <- as.numeric(coord_df$x)
  coord_df$y <- as.numeric(coord_df$y)
  coord_df$x <- coord_df$x - min(coord_df$x)
  coord_df$y <- coord_df$y - min(coord_df$y)
  scaleFactor <- max(coord_df$x, coord_df$y)
  coord_df$x <- coord_df$x / scaleFactor
  coord_df$y <- coord_df$y / scaleFactor
  coord_df <- as.matrix(coord_df[, c("x", "y")])
  zero_thresh <- 1
  expr_matrix <- expr_matrix[rowMeans(expr_matrix == 0) < zero_thresh, ]
  expr_list[[sample]] <- expr_matrix
  coord_list[[sample]] <- coord_df
  label_list[[sample]] <- true_labels
}

#### 第三步：处理空间数据 - 使用rliger筛选高变异基因 ####
ifnb_liger <- createLiger(expr_list, remove.missing = F, take.gene.union = T, verbose = F)
ifnb_liger <- rliger::normalize(ifnb_liger)
ifnb_liger <- selectGenes(ifnb_liger)  # 筛选空间数据的高变异基因
ifnb_liger <- scaleNotCenter(ifnb_liger)

# 获取空间数据的高变异基因
st_hvg <- ifnb_liger@varFeatures
cat("Number of HVG in spatial data:", length(st_hvg), "\n")

#### 第四步：获取共同基因（空间HVG与单细胞所有基因的交集） ####
common_genes <- intersect(sc_genes, st_hvg)
cat("Number of common genes:", length(common_genes), "\n")

# 检查共同基因数量是否足够
if (length(common_genes) < 100) {
  warning("Too few common genes! Consider increasing HVG number in spatial data.")
}

#### 第五步：基于共同基因生成细胞类型表达矩阵 ####
# 使用单细胞数据中的共同基因
sc_count_common <- sc_count[common_genes, ]

# 生成细胞类型平均表达矩阵
cell_types <- unique(sc_meta[[ct.varname]])
celltype_matrix <- matrix(0, nrow = length(common_genes), ncol = length(cell_types))
rownames(celltype_matrix) <- common_genes
colnames(celltype_matrix) <- cell_types

for (ct in cell_types) {
  ct_cells <- which(sc_meta[[ct.varname]] == ct)
  if (length(ct_cells) > 0) {
    celltype_matrix[, ct] <- rowMeans(sc_count_common[, ct_cells, drop = FALSE])
  }
}

cat("Cell type matrix dimensions:", dim(celltype_matrix), "\n")

Y_list <- lapply(sample_ids, function(sample) {
  # 提取 scaleData 并转换为普通矩阵（如果不是 matrix 格式）
  as.matrix(ifnb_liger@datasets[[sample]]@scaleData)
})
# 添加切片名称作为列表元素名
names(Y_list) <- sample_ids
#### 第六步：过滤空间数据只保留共同基因 ####
for (sample in sample_ids) {
  # 只保留共同基因
  Y_list[[sample]] <- Y_list[[sample]][common_genes, ]
}





library(RANN)
library(Matrix)
createA <- function(locationList) {
  nSlices = length(locationList)
  AList = list()
  for (islice in 1:nSlices) {
    location = as.data.frame(locationList[[islice]])
    norm_cords = location[, c("x", "y")]
    rownames(norm_cords) <- rownames(location)
    ineibor = 11
    near_data = nn2(norm_cords[, 1:2], k = ineibor)
    neibors = near_data$nn.idx
    neibors = neibors[, -1]
    Nmat = Matrix(0, nrow = nrow(neibors), ncol = nrow(neibors), sparse = TRUE)
    for (icol in 1:ncol(neibors)) {
      edges = data.frame(i = 1:nrow(neibors), j = neibors[, icol])
      adjacency = sparseMatrix(i = as.integer(edges$i),
                               j = as.integer(edges$j), x = 1, dims = rep(nrow(neibors), 2), use.last.ij = TRUE)
      Nmat = Nmat + adjacency
    }
    Nmat = Nmat * t(Nmat)
    rownames(Nmat) = colnames(Nmat) = rownames(norm_cords)
    AList[[islice]] = Nmat
  }
  return(AList)
}

adj_list <- createA(coord_list)

K <- 3  # 高阶扩散最大阶数
alpha_func <- function(k) { 1 / k }  # 衰减策略，可改为 exp(-k)
#alpha_func <- function(k) {exp(-k)}
# 初始化输出列表
PP_list <- list()         # 存储转移概率矩阵
A_final_list <- list()   # 存储多阶融合邻接矩阵
D_final_list <- list()         # 存储拉普拉斯矩阵（非归一化）

# 遍历每个切片的邻接矩阵
for (i in seq_along(adj_list)) {
  A <- adj_list[[i]]
  n <- nrow(A)
  
  # Step 1: 构造转移概率矩阵 P = D^{-1} A
  row_sums <- rowSums(A)
  PP <- A / row_sums
  PP[is.na(PP)] <- 0  # 处理孤立节点
  
  PP_list[[i]] <- PP
  
  # Step 2: 多阶扩散并加权
  A_final <- matrix(0, n, n)
  PP_k <- PP  # 初始为 P^1
  
  for (k in 1:K) {
    alpha_k <- alpha_func(k)
    A_final <- A_final + alpha_k * PP_k
    PP_k <- PP_k %*% PP
  }
  
  A_final_list[[i]] <- A_final
  
  # Step 3: 构造拉普拉斯矩阵（非归一化）
  D_final <- diag(rowSums(A_final))
  #L <- D_final - A_final
  
  D_final_list[[i]] <- D_final
}


kmeansFunc<-function (data, k) 
{
  set.seed(12345678)
  if (nrow(data) < 3e+05) {
    numStart = 100
  } else {
    numStart = 1
  }
  cl <- suppressWarnings(try(kmeans(data, k, nstart = numStart, 
                                    iter.max = 100), silent = TRUE))
  
  # 返回聚类标签和聚类质心
  return(list(cluster = cl$cluster, centers = cl$centers, size = cl$size))
}
infoNCE <- function(P_list, Kcu, tau) {
  Z <- list()
  cluster_size_list <- list()
  kmeans_result_list <- list()
  
  for (i in seq_along(P_list)) {
    sample_id <- names(P_list)[i]
    P <- P_list[[i]]
    kmeans_result <- kmeansFunc(P, Kcu)
    centers <- kmeans_result[["centers"]]
    sizes <- kmeans_result[["size"]]
    
    Z[[sample_id]] <- centers
    cluster_size_list[[sample_id]] <- sizes
    kmeans_result_list[[sample_id]] <- kmeans_result
  }
  
  library(Matrix)
  library(proxy)
  
  d <- ncol(Z[[1]])  # 特征维度
  Q_list <- list()
  W_list <- list()
  
  for (s_query in seq_along(Z)) {
    sample_name <- names(Z)[s_query]
    Z_s <- Z[[sample_name]]  # 聚类中心矩阵
    Q_s <- matrix(0, nrow = nrow(Z_s), ncol = d)
    W_s <- matrix(0, nrow = nrow(Z_s), ncol = d)
    
    for (k_query in 1:nrow(Z_s)) {
      query_centroid <- Z_s[k_query, , drop = FALSE]
      
      # ------- 正样本 -------
      pos_samples <- list()
      pos_indices <- list()
      for (s_other in seq_along(Z)) {
        if (s_other == s_query) next
        Z_other <- Z[[s_other]]
        sims <- proxy::simil(query_centroid, Z_other, method = "cosine")
        sims <- as.vector(sims)
        best_idx <- which.max(sims)
        pos_samples[[length(pos_samples) + 1]] <- Z_other[best_idx, , drop = FALSE]
        pos_indices[[length(pos_indices) + 1]] <- list(slice = s_other, region = best_idx)
      }
      
      pos_matrix <- do.call(rbind, pos_samples)
      euclidean_distances <- apply(pos_matrix, 1, function(x) exp(sqrt(sum((x - query_centroid)^2))))
      w <- exp(1 / euclidean_distances) / sum(exp(1 / euclidean_distances))
      weighted_pos_samples <- matrix(w, nrow = 1) %*% pos_matrix
      
      # ------- 负样本 -------
      neg_samples <- list()
      for (s in seq_along(Z)) {
        Z_s_other <- Z[[s]]
        keep_rows <- rep(TRUE, nrow(Z_s_other))
        if (s == s_query) {
          keep_rows[k_query] <- FALSE
        } else {
          for (p in pos_indices) {
            if (p$slice == s) keep_rows[p$region] <- FALSE
          }
        }
        if (any(keep_rows)) {
          neg_samples[[length(neg_samples) + 1]] <- Z_s_other[keep_rows, , drop = FALSE]
        }
      }
      neg_matrix <- do.call(rbind, neg_samples)
      
      # ------- 相似度与导数函数 -------
      cosine_similarity <- function(x, y) {
        sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
      }
      custom_function1 <- function(x, y) {
        norm_x <- sqrt(sum(x^2))
        norm_y <- sqrt(sum(y^2))
        numerator <- sum(y) * (norm_x^2)
        denominator <- norm_y * (norm_x^3)
        return(numerator / denominator)
      }
      custom_function2 <- function(x, y) {
        norm_x <- sqrt(sum(x^2))
        norm_y <- sqrt(sum(y^2))
        numerator <- sum(x %*% t(y)) * x
        denominator <- norm_y * (norm_x^3)
        return(numerator / denominator)
      }
      
      cos_sim_positive <- cosine_similarity(query_centroid, weighted_pos_samples)
      loss_positive <- exp(cos_sim_positive / tau)
      
      sum_neg_sim <- 0
      sum_neg_dao1 <- 0
      sum_neg_dao2 <- rep(0, d)
      for (j in 1:nrow(neg_matrix)) {
        neg_vec <- neg_matrix[j, ]
        cos_sim_negative <- cosine_similarity(query_centroid, neg_vec)
        sim_term <- exp(cos_sim_negative / tau)
        sum_neg_sim <- sum_neg_sim + sim_term
        sum_neg_dao1 <- sum_neg_dao1 + sim_term * custom_function1(query_centroid, neg_vec)
        sum_neg_dao2 <- sum_neg_dao2 + sim_term * custom_function2(query_centroid[1, ], neg_vec)
      }
      
      cluster_size <- cluster_size_list[[sample_name]][k_query]
      weight_factor <- 1 / cluster_size
      N <- sum(sapply(Z, nrow))  # 所有聚类中心总数
      
      Q_s[k_query, ] <- weight_factor * (1 / (N * tau)) * (
        custom_function2(query_centroid, weighted_pos_samples) +
          (1 / (loss_positive + sum_neg_sim)) *
          (loss_positive * custom_function1(query_centroid, weighted_pos_samples) + sum_neg_sim * sum_neg_dao1)
      )
      
      W_s[k_query, ] <- weight_factor * (1 / (N * tau)) * (
        custom_function1(query_centroid, weighted_pos_samples) +
          (1 / (loss_positive + sum_neg_sim)) *
          (loss_positive * custom_function2(query_centroid, weighted_pos_samples) + sum_neg_sim * sum_neg_dao2)
      )
    }
    
    # ------ 映射到每个 spot 的 Q/W ------
    cluster_labels <- kmeans_result_list[[sample_name]]$cluster
    if (min(cluster_labels) == 0) {
      cluster_labels <- cluster_labels + 1  # 若从0开始则加1
    }
    
    Q_spot <- Q_s[cluster_labels, , drop = FALSE]
    W_spot <- W_s[cluster_labels, , drop = FALSE]
    
    Q_list[[sample_name]] <- Q_spot
    W_list[[sample_name]] <- W_spot
  }
  
  return(list(Q_list = Q_list, W_list = W_list))
}



multiNMF_graph_v2 <- function(B,Y_list, A_list, D_list, rank, lambda1,lambda2,Kcu=7,tau=0.3,max_iter = 500, tol = 1e-5, verbose = TRUE) {
  S <- length(Y_list)
  G <- nrow(Y_list[[1]])
  P_list <- lapply(Y_list, function(Y_s) {
    N_s <- ncol(Y_s)
    set.seed(123)
    abs(matrix(runif(N_s * rank), N_s, rank))
  })
  infoNCE<-infoNCE(P_list,Kcu,tau)
  Q_list<-infoNCE$Q_list
  W_list<-infoNCE$W_list
  loss_history <- numeric(max_iter)
  for (iter in 1:max_iter) {
    
    # === 更新每个 P_s（图正则） ===
    for (s in 1:S) {
      Y_s <- Y_list[[s]]         # G x N_s
      P_s <- P_list[[s]]         # N_s x d
      A_s <-  A_final_list[[s]]         # N_s x N_s
      D_s <-  D_final_list[[s]]         # N_s x N_s
      Q_s <- Q_list[[s]]
      W_s <- W_list[[s]]
      numerator_P <- t(Y_s) %*% B + lambda1 * (A_s %*% P_s)+(1/2)*lambda2*W_s
      denominator_P <- P_s %*% (t(B) %*% B) + lambda1 * (D_s %*% P_s)+(1/2)*lambda2*Q_s
      
      P_s <- P_s * (numerator_P / (denominator_P + 1e-10))
      P_s[P_s < 1e-10] <- 1e-10
      P_list[[s]] <- P_s
    }
    # === 计算重构损失（可选） ===
    loss <- 0
    for (s in 1:S) {
      Y_s <- Y_list[[s]]
      P_s <- P_list[[s]]
      recon <- B %*% t(P_s)
      loss <- loss + sum((Y_s - recon)^2)
    }
    loss_history[iter] <- loss
    
    if (iter > 1) {
      rel_change <- abs(loss_history[iter] - loss_history[iter - 1]) / loss_history[iter - 1]
      if (rel_change < tol) {
        if (verbose) cat("Converged at iteration", iter, "with loss:", loss, "\n")
        loss_history <- loss_history[1:iter]
        break
      }
    }
  }
  
  return(list(
    #B = B,
    P_list = P_list,
    loss_history = loss_history,
    n_iter = iter
  ))
}


library(mclust)

lambda1_values <- 5280  # 从1000到5000，步长为100
lambda2_value <- 151  # 固定lambda2
rank_value <- 44
k <- length(unique(label_list[[1]]))
B <- celltype_matrix

# 存储结果
results <- data.frame(lambda1 = numeric(), avg_ari = numeric())

for (lambda1 in lambda1_values) {
  cat(sprintf("\nlambda1 = %d:\n", lambda1))
  
  result_graph <- multiNMF_graph_v2(B, Y_list, A_final_list, D_final_list, 
                                    rank = rank_value, 
                                    lambda1 = lambda1, 
                                    lambda2 = lambda2_value)
  
  P_matrices <- result_graph$P_list
  ari_values <- numeric()
  
  for (i in seq_along(P_matrices)) {
    sample_id <- names(Y_list)[i]
    true_labels <- label_list[[sample_id]]
    P <- P_matrices[[i]]
    
    set.seed(123)
    kmeans_result <- kmeans(P, centers = k, nstart = 25)
    ari_value <- adjustedRandIndex(kmeans_result$cluster, true_labels)
    ari_values[i] <- ari_value
    
    cat(sprintf("  %s: %.4f\n", sample_id, ari_value))
  }
  
  avg_ari <- mean(ari_values)
  cat(sprintf("Average: %.4f\n", avg_ari))
  
  results <- rbind(results, data.frame(lambda1 = lambda1, avg_ari = avg_ari))
}

# 输出最佳结果
best_idx <- which.max(results$avg_ari)
cat(sprintf("\nBest: lambda1 = %d, ARI = %.4f\n", 
            results$lambda1[best_idx], results$avg_ari[best_idx]))




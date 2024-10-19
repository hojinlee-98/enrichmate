##############
# enrichmate
# hojinlee
# 2024.10.19
##############

library(dplyr)
library(plyr)
library(fgsea)
library(tibble)
library(ggplot2)
library(GO.db)

# 1. Setting environment for analysis ------
## get database ----
gobp <- gmtPathways("c5.go.bp.v2023.2.Hs.symbols.gmt")
gobp_hj <- read.table("20240102_c5.go.bp.v2023.2.Hs.term.id_mapping_hj.txt",header = T, sep = "\t") # this file is made from json file

gobp_ancestor <- as.list(GOBPANCESTOR)
names(gobp) %>% length() # 7647
which(gobp_hj$GOID %in% names(gobp_ancestor)) %>% length() # 7636

library(GOxploreR)
df <- data.frame()
for (level_tmp in seq(1,18)) {
  goid_tmp <- GOxploreR::Level2GOTermBP(level = level_tmp) # Level2GOTermBP
  df_tmp <- data.frame(LEVEL = rep(level_tmp, length(goid_tmp)),
                       GOID = goid_tmp)
  df <- rbind(df, df_tmp)
}

goid_level_df <- df

## do not run this step ----
# this step is time consuming, so the results object was already saved as object.

# df <- data.frame()
# for (gobp_tmp in names(gobp_ancestor)) {
#   ancestor_tmp <- gobp_ancestor[[gobp_tmp]]
#   df_tmp <- data.frame(GOID = rep(gobp_tmp, length(ancestor_tmp)+1),
#                        GOBP_ANCESTOR = c(ancestor_tmp, gobp_tmp)) # ANCESTOR + TARGET GOID
#   df <- rbind(df, df_tmp)
# } 
# save(df, file = "20240122_GOBP_ANCESTOR_df_hj.rda")

## run here ----
load("20240122_GOBP_ANCESTOR_df_hj.rda")

head(df)
goid_level_df %>% head()
colnames(goid_level_df) <- c("LEVEL", "GOBP_ANCESTOR")
colnames(df) <- c("TARGET_GOID", "GOBP_ANCESTOR")
goid_ancestor_level_df <- merge(goid_level_df, df, by = "GOBP_ANCESTOR", all.y = T)
goid_ancestor_level_df <- goid_ancestor_level_df %>% dplyr::filter(LEVEL != "all")
goid_ancestor_level_df$LEVEL %>% is.na() %>% which() # there are no NA data. 

goid_ancestor_level_df # TARGET_GOID is base and GOBP_ANCESTOR is the ancestor of the base.
# and the level is for GOBP_ANCESTOR.

gobp_with_ancestor_info <- gobp_hj$GOID[which(gobp_hj$GOID %in% names(gobp_ancestor))] # this gobp is used for GSEA.
gobp_hj <- gobp_hj %>% dplyr::filter(GOID %in% gobp_with_ancestor_info)
gobp <- gobp[gobp_hj$GOTERM] # 7636 GOBPs

# 2. Get GOBP results ------
mydf <- read.table("20240122_cr_base_y2_gsea_results_hj.txt", header = T, sep = "\t")

# 3. clustering -------
## input GOIDs ----------
mydf <- merge(mydf, gobp_hj, by.x = "pathway", by.y = "GOTERM", all.x = T)
goids <- mydf$GOID

## compute IC score ------
library(GOSemSim)
library(simplifyEnrichment)
if (!("hsGO" %in% ls())) {
  set.seed(1234); hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC = T) # make IC score using GOsemsim package 
}

sim_method <- "Rel"
set.seed(1234); w <- GOSemSim::mgoSim(goids, goids, semData = hsGO, measure = sim_method, combine = NULL) # w is similarity matrix. 
na_goid <- names(which(is.na(w[1,]))) # select NA goid.
goids <- goids[which(!goids %in% na_goid)]
set.seed(1234); w <- GOSemSim::mgoSim(goids, goids, semData = hsGO, measure = sim_method, combine = NULL) # w is similarity

## setting for binary cut clustering ------
library(ggplot2)
cutoff <- seq(0.6, 0.9, by = 0.01)
cutoff = cutoff[cutoff >= 0.5 & cutoff <= 1]
s1 = s2 = s3 = NULL

min_cluster <- 2

for (i in seq_along(cutoff)) {
  set.seed(1234)
  cl = binary_cut(w, cutoff = cutoff[i])
  s1[i] = difference_score(w, cl)
  tb = table(cl)
  s2[i] = length(tb)
  s3[i] = sum(tb < min_cluster)
}

df1 = data.frame(cutoff, s2)
colnames(df1) = c("method", "value")
df2 = data.frame(cutoff, s3)
colnames(df2) = c("method", "value")
df1$type = "All sizes"
df2$type = paste("Size", "<", min_cluster)
df = rbind(df1, df2)

p <- df %>% ggplot(aes(x = method, y = value, color = type)) +
  geom_point() +
  theme_classic() +
  xlab("Cutoff") +
  ylab("# of clusters") +
  ggtitle("Outlier detection") +
  scale_color_discrete("Type of clusters") +
  theme(aspect.ratio = 0.5,
        plot.title = element_text(hjust = 0.5, size = 15))

print(p) # select outlier cutoff using this plot.

## select cutoff -----
cutoff <- 0.84 # change this cutoff using the above plot.

## binary clustering -----
set.seed(1234); cl <- binary_cut(mat = w, partition_fun = partition_by_pam, cutoff = cutoff, partial = F)
min_term_numb <- 1
set.seed(1234); cluster_plot <- ht_clusters(mat = w, cl = cl, draw_word_cloud = F, min_term = min_term_numb) # in this heatmap, small clusters will be considered as outliers.
bc_res <- data.frame(GOID = colnames(w), cluster = cl) 
outlier_clusters <- bc_res %>% dplyr::group_by(cluster) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n < min_cluster) %>% dplyr::pull(cluster)

## assign pseudo-name for outliers ----
# if you want to not assign outliers, change the outlier_detect to FALSE.
outlier_detect <- T
if (outlier_detect == T) {
  set.seed(1234); pseudoname_outliers <- sample(10000:20000, size = length(outlier_clusters), replace = F)
  pseudoname_outliers_df <- data.frame(pseudoname = pseudoname_outliers, cluster = outlier_clusters) 
  pseudoname_outliers_df <- merge(pseudoname_outliers_df, bc_res, by = "cluster", all.x = T)
} 

## hclust ward.D using the remaining ----
set.seed(1234); dmat <- stats::as.dist(1 - w) # distance matrix.
set.seed(1234); hc <- stats::hclust(dmat, method = "ward.D") # h.clust using ward.D.

## clustering results for various k_val ----
test_range = seq(1,30) # if you want to increase or decrease the number of clusters to be tested, change this range. 

hc_res_list <- list()
for (k_val in seq(1,30)) {
  
  set.seed(1234); hclust_wardD <- stats::cutree(hc, k = k_val)
  hc_res_tmp <- data.frame(k = rep(k_val, length(hclust_wardD)),
                           GOID = names(hclust_wardD),
                           cluster = hclust_wardD,
                           row.names = NULL)
  
  hc_res_tmp <- hc_res_tmp[hc$order, ] # ordering
  
  hc_res_list[[paste0("hclust_", k_val)]] <- hc_res_tmp
}

hc_res_df <- do.call(rbind, hc_res_list)
rownames(hc_res_df) <- NULL

## binary cut and hclust ward.D --------
# if you do not want to assign outliers, change outlier_detect.
outlier_detect <- T
if (outlier_detect == T) {
  pseudoname_outliers_df <- pseudoname_outliers_df %>% dplyr::select(-c("cluster")) %>% dplyr::rename("cluster" = pseudoname) 
}

## get GOTERM for GOID ----
library("GO.db")
go_term <- as.list(GOTERM)

goids <- colnames(w)
goids_term_list <- list()

for (i in 1:length(goids)) {
  goid_tmp <- goids[i]
  detector <- go_term[[goid_tmp]]
  
  if (!is.null(detector)) {
    goterm_tmp <- Term(go_term[[goid_tmp]])
    
    goids_term_list[[goid_tmp]] <- data.frame(GOID = goid_tmp, GOTERM = goterm_tmp)
  }
}

goids_term_df <- do.call(rbind, goids_term_list)

## GOID GOTERM df for each h.clust k-value  ----
final_res_df <- merge(hc_res_df, goids_term_df, by = "GOID", all.x = T)
goid_level_df <- GOxploreR::GOTermBPOnLevel(goterm = goids)
final_res_df <- merge(final_res_df, goid_level_df, by.x = "GOID", by.y = "Term")
final_res_df <- final_res_df %>% dplyr::arrange(cluster)

mydf_tmp <- mydf %>% dplyr::select(c("GOID", "NES")) # if the results of GSEA are used, selecting GOID and NES were recommended. 
final_res_df <- merge(final_res_df, mydf_tmp, by.x = "GOID", by.y = "GOID") %>% dplyr::arrange(cluster)
#write.table(final_res_df, "filename", col.names = T, sep = "\t", row.names = F, quote = F) 

# 4. representative terms for cluster (common ancestor) --------
#ignore_go_level <- 3 # set ignore GOBP levels.
# k_val <- seq(15, 20) 
k_val <- 21 

df <- data.frame()
ans_df <- data.frame()
first_rep_c_goid_df <- data.frame()
for (k_tmp in k_val) {
  final_res_df_tmp <- final_res_df %>% dplyr::filter(k == k_tmp)
  
  for (c_tmp in seq(1:max(final_res_df_tmp$cluster))) {
    #c_tmp <- 7
    final_res_df_c_tmp <- final_res_df_tmp %>% dplyr::filter(cluster == c_tmp)
    ancestor_df_tmp <- goid_ancestor_level_df %>% dplyr::filter(TARGET_GOID %in% final_res_df_c_tmp$GOID)
    ancestor_df_tmp$k <- k_tmp
    ancestor_df_tmp$cluster <- c_tmp
    ans_df <- rbind(ans_df, ancestor_df_tmp)
    
    seed_c_tmp <- final_res_df_c_tmp
    remained_goid <- final_res_df_c_tmp$GOID
    
    while (TRUE) {
      
      explained_go <- floor(nrow(seed_c_tmp)/2)
      ancestor_sub_df_tmp <- goid_ancestor_level_df %>% dplyr::filter(TARGET_GOID %in% seed_c_tmp$GOID)
      
      ancestor_sub_df_filtered <- ancestor_sub_df_tmp %>%
        dplyr::group_by(GOBP_ANCESTOR, LEVEL) %>%
        dplyr::summarise(n=n()) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::filter(n > explained_go)
      
      if (nrow(ancestor_sub_df_filtered) != 0) {
        
        max_ancestor_tmp <- ancestor_sub_df_filtered %>%
          dplyr::ungroup(GOBP_ANCESTOR, LEVEL) %>%
          dplyr::filter(LEVEL == max(LEVEL)) %>%
          dplyr::filter(n == max(n))
        
        #print(max_ancestor_tmp)
        
        first_rep_c_goid <- ancestor_sub_df_tmp %>%
          dplyr::filter(GOBP_ANCESTOR %in% max_ancestor_tmp$GOBP_ANCESTOR) %>%
          dplyr::pull(TARGET_GOID) %>% unique()
        
        # get the common ancestor terms and child term 
        first_rep_c_goid_df_tmp <- ancestor_sub_df_tmp %>% dplyr::filter(TARGET_GOID %in% first_rep_c_goid &
                                                                           GOBP_ANCESTOR %in% max_ancestor_tmp$GOBP_ANCESTOR)
        first_rep_c_goid_df_tmp$cluster <- c_tmp
        first_rep_c_goid_df <- rbind(first_rep_c_goid_df, first_rep_c_goid_df_tmp) # ancestor and child dataframe
        
        remained_goid <- seed_c_tmp %>%
          dplyr::filter(!GOID %in% first_rep_c_goid) %>%
          dplyr::pull(GOID)
        
        seed_c_tmp <- seed_c_tmp %>% dplyr::filter(GOID %in% remained_goid)
        #print(seed_c_tmp)
        
        df_tmp <- data.frame(k = rep(k_tmp, nrow(max_ancestor_tmp)),
                             cluster = rep(c_tmp, nrow(max_ancestor_tmp)),
                             max_level = rep(unique(max_ancestor_tmp$LEVEL), nrow(max_ancestor_tmp)),
                             n = rep(unique(max_ancestor_tmp$n), nrow(max_ancestor_tmp)),
                             total_gobp = rep(nrow(final_res_df_c_tmp), nrow(max_ancestor_tmp)),
                             n_total_gobp = rep(paste(c(unique(max_ancestor_tmp$n), nrow(final_res_df_c_tmp)), collapse = "/"), nrow(max_ancestor_tmp)),
                             GOBP_ANCESTOR = max_ancestor_tmp$GOBP_ANCESTOR)
        
        df <- rbind(df, df_tmp)
        #print(1)
        #print(df_tmp)
        
        if (length(remained_goid) == 1) {
          df_tmp <- data.frame(k = rep(k_tmp, nrow(seed_c_tmp)),
                               cluster = rep(c_tmp, nrow(seed_c_tmp)),
                               max_level = rep(".", nrow(seed_c_tmp)),
                               n = rep(1, nrow(max_ancestor_tmp)),
                               total_gobp = rep(nrow(final_res_df_c_tmp), nrow(seed_c_tmp)),
                               n_total_gobp = rep(paste(c(1, nrow(final_res_df_c_tmp)), collapse = "/"), nrow(seed_c_tmp)),
                               GOBP_ANCESTOR = seed_c_tmp$GOID)
          df <- rbind(df, df_tmp)
          #print(2)
          #print(df_tmp)
          break
        } else if (length(remained_goid) == 0) {
          break
        }
        
      } else {
        
        df_tmp <- data.frame(k = rep(k_tmp, nrow(seed_c_tmp)),
                             cluster = rep(c_tmp, nrow(seed_c_tmp)),
                             max_level = rep(".", nrow(seed_c_tmp)),
                             n = rep(1, nrow(max_ancestor_tmp)),
                             total_gobp = rep(nrow(final_res_df_c_tmp), nrow(seed_c_tmp)),
                             n_total_gobp = rep(paste(c(1, nrow(final_res_df_c_tmp)), collapse = "/"), nrow(seed_c_tmp)),
                             GOBP_ANCESTOR = seed_c_tmp$GOID)
        df <- rbind(df, df_tmp)
        #print(3)
        #print(df_tmp)
        break
      }
      
    }
    
    # df <- rbind(df, df_tmp)
  }
}


rep_df <- df

write.table(ans_df, "filename", col.names = T, row.names = F, sep = "\t", quote = F) # all ancestors for each GO terms are included. (tmp file)
write.table(rep_df, "filename", col.names = T, row.names = F, sep = "\t", quote = F) # representative terms are included. (tmp file)
write.table(first_rep_c_goid_df, "filename", col.names = T, row.names = F, sep = "\t", quote = F) # representative terms and the GO terms that matched with the representative terms are included. (tmp file)

# 5. visualization --------
library(ComplexHeatmap)
library(colorRamp2)
col <- my_palette <- colorRampPalette(c("white","red"))(n = 300)
col_fun <- colorRamp2(seq(0, quantile(w[w > 0], 0.975), length = length(c("white","red"))), c("white","red"))
set.seed(1234); dmat <- stats::as.dist(1 - w) # distance matrix.
set.seed(1234); hc <- stats::hclust(dmat, method = "ward.D") # h.clust using ward.D.
hclust_wardD <- stats::cutree(hc, k = k_val)

hc_cluster_df <- data.frame(GOID = names(hclust_wardD), cluster = hclust_wardD, row.names = NULL)
hc_cluster_df <- hc_cluster_df[hc$order, ] # ordering 
rownames(hc_cluster_df) <- NULL

# if you do not want to assign outliers, change outlier_detect
outlier_detect <- T
if (outlier_detect == T) {
  hc_cluster_df <- hc_cluster_df %>% dplyr::filter(!(GOID %in% pseudoname_outliers_df$GOID))
  hc_cluster_df$order <- as.numeric(rownames(hc_cluster_df))
  bc_cluster_df <- pseudoname_outliers_df %>% dplyr::arrange(cluster) %>% dplyr::select(c("GOID"))
  bc_cluster_df <- pseudoname_outliers_df %>% dplyr::arrange(cluster)
  
  bc_cluster_tb <- table(pseudoname_outliers_df$cluster)
  bc_cluster_numb <- length(bc_cluster_tb)
  bc_cluster_orig_name <- as.numeric(names(bc_cluster_tb))
  
  max_cl <- max(hc_cluster_df$cluster)
  cl_mold <- data.frame(new_cluster = seq(max_cl+1, max_cl+bc_cluster_numb),
                        cluster = bc_cluster_orig_name)
  
  bc_cluster_df <- merge(pseudoname_outliers_df, cl_mold, by = "cluster", all.x = T) %>% dplyr::select(c("GOID", "new_cluster")) %>% dplyr::rename("cluster" = new_cluster)
  bc_cluster_df$order <- max(as.numeric(rownames(hc_cluster_df))) + as.numeric(rownames(bc_cluster_df))
  
  final_order <- rbind(hc_cluster_df, bc_cluster_df)
  rownames(final_order) <- final_order$GOID
  
} else if (outlier_detect == F) {
  hc_cluster_df$order <- as.numeric(rownames(hc_cluster_df))
  
  max_cl <- max(hc_cluster_df$cluster)
  final_order <- hc_cluster_df
  rownames(final_order) <- final_order$GOID
}

final_order <- merge(final_order, mydf_tmp, by.x = "GOID", by.y = "GOID") %>% dplyr::arrange(order)
total_clusters = final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(n=n()) %>% dplyr::select(cluster) %>% table() %>% names()
orig_order <- data.frame(GOID = colnames(w), orig_order = seq(1, nrow(w)), row.names = colnames(w))
orig_order <- orig_order[final_order$GOID, ]

# frequency for clusters
cl_vec <- c()
freq_cl_vec <- c()
for (anode in unique(final_order$cluster)) {
  cl_tmp <- as.numeric(table(final_order$cluster)[as.character(anode)])
  cl_vec <- append(cl_vec, cl_tmp)
  term_sum <- length(final_order$cluster)
  freq_cl_tmp <- cl_tmp/term_sum
  freq_cl_vec <- append(freq_cl_vec, freq_cl_tmp)
}

library(ComplexHeatmap)
pdf("filename", width = 20, height = 20)
Heatmap(mat = w, col = col_fun, name = "Similarity",
        row_order = orig_order$orig_order,
        column_order = orig_order$orig_order,
        border = "#404040", row_title = NULL,
        show_column_dend = F,
        show_row_dend = F,
        show_row_names = F,
        show_column_names = F,
        row_dend_width = unit(2, "cm"),
        width = unit(20, "cm"),
        height = unit(20, "cm"),
        left_annotation = rowAnnotation(ggplot = anno_empty(height = unit(20, "cm"), width = unit(1, "cm")), show_annotation_name = FALSE, gp = gpar(col = "white")))

freq_cl_cumsum <- cumsum(freq_cl_vec)
freq_cl_cumsum <- freq_cl_cumsum[-length(freq_cl_cumsum)]

gap_filler <- 0.0001
decorate_heatmap_body("Similarity", {
  for (freq_cl in freq_cl_cumsum) {
    grid.lines(c(freq_cl + gap_filler, freq_cl + gap_filler), c(0, 1), gp = gpar(lty = 1, lwd = 1.5))
  }
})

freq_cl_cumsum_rev <- cumsum(rev(freq_cl_vec))
freq_cl_cumsum_rev <- freq_cl_cumsum_rev[-length(freq_cl_cumsum_rev)]
decorate_heatmap_body("Similarity", {
  for (freq_cl in freq_cl_cumsum_rev) {
    grid.lines(c(1, 0), c(freq_cl + gap_filler, freq_cl + gap_filler), gp = gpar(lty = 1, lwd = 1.5))
  }
})
coor_df <- data.frame(x1 = 0, x2 = 5,
                      y1 = c(0,freq_cl_cumsum_rev), y2 = c(freq_cl_cumsum_rev, 1))

final_order$cluster %>% table()

p <- ggplot() + 
  geom_rect(data=coor_df,
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=as.factor(y1)),
            color= "black",
            alpha=1) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_void() +
  theme(panel.background = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.background=element_blank()) + 
  labs(x=NULL, y=NULL) +
  geom_text(data = coor_df, # only 17 clusters are annotated. 
            aes(x=(x1+x2)/2, y=y1+((y2-y1)/2), label = rev(c(paste0("c", total_clusters[1:k_val]),
                                                             rep("", length(total_clusters[(k_val+1):length(total_clusters)])))
                                                           ), size=3))

decorate_annotation("ggplot", {
  vp = current.viewport()$name
  print(p, vp = vp)
})

dev.off()

# 6. Summary results -----
final_order <- merge(final_order, gobp_hj, by = "GOID", all.x = T) %>% dplyr::arrange(order)
clusters_orig <- final_order$cluster[(!duplicated(final_order$cluster))]
c_df <- data.frame(clusters_orig = clusters_orig, clustsers_anno = seq(1, length(clusters_orig)))
final_order <- merge(final_order, c_df, by.x = "cluster", by.y = "clusters_orig", all.x = T)
c_nes_mean <- final_order %>% dplyr::group_by(cluster) %>% dplyr::summarise(cluster_nes = mean(abs(NES)))
final_order <- merge(final_order, c_nes_mean, by = "cluster", all.x = T)
final_order <- final_order %>% dplyr::arrange(desc(cluster_nes), desc(NES))
final_order$GOTERM <- gsub(x = final_order$GOTERM, pattern = "GOBP_", replacement = "") %>%
  gsub(x = ., pattern = "_", replacement = " ") %>% tolower()

# total GO terms with the information
write.table(final_order, "filename", col.names = T, row.names = F, sep = "\t", quote = F) # this dataframe includes GO terms and cluster information.

final_order <- merge(final_order, rep_df, by.x = "cluster", by.y = "cluster", all.x = T)
final_order <- final_order %>% dplyr::arrange(desc(cluster_nes), desc(NES))
rep_results_anno <- final_order %>% dplyr::select(c("clustsers_anno", "cluster_nes", "max_level", "GOBP_ANCESTOR", "n_total_gobp"))
rep_results_anno <- rep_results_anno[!duplicated(rep_results_anno[c("clustsers_anno", "cluster_nes", "max_level", "GOBP_ANCESTOR", "n_total_gobp")]),]
rep_results_anno <- rep_results_anno %>% dplyr::filter(!is.na(GOBP_ANCESTOR)) # filtering out outliers
rep_results_anno_final <- data.frame()
for (i in seq(1:nrow(rep_results_anno))) {
  rep_results_anno_tmp <- rep_results_anno[i,]
  goterm_matched <- GOTERM[[rep_results_anno_tmp$GOBP_ANCESTOR]]
  GOTERM_tmp <- goterm_matched@Term
  rep_results_anno_tmp$GOBP_ANCESTOR_GOTERM <- GOTERM_tmp
  rep_results_anno_final <- rbind(rep_results_anno_final, rep_results_anno_tmp)
}

# final results
write.table(rep_results_anno_final, "filename", col.names = T, row.names = F, sep = "\t", quote = F) # this dataframe is the final results. put the table with final heatmap.

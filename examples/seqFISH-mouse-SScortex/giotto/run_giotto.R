# http://spatialgiotto.rc.fas.harvard.edu/giotto.seqfish.html

library(Giotto)

setwd('/home/cang/Dropbox/Projects_UCI/DeepLearningSpatialTranscriptomics/dev/package/examples/seqFISH-mouse-SScortex/giotto')
my_working_dir = "./data"
python_path = "/home/cang/anaconda3/envs/stdl_env/bin/python3"

getSpatialDataset(dataset = 'seqfish_SS_cortex', directory = my_working_dir, method = 'wget')



# 1. (optional) set Giotto instructions
instrs = createGiottoInstructions(save_plot = TRUE,show_plot = FALSE,save_dir = my_working_dir,python_path = python_path)
# 2. create giotto object from provided paths ####
expr_path = fs::path(my_working_dir, "cortex_svz_expression.txt")
loc_path = fs::path(my_working_dir, "cortex_svz_centroids_coord.txt")
meta_path = fs::path(my_working_dir, "cortex_svz_centroids_annot.txt")
# 3. This dataset contains multiple field of views which need to be stitched together
## first merge location and additional metadata
SS_locations = data.table::fread(loc_path)
cortex_fields = data.table::fread(meta_path)
SS_loc_annot = data.table::merge.data.table(SS_locations, cortex_fields, by = 'ID')
SS_loc_annot[, ID := factor(ID, levels = paste0('cell_',1:913))]
data.table::setorder(SS_loc_annot, ID)
## create file with offset information
my_offset_file = data.table::data.table(field = c(0, 1, 2, 3, 4, 5, 6), x_offset = c(0, 1654.97, 1750.75, 1674.35, 675.5, 2048, 675), y_offset = c(0, 0, 0, 0, -1438.02, -1438.02, 0))
## create a stitch file
stitch_file = stitchFieldCoordinates(location_file = SS_loc_annot,offset_file = my_offset_file,cumulate_offset_x = T,cumulate_offset_y = F,field_col = 'FOV', reverse_final_x = F, reverse_final_y = T)
stitch_file = stitch_file[,.(ID, X_final, Y_final)]
my_offset_file = my_offset_file[,.(field, x_offset_final, y_offset_final)]



SS_seqfish <- createGiottoObject(raw_exprs = expr_path,spatial_locs = stitch_file,offset_file = my_offset_file, instructions = instrs)
SS_seqfish = addCellMetadata(SS_seqfish,new_metadata = cortex_fields,by_column = T,column_cell_ID = 'ID')
cell_metadata = pDataDT(SS_seqfish)
cortex_cell_ids = cell_metadata[FOV %in% 0:4]$cell_ID
SS_seqfish = subsetGiotto(SS_seqfish, cell_ids = cortex_cell_ids)
SS_seqfish <- filterGiotto(gobject = SS_seqfish,expression_threshold = 1,gene_det_in_min_cells = 10,min_det_genes_per_cell = 10, expression_values = c('raw'),verbose = T)
## normalize
SS_seqfish <- normalizeGiotto(gobject = SS_seqfish, scalefactor = 6000, verbose = T)
## add gene & cell statistics
SS_seqfish <- addStatistics(gobject = SS_seqfish)
## adjust expression matrix for technical or known variables
SS_seqfish <- adjustGiottoMatrix(gobject = SS_seqfish, expression_values = c('normalized'),batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),return_gobject = TRUE,update_slot = c('custom'))




SS_seqfish <- calculateHVG(gobject = SS_seqfish, method = 'cov_loess', difference_in_cov = 0.1, save_param = list(save_name = '3_a_HVGplot', base_height = 5, base_width = 5))
## select genes based on HVG and gene statistics, both found in gene metadata
gene_metadata = fDataDT(SS_seqfish)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 4 & mean_expr_det > 0.5]$gene_ID




# runPCA: normal and recommended usage, set center = T:
# SS_seqfish <- runPCA(gobject = SS_seqfish, genes_to_use = featgenes, scale_unit = F, center = T)
# center=F for compatibility reason (with paper and previous Giotto version):
SS_seqfish <- runPCA(gobject = SS_seqfish, genes_to_use = featgenes, scale_unit = F, center = F)
screePlot(SS_seqfish, save_param = list(save_name = '3_b_screeplot'))
plotPCA(gobject = SS_seqfish,save_param = list(save_name = '3_c_PCA_reduction'))
## run UMAP and tSNE on PCA space (default)
SS_seqfish <- runUMAP(SS_seqfish, dimensions_to_use = 1:15, n_threads = 10)
plotUMAP(gobject = SS_seqfish,save_param = list(save_name = '3_d_UMAP_reduction'))
SS_seqfish <- runtSNE(SS_seqfish, dimensions_to_use = 1:15)
plotTSNE(gobject = SS_seqfish,save_param = list(save_name = '3_e_tSNE_reduction'))




SS_seqfish <- createNearestNetwork(gobject = SS_seqfish, dimensions_to_use = 1:15, k = 15)
SS_seqfish <- doLeidenCluster(gobject = SS_seqfish, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = SS_seqfish,cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,save_param = list(save_name = '4_a_UMAP_leiden'))
## Leiden subclustering for specified clusters
SS_seqfish = doLeidenSubCluster(gobject = SS_seqfish, cluster_column = 'leiden_clus',resolution = 0.2, k_neighbors = 10,hvg_param = list(method = 'cov_loess', difference_in_cov = 0.1),pca_param = list(expression_values = 'normalized', scale_unit = F),nn_param = list(dimensions_to_use = 1:5), selected_clusters = c(5, 6, 7),name = 'sub_leiden_clus_select')
## set colors for clusters
subleiden_order = c( 1.1, 5.1, 5.2,  2.1, 3.1,4.1,  4.2, 4.3, 6.2, 6.1,7.1, 7.2, 9.1, 8.1)
subleiden_colors = Giotto:::getDistinctColors(length(subleiden_order)) 
names(subleiden_colors) = subleiden_order
plotUMAP(gobject = SS_seqfish,cell_color = 'sub_leiden_clus_select', cell_color_code = subleiden_colors,show_NN_network = T, point_size = 2.5, show_center_label = F, legend_text = 12, legend_symbol_size = 3,save_param = list(save_name = '4_b_UMAP_leiden_subcluster'))
showClusterHeatmap(gobject = SS_seqfish, cluster_column = 'sub_leiden_clus_select',save_param = list(save_name = '4_c_heatmap', units = 'cm'),row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9))
showClusterDendrogram(SS_seqfish, h = 0.5, rotate = T, cluster_column = 'sub_leiden_clus_select', save_param = list(save_name = '4_d_dendro', units = 'cm'))




spatDimPlot(gobject = SS_seqfish, cell_color = 'sub_leiden_clus_select', cell_color_code = subleiden_colors, dim_point_size = 2, spat_point_size = 2,save_param = list(save_name = '5_a_covis_leiden'))
# selected groups and provide new colors
groups_of_interest = c(6.1, 6.2, 7.1, 7.2)
group_colors = c('red', 'green', 'blue', 'purple'); names(group_colors) = groups_of_interest
spatDimPlot(gobject = SS_seqfish, cell_color = 'sub_leiden_clus_select', dim_point_size = 2, spat_point_size = 2,select_cell_groups = groups_of_interest, cell_color_code = group_colors,save_param = list(save_name = '5_b_covis_leiden_selected'))




gini_markers_subclusters = findMarkers_one_vs_all(gobject = SS_seqfish,method = 'gini', expression_values = 'normalized',cluster_column = 'sub_leiden_clus_select', min_genes = 20, min_expr_gini_score = 0.5, min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']
# violinplot
violinPlot(SS_seqfish, genes = unique(topgenes_gini$genes), cluster_column = 'sub_leiden_clus_select',strip_text = 8, strip_position = 'right', cluster_custom_order = unique(topgenes_gini$cluster),save_param = c(save_name = '6_a_violinplot_gini', base_width = 5, base_height = 10))
# cluster heatmap
topgenes_gini2 = gini_markers_subclusters[, head(.SD, 6), by = 'cluster']
plotMetaDataHeatmap(SS_seqfish, selected_genes = unique(topgenes_gini2$genes), custom_gene_order = unique(topgenes_gini2$genes),custom_cluster_order = unique(topgenes_gini2$cluster),metadata_cols = c('sub_leiden_clus_select'), x_text_size = 10, y_text_size = 10,save_param = c(save_name = '6_b_metaheatmap_gini'))




#create vector with names
clusters_cell_types_cortex = c('L6 eNeuron', 'L4 eNeuron', 'L2/3 eNeuron', 'L5 eNeuron', 'Lhx6 iNeuron', 'Adarb2 iNeuron', 'endothelial', 'mural','OPC','Olig','astrocytes', 'microglia')
names(clusters_cell_types_cortex) = c(1.1, 2.1, 3.1, 4.1,5.1, 5.2,6.1, 6.2, 7.1, 7.2,8.1, 9.1)
SS_seqfish = annotateGiotto(gobject = SS_seqfish, annotation_vector = clusters_cell_types_cortex,
                            cluster_column = 'sub_leiden_clus_select', name = 'cell_types')
# cell type order and colors
cell_type_order = c('L6 eNeuron', 'L5 eNeuron', 'L4 eNeuron', 'L2/3 eNeuron','astrocytes', 'Olig', 'OPC','Adarb2 iNeuron', 'Lhx6 iNeuron','endothelial', 'mural', 'microglia')
cell_type_colors = subleiden_colors
names(cell_type_colors) = clusters_cell_types_cortex[names(subleiden_colors)]
cell_type_colors = cell_type_colors[cell_type_order]
## violinplot
violinPlot(gobject = SS_seqfish, genes = unique(topgenes_gini$genes),strip_text = 7, strip_position = 'right', cluster_custom_order = cell_type_order,cluster_column = 'cell_types', color_violin = 'cluster',save_param = c(save_name = '7_a_violinplot', base_width = 5))
## co-visualization
spatDimPlot(gobject = SS_seqfish, cell_color = 'cell_types',dim_point_size = 2, spat_point_size = 2, dim_show_cluster_center = F, dim_show_center_label = T,save_param = c(save_name = '7_b_covisualization'))





# Output expression matrix, locations, and assigned cell types
write.csv(SS_seqfish@cell_metadata$cell_types, file="./data/celltype_annotation.csv")
write.csv(t(as.matrix(SS_seqfish@raw_exprs)), file="./data/raw_exprs.csv")
write.csv(SS_seqfish@gene_metadata$gene_ID, file="./data/gene_id.csv")
write.csv(SS_seqfish@spatial_locs, file="./data/spatial_locs.csv")

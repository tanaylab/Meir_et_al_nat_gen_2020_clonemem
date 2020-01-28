##
##
## expression_analysis.r
## Date: 18.01.2020
## This file summarizes main steps use for transcription analysis of single-cells and single-cell-derived clonal populations in Nat. gen. Meir et al. 2020
##
##


#load dependency packages and project functions
source("Meir_et_al_2020_nat_gen_functions.r")

#load M scRNA MARS-seq umi_tables (cmem_sumi_PAR: HCT116 WT cells, cmem_sumi_DKO: HCT116 DKO (DNMT1;DNMT3B) cells, cmem_sumi_WI38: WI38 cells)
cmem_load_singles_marsSeq()

#downolad and load 10x scRNA-seq umi_tables
if(!dir.exists("expression_data/10x_umi_tables/")) { dir.create(paste0(base_dir,"expression_data/10x_umi_tables/")) }
download.file(url = "http://www.wisdom.weizmann.ac.il/~zoharme/Meir_et_al_NatGen_2020/expression_data/10x_umi_tables/scRNA_H1299_10x_umiTable.txt", 
			  destfile = "expression_data/10x_umi_tables/scRNA_H1299_10x_umiTable.txt")
download.file(url = "http://www.wisdom.weizmann.ac.il/~zoharme/Meir_et_al_NatGen_2020/expression_data/10x_umi_tables/scRNA_HCT116_DKO_10x_umiTable.txt", 
			  destfile = "expression_data/10x_umi_tables/scRNA_HCT116_DKO_10x_umiTable.txt")
download.file(url = "http://www.wisdom.weizmann.ac.il/~zoharme/Meir_et_al_NatGen_2020/expression_data/10x_umi_tables/scRNA_HCT116_WT_10x_umiTable.txt", 
			  destfile = "expression_data/10x_umi_tables/scRNA_HCT116_WT_10x_umiTable.txt")
cmem_sumi_10x_PAR = read.table("expression_data/10x_umi_tables/scRNA_HCT116_WT_10x_umiTable.txt",sep="\t", check.names=FALSE)
cmem_sumi_10x_DKO = read.table("expression_data/10x_umi_tables/scRNA_HCT116_DKO_10x_umiTable.txt",sep="\t", check.names=FALSE)
cmem_sumi_10x_H1299 = read.table("expression_data/10x_umi_tables/scRNA_H1299_10x_umiTable.txt",sep="\t", check.names=FALSE)

#load clones bulk MARS-seq umi_tables
cmem_load_clones()

#load cellCycle_genes table (Sup. Table 1)
cc_genes_table = read.csv ("expression_data/Supplementary_Table1_ccGenes.csv",row.names = 1)





################################################################################################
######
###### part 1: Analysis of transient and clonally-stable transcription in the Luria-Delbruck scheme.
######
###### 		module 1.1 - identification of cell-cycle (cc) independent expression variation
###### 		module 1.2 - Longitudinal analysis of cc-independent programs
################################################################################################




################################################################################################
######
###### module 1.1: identification of cell-cycle (cc) independent expression variation
######
######	step 1.1.1 - Create cc-model for each cell line and generate cell-cell similarity graph, based on this cc-model alone.
######  step 1.1.2 - Randomize expression of each cell, according to K most adjacent cells to it, according to the cc-model.
######  step 1.1.3 - Find co-expression patterns that survived the randmoization process, and mark them as cc-independent.
######
################################################################################################


######	step 1.1.1_a - Create cc-model for each cell line and generate cell-cell similarity graph, based on this cc-model alone.
#
# cc-model of each cell-line generated with the MetaCell package (Baran et al. 2019) by:
# a) finding genes with correlation to strong cell cycle markers (such as MKi67, PTTG1, GINS2 or others).
# b) generate balanced cell-cell similarity graph and assign cells to Metacells based on expression of these genes only.
# *) A detailed example for this process attached in a separate .html vignette in this gitHub directory.

# we used MetaCell to compute expression foot-print (fp) of each gene in each metacell, and to generate 2D-projection of the expression manifold according to Baran et al. 2019.
hct116_wt_10x_metacells = read.table("expression_data/MetaCell_objects/HCT116_wt_10x_ccModel_mcCellAssignment.txt",sep="\t")
hct116_wt_10x_metacells_fp = read.table("expression_data/MetaCell_objects/HCT116_wt_10x_ccModel_mcFP.txt",sep="\t",check.names=FALSE)
hct116_wt_10x_metacells_2d = read.table("expression_data/MetaCell_objects/HCT116_wt_10x_ccModel_mc2D.txt",sep="\t")

h1299_10x_metacells = read.table("expression_data/MetaCell_objects/H1299_10x_ccModel_mcCellAssignment.txt",sep="\t")
h1299_10x_metacells_fp = read.table("expression_data/MetaCell_objects/H1299_10x_ccModel_mcFP.txt",sep="\t", check.names=FALSE)
h1299_10x_metacells_2d = read.table("expression_data/MetaCell_objects/H1299_10x_ccModel_mc2D.txt",sep="\t")

wi38_mars_metacells = read.table("expression_data/MetaCell_objects/WI38_scMARSseq_ccModel_mcCellAssignment.txt",sep="\t")
wi38_mars_metacells_fp = read.table("expression_data/MetaCell_objects/WI38_scMARSseq_ccModel_mcFP.txt",sep="\t", check.names=FALSE)
wi38_mars_metacells_2d = read.table("expression_data/MetaCell_objects/WI38_scMARSseq_ccModel_mc2D.txt",sep="\t")



######	step 1.1.1_b - plot the cc-model in each cell line
# 
# a) expression heatmap
# b) 2-D projection, colored by expression of specific gene
# c) scatter of total output of S- and M-phase genes.


### HCT116 cells
#
# expression heatmap
hct116_wt_10x_mc_assign_vector = hct116_wt_10x_metacells[,"metacell"]
names(hct116_wt_10x_mc_assign_vector) = rownames(hct116_wt_10x_metacells)
cmem_mcell_mc_plot_marks (metacell_assignment = hct116_wt_10x_mc_assign_vector,
						  metacell_fp = hct116_wt_10x_metacells_fp,
						  genes_to_plot = rownames(cc_genes_table)[!is.na(cc_genes_table[,"HCT116_plot"])][order(cc_genes_table[,"HCT116_plot"][!is.na(cc_genes_table[,"HCT116_plot"])])], 
						  cells_umi_tab = cmem_sumi_10x_PAR,				  
						  height = 440,
						  width  = 1139.25,
						  fig_fn = "./HCT116_ccMarkers_heatmap.png",
						  k_log=7,
						  text_cex = 2.85,
						  lateral_genes = NULL,
						  lateral_genes_color = "red",
                          plot_cells = TRUE,
						  ideal_umi = NULL,
						  plot_top_marg = c(0,13,5,20),
						  pminpmax_cells = FALSE,
						  heatmap_lwd=0.5,
						  alternate_side_text = TRUE,
						  min_quantile_cells = 0,
						  max_quantile_cells = 1
						  )						 
#
# 2-D projection, colored by expression of specific gene
cmem_plot_gene_fp_over_mc2d (gene_to_plot = "MKI67",
							 z_score_range = 2,
							 metacells_fp = hct116_wt_10x_metacells_fp,
							 metacells_2d = hct116_wt_10x_metacells_2d,
							 cells_2d = hct116_wt_10x_metacells,
 							 my_xlim = c(-50,500),
							 my_ylim = c(-50,450),
							 my_text_cex = 3.6,
							 plot_w = 450,
							 plot_h = 450,
							 cells_cex = 0.5,
							 cells_color = "grey",
							 metacells_cex = 8,
							 plot_margins = c(4,4,4,4),
							 fig_fn = "./HCT116_2d_cellCycle_model_MKI67_lfp.png")
#							 
# scatter of total output of S- and M-phase genes
hct116_s_genes = unlist(lapply(na.omit(rownames(cc_genes_table)[cc_genes_table[,"HCT116"]=="S"]),"[[",1))
hct116_m_genes = unlist(lapply(na.omit(rownames(cc_genes_table)[cc_genes_table[,"HCT116"]=="M"]),"[[",1))
cmem_plot_gene_modules_scatter (genes_x = hct116_s_genes,
								genes_y = hct116_m_genes,
								x_label = "S-phase (log)",
								y_label = "M-phase (log)",
								cells_umi_tab = cmem_sumi_10x_PAR,
								metacell_assignment = hct116_wt_10x_mc_assign_vector,
								metacells_fp = hct116_wt_10x_metacells_fp,
								lfp_color_gene = "MKI67",
								z_score_range = 2,
								my_xlim = c(3.5,7),
								my_ylim = c(3,8),
								log2_exp = TRUE,
								eps = 1,
								my_text_cex = 3.6,
								cells_color = "#808080",
							    plot_w = 450,
							    plot_h = 450,
								cex_cells = 0.5,
								cex_metacells = 8,
								plot_margins = c(4,4,4,4),
								fig_fn = "./HCT116_S_M_phases_scatter_coloredBy_MKi67.png"
								)


								
### H1299 cells
#
# expression heatmap
h1299_10x_mc_assign_vector = h1299_10x_metacells[,"metacell"]
names(h1299_10x_mc_assign_vector) = rownames(h1299_10x_metacells)
cmem_mcell_mc_plot_marks (metacell_assignment = h1299_10x_mc_assign_vector,
						  metacell_fp = h1299_10x_metacells_fp,
						  genes_to_plot = rownames(cc_genes_table)[!is.na(cc_genes_table[,"H1299_plot"])][order(cc_genes_table[,"H1299_plot"][!is.na(cc_genes_table[,"H1299_plot"])])], 
						  cells_umi_tab = cmem_sumi_10x_H1299,				  
						  height = 484,
						  width  = 1254.4,
						  fig_fn = "./H1299_ccMarkers_heatmap.png",
						  k_log=7,
						  text_cex = 3.135,
						  lateral_genes = NULL,
						  lateral_genes_color = "red",
                          plot_cells = TRUE,
						  ideal_umi = NULL,
						  plot_top_marg = c(0,13,5,20),
						  pminpmax_cells = FALSE,
						  heatmap_lwd=0.5,
						  alternate_side_text = TRUE,
						  min_quantile_cells = 0,
						  max_quantile_cells = 0.8
						  )						 
#
# 2-D projection, colored by expression of specific gene
cmem_plot_gene_fp_over_mc2d (gene_to_plot = "MKI67",
							 z_score_range = 2,
							 metacells_fp = h1299_10x_metacells_fp,
							 metacells_2d = h1299_10x_metacells_2d,
							 cells_2d = h1299_10x_metacells,
 							 my_xlim = c(0,460),
							 my_ylim = c(0,430),
							 my_text_cex = 3.6,
							 plot_w = 450,
							 plot_h = 450,
							 cells_cex = 0.5,
							 cells_color = "grey",
							 metacells_cex = 8,
							 plot_margins = c(4,4,4,4),
							 fig_fn = "./H1299_2d_cellCycle_model_MKI67_lfp.png")
#							 
# scatter of total output of S- and M-phase genes
h1299_s_genes = unlist(lapply(na.omit(rownames(cc_genes_table)[cc_genes_table[,"H1299"]=="S"]),"[[",1))
h1299_m_genes = unlist(lapply(na.omit(rownames(cc_genes_table)[cc_genes_table[,"H1299"]=="M"]),"[[",1))
cmem_plot_gene_modules_scatter (genes_x = h1299_s_genes,
								genes_y = h1299_m_genes,
								x_label = "S-phase (log)",
								y_label = "M-phase (log)",
								cells_umi_tab = cmem_sumi_10x_H1299,
								metacell_assignment = h1299_10x_mc_assign_vector,
								metacells_fp = h1299_10x_metacells_fp,
								lfp_color_gene = "MKI67",
								z_score_range = 2,
								my_xlim = c(3.5,7),
								my_ylim = c(2.5,7.5),
								log2_exp = TRUE,
								eps = 1,
								my_text_cex = 3.6,
								cells_color = "#808080",
							    plot_w = 450,
							    plot_h = 450,
								cex_cells = 0.5,
								cex_metacells = 8,
								plot_margins = c(4,4,4,4),
								fig_fn = "./H1299_S_M_phases_scatter_coloredBy_MKi67.png"
								)


### WI38 cells
#
# expression heatmap
wi38_mars_mc_assign_vector = wi38_mars_metacells[,"metacell"]
names(wi38_mars_mc_assign_vector) = rownames(wi38_mars_metacells)
cmem_mcell_mc_plot_marks (metacell_assignment = wi38_mars_mc_assign_vector,
						  metacell_fp = wi38_mars_metacells_fp,
						  genes_to_plot = rownames(cc_genes_table)[!is.na(cc_genes_table[,"WI38_plot"])][order(cc_genes_table[,"WI38_plot"][!is.na(cc_genes_table[,"WI38_plot"])])], 
						  cells_umi_tab = cmem_sumi_WI38,				  
						  height = 418,
						  width  = 1000,
						  fig_fn = "./WI38_ccMarkers_heatmap.png",
						  k_log=7,
						  text_cex = 2.7075,
						  lateral_genes = NULL,
						  lateral_genes_color = "red",
                          plot_cells = TRUE,
						  ideal_umi = NULL,
						  plot_top_marg = c(0,13,5,20),
						  pminpmax_cells = FALSE,
						  heatmap_lwd=0.5,
						  alternate_side_text = TRUE,
						  min_quantile_cells = 0,
						  max_quantile_cells = 0.8
						  )						 
#
# 2-D projection, colored by expression of specific gene
cmem_plot_gene_fp_over_mc2d (gene_to_plot = "MKI67",
							 z_score_range = 1.5,
							 metacells_fp = wi38_mars_metacells_fp,
							 metacells_2d = wi38_mars_metacells_2d,
							 cells_2d = wi38_mars_metacells,
 							 my_xlim = c(0,550),
							 my_ylim = c(0,550),
							 my_text_cex = 3.6,
							 plot_w = 450,
							 plot_h = 450,
							 cells_cex = 0.5,
							 cells_color = "grey",
							 metacells_cex = 8,
							 plot_margins = c(4,4,4,4),
							 fig_fn = "./WI38_2d_cellCycle_model_MKI67_lfp.png")
#							 
# scatter of total output of S- and M-phase genes
wi38_s_genes = unlist(lapply(na.omit(rownames(cc_genes_table)[cc_genes_table[,"WI38"]=="S"]),"[[",1))
wi38_m_genes = unlist(lapply(na.omit(rownames(cc_genes_table)[cc_genes_table[,"WI38"]=="M"]),"[[",1))
cmem_plot_gene_modules_scatter (genes_x = wi38_s_genes,
								genes_y = wi38_m_genes,
								x_label = "S-phase (log)",
								y_label = "M-phase (log)",
								cells_umi_tab = cmem_sumi_WI38,
								metacell_assignment = wi38_mars_mc_assign_vector,
								metacells_fp = wi38_mars_metacells_fp,
								lfp_color_gene = "MKI67",
								z_score_range = 1.5,
								my_xlim = c(3,8),
								my_ylim = c(3.5,7.5),
								log2_exp = TRUE,
								eps = 1,
								my_text_cex = 3.6,
								cells_color = "#808080",
							    plot_w = 450,
							    plot_h = 450,
								cex_cells = 0.5,
								cex_metacells = 8,
								plot_margins = c(4,4,4,4),
								fig_fn = "./WI38_S_M_phases_scatter_coloredBy_MKi67.png"
								)





##  step 1.1.2 - Randomize expression of each cell, according to K most adjacent cells to it, according to the cc-model.
#
#	a) Load cell-cell balanced KNN graph generated with MetaCell based on the cell-cycle model
#   b) Randomize expression of cells by their adjacent K cells based on this graph.
#	c) Define genes that maintain correlation with another gene despite the randomization process as cell-cycle independent.
#	*) Adjust correlation of cell-cycle independent genes by sampling depth.

#HCT116
#
# load graph
hct116_cc_cgraph = read.table("expression_data/MetaCell_objects/HCT116_wt_10x_ccModel_cgraph.txt",sep="\t",check.names=FALSE,stringsAsFactors=TRUE)
#
# Randomize expression by CC model
hct116_randomized_geneCorrs = cmem_mcell_cgraph_norm_gcor (cgraph_edges = hct116_cc_cgraph,
														   cells_umi_tab = cmem_sumi_10x_PAR,
														   K = 20,
														   dowsamp_value = 6321,
														   min_gtot = 100)
#
# Define cc independent genes
hct116_cc_independet_genes = cmem_plot_ccNorm (cc_cors = hct116_randomized_geneCorrs,
														 fig_nm = "hct116_cc_maxCorrs.png",
														 x_scatter_lim=c(0,.5),
														 y_scatter_lim=c(0,.5),
														 min_cor = 0.1,
														 max_r_cor = 0.15,
														 min_dif = 0.04)	


														   
#NCI-H1299
#
# load graph
h1299_cc_cgraph = read.table("expression_data/MetaCell_objects/H1299_10x_ccModel_cgraph.txt",sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
# Randomize expression by CC model
h1299_randomized_geneCorrs = cmem_mcell_cgraph_norm_gcor (cgraph_edges = h1299_cc_cgraph,
														   cells_umi_tab = cmem_sumi_10x_H1299,
														   K = 20,
														   dowsamp_value = 7373,
														   min_gtot=50)
#
# Define cc independent genes
h1299_cc_independet_genes = cmem_plot_ccNorm (cc_cors = h1299_randomized_geneCorrs,
														 fig_nm = "h1299_cc_maxCorrs.png",
														 x_scatter_lim=c(0,.5),
														 y_scatter_lim=c(0,.5),
														 min_cor = 0.16,
														 max_r_cor = 0.16,
														 min_dif = 0.1)

#WI38
#
# load graph
wi38_cc_cgraph = read.table("expression_data/MetaCell_objects/WI38_scMARSseq_ccModel_cgraph.txt",sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
# Randomize expression by CC model
wi38_randomized_geneCorrs = cmem_mcell_cgraph_norm_gcor (cgraph_edges = wi38_cc_cgraph,
														   cells_umi_tab = cmem_sumi_WI38,
														   K = 20,
														   dowsamp_value = 2844,
														   min_gtot = 50)

wi38_cc_independet_genes = cmem_plot_ccNorm (cc_cors = wi38_randomized_geneCorrs,
														 fig_nm = "wi38_cc_maxCorrs.png",
														 x_scatter_lim=c(0,.5),
														 y_scatter_lim=c(0,.5),
														 min_cor = 0.12,
														 max_r_cor = 0.14,
														 min_dif = 0.02)															   
														   
														   
														   


# *) Adjusting correlations by sampling depth, and plotting one Example										
#
#	a) adjust correaltions - by subtracting the sampling-depth of both genes from the correlation between them.
#		(so we will find surprising correlations between genes, given their overall coverage in a sparse scRNA-seq dataset)
#	b) find co-expressed partners of cc-independent genes


### HCT116
#
# adjust correlation and plot examples of some genes
gcor_hct116 = cmem_adjust_gcor_by_samplingDepth (cells_umi_tab = cmem_sumi_PAR,
											     min_gtot = 50,
											     min_cell = 1000,
											     max_cell = 40000,
											     rollmean_n_genes = 101,
											     cor_method_spearman = TRUE,
											     symmertrize_adj_corr = TRUE,
											     plot_example = TRUE,
											     gene_to_plot = c("EPCAM","VIM","IDI1"),
											     base_dir = "./",
											     cex_text = 1.8,
											     fig_w = 400,
											     fig_h = 400,
											     x_text_label_factor = 0.018,
											     y_text_label_factor = 0.015,
											     low_corr_t = (-0.16),
											     high_corr_t = 0.17,
											     trend_col = "indianred4",
											     positive_corr_col = "blue",
											     negative_corr_col = "red",
											     ylim_example1 = c(-0.15,0.4),
											     ylim_example2 = c(-0.3,0.23),
											     xlim_example1 = c(7.5,13.5)
											     )							
#												 
# plot the most correlated genes to EpCAM, over shuffled control
gcor_hct116_perm_EpCAM = cmem_permute_gene_adjust_gcor (cells_umi_tab = cmem_sumi_PAR,
														min_gtot = 50,
														min_cell = 1000,
														max_cell = 40000,
														rollmean_n_genes = 101,
														cor_method_spearman = TRUE,
														symmertrize_adj_corr = FALSE,
														gene_to_permute = "EPCAM")
#							   
cmem_barplot_top_gcor (adj_gcor_mat = gcor_hct116$gcor_n,
					   perm_adj_gcor_mat = gcor_hct116_perm_EpCAM$perm_gcor_n,
					   gene_name = "EPCAM",
					   n_genes = 10,
					   positive_corr_col = "blue",
					   negative_corr_col = "red",
					   positive_permuted_corr_col = "lightblue",
					   negative_permuted_corr_col = "tomato",
					   positive_corr_margins = c(14,4.5,7,4),
					   negative_corr_margins = c(7,4.5,14,4),
					   ylim_positive_corr = c(-0.075,0.25),
					   ylim_negative_corr = c(-0.25,0.075),
					   ylim_density = c(0,15),
					   xlim_ecdf_positive = c(0.075,0.23),
					   ylim_ecdf_positive = c(0,100),
					   xlim_ecdf_negative = c(-0.25,-0.075),
					   ylim_ecdf_negative = c(0,100),
					   text_cex=1.2,
					   barplot_h = 400,
					   barplot_w = 800,
					   density_plot_w = 250,
					   density_plot_h = 250,
					   ecdf_plot_w = 220,
					   ecdf_plot_h = 220,
					   base_dir = "./"
					   )


								   
								   
								   



														   
######  step 1.1.3 - compare expression distributions of cell-cycle associated (transient), and cell-cycle independent (testing if stable..) programs in single cells and in clones. 
#
#
# load cc independent gene modules
#
cc_independent_geneModules = read.csv("expression_data/Supplementary_Table3_ccIndependet_geneModules.csv",row.names="gene")

### HCT116 
#
# compare cells and clones total output of programs
cmem_compare_cells_clones (mod_to_plot = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$HCT116_WT=="Epithelial"]),"[[",1)),
						mod_nm = "EpCAM module",
						sc_mat = cmem_sumi_PAR,
						c_mat = cmem_umi_PAR,
						minimal_cell = 6e3,
						minimal_clone = 1e4,
						base_dir = "./",
						legend_x = 270,
						legend_y = 0.0058,
						density_x_lim = c(-50,650),
						density_y_lim = c(0,0.00525),
						ecdf_x_lim = c(-50,650),
						plot_ecdf=FALSE,
						density_mar = c(4.5,4,0.5,0.5),
						show_legend_density=TRUE)
#						
cmem_compare_cells_clones (mod_to_plot = hct116_s_genes,
						mod_nm = "S module",
						sc_mat = cmem_sumi_PAR,
						c_mat = cmem_umi_PAR,
						minimal_cell = 6e3,
						minimal_clone = 1e4,
						base_dir = "./",
						density_x_lim = c(-100,3000),
						density_y_lim = c(0,0.004),
						ecdf_x_lim = c(-100,3000),
						plot_ecdf=FALSE,
						density_mar = c(4.5,4,0.5,0.5),
						show_legend_density=FALSE)

# plot distinct VIM-high subpopulation of cells and clones
cmem_hist_gene_cells_clones (gene_to_plot = "VIM",
							 min_cell = 6e3,
							 max_cell = 4e4,
							 min_clone = 1e4,
							 cells_umi_table = cmem_sumi_PAR,
							 clones_umi_table = cmem_umi_PAR,
							 cells_umi_table_norm = cmem_sumi_PAR_norm,
							 clones_umi_table_norm = cmem_umi_PAR_norm,
							 log2_exp = TRUE,
							 eps = 1,
							 cex_text=1.8,
							 base_dir = "./",
							 x_label = "VIM (log)",
							 plot_h = 273,
							 plot_w = 254.8,
							 plot_margins = c(5.5,4.6,1,2),
							 plot_legend_cells = TRUE,
							 plot_legend_clones = TRUE,
							 population_legend = c("normal","VIM-high"),
							 legend_x_cells = 1.6,
							 legend_y_cells = 400,
							 legend_x_clones = 1.6,
							 legend_y_clones = 150,
							 hist_breaks = seq(0,6,0.75),
							 hist_colors = c(rep("snow2",5),rep("indianred4",17)))

							 
### H1299
#
# compare cells and clones total output of programs
cmem_compare_cells_clones (mod_to_plot =  unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$H1299=="ID"]),"[[",1)),
						mod_nm = "ID module",
						sc_mat = cmem_sumi_10x_H1299,
						c_mat = cmem_umi_H1299,
						minimal_cell = 6e3,
						minimal_clone = 1e4,
						base_dir = "./",
						legend_x = 475,
						legend_y = 0.0036,
						density_x_lim = c(-150,1200),
						density_y_lim = c(0,0.00325),
						ecdf_x_lim = c(-150,1200),
						plot_ecdf=FALSE,
						density_mar = c(4.5,4,0.5,0.5),
						show_legend_density=TRUE)
#
cmem_compare_cells_clones (mod_to_plot = h1299_m_genes,
						mod_nm = "M module",
						sc_mat = cmem_sumi_10x_H1299,
						c_mat = cmem_umi_H1299,
						minimal_cell = 1e3,
						minimal_clone = 3e3,
						base_dir = ".",
						legend_x = 900,
						legend_y = 0.45,
						density_x_lim = c(-200,4000),
						density_y_lim = c(0,0.003),
						ecdf_x_lim = c(-200,3000),
						plot_ecdf=FALSE,
						density_mar=c(4.5,3,0.5,1.5))					

### WI38 
#
# compare cells and clones total output of programs
cmem_compare_cells_clones (mod_to_plot = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$WI38=="Collagen"]),"[[",1)),
						mod_nm = "Collagen module",
						sc_mat = cmem_ds(cmem_sumi_WI38,3e3),
						c_mat = cmem_ds(cmem_umi_WI38,3e4),
						minimal_cell = 1e3,
						minimal_clone = 1e4,
						base_dir = "./",
						legend_x = 1875,
						legend_y = 0.00125,
						density_x_lim = c(0,4000),
						density_y_lim = c(0,0.00115),
						ecdf_x_lim = c(0,4000),
						plot_ecdf=FALSE,
						density_mar = c(4.5,4,0.5,0.5),
						show_legend_density=TRUE)
#
cmem_compare_cells_clones (mod_to_plot = wi38_m_genes,
						mod_nm = "M module",
						sc_mat = cmem_ds(cmem_sumi_WI38,3e3),
						c_mat = cmem_ds(cmem_umi_WI38,3e4),
						minimal_cell = 1e3,
						minimal_clone = 3e4,
						base_dir = "./",
						legend_x = 640,
						legend_y = 0.004,
						density_x_lim = c(-250,3050),
						density_y_lim = c(0,0.0045),
						ecdf_x_lim = c(-250,3050),
						plot_ecdf=FALSE,
						density_mar = c(4,4,1,0.5),
						show_legend_density=FALSE)
						

						
################################################################################################
######
###### module 1.2: Longitudinal analysis of cc-independent programs
######
######	step 1.2.1 - map heterogeneity in bulk-RNA (days 10, 18) and scRNA (days 33, 62, 98 and 148) longitudinal tracking of single cell-derived HCT116 clones
######			
######  step 1.2.2 - Follow Epithelial and EMT-like clonal expression over time 
######
################################################################################################

# load long-term RNA objects
#
# day 10: low-depth bulk RNA sampling
# day 18: low-depth bulk RNA sampling
# day 33: scRNA seq
# day 62: scRNA seq
# day 98: scRNA seq
# day 148: scRNA seq
# *Exome data on these clones from two time-points (d=78, d=168) is also available 
cmem_load_longitudinal_HCT116 (min_cell=1e3,
							   min_clone_day10 = 3e3,
							   min_clone_day18 = 8e3
							   )
							   

######	step 1.2.1 - map heterogeneity in bulk-RNA (days 10, 18) and scRNA (days 33, 62, 98 and 148) longitudinal tracking of single cell-derived HCT116 clones
#
# we used MetaCell to compute expression foot-print (fp) of each gene in each metacell, and to generate 2D-projection of the expression manifold according to Baran et al. 2019.
scLongitudinal_metacells = read.csv("expression_data/MetaCell_objects/HCT116_WT_scMARSseq_longitudinal_mcCellAssignment.csv",row.names=1)
scLongitudinal_metacells_fp = read.csv("expression_data/MetaCell_objects/HCT116_WT_scMARSseq_longitudinal_mcFP.csv",row.names=1,check.names=FALSE)
scLongitudinal_metacells_2d = read.csv("expression_data/MetaCell_objects/HCT116_WT_scMARSseq_longitudinal_mc2D.csv",row.names=1,check.names=FALSE)

#
# a) plot heatmap of epithelial markers over longitudinal metacells
#
scLongitudinal_assign_vector = scLongitudinal_metacells[,"metacell"]
names(scLongitudinal_assign_vector) = rownames(scLongitudinal_metacells)
epithelial_genes = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$HCT116_WT=="Epithelial"]),"[[",1))
emt_mrks = c("ZEB1","VIM","AK311497","INS-IGF2;INS")
cmem_mcell_mc_plot_marks (metacell_assignment = scLongitudinal_assign_vector,
						  metacell_fp = scLongitudinal_metacells_fp,
						  genes_to_plot = c(epithelial_genes,emt_mrks), 
						  cells_umi_tab = cmem_sumi_lt_PAR,				  
						  height = 440*1.6,
						  width  = 1139.25*1.85,
						  fig_fn = "./HCT116_scLongitudinal_heatmap.png",
						  k_log=7,
						  text_cex = 2.65,
						  lateral_genes = NULL,
						  lateral_genes_color = "red",
                          plot_cells = TRUE,
						  ideal_umi = NULL,
						  plot_top_marg = c(0,13,5,20),
						  pminpmax_cells = FALSE,
						  heatmap_lwd=0.5,
						  alternate_side_text = TRUE,
						  min_quantile_cells = 0.,
						  max_quantile_cells = 1,
						  las_metacellID = 2
						  )	
#						  
# b) project EpCAM gene expression over 2-D projection
#
cmem_plot_gene_fp_over_mc2d (gene_to_plot = "EPCAM",
							 z_score_range = 2,
							 z_colors = rev(RColorBrewer::brewer.pal(n=11, 'Spectral')),
							 metacells_fp = scLongitudinal_metacells_fp,
							 metacells_2d = scLongitudinal_metacells_2d,
							 cells_2d = scLongitudinal_metacells,
							 my_text_cex = 3.6,
							 plot_w = 900,
							 plot_h = 900,
							 cells_cex = 0.5,
							 cells_color = "grey",
							 metacells_cex = 8,
							 plot_margins = c(4,4,4,4),
							 fig_fn = "./HCT116_WT_longitudinal_EpCAM_lfp.png")							   
#
######  step 1.2.2 - Follow Epithelial and EMT-like clonal expression over time 
#
# a) project clonal distributions at different times by highlighting it over the metacell model 2D-projection
#
scLongitudinal_metacells$clones_over_time = paste0(scLongitudinal_metacells$clone,"_",scLongitudinal_metacells$day)
scLongitudinal_metacells$clones_over_time = factor (scLongitudinal_metacells$clones_over_time,
											levels = c("1D12_d_33","1D12_d_62","1D12_d_98","1D12_d_148",
													   "4E1_d_33","4E1_d_62","4E1_d_98","4E1_d_148",
													   "7A2_d_33","7A2_d_62","7A2_d_98","7A2_d_148",
													   "3B3_d_33","3B3_d_62","3B3_d_98","3B3_d_148",
													   "4B10_d_33","4B10_d_62","4B10_d_98","4B10_d_148",
													   "7B11_d_33","7B11_d_62","7B11_d_98","7B11_d_148"))
clones_over_time_colors = c(rep("#B79F00",4),rep("#F564E3",4),rep("#00BFC4",4),rep("#F8766D",4),rep("#619CFF",4),rep("#00BA38",4))
names(clones_over_time_colors) = levels(scLongitudinal_metacells$clones_over_time)
#
cmem_mcell_mc2d_plot_by_factor(				metacell_assignment = scLongitudinal_assign_vector,
											cells_metadata = scLongitudinal_metacells,
											meta_field = "clones_over_time",
											single_plot = TRUE,
											filter_values = NULL, 
											filter_name = NULL,
											ncols=4,
											show_titles = TRUE,
											title_colors = clones_over_time_colors,
											cex_main = 4,
											fig_margins = c(1,1,4.5,1),
											fig_nm = "./6_clones_over_time.png",
											metacell_colors = cmem_lfp_col_z(nm = "EPCAM", 
																			 lfp =  log2(scLongitudinal_metacells_fp),
																			 z = 2,
																			 col_spectrum = rev(RColorBrewer::brewer.pal(n=11, 'Spectral'))),
											plot_h = 1500,
											plot_w = 1500,
											mcp_2d_cex = 3,
											mcp_2d_legend_cex = 2,
											color_bg_cells = "grey90")											
#
#
# b_1) plot Epithelial program and Vimentin expression over time (scatter showing single cells)
cmem_plot_longitudinal_cells_xy (cells_metadata = cmem_sumi_lt_PAR_metadata,
								 cells_umi_tab_norm = cmem_sumi_lt_PAR_norm,
								 genes_a = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$HCT116_WT=="Epithelial"]),"[[",1)),
								 genes_b = "VIM",
								 dot_cex = 0.3,
								 text_cex = 12,
								 lab_a = "EpCAM module (log)",
								 lab_b = "Vimentin",
								 log2_exp_a = TRUE,
								 log2_exp_b = FALSE,
								 eps = 1)

							 
#								 
# b_1) plot Epithelial program and Vimentin expression over time (averaged output per clone)
cmem_plot_longitudinal_lines (cells_metadata = cmem_sumi_lt_PAR_metadata,
							  cells_umi_tab = cmem_sumi_lt_PAR,
							  clones_umi_tab1 = cmem_umi_lt_PAR_day10,
							  clones_umi_tab2 = cmem_umi_lt_PAR_day18,
							  genes_to_plot = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$HCT116_WT=="Epithelial"]),"[[",1)),
							  line_lwd = 2.5,
							  lab_y = "Epithelial module",
							  base_dir = "./",
							  gg_text_sz = 22,
							  clones_line_colors = c("#B79F00","#F564E3","#00BFC4","#F8766D","#619CFF","#00BA38"),
							  error_bar_width = .4,
							  erro_bar_lwd = 1.5,
							  ggplot_w = 6,
							  ggplot_h = 3,
							  gg_legend_position="right",
							  gg_margins = c(0,0,0,2)
							  )		
#							  
cmem_plot_longitudinal_lines (cells_metadata = cmem_sumi_lt_PAR_metadata,
							  cells_umi_tab = cmem_sumi_lt_PAR,
							  clones_umi_tab1 = cmem_umi_lt_PAR_day10,
							  clones_umi_tab2 = cmem_umi_lt_PAR_day18,
							  genes_to_plot = "VIM",
							  line_lwd = 2.5,
							  lab_y = "Vimentin",
							  base_dir = "./",
							  gg_text_sz = 22,
							  clones_line_colors = c("#B79F00","#F564E3","#00BFC4","#F8766D","#619CFF","#00BA38"),
							  error_bar_width = .4,
							  erro_bar_lwd = 1.5,
							  ggplot_w = 6,
							  ggplot_h = 3,
							  gg_legend_position="right",
							  gg_margins = c(0,0,0,2)
							  )		
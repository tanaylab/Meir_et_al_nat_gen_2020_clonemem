##
##
## Meir_et_al_2020_nat_gen_functions.r
## Date: 18.01.2020
## This file gather functions used in main analysis steps of Meir et al. 2020
##



## Load Dependency packages

# Tanay group packages
require (tgstat)

# Other packages
require(zoo)
require(Matrix)
require(RColorBrewer)
require(pheatmap)
require(cowplot)
require(dplyr)
require(reshape2)
require(glue)

### list of functions in this file:
#
# 1. cmem_load_singles_marsSeq()
#
# 2. cmem_load_clones ()
#
# 3. cmem_load_longitudinal_HCT116 <- function (min_cell = 1e3,
#										        max_cell = 2.5e4,
#										        min_clone_day10 = 8e3,
#										        min_clone_day18 = 8e3,
#										        sc_umi_dir = "expression_data/scMARSseq_umi_tables/",
#										        clones_umi_dir = "expression_data/clones_MARSseq_umi_tables/",
#										        clones_pd_dir = "expression_data/clones_MARSseq_plates_design/"
#												)
#
# 4. cmem_agg_UMIs <- function (object_name, 
#						        umiTab_dir = "expression_data/clones_MARSseq_umi_tables/",
#						        pd_dir = "expression_data/clones_MARSseq_plates_design/",
#						        b1_umiTab, b2_umiTab)
#
# 5. cmem_mcell_mc_plot_marks <- function(metacell_assignment = NULL,
#									      metacell_fp = NULL,
#									      genes_to_plot = NULL, 
#									      cells_umi_tab = NULL,
#									      height,
#									      width,
#									      fig_fn = ".",
#									      k_log=7,
#									      text_cex = 2.85,
#									      lateral_genes = NULL,
#									      lateral_genes_color = "red",
#                                          plot_cells = TRUE,
#									      ideal_umi = NULL,
#									      plot_top_marg = c(0,13,5,20),
#									      pminpmax_cells = FALSE,
#									      heatmap_lwd=0.5,
#									      alternate_side_text = TRUE,
#									      min_quantile_cells = 0,
#									      max_quantile_cells = 1
#									      )
#
# 6. cmem_plot_gene_fp_over_mc2d <- function (gene_to_plot = NULL,
#										      z_score_range = 2,
#										      metacells_fp = NULL,
#										      metacells_2d = NULL,
#										      cells_2d = NULL,
#										      my_ylim = NULL,
#										      my_xlim = NULL,
#										      my_text_cex = 3.6,
#										      plot_w = 450,
#										      plot_h = 450,
#										      cells_cex=0.5,
#										      cells_color = "grey",
#										      metacells_cex = 8,
#										      plot_margins = c(4,4,4,4),
#										      fig_fn = "./gene_over_mc2d.png"
#										      )
#
# 7. cmem_plot_gene_modules_scatter <- function (genes_x = NULL,
# 								                 genes_y = NULL,
# 								                 x_label = "lab1",
# 								                 y_label = "lab2",
# 								                 cells_umi_tab = NULL,
# 								                 metacell_assignment = NULL,
# 								                 metacells_fp = NULL,
# 								                 lfp_color_gene = NULL,
# 								                 z_score_range = 2,
# 								                 my_xlim = NULL,
# 								                 my_ylim = NULL,
# 								                 log2_exp = TRUE,
# 								                 eps=1,
# 								                 my_text_cex = 3.6,
# 								                 cells_color = "#808080",
# 							                     plot_w = 450,
# 							                     plot_h = 450,
# 								                 cex_cells = 0.5,
# 								                 cex_metacells = 8,
# 								                 plot_margins = c(4,4,4,4),
# 								                 fig_fn = "./lfp_byGene_xy_scatter.png"
# 								                 )
#
# 8. cmem_lfp_col_z = function (nm, lfp, z=4, col_spectrum = rev(RColorBrewer::brewer.pal(n=11, 'RdYlBu'))) 
#
# 9. cmem_mcell_cgraph_norm_gcor <- function (cgraph_edges = NULL,
# 										      cells_umi_tab = NULL,
# 										      K=-1,
# 										      dowsamp_value = 5000,
# 										      min_gtot=1000
#										      )
#
# 10. cmem_plot_ccNorm <- function (cc_cors,
# 							       base_dir = "./",
# 							       fig_nm,
# 							       show_legend_text=TRUE,
# 							       x_scatter_lim,
# 							       y_scatter_lim,
# 							       min_cor,
# 							       max_r_cor,
# 							       min_dif
#								   )
#
#
# 11. cmem_compare_cells_clones <- function (mod_to_plot,
#								   			 mod_nm,
#								   			 sc_mat,
#								   			 c_mat,
#								   			 minimal_cell,
#								   			 minimal_clone,
#								   			 base_dir = "./",
#								   			 legend_x,
#								   			 legend_y,
#								   			 density_x_lim,
#								   			 density_y_lim,
#								   			 ecdf_x_lim,
#								   			 density_mar = c(4.5,4,0.5,1.5),
#								   			 plot_ecdf=TRUE,
#								   			 show_legend_density=FALSE
#											 )
#
# 12. cmem_hist_gene_cells_clones <- function (gene_to_plot = NULL,
#									       min_cell = 6e3,
#									       max_cell = 4e4,
#									       min_clone = 1e4,
#									       cells_umi_table = NULL,
#									       clones_umi_table = NULL,
#										   cells_umi_table_norm = NULL,
#									       clones_umi_table_norm = NULL,
#										   log2_exp = FALSE,
#										   eps = 1,
#										   cex_text=1.8,
#									       base_dir = "./",
#										   x_label = "lab1",
#										   plot_h = 273,
#										   plot_w = 254.8,
#										   plot_margins = c(5.5,4.6,1,2),
#										   plot_legend_cells = FALSE,
#										   plot_legend_clones = FALSE,
#										   population_legend = NULL,
#										   legend_x_cells = 1.6,
#										   legend_y_cells = 400,
#										   legend_x_clones = 1.6,
#										   legend_y_clones = 400,
#								           hist_breaks = NULL,
#									       hist_colors = NULL
#										   )
#
# 13. cmem_adjust_gcor_by_samplingDepth <- function (cells_umi_tab = NULL,
#											         min_gtot = 50,
#											         min_cell = 1000,
#											         max_cell = 40000,
#											         rollmean_n_genes = 101,
#											         cor_method_spearman = TRUE,
#											         symmertrize_adj_corr = FALSE
#												     )
#
# 14. cmem_norm_gene = function(x, w, gtot_ord) 
#
# 15. cmem_permute_gene_adjust_gcor <- function (cells_umi_tab = NULL,
#										   min_gtot = 50,
#										   min_cell = 1000,
#										   max_cell = 40000,
#										   rollmean_n_genes = 101,
#										   cor_method_spearman = TRUE,
#										   symmertrize_adj_corr = FALSE,
#										   gene_to_permute=NULL)
#
# 16. cmem_barplot_top_gcor <- function (adj_gcor_mat = NULL,
#								   perm_adj_gcor_mat = NULL,
#								   gene_name = NULL,
#								   n_genes = 10,
#								   positive_corr_col = "blue",
#								   negative_corr_col = "red",
#								   positive_permuted_corr_col = "lightblue",
#								   negative_permuted_corr_col = "tomato",
#								   positive_corr_margins = c(14,4.5,7,4),
#								   negative_corr_margins = c(7,4.5,14,4),
#								   ylim_positive_corr = NULL,
#								   ylim_negative_corr = NULL,
#								   ylim_density = NULL,
#								   xlim_ecdf_positive = NULL,
#								   ylim_ecdf_positive = NULL,
#								   xlim_ecdf_negative = NULL,
#								   ylim_ecdf_negative = NULL,
#								   text_cex=1.2,
#								   barplot_w = 800,
#								   barplot_h = 400,
#								   density_plot_w = 250,
#								   density_plot_h = 250,
#								   ecdf_plot_w = 220,
#								   ecdf_plot_h = 220,
#								   base_dir = "./"
#								   )
#
# 17. cmem_ds <- function (umi_tab, ds_val)
#
# 18. cmem_scm_downsamp <- function (umis, n)
#
# cmem_plot_longitudinal_cells_xy <- function (cells_metadata = NULL,
#								 			   cells_umi_tab_norm = NULL,
#								 			   genes_a = NULL,
#								 			   genes_b = NULL,
#								 			   dot_cex = 0.3,
#								 			   text_cex = 12,
#								 			   lab_a = "lab1",
#								 			   lab_b = "lab2",
#								 			   base_dir = "./",
#								 			   log2_exp_a = FALSE,
#								 			   log2_exp_b = FALSE,
#								 			   ggplot_w = 4.123,
#								 			   ggplot_h = 5.15375,
#								 			   eps = 1
#								 			   )
##################################################################################################################



# this function loads MARS-umiTables, remove spike-INs, combines and saves them into one matrix for each cell-line
##################################################################################################################
cmem_load_singles_marsSeq <- function (umi_dir = "expression_data/scMARSseq_umi_tables/",
									   minUMI_cell=1000, 
									   maxUMI_cell=25000
									  )
{
#HCT116 parental singles
sc1=read.table(paste0(umi_dir,"ABZM0005.txt"),sep="\t")
sc2=read.table(paste0(umi_dir,"ABZM0006.txt"),sep="\t")
sc3=read.table(paste0(umi_dir,"ABZM0007.txt"),sep="\t")
sc4=read.table(paste0(umi_dir,"ABZM0008.txt"),sep="\t")
sc5=read.table(paste0(umi_dir,"ABZM0009.txt"),sep="\t")
sc6=read.table(paste0(umi_dir,"ABZM0010.txt"),sep="\t")
scTable=cbind(sc1,sc2,sc3,sc4,sc5,sc6)

#HCT116 DKO singles
sc1_DKO=read.table(paste0(umi_dir,"ABZM0019.txt"),sep="\t")
sc2_DKO=read.table(paste0(umi_dir,"ABZM0020.txt"),sep="\t")
sc3_DKO=read.table(paste0(umi_dir,"ABZM0021.txt"),sep="\t")
sc4_DKO=read.table(paste0(umi_dir,"ABZM0022.txt"),sep="\t")
sc5_DKO=read.table(paste0(umi_dir,"ABZM0023.txt"),sep="\t")
sc6_DKO=read.table(paste0(umi_dir,"ABZM0024.txt"),sep="\t")
sc7_DKO=read.table(paste0(umi_dir,"ABZM0055.txt"),sep="\t")
sc8_DKO=read.table(paste0(umi_dir,"ABZM0056.txt"),sep="\t")
sc9_DKO=read.table(paste0(umi_dir,"ABZM0057.txt"),sep="\t")
sc10_DKO=read.table(paste0(umi_dir,"ABZM0058.txt"),sep="\t")
sc11_DKO=read.table(paste0(umi_dir,"ABZM0059.txt"),sep="\t")
sc12_DKO=read.table(paste0(umi_dir,"ABZM0060.txt"),sep="\t")
sc13_DKO=read.table(paste0(umi_dir,"ABZM0039.txt"),sep="\t")
sc14_DKO=read.table(paste0(umi_dir,"ABZM0040.txt"),sep="\t")
scTable_DKO=cbind(sc1_DKO,sc2_DKO,sc3_DKO,sc4_DKO,sc5_DKO,sc6_DKO,sc7_DKO,sc8_DKO,sc9_DKO,sc10_DKO,sc11_DKO,sc12_DKO,sc13_DKO,sc14_DKO)

## WI38 singles
sc1_WI38=read.table(paste0(umi_dir,"ABZM0337.txt"),sep="\t")
sc2_WI38=read.table(paste0(umi_dir,"ABZM0338.txt"),sep="\t")
sc3_WI38=read.table(paste0(umi_dir,"ABZM0339.txt"),sep="\t")
sc4_WI38=read.table(paste0(umi_dir,"ABZM0340.txt"),sep="\t")
scTable_WI38=cbind(sc1_WI38,sc2_WI38,sc3_WI38,sc4_WI38)

#clear spikes
scTable=scTable[-c(grep("ERCC-00",rownames(scTable))),]
scTable_DKO=scTable_DKO[-c(grep("ERCC-00",rownames(scTable_DKO))),]
scTable_WI38=scTable_WI38[-c(grep("ERCC-00",rownames(scTable_WI38))),]

#filter by minimal and maximal UMI counts
filt_scTable=scTable[,( colSums(scTable)>=minUMI_cell & colSums(scTable)<=maxUMI_cell )]
filt_scTable_DKO=scTable_DKO[,( colSums(scTable_DKO)>=minUMI_cell & colSums(scTable_DKO)<=maxUMI_cell )]
filt_scTable_WI38=scTable_WI38[,( colSums(scTable_WI38)>=minUMI_cell )]

#HCT116 WT
cmem_sumi_PAR <<- filt_scTable;
cmem_sumi_PAR_norm <<- apply(cmem_sumi_PAR,2,function(x){(x/sum(x))*1e4}) 

#HCT116 DKO
cmem_sumi_DKO <<- filt_scTable_DKO;
cmem_sumi_DKO_norm <<- apply(cmem_sumi_DKO,2,function(x){(x/sum(x))*1e4}) 

#WI38
cmem_sumi_WI38 <<- filt_scTable_WI38;
cmem_sumi_WI38_norm <<- apply(cmem_sumi_WI38,2,function(x){(x/sum(x))*1e4}) 
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# this function load aggregates umi tables of clones into global variables
##################################################################################################################
cmem_load_clones = function (parental=TRUE,
							 DKO=TRUE,
							 KO=FALSE,
							 A549=FALSE,
							 H1299=TRUE,
							 WI38=TRUE,
							 minimal_UMI_per_clone = 1e4,
						 	 umi_dir = "expression_data/clones_MARSseq_umi_tables/",
							 pd_dir = "expression_data/clones_MARSseq_plates_design/"
							)
{
if (parental){
par1 = cmem_agg_UMIs (object_name="parental_clones1",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0013", b2_umiTab="ABZM0014")
par2 = cmem_agg_UMIs (object_name="parental_clones2",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0017", b2_umiTab="ABZM0018")
cmem_umi_PAR_wSpike = cbind(par1,par2)
cmem_umi_PAR = cmem_umi_PAR_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_PAR_wSpike))),]
cmem_umi_PAR = cmem_umi_PAR[,colSums(cmem_umi_PAR) >= minimal_UMI_per_clone]
cmem_umi_PAR <<- cmem_umi_PAR
cmem_umi_PAR_norm <<- apply(cmem_umi_PAR,2,function(x){(x/sum(x))*1e4})}


if (KO){
ko1 = cmem_agg_UMIs (object_name="KO_clones1",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0113", b2_umiTab="ABZM0114")
ko2 = cmem_agg_UMIs (object_name="KO_clones2",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0115", b2_umiTab="ABZM0116")

cmem_umi_KO_wSpike = cbind(ko1,ko2)
cmem_umi_KO = cmem_umi_KO_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_KO_wSpike))),]
cmem_umi_KO = cmem_umi_KO[,colSums(cmem_umi_KO) >= minimal_UMI_per_clone]
cmem_umi_KO <<- cmem_umi_KO
cmem_umi_KO_norm <<- apply(cmem_umi_KO,2,function(x){(x/sum(x))*1e4})}

if (DKO){
DKO1 = cmem_agg_UMIs (object_name="DKO_clones1",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0047", b2_umiTab="ABZM0048")
DKO2 = cmem_agg_UMIs (object_name="DKO_clones2",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0041", b2_umiTab="ABZM0042")
DKO3 = cmem_agg_UMIs (object_name="DKO_clones3",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0043", b2_umiTab="ABZM0044")
DKO4 = cmem_agg_UMIs (object_name="DKO_clones4",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0045", b2_umiTab="ABZM0046")
cmem_umi_DKO_wSpike = cbind(DKO1,DKO2,DKO3,DKO4)
cmem_umi_DKO = cmem_umi_DKO_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_DKO_wSpike))),]
cmem_umi_DKO = cmem_umi_DKO[,colSums(cmem_umi_DKO) >= minimal_UMI_per_clone]
cmem_umi_DKO <<- cmem_umi_DKO
cmem_umi_DKO_norm <<- apply(cmem_umi_DKO,2,function(x){(x/sum(x))*1e4}) }


if (A549){
a549_1 = cmem_agg_UMIs (object_name="A549_clones1",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0327", b2_umiTab="ABZM0328")
a549_2 = cmem_agg_UMIs (object_name="A549_clones2",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0329", b2_umiTab="ABZM0330")
a549_3 = cmem_agg_UMIs (object_name="A549_clones3",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0331", b2_umiTab="ABZM0332")
cmem_umi_A549_wSpike = cbind(a549_1, a549_2, a549_3)
cmem_umi_A549 = cmem_umi_A549_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_A549_wSpike))),]
cmem_umi_A549 = cmem_umi_A549[,colSums(cmem_umi_A549) >= minimal_UMI_per_clone]
cmem_umi_A549 <<- cmem_umi_A549
cmem_umi_A549_norm <<- apply(cmem_umi_A549,2,function(x){(x/sum(x))*1e4})}

if (H1299){
h1299_1 = cmem_agg_UMIs (object_name="H1299_clones1",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0317", b2_umiTab="ABZM0318")
h1299_2 = cmem_agg_UMIs (object_name="H1299_clones2",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0319", b2_umiTab="ABZM0320")
h1299_3 = cmem_agg_UMIs (object_name="H1299_clones3",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0321", b2_umiTab="ABZM0322")
h1299_4 = cmem_agg_UMIs (object_name="H1299_clones4",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0323", b2_umiTab="ABZM0324")
cmem_umi_H1299_wSpike = cbind(h1299_1, h1299_2, h1299_3, h1299_4)
cmem_umi_H1299 = cmem_umi_H1299_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_H1299_wSpike))),]
cmem_umi_H1299 = cmem_umi_H1299[,colSums(cmem_umi_H1299) >= minimal_UMI_per_clone]
cmem_umi_H1299 <<- cmem_umi_H1299
cmem_umi_H1299_norm <<- apply(cmem_umi_H1299,2,function(x){(x/sum(x))*1e4})}

if (WI38){
wi38_1 = cmem_agg_UMIs (object_name="WI38_clones1",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0333", b2_umiTab="ABZM0334")
wi38_2 = cmem_agg_UMIs (object_name="WI38_clones2",
					  umiTab_dir=umi_dir,
					  pd_dir=pd_dir,
					  b1_umiTab="ABZM0335", b2_umiTab="ABZM0336")

cmem_umi_WI38_wSpike = cbind(wi38_1, wi38_2)
cmem_umi_WI38 = cmem_umi_WI38_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_WI38_wSpike))),]
cmem_umi_WI38 = cmem_umi_WI38[,colSums(cmem_umi_WI38) >= minimal_UMI_per_clone]
cmem_umi_WI38 <<- cmem_umi_WI38
cmem_umi_WI38_norm <<- apply(cmem_umi_WI38,2,function(x){(x/sum(x))*1e4})}
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# this function loads the longitudinal transcriptome trancking of 6 selected HCT116 clones over 6 time points
# day 10: low-depth bulk RNA sampling
# day 18: low-depth bulk RNA sampling
# day 33: scRNA seq
# day 62: scRNA seq
# day 98: scRNA seq
# day 148: scRNA seq
# *Exome data on these clones from two time-points (d=78, d=168) is also available
##################################################################################################################
cmem_load_longitudinal_HCT116 <- function (min_cell = 1e3,
										   max_cell = 2.5e4,
										   min_clone_day10 = 8e3,
										   min_clone_day18 = 8e3,
										   sc_umi_dir = "expression_data/scMARSseq_umi_tables/",
										   clones_umi_dir = "expression_data/clones_MARSseq_umi_tables/",
										   clones_pd_dir = "expression_data/clones_MARSseq_plates_design/")
{										   
#load long term singles
a_3b3=read.table(paste0(sc_umi_dir,"ABZM0141.txt"),sep="\t")
b_3b3=read.table(paste0(sc_umi_dir,"ABZM0142.txt"),sep="\t")
a_1d12=read.table(paste0(sc_umi_dir,"ABZM0143.txt"),sep="\t")
b_1d12=read.table(paste0(sc_umi_dir,"ABZM0144.txt"),sep="\t")
a_7b11=read.table(paste0(sc_umi_dir,"ABZM0145.txt"),sep="\t")
b_7b11=read.table(paste0(sc_umi_dir,"ABZM0146.txt"),sep="\t")
a_7a2=read.table(paste0(sc_umi_dir,"ABZM0147.txt"),sep="\t")
b_7a2=read.table(paste0(sc_umi_dir,"ABZM0148.txt"),sep="\t")
a_4b10=read.table(paste0(sc_umi_dir,"ABZM0149.txt"),sep="\t")
b_4b10=read.table(paste0(sc_umi_dir,"ABZM0150.txt"),sep="\t")
a_4e1=read.table(paste0(sc_umi_dir,"ABZM0151.txt"),sep="\t")
b_4e1=read.table(paste0(sc_umi_dir,"ABZM0152.txt"),sep="\t")

scTable_lt3_3b3=cbind(a_3b3,b_3b3)
scTable_lt3_1d12=cbind(a_1d12,b_1d12)
scTable_lt3_7b11=cbind(a_7b11,b_7b11)
scTable_lt3_7a2=cbind(a_7a2,b_7a2)
scTable_lt3_4b10=cbind(a_4b10,b_4b10)
scTable_lt3_4e1=cbind(a_4e1,b_4e1)

lt4_a_3b3=read.table(paste0(sc_umi_dir,"ABZM0163.txt"),sep="\t")
lt4_b_3b3=read.table(paste0(sc_umi_dir,"ABZM0164.txt"),sep="\t")
lt4_a_1d12=read.table(paste0(sc_umi_dir,"ABZM0165.txt"),sep="\t")
lt4_b_1d12=read.table(paste0(sc_umi_dir,"ABZM0166.txt"),sep="\t")
lt4_a_7b11=read.table(paste0(sc_umi_dir,"ABZM0167.txt"),sep="\t")
lt4_b_7b11=read.table(paste0(sc_umi_dir,"ABZM0168.txt"),sep="\t")
lt4_a_7a2=read.table(paste0(sc_umi_dir,"ABZM0169.txt"),sep="\t")
lt4_b_7a2=read.table(paste0(sc_umi_dir,"ABZM0170.txt"),sep="\t")
lt4_a_4b10=read.table(paste0(sc_umi_dir,"ABZM0171.txt"),sep="\t")
lt4_b_4b10=read.table(paste0(sc_umi_dir,"ABZM0172.txt"),sep="\t")
lt4_a_4e1=read.table(paste0(sc_umi_dir,"ABZM0173.txt"),sep="\t")
lt4_b_4e1=read.table(paste0(sc_umi_dir,"ABZM0174.txt"),sep="\t")
scTable_lt4_3b3=cbind(lt4_a_3b3,lt4_b_3b3)
scTable_lt4_1d12=cbind(lt4_a_1d12,lt4_b_1d12)
scTable_lt4_7b11=cbind(lt4_a_7b11,lt4_b_7b11)
scTable_lt4_7a2=cbind(lt4_a_7a2,lt4_b_7a2)
scTable_lt4_4b10=cbind(lt4_a_4b10,lt4_b_4b10)
scTable_lt4_4e1=cbind(lt4_a_4e1,lt4_b_4e1)

lt5_a_3b3=read.table(paste0(sc_umi_dir,"ABZM0189.txt"),sep="\t")
lt5_b_3b3=read.table(paste0(sc_umi_dir,"ABZM0190.txt"),sep="\t")
lt5_a_1d12=read.table(paste0(sc_umi_dir,"ABZM0191.txt"),sep="\t")
lt5_b_1d12=read.table(paste0(sc_umi_dir,"ABZM0192.txt"),sep="\t")
lt5_a_7b11=read.table(paste0(sc_umi_dir,"ABZM0193.txt"),sep="\t")
lt5_b_7b11=read.table(paste0(sc_umi_dir,"ABZM0194.txt"),sep="\t")
lt5_a_7a2=read.table(paste0(sc_umi_dir,"ABZM0195.txt"),sep="\t")
lt5_b_7a2=read.table(paste0(sc_umi_dir,"ABZM0196.txt"),sep="\t")
lt5_a_4b10=read.table(paste0(sc_umi_dir,"ABZM0197.txt"),sep="\t")
lt5_b_4b10=read.table(paste0(sc_umi_dir,"ABZM0198.txt"),sep="\t")
lt5_a_4e1=read.table(paste0(sc_umi_dir,"ABZM0199.txt"),sep="\t")
lt5_b_4e1=read.table(paste0(sc_umi_dir,"ABZM0200.txt"),sep="\t")
scTable_lt5_3b3=cbind(lt5_a_3b3,lt5_b_3b3)
scTable_lt5_1d12=cbind(lt5_a_1d12,lt5_b_1d12)
scTable_lt5_7b11=cbind(lt5_a_7b11,lt5_b_7b11)
scTable_lt5_7a2=cbind(lt5_a_7a2,lt5_b_7a2)
scTable_lt5_4b10=cbind(lt5_a_4b10,lt5_b_4b10)
scTable_lt5_4e1=cbind(lt5_a_4e1,lt5_b_4e1)

lt6_a_3b3=read.table(paste0(sc_umi_dir,"ABZM0211.txt"),sep="\t")
lt6_b_3b3=read.table(paste0(sc_umi_dir,"ABZM0212.txt"),sep="\t")
lt6_a_1d12=read.table(paste0(sc_umi_dir,"ABZM0213.txt"),sep="\t")
lt6_b_1d12=read.table(paste0(sc_umi_dir,"ABZM0214.txt"),sep="\t")
lt6_a_7b11=read.table(paste0(sc_umi_dir,"ABZM0215.txt"),sep="\t")
lt6_b_7b11=read.table(paste0(sc_umi_dir,"ABZM0216.txt"),sep="\t")
lt6_a_7a2=read.table(paste0(sc_umi_dir,"ABZM0217.txt"),sep="\t")
lt6_b_7a2=read.table(paste0(sc_umi_dir,"ABZM0218.txt"),sep="\t")
lt6_a_4b10=read.table(paste0(sc_umi_dir,"ABZM0219.txt"),sep="\t")
lt6_b_4b10=read.table(paste0(sc_umi_dir,"ABZM0220.txt"),sep="\t")
lt6_a_4e1=read.table(paste0(sc_umi_dir,"ABZM0221.txt"),sep="\t")
lt6_b_4e1=read.table(paste0(sc_umi_dir,"ABZM0222.txt"),sep="\t")
scTable_lt6_3b3=cbind(lt6_a_3b3,lt6_b_3b3)
scTable_lt6_1d12=cbind(lt6_a_1d12,lt6_b_1d12)
scTable_lt6_7b11=cbind(lt6_a_7b11,lt6_b_7b11)
scTable_lt6_7a2=cbind(lt6_a_7a2,lt6_b_7a2)
scTable_lt6_4b10=cbind(lt6_a_4b10,lt6_b_4b10)
scTable_lt6_4e1=cbind(lt6_a_4e1,lt6_b_4e1)

#clear spikes
scTable_lt3_3b3=scTable_lt3_3b3[-c(grep("ERCC-00",rownames(scTable_lt3_3b3))),]
scTable_lt3_1d12=scTable_lt3_1d12[-c(grep("ERCC-00",rownames(scTable_lt3_1d12))),]
scTable_lt3_7b11=scTable_lt3_7b11[-c(grep("ERCC-00",rownames(scTable_lt3_7b11))),]
scTable_lt3_7a2=scTable_lt3_7a2[-c(grep("ERCC-00",rownames(scTable_lt3_7a2))),]
scTable_lt3_4b10=scTable_lt3_4b10[-c(grep("ERCC-00",rownames(scTable_lt3_4b10))),]
scTable_lt3_4e1=scTable_lt3_4e1[-c(grep("ERCC-00",rownames(scTable_lt3_4e1))),]

scTable_lt4_3b3=scTable_lt4_3b3[-c(grep("ERCC-00",rownames(scTable_lt4_3b3))),]
scTable_lt4_1d12=scTable_lt4_1d12[-c(grep("ERCC-00",rownames(scTable_lt4_1d12))),]
scTable_lt4_7b11=scTable_lt4_7b11[-c(grep("ERCC-00",rownames(scTable_lt4_7b11))),]
scTable_lt4_7a2=scTable_lt4_7a2[-c(grep("ERCC-00",rownames(scTable_lt4_7a2))),]
scTable_lt4_4b10=scTable_lt4_4b10[-c(grep("ERCC-00",rownames(scTable_lt4_4b10))),]
scTable_lt4_4e1=scTable_lt4_4e1[-c(grep("ERCC-00",rownames(scTable_lt4_4e1))),]

scTable_lt5_3b3=scTable_lt5_3b3[-c(grep("ERCC-00",rownames(scTable_lt5_3b3))),]
scTable_lt5_1d12=scTable_lt5_1d12[-c(grep("ERCC-00",rownames(scTable_lt5_1d12))),]
scTable_lt5_7b11=scTable_lt5_7b11[-c(grep("ERCC-00",rownames(scTable_lt5_7b11))),]
scTable_lt5_7a2=scTable_lt5_7a2[-c(grep("ERCC-00",rownames(scTable_lt5_7a2))),]
scTable_lt5_4b10=scTable_lt5_4b10[-c(grep("ERCC-00",rownames(scTable_lt5_4b10))),]
scTable_lt5_4e1=scTable_lt5_4e1[-c(grep("ERCC-00",rownames(scTable_lt5_4e1))),]

scTable_lt6_3b3=scTable_lt6_3b3[-c(grep("ERCC-00",rownames(scTable_lt6_3b3))),]
scTable_lt6_1d12=scTable_lt6_1d12[-c(grep("ERCC-00",rownames(scTable_lt6_1d12))),]
scTable_lt6_7b11=scTable_lt6_7b11[-c(grep("ERCC-00",rownames(scTable_lt6_7b11))),]
scTable_lt6_7a2=scTable_lt6_7a2[-c(grep("ERCC-00",rownames(scTable_lt6_7a2))),]
scTable_lt6_4b10=scTable_lt6_4b10[-c(grep("ERCC-00",rownames(scTable_lt6_4b10))),]
scTable_lt6_4e1=scTable_lt6_4e1[-c(grep("ERCC-00",rownames(scTable_lt6_4e1))),]

#filter by minimal and maximal UMI counts
filt_scTable_lt3_3b3=scTable_lt3_3b3[,( colSums(scTable_lt3_3b3)>=min_cell & colSums(scTable_lt3_3b3)<=max_cell )]
filt_scTable_lt3_1d12=scTable_lt3_1d12[,( colSums(scTable_lt3_1d12)>=min_cell & colSums(scTable_lt3_1d12)<=max_cell )]
filt_scTable_lt3_7b11=scTable_lt3_7b11[,( colSums(scTable_lt3_7b11)>=min_cell & colSums(scTable_lt3_7b11)<=max_cell )]
filt_scTable_lt3_7a2=scTable_lt3_7a2[,( colSums(scTable_lt3_7a2)>=min_cell & colSums(scTable_lt3_7a2)<=max_cell )]
filt_scTable_lt3_4b10=scTable_lt3_4b10[,( colSums(scTable_lt3_4b10)>=min_cell & colSums(scTable_lt3_4b10)<=max_cell )]
filt_scTable_lt3_4e1=scTable_lt3_4e1[,( colSums(scTable_lt3_4e1)>=min_cell & colSums(scTable_lt3_4e1)<=max_cell )]

filt_scTable_lt4_3b3=scTable_lt4_3b3[,( colSums(scTable_lt4_3b3)>=min_cell & colSums(scTable_lt4_3b3)<=max_cell )]
filt_scTable_lt4_1d12=scTable_lt4_1d12[,( colSums(scTable_lt4_1d12)>=min_cell & colSums(scTable_lt4_1d12)<=max_cell )]
filt_scTable_lt4_7b11=scTable_lt4_7b11[,( colSums(scTable_lt4_7b11)>=min_cell & colSums(scTable_lt4_7b11)<=max_cell )]
filt_scTable_lt4_7a2=scTable_lt4_7a2[,( colSums(scTable_lt4_7a2)>=min_cell & colSums(scTable_lt4_7a2)<=max_cell )]
filt_scTable_lt4_4b10=scTable_lt4_4b10[,( colSums(scTable_lt4_4b10)>=min_cell & colSums(scTable_lt4_4b10)<=max_cell )]
filt_scTable_lt4_4e1=scTable_lt4_4e1[,( colSums(scTable_lt4_4e1)>=min_cell & colSums(scTable_lt4_4e1)<=max_cell )]

filt_scTable_lt5_3b3=scTable_lt5_3b3[,( colSums(scTable_lt5_3b3)>=min_cell & colSums(scTable_lt5_3b3)<=max_cell )]
filt_scTable_lt5_1d12=scTable_lt5_1d12[,( colSums(scTable_lt5_1d12)>=min_cell & colSums(scTable_lt5_1d12)<=max_cell )]
filt_scTable_lt5_7b11=scTable_lt5_7b11[,( colSums(scTable_lt5_7b11)>=min_cell & colSums(scTable_lt5_7b11)<=max_cell )]
filt_scTable_lt5_7a2=scTable_lt5_7a2[,( colSums(scTable_lt5_7a2)>=min_cell & colSums(scTable_lt5_7a2)<=max_cell )]
filt_scTable_lt5_4b10=scTable_lt5_4b10[,( colSums(scTable_lt5_4b10)>=min_cell & colSums(scTable_lt5_4b10)<=max_cell )]
filt_scTable_lt5_4e1=scTable_lt5_4e1[,( colSums(scTable_lt5_4e1)>=min_cell & colSums(scTable_lt5_4e1)<=max_cell )]

filt_scTable_lt6_3b3=scTable_lt6_3b3[,( colSums(scTable_lt6_3b3)>=min_cell & colSums(scTable_lt6_3b3)<=max_cell )]
filt_scTable_lt6_1d12=scTable_lt6_1d12[,( colSums(scTable_lt6_1d12)>=min_cell & colSums(scTable_lt6_1d12)<=max_cell )]
filt_scTable_lt6_7b11=scTable_lt6_7b11[,( colSums(scTable_lt6_7b11)>=min_cell & colSums(scTable_lt6_7b11)<=max_cell )]
filt_scTable_lt6_7a2=scTable_lt6_7a2[,( colSums(scTable_lt6_7a2)>=min_cell & colSums(scTable_lt6_7a2)<=max_cell )]
filt_scTable_lt6_4b10=scTable_lt6_4b10[,( colSums(scTable_lt6_4b10)>=min_cell & colSums(scTable_lt6_4b10)<=max_cell )]
filt_scTable_lt6_4e1=scTable_lt6_4e1[,( colSums(scTable_lt6_4e1)>=min_cell & colSums(scTable_lt6_4e1)<=max_cell )]

cells_3b3 = c(colnames(filt_scTable_lt3_3b3),colnames(filt_scTable_lt4_3b3),colnames(filt_scTable_lt5_3b3),colnames(filt_scTable_lt6_3b3))
cells_1d12 = c(colnames(filt_scTable_lt3_1d12),colnames(filt_scTable_lt4_1d12),colnames(filt_scTable_lt5_1d12),colnames(filt_scTable_lt6_1d12))
cells_7b11 = c(colnames(filt_scTable_lt3_7b11),colnames(filt_scTable_lt4_7b11),colnames(filt_scTable_lt5_7b11),colnames(filt_scTable_lt6_7b11))
cells_7a2 = c(colnames(filt_scTable_lt3_7a2),colnames(filt_scTable_lt4_7a2),colnames(filt_scTable_lt5_7a2),colnames(filt_scTable_lt6_7a2))
cells_4b10 = c(colnames(filt_scTable_lt3_4b10),colnames(filt_scTable_lt4_4b10),colnames(filt_scTable_lt5_4b10),colnames(filt_scTable_lt6_4b10))
cells_4e1 = c(colnames(filt_scTable_lt3_4e1),colnames(filt_scTable_lt4_4e1),colnames(filt_scTable_lt5_4e1),colnames(filt_scTable_lt6_4e1))
cells_lt3 = c(colnames(filt_scTable_lt3_3b3),colnames(filt_scTable_lt3_1d12),colnames(filt_scTable_lt3_7b11),colnames(filt_scTable_lt3_7a2),colnames(filt_scTable_lt3_4b10),colnames(filt_scTable_lt3_4e1))
cells_lt4 = c(colnames(filt_scTable_lt4_3b3),colnames(filt_scTable_lt4_1d12),colnames(filt_scTable_lt4_7b11),colnames(filt_scTable_lt4_7a2),colnames(filt_scTable_lt4_4b10),colnames(filt_scTable_lt4_4e1))
cells_lt5 = c(colnames(filt_scTable_lt5_3b3),colnames(filt_scTable_lt5_1d12),colnames(filt_scTable_lt5_7b11),colnames(filt_scTable_lt5_7a2),colnames(filt_scTable_lt5_4b10),colnames(filt_scTable_lt5_4e1))
cells_lt6 = c(colnames(filt_scTable_lt6_3b3),colnames(filt_scTable_lt6_1d12),colnames(filt_scTable_lt6_7b11),colnames(filt_scTable_lt6_7a2),colnames(filt_scTable_lt6_4b10),colnames(filt_scTable_lt6_4e1))

cmem_sumi_lt_PAR = cbind(filt_scTable_lt3_3b3,filt_scTable_lt3_1d12,filt_scTable_lt3_7b11,filt_scTable_lt3_7a2,filt_scTable_lt3_4b10,filt_scTable_lt3_4e1,
					     filt_scTable_lt4_3b3,filt_scTable_lt4_1d12,filt_scTable_lt4_7b11,filt_scTable_lt4_7a2,filt_scTable_lt4_4b10,filt_scTable_lt4_4e1,
					     filt_scTable_lt5_3b3,filt_scTable_lt5_1d12,filt_scTable_lt5_7b11,filt_scTable_lt5_7a2,filt_scTable_lt5_4b10,filt_scTable_lt5_4e1,
					     filt_scTable_lt6_3b3,filt_scTable_lt6_1d12,filt_scTable_lt6_7b11,filt_scTable_lt6_7a2,filt_scTable_lt6_4b10,filt_scTable_lt6_4e1)

cmem_sumi_lt_PAR_norm = apply(cmem_sumi_lt_PAR,2,function(x){(x/sum(x))*1e4}) 
cmem_sumi_lt_PAR_metadata = data.frame(clone = rep(NA,ncol(cmem_sumi_lt_PAR)), day = rep(NA,ncol(cmem_sumi_lt_PAR)), row.names = colnames(cmem_sumi_lt_PAR))
cmem_sumi_lt_PAR_metadata[cells_3b3,"clone"] = "3B3"
cmem_sumi_lt_PAR_metadata[cells_1d12,"clone"] = "1D12"
cmem_sumi_lt_PAR_metadata[cells_7b11,"clone"] = "7B11"
cmem_sumi_lt_PAR_metadata[cells_7a2,"clone"] = "7A2"
cmem_sumi_lt_PAR_metadata[cells_4b10,"clone"] = "4B10"
cmem_sumi_lt_PAR_metadata[cells_4e1,"clone"] = "4E1"
cmem_sumi_lt_PAR_metadata[cells_lt3,"day"] = "d_33"
cmem_sumi_lt_PAR_metadata[cells_lt4,"day"] = "d_62"
cmem_sumi_lt_PAR_metadata[cells_lt5,"day"] = "d_98"
cmem_sumi_lt_PAR_metadata[cells_lt6,"day"] = "d_148"

cmem_sumi_lt_PAR_metadata[,"day"] = factor(cmem_sumi_lt_PAR_metadata[,"day"], levels = c("d_33","d_62","d_98","d_148"))
cmem_sumi_lt_PAR_metadata[,"clone"] = factor(cmem_sumi_lt_PAR_metadata[,"clone"], levels=rev(c("7B11","4B10","3B3","7A2","4E1","1D12")))
cmem_sumi_lt_PAR_metadata[,"color"] = c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")[ifelse(rownames(cmem_sumi_lt_PAR_metadata) %in% cells_3b3,1,
																									 ifelse(rownames(cmem_sumi_lt_PAR_metadata) %in% cells_1d12,2,
																									 ifelse(rownames(cmem_sumi_lt_PAR_metadata) %in% cells_7b11,3,
																									 ifelse(rownames(cmem_sumi_lt_PAR_metadata) %in% cells_7a2,4,
																									 ifelse(rownames(cmem_sumi_lt_PAR_metadata) %in% cells_4b10,5,
																									 6)))))]


cmem_sumi_lt_PAR <<- cmem_sumi_lt_PAR
cmem_sumi_lt_PAR_norm <<- cmem_sumi_lt_PAR_norm
cmem_sumi_lt_PAR_metadata <<- cmem_sumi_lt_PAR_metadata

# load low depth long term clones (sampled after 10 and 18 days)
lt1_1 = cmem_agg_UMIs (object_name="parental_longTerm_time1_1",
					   umiTab_dir = clones_umi_dir,
					   pd_dir = clones_pd_dir,
					   b1_umiTab="ABZM0117", b2_umiTab="ABZM0118")
lt1_2 = cmem_agg_UMIs (object_name="parental_longTerm_time1_2",
					   umiTab_dir = clones_umi_dir,
					   pd_dir = clones_pd_dir,
					   b1_umiTab="ABZM0119", b2_umiTab="ABZM0120")
cmem_umi_lt1_wSpike = cbind(lt1_1,lt1_2)
cmem_umi_lt1  = cmem_umi_lt1_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_lt1_wSpike))),]
cmem_umi_lt_PAR_day10 <<- cmem_umi_lt1[,colSums(cmem_umi_lt1) >= min_clone_day10]
cmem_umi_lt_PAR_day10_norm <<- apply(cmem_umi_lt_PAR_day10,2,function(x){(x/sum(x))*1e4})

lt2_1 = cmem_agg_UMIs (object_name="parental_longTerm_time2_1",
					   umiTab_dir = clones_umi_dir,
					   pd_dir = clones_pd_dir,
					   b1_umiTab="ABZM0223", b2_umiTab="ABZM0224")
lt2_2 = cmem_agg_UMIs (object_name="parental_longTerm_time2_2",
					   umiTab_dir = clones_umi_dir,
					   pd_dir = clones_pd_dir,
					   b1_umiTab="ABZM0225", b2_umiTab="ABZM0226")
lt2_3 = cmem_agg_UMIs (object_name="parental_longTerm_time2_3",
					   umiTab_dir = clones_umi_dir,
					   pd_dir = clones_pd_dir,
					   b1_umiTab="ABZM0227", b2_umiTab="ABZM0228")
					   
cmem_umi_lt2_wSpike = cbind(lt2_1,lt2_2,lt2_3)
cmem_umi_lt2 = cmem_umi_lt2_wSpike[-c(grep("ERCC-00",rownames(cmem_umi_lt2_wSpike))),]
cmem_umi_lt_PAR_day18 <<- cmem_umi_lt2[,colSums(cmem_umi_lt2) >= min_clone_day18]
cmem_umi_lt_PAR_day18_norm <<- apply(cmem_umi_lt_PAR_day18,2,function(x){(x/sum(x))*1e4})
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# This function sums up UMI counts for custom 384 MARS-Seq plate (pooled in 2 batches), where multiple wells are assigned to a single sample.
# In addition to the 2 umi tables, the function needs 2 tables
#	• wellsCells_<object_name>.txt file, with conversion key of well names as appear umi table ("PZM0010382","PZM0021007"..) to well coordinates ("A1,O22..")
#	• plateDesign_<object_name>.txt conversion key of well_coordinate ("A1,O22..") to sample name("cloneDKO246","clone115"...)
#
# note: this function should be used for pooled 384-wells plates in 2 sets of 192 barcodes (therefore having 2 amplification batches per plate)
##################################################################################################################
cmem_agg_UMIs <- function (object_name, 
						   umiTab_dir = "expression_data/clones_MARSseq_umi_tables/",
						   pd_dir = "expression_data/clones_MARSseq_plates_design/",
						   b1_umiTab, b2_umiTab)
{
wc_obj=read.table(paste0(pd_dir,"wellsCells_",object_name,".txt"),
				  sep="\t",header=T,colClasses=c(rep("character",2),rep("NULL",11)))		  
pd_obj=read.table(paste0(pd_dir,"plateDesign_",object_name,".txt"),sep="\t",header=TRUE)
b1_obj=read.table(paste0(umiTab_dir,b1_umiTab,".txt"),sep="\t")
b2_obj=read.table(paste0(umiTab_dir,b2_umiTab,".txt"),sep="\t")
umis_obj=cbind(b1_obj,b2_obj)
wc_pd=data.frame(well=wc_obj$Well_ID,
				 clone=pd_obj$ID[match(wc_obj$well_coordinates,pd_obj$well)])
order_clones=match(colnames(umis_obj),wc_pd$well)
return(t(apply(umis_obj,1,function(x){tapply(x,wc_pd$clone[order_clones],sum)})))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################





# This function plots heatmap of expression foot-print (fp) of marker genes in cells or in metacells
##################################################################################################################
cmem_mcell_mc_plot_marks <- function(metacell_assignment = NULL,
									metacell_fp = NULL,
									genes_to_plot = NULL, 
									cells_umi_tab = NULL,
									height,
									width,
									fig_fn = "./cells_heatmap_mrks.png",
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
									max_quantile_cells = 1,
									las_metacellID = 1
									)
{
        mcp_heatmap_height = height
        mcp_heatmap_width = width
		mcp_heatmap_text_cex = text_cex
        mcp_heatmap_ideal_umi = ideal_umi
		mcp_heatmap_fp_shades = colorRampPalette(list(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"),
													  c("white","#FFFFCC","#FFEDA0","#FED976","#FEB24C","#FD8D3C","#FC4E2A","#E31A1C","#BD0026","#800026","black"))[[ifelse(plot_cells,2,1)]])(1000)

						
        if(is.null(metacell_assignment)) {
                stop("undefined meta cell assignment ")
        }
		  if(is.null(metacell_fp) & (!plot_cells)) {
                stop("undefined meta cell foot-print expression ")
        }
		
        n_mc = ncol(metacell_fp)
        if(is.null(genes_to_plot)) {
                stop("undefined genes to when trying to plot markers")
        }
        scmat = cells_umi_tab
        if(is.null(scmat)) {
                stop("undefined mat object when trying to plot markers")
        }
		
        if(length(intersect(names(metacell_assignment), colnames(scmat))) != length(names(metacell_assignment))) {
                stop("cells in meta cell are missing from provided matrix in mc_plot_marks")
        }
        if(is.null(mcp_heatmap_ideal_umi)) {
                mcp_heatmap_ideal_umi = quantile(colSums(scmat), 0.25)
        }

        cell_ord = names(metacell_assignment)[order(metacell_assignment)]

        good_marks = intersect(genes_to_plot, rownames(metacell_fp))
        lateral_marks = c()
        if(!is.null(lateral_genes)) {
                lateral_marks = intersect(lateral_genes, rownames(metacell_fp))
                good_marks = setdiff(good_marks, lateral_marks)
        }

        if(length(good_marks) < 2) {
                stop("Could not get >=2 markers to plot marker matrix")
        }

        gene_folds = metacell_fp

        if(!is.null(fig_fn)) {
          if(is.null(mcp_heatmap_height)){
            mcp_heatmap_height = 16*length(c(good_marks,lateral_marks)) + 250
          }
          if(is.null(mcp_heatmap_width)){
            mcp_heatmap_width= max(min(3000,length(cell_ord)+200),800)
          }
          png(fig_fn, w = mcp_heatmap_width, h = mcp_heatmap_height);
        }
        layout(matrix(c(1,2),nrow=2),heights=c(mcp_heatmap_height, 100))

        par (mar = plot_top_marg)

        if(plot_cells) {
                mat = as.matrix(scmat[c(good_marks, lateral_marks),names(metacell_assignment)])
                mat = mat[,cell_ord]
                totu = colSums(scmat[,names(metacell_assignment)])
                mat = t(t(mat)/totu)*mcp_heatmap_ideal_umi
					if (pminpmax_cells){
					mat= pmin(pmax(mat,quantile(mat,min_quantile_cells)),quantile(mat,max_quantile_cells))
					}
                lus_1 = log2(1+k_log*mat)
                lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
                if (length(cell_ord) < mcp_heatmap_width) {
                        smooth_n = 1
                }
                smooth_n = max(2,ceiling(2*length(cell_ord)/mcp_heatmap_width))
                lus_smoo = t(apply(lus, 1, function(x) rollmean(x,smooth_n, fill=0)))

                image(t(lus_smoo), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n')

                cell_x = rep(NA, time=length(cell_ord))
                names(cell_x) = names(metacell_assignment)
                cell_x[cell_ord] = 1:length(cell_ord)
                cl_x = tapply(cell_x, metacell_assignment, mean)/length(cell_ord)
				mtext(1:n_mc, side = 3, at=cl_x, las=las_metacellID, line = 0.2, cex=mcp_heatmap_text_cex,font=2)
				cl_x_b = tapply(cell_x, metacell_assignment, max)/length(cell_ord)
                abline(v=cl_x_b, lwd=heatmap_lwd)
        } else {
                mat = log2(gene_folds[c(good_marks, lateral_marks), ])
                mat = pmax(pmin(mat,3),-3)
                image(t(mat), col=mcp_heatmap_fp_shades, xaxt='n', yaxt='n', zlim=c(-3,3))
                mtext(1:n_mc, side = 3, at=seq(0,1,l=n_mc), las=2, line = 2, cex=mcp_heatmap_text_cex)
        }

        all_marks = c(good_marks, lateral_marks)
        g_n = length(all_marks)
        gene_cols = rep("black", g_n)
        if(length(lateral_marks) != 0) {
                gene_cols[(length(good_marks)+1):length(all_marks)] = lateral_genes_color
        }

        if(alternate_side_text) {
                odd = seq(1,g_n,2)
                even = seq(2,g_n,2)
                mtext(substr(all_marks[odd],1,8),
                        at=seq(0,1,length.out=g_n)[odd],
                        side=2, las=2, cex=mcp_heatmap_text_cex, col=gene_cols[odd])
                mtext(substr(all_marks[even],1,8),
                        at=seq(0,1,length.out=g_n)[even],
                        side=4, las=2, cex=mcp_heatmap_text_cex, col=gene_cols[even])
        } else {
                mtext(substr(all_marks,1,8),
                        at=seq(0,1,length.out=g_n),
                        side=2, las=2, cex=mcp_heatmap_text_cex, col=gene_cols)
                mtext(substr(all_marks,1,8),
                        at=seq(0,1,length.out=g_n),
                        side=4, las=2, cex=mcp_heatmap_text_cex, col=gene_cols)
        }

        dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# this function gets metacell expression footprint, as well as 2D-projection coordinated of cells and metacells,
# and plots fp of specific gene over the 2D-projection
##################################################################################################################
cmem_plot_gene_fp_over_mc2d <- function (gene_to_plot = NULL,
										 z_score_range = 2,
										 z_colors = rev(RColorBrewer::brewer.pal(n=11, 'RdYlBu')),
										 metacells_fp = NULL,
										 metacells_2d = NULL,
										 cells_2d = NULL,
										 my_ylim = NULL,
										 my_xlim = NULL,									 
										 my_text_cex = 3.6,
										 plot_w = 450,
										 plot_h = 450,
										 cells_cex=0.5,
										 cells_color = "grey",
										 metacells_cex = 8,
										 plot_margins = c(4,4,4,4),
										 fig_fn = "./gene_over_mc2d.png")
{
        if(is.null(gene_to_plot) | (!gene_to_plot %in% rownames(metacells_fp)) ) {
                stop("undefined gene to plot ")
        }
		  if(is.null(metacells_fp) ) {
                stop("undefined meta cell foot-print expression ")
        }
		  if(is.null(cells_2d) | is.null(metacells_2d) ) {
                stop("undefined 2d-projections of cells or metacells")
        }		
lfp = log2(metacells_fp)
png(fig_fn,height=plot_h,width=plot_w)
par(mar=plot_margins)
plot(cells_2d$x_proj, cells_2d$y_proj,
	 pch = 19,
	 cex = cells_cex, col = cells_color,
	 xaxt = "n",
	 yaxt = "n",
	 xlim = list(range(cells_2d$x_proj,na.rm=T),my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
	 ylim = list(range(cells_2d$y_proj,na.rm=T),my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 cex.axis = my_text_cex,
	 cex.lab = my_text_cex,
	 xlab="",
	 ylab="") 
title(xlab="X projection", cex.lab=my_text_cex, line=2)
title(ylab="Y projection", cex.lab=my_text_cex, line=1)
points(metacells_2d$x_proj, metacells_2d$y_proj, pch=21, cex=8, bg = cmem_lfp_col_z(nm = gene_to_plot, 
																					lfp =  lfp,
																					z = z_score_range,
																					col_spectrum = z_colors))
text(metacells_2d$x_proj, metacells_2d$y_proj, 1:ncol(lfp), cex=my_text_cex)
dev.off()

}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# this function should plot something over the 2d space by factor. colors of cells defined by their metacells's lfp of any gene 
##################################################################################################################
cmem_mcell_mc2d_plot_by_factor <- function (metacell_assignment = NULL,
											genes_to_plot = NULL, 
											cells_metadata = NULL,
											meta_field,
											show_titles = TRUE,
											title_colors = NULL,
											single_plot = TRUE,
											filter_values = NULL, 
											filter_name = NULL,
											ncols=NULL,
											fig_nm = NULL,
											cex_main = 3,
											metacell_colors = NULL,
											plot_h = 2000,
											plot_w = 2000,
											fig_margins = c(0.5,0.5,3,0.5),
											mcp_2d_cex = 3,
											mcp_2d_legend_cex = 2,
											color_bg_cells = "grey90"
											)
{
  bg_col="grey90"

  if(is.null(metacell_assignment)) {
    stop("missing mc assignment")
  }
  
  if(is.null(metacell_colors)) {
    warning("missing mc colors to highlight cells coordinates, coloring all black")
  }


  if (is.null(cells_metadata)) {
    stop("missing metadata on cells")
  }

  if (any(rownames(cells_metadata) != names(metacell_assignment))) {
    stop("cells mismatch between metacells and umi table")
  }

  c_by_f = split(names(metacell_assignment), cells_metadata[names(metacell_assignment), meta_field])
  if (is.null(filter_values)) {
    filter_values = names(c_by_f)
  }
  else {
    c_by_f = c_by_f[names(c_by_f) %in% filter_values]
  }

 cols=list(metacell_colors,rep("black",length(table(metacell_assignment))))[[ifelse(is.null(metacell_colors),2,1)]]
 
  if (single_plot) {
    ny = ifelse(is.null(ncols), floor(sqrt(length(c_by_f))), ncols)
    nx = ceiling((length(c_by_f)/ny))

    if (is.null(fig_nm)){ stop ("insert file name to save the plot to") }
	
    png(fig_nm, width = plot_w, height = plot_w / ny * nx)

    layout(matrix(1:(nx*ny), nx, ny, byrow=T))
    par(mar = fig_margins)
  }
  for (meta_field_v in filter_values) {
    ccells = c_by_f[[meta_field_v]]

    if (!single_plot) {
      fig_nm = paste0("mc2d", sprintf("2d_proj_%s", meta_field_v, ifelse(is.null(filter_name), "", paste0(filter_name, "_"))), sprintf("%s/%s.by_%s", base_dir, "mc2D", meta_field))
      png(fig_nm, width = plot_w, height = plot_h)
      par(mar = fig_margins)
    }

    plot(cells_metadata$x_proj, cells_metadata$y_proj, cex=mcp_2d_cex, pch=21, col=bg_col, bg=bg_col, xlab="", ylab="", xaxt='n', yaxt='n')
    points(cells_metadata[ccells,"x_proj"], cells_metadata[ccells,"y_proj"], cex= mcp_2d_cex, col="black", lwd=0.5, pch=21, bg=cols[metacell_assignment[ccells]])

	if(show_titles)
	{ title(main=sprintf("%s (%d)", meta_field_v, length(ccells)), cex.main = cex_main,
			col.main = c(title_colors[meta_field_v],"black")[ifelse(is.null(metacell_colors),2,1)] ) }
	
    if (!single_plot) {
      dev.off()
    }
  }

  if (single_plot) {
    dev.off()
  }

}

#
#
#
#


# this function gets lfp colors of specific gene, and plot them over scatter of total output coming of 2 gene modules
##################################################################################################################
cmem_plot_gene_modules_scatter <- function (genes_x = NULL,
								            genes_y = NULL,
								            x_label = "lab1",
								            y_label = "lab2",
								            cells_umi_tab = NULL,
								            metacell_assignment = NULL,
								            metacells_fp = NULL,
								            lfp_color_gene = NULL,
								            z_score_range = 2,
											z_colors = rev(RColorBrewer::brewer.pal(n=11, 'RdYlBu')),
								            my_xlim = NULL,
								            my_ylim = NULL,
								            log2_exp = TRUE,
								            eps=1,
								            my_text_cex = 3.6,
								            cells_color = "#808080",
							                plot_w = 450,
							                plot_h = 450,
								            cex_cells = 0.5,
								            cex_metacells = 8,
								            plot_margins = c(4,4,4,4),
								            fig_fn = "./lfp_byGene_xy_scatter.png"
								            )
{

      if(is.null(genes_x) | is.null(genes_y) | (sum(! c(genes_x,genes_y) %in% rownames(cells_umi_tab))>0) ) {
               stop("undefined genes to plot ")
       }
	  if(is.null(metacells_fp) ) {
               stop("undefined meta cell foot-print expression ")
       }
	  if(is.null (lfp_color_gene) | (!lfp_color_gene %in% rownames(cells_umi_tab)) ) {
               stop("undefined gene to color by its lfp in metacells")
       }
	  if(is.null(cells_umi_tab) ) {
               stop("undefined expression umi-table")
       }	   

good_cells = intersect(names(metacell_assignment), colnames(cells_umi_tab))
sc_n = apply(cells_umi_tab[,good_cells], 2, function(x) {(x/sum(x))*1e4} )
sc_x =  colSums(sc_n[genes_x,])
sc_y =  colSums(sc_n[genes_y,])
if (log2_exp) { sc_x = log2(sc_x + eps); sc_y = log2(sc_y + eps) }

mc_total_exp = tapply(  colSums(cells_umi_tab[,good_cells]),metacell_assignment[good_cells],sum)
mc_x_total_exp = tapply(colSums(cells_umi_tab[genes_x,good_cells]),metacell_assignment[good_cells],sum)
mc_y_total_exp = tapply(colSums(cells_umi_tab[genes_y,good_cells]),metacell_assignment[good_cells],sum)
mc_x = (mc_x_total_exp/mc_total_exp)*1e4
mc_y = (mc_y_total_exp/mc_total_exp)*1e4
if (log2_exp) { mc_x = log2(mc_x + eps); mc_y = log2(mc_y + eps) }

png(fig_fn, height=plot_h, width=plot_w)
par(mar = plot_margins)
plot(sc_x,sc_y,
	 cex = cex_cells,
	 col = cells_color,
	 xlim = list(c(min(sc_x),max(sc_x)),my_xlim)[[ifelse(is.null(my_xlim),1,2)]],
	 ylim = list(c(min(sc_y),max(sc_y)),my_ylim)[[ifelse(is.null(my_ylim),1,2)]],
	 cex.axis=my_text_cex,
	 cex.lab=my_text_cex,
	 xlab="",
	 ylab="",
	 xaxt="n",
	 yaxt="n",) 
axis(side = 3,cex.axis = my_text_cex, cex = my_text_cex, line = -0.5, tick = FALSE)
axis(side = 4,cex.axis = my_text_cex, cex = my_text_cex, las = 3,line = 0.7, tick = FALSE)
title(xlab = x_label, cex.lab = my_text_cex, line=2)
title(ylab = y_label, cex.lab = my_text_cex, line=1)
points(mc_x, mc_y,
	   pch=21,
	   cex=cex_metacells,
	   bg = cmem_lfp_col_z(nm = lfp_color_gene, lfp =  log2(metacells_fp), z = z_score_range, col_spectrum = z_colors))
text(mc_x, mc_y, 1:ncol(metacells_fp), cex = my_text_cex)
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This function is color scaling log-footprint (lfp) expression, ranged by z parameter.
##################################################################################################################
cmem_lfp_col_z = function(nm, lfp, z = 4, col_spectrum = rev(RColorBrewer::brewer.pal(n=11, 'RdYlBu'))) 
{ 
x = pmin(pmax(lfp[nm, ], -z), z) + z;
return(colorRampPalette(col_spectrum)(200 * z + 1)[round(100 * as.numeric(x)) + 1]) 
}
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This gets a balanced cell-cell graph, and randomized expression in each cell based on its K neighbors
##################################################################################################################
cmem_mcell_cgraph_norm_gcor <- function(cgraph_edges = NULL,
										cells_umi_tab = NULL,
										K=-1,
										dowsamp_value = 5000,
										min_gtot=1000)
{
        if(is.null(cgraph_edges)) {
                stop("missing cgraph for gene cor normalization, id = ", cgraph_id)
        }
        if(is.null(cells_umi_tab)) {
                stop("missing mat for gene cor normalization, id = ", mat_id)
        }

        mat_ds = cmem_scm_downsamp(umis = cells_umi_tab, 
								   n = dowsamp_value)
        c_nms = colnames(mat_ds)

        adjs = split(cgraph_edges$mc2, cgraph_edges$mc1)
		c_nms = intersect(c_nms,names(adjs))
		
        adjs = adjs[c_nms]
        if(K==-1) {
          K = max(unlist(lapply(adjs, length)))
        }
        m_adjs = t(sapply(adjs, function(x) {
                                                x1 = intersect(x, c_nms);
                                                x[1+(0:(K-1))%%length(x)]
                                                }))
        C = nrow(m_adjs)
        f_gene = rowSums(mat_ds) > min_gtot
        r_mat = t(apply(as.matrix(mat_ds[f_gene,]), 1, function(x) {
                                shuf = m_adjs[(1:C)+C*(sample(1:K,s=C, r=T)-1)];
								x[shuf]
                                }))

        gcor = tgs_cor(as.matrix(t(mat_ds[f_gene,])), spearman = TRUE)
        r_mat[is.na(r_mat)] = 0
        r_gcor = tgs_cor(t(r_mat), spearman = TRUE)
        return(list(gcor=gcor, r_gcor=r_gcor))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



### This is function is plotting max gene corrs after randomization, and used to define cc-independent genes
cmem_plot_ccNorm <- function (cc_cors,
							  base_dir = "./",
							  fig_nm,
							  show_legend_text = TRUE,
							  x_scatter_lim,
							  y_scatter_lim,
							  min_cor,
							  max_r_cor,
							  min_dif)
{
diag(cc_cors$gcor)=0
diag(cc_cors$r_gcor)=0

max_gc = apply(cc_cors$gcor,1,max)
max_rgc = apply(cc_cors$r_gcor,1,max)

interesting_genes = names(max_gc) [max_gc>=min_cor & max_rgc<=max_r_cor & (max_gc-max_rgc) >= min_dif]
cat(paste0(length(interesting_genes)," non-cc genes\n\nplotting..\n"))
#plot
png(fig_nm,width=280,height=280)
par(mar=c(5,5,1,1))
plot(
max_gc,
max_rgc,

cex=.4, cex.lab=1.5,cex.axis=1.5,
xlab = "original max corr.",
ylab = "randomized max corr.",
pch=19,
col="#00000030")
points(max_gc[interesting_genes],
		max_rgc[interesting_genes],
		cex=.8, 
		pch=19,
		col="#0000ff50")
if (show_legend_text){
legend(x=par("usr")[1]*-0.01,
	   y=par("usr")[4]*1.02,
	   xjust=0,
	   bty="n", x.intersp=0.2, y.intersp=0.8,
	   cex=1.5,
	   text.col="blue",
	   legend = "cc-independent\ngenes")
	   }
grid()
dev.off()
return(interesting_genes);
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




### This function compares expression distributions of gene modules in cells and clones
cmem_compare_cells_clones <- function(mod_to_plot,
								      mod_nm,
								      sc_mat,
								      c_mat,
								      minimal_cell,
								      minimal_clone,
								      base_dir = "./",
								      legend_x,
								      legend_y,
								      density_x_lim,
								      density_y_lim,
								      ecdf_x_lim,
								      density_mar = c(4.5,4,0.5,1.5),
								      plot_ecdf=TRUE,
								      show_legend_density=FALSE)
{
cmn_g = rownames(c_mat)[rownames(c_mat) %in% rownames(sc_mat)]
mod_to_plot=mod_to_plot[mod_to_plot %in% cmn_g]
c_mat = c_mat[cmn_g,]
sc_mat = sc_mat[cmn_g,]
c_n = apply(c_mat[,colSums(c_mat)>=minimal_clone],2,function(x){(x/sum(x))*1e5})
sc_n = apply(sc_mat[,colSums(c_mat)>=minimal_clone],2,function(x){(x/sum(x))*1e5})
tot_mod_sc = (sum(sc_mat[mod_to_plot,]) / sum(sc_mat)) * 1e5
tot_mod_c = (sum(c_mat[mod_to_plot,]) / sum(c_mat)) * 1e5

if (plot_ecdf){
png(paste0(base_dir,"density_",mod_nm,".png"),width=270*0.85,height=500*0.85)
par(mar=c(0.5,4,4.5,1.5),mfcol=c(2,1))
plot(ecdf(colSums(c_n[mod_to_plot,])),
	lwd=3.5,col="black",cex.axis=1.5,
	xlim=c(ecdf_x_lim[1],ecdf_x_lim[2]),
	cex.lab=1.5,cex.main=2,
	xaxt="n",
	ylab="",
	main="",
	xlab="")
mtext(side=2,text="ecdf",line=2.4,cex=1.5)
axis(side=1,labels=FALSE)
lines(ecdf(colSums(sc_n[mod_to_plot,])),lwd=2.5,col="red")
legend(x=legend_x,y=legend_y,
		fill=c("black","red"),legend=c("clones","singles"),bty="n",cex=1.5,x.intersp=0.2);

par(mar=density_mar)
plot(density(colSums(c_n[mod_to_plot,])),lwd=3.5,col="black",main="",cex.axis=1.5,
				xlim=c(density_x_lim[1],density_x_lim[2]),
				ylim=c(density_y_lim[1],density_y_lim[2]),
				cex.lab=1.5,cex.main=2,
				ylab="",
				xlab=mod_nm)
mtext(side=2,text="density",line=2.4,cex=1.5)
lines(density(colSums(sc_n[mod_to_plot,])),lwd=2.5,col="red")
abline(v=tot_mod_sc,lty=2,lwd=1.5,col="red")
abline(v=tot_mod_c,lty=2,lwd=1.5,col="black")
dev.off()
} else {
png(paste0(base_dir,"density_",mod_nm,".png"),width=270*0.85,height=500*0.85/2)


par(mar=density_mar)
plot(density(colSums(c_n[mod_to_plot,])),lwd=3.5,col="black",main="",cex.axis=1.5,
				xlim=c(density_x_lim[1],density_x_lim[2]),
				ylim=c(density_y_lim[1],density_y_lim[2]),
				cex.lab=1.5,cex.main=2,
				ylab="",
				xlab="")
title(xlab=mod_nm,cex.lab=1.5,line=2.25)
mtext(side=2,text="density",line=2.4,cex=1.5)
lines(density(colSums(sc_n[mod_to_plot,])),lwd=2.5,col="red")
abline(v=tot_mod_sc,lty=2,lwd=1.5,col="red")
abline(v=tot_mod_c,lty=2,lwd=1.5,col="black")
	if(show_legend_density){
		legend(#x="topright",
			x=legend_x,
			y=legend_y,
			cex=1.5,bty="n",
			x.intersp=0.2,
			y.intersp=0.8,
			legend=c("clones","cells"),fill=c("black","red"))
	}
dev.off()
}
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



### This function compares expression histograms of specific genes in cells and clones
cmem_hist_gene_cells_clones <- function(gene_to_plot = NULL,
									    min_cell = 6e3,
									    max_cell = 4e4,
									    min_clone = 1e4,
									    cells_umi_table = NULL,
									    clones_umi_table = NULL,
										cells_umi_table_norm = NULL,
									    clones_umi_table_norm = NULL,
										log2_exp = FALSE,
										eps = 1,
										cex_text=1.8,
									    base_dir = "./",
										x_label = "lab1",
										plot_h = 273,
										plot_w = 254.8,
										plot_margins = c(5.5,4.6,1,2),
										plot_legend_cells = FALSE,
										plot_legend_clones = FALSE,
										population_legend = NULL,
										legend_x_cells = 1.6,
										legend_y_cells = 400,
										legend_x_clones = 1.6,
										legend_y_clones = 400,
								        hist_breaks = NULL,
									    hist_colors = NULL)
{
sc_exp = cells_umi_table_norm[gene_to_plot,colSums(cells_umi_table) >= min_cell & colSums(cells_umi_table) <= max_cell]
c_exp = clones_umi_table_norm[gene_to_plot,colSums(clones_umi_table) >= min_clone]
if (log2_exp) {sc_exp = log2(sc_exp+eps); c_exp = log2(c_exp+eps)}

png(file = paste0(base_dir,gene_to_plot,"_cells_hist.png"), width=plot_w, height=plot_h)
par(mar = plot_margins)
hist(sc_exp,
	 col=list(rep("grey",100),hist_colors)[[ifelse(is.null(hist_colors),1,2)]],
	 cex.lab=cex_text, cex.axis=cex_text,
	 breaks=list(20,hist_breaks)[[ifelse(is.null(hist_breaks),1,2)]],	 
	 xlab="",
	 ylab="",
	 main="")
title(ylab="# cells",cex.axis=cex_text,cex.lab=cex_text,cex=cex_text,line=3.2)
title(xlab=x_label,cex.axis=cex_text,cex.lab=cex_text,cex=cex_text,line=3.2,las=3)
if (plot_legend_cells)
	{
	legend(x=legend_x_cells,
			y=legend_y_cells,
			x.intersp=0.2,
			y.intersp=0.8,
			fill=list(rep("grey",2),unique(hist_colors))[[ifelse(is.null(hist_colors),1,2)]],
			legend=list(c("pop1","pop2"),population_legend)[[ifelse(is.null(population_legend),1,2)]],
			cex=cex_text,
			bty="n")
	}
dev.off()

png(file = paste0(base_dir,gene_to_plot,"_clones_hist.png"),width = plot_w, height=plot_h)
par(mar=c(5.5,4.6,1,2))
hist(c_exp,
	 col=list(rep("grey",100),hist_colors)[[ifelse(is.null(hist_colors),1,2)]],
	 cex.lab=cex_text, cex.axis=cex_text,
	 breaks=list(20,hist_breaks)[[ifelse(is.null(hist_breaks),1,2)]],	 
	 xlab="",
	 ylab="",
	 main="")
title(ylab = "# clones", cex.axis=cex_text, cex.lab=cex_text, cex=cex_text, line=3.2)
title(xlab = x_label, cex.axis=cex_text, cex.lab=cex_text, cex=cex_text, line=3.2)
if (plot_legend_clones)
	{
	legend(x=legend_x_clones,
			y=legend_y_clones,
			x.intersp=0.2,
			y.intersp=0.8,
			fill=list(rep("grey",2),unique(hist_colors))[[ifelse(is.null(hist_colors),1,2)]],
			legend=list(c("pop1","pop2"),population_legend)[[ifelse(is.null(population_legend),1,2)]],
			cex=cex_text,
			bty="n")
	}	 
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This function get umi_matrix and retrieve gene-gene correlation matrix and gene-gene correlation matrix adjusted
# to sampling depth of each pair of genes
#################################################################################################################
cmem_adjust_gcor_by_samplingDepth <- function (cells_umi_tab = NULL,
											   min_gtot = 50,
											   min_cell = 1000,
											   max_cell = 40000,
											   rollmean_n_genes = 101,
											   cor_method_spearman = TRUE,
											   symmertrize_adj_corr = FALSE,
											   plot_example = FALSE,
											   gene_to_plot = NULL,
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
											   ylim_example1 = NULL,
											   ylim_example2 = NULL,
											   xlim_example1 = NULL
											   )
{

mat_f = cells_umi_tab[rowSums(cells_umi_tab) >= min_gtot, colSums(cells_umi_tab) <= max_cell & colSums(cells_umi_tab) >= min_cell]
g_cor = tgs_cor(t(as.matrix(mat_f)), spearman = cor_method_spearman)
gsum = rowSums(mat_f)
gtot_ord = order( gsum)

gcor_n = apply(g_cor, 1, cmem_norm_gene, rollmean_n_genes, gtot_ord)
if (symmertrize_adj_corr) { gcor_n = (gcor_n + t(gcor_n)) / 2 }


if (plot_example) {
 if (is.null(gene_to_plot)) {gene_to_plot = sample(rownames(mat_f),size=1)}
 
# plot examples for the effecto of gene corr adjustment for specific genes
for (i in gene_to_plot){
	x_gcorn = gcor_n[,i]
	int_genes_down = names(x_gcorn)[x_gcorn < low_corr_t]
	int_genes_up = names(x_gcorn)[x_gcorn > high_corr_t]
	
	
	png(paste0(base_dir,i,"_corTrend_example1.png"), width=fig_w, height=fig_h)                                                                                                                          
	par(mar=c(5,5,1,1))
	plot(log2(gsum), g_cor[,i],
		cex=0.3,pch=19,
		ylim = list(range(log2(gsum)),ylim_example1)[[ifelse(is.null(ylim_example1),1,2)]],
		xlim = list(range(log2(gsum)),xlim_example1)[[ifelse(is.null(xlim_example1),1,2)]],
		cex.lab=cex_text,cex.axis=cex_text,
		col="#00000050",
		ylab=paste0("corr. to ",i),
		xlab="",
		xaxt="n")     	
	points(log2(gsum), g_cor[,i]-x_gcorn,pch=19,cex=1, col=trend_col)                                                                          
	if (length(int_genes_down)>0) {points(log2(gsum)[c(int_genes_down)], g_cor[,i][c(int_genes_down)],
		pch=19, cex=1, col=negative_corr_col)}
	if (length(int_genes_up)>0) { points(log2(gsum)[c(int_genes_up)], g_cor[,i][c(int_genes_up)],
		pch=19, cex=1, col=positive_corr_col)}
	axis(1,labels=FALSE)
	if (length(int_genes_down)>0) {text(log2(gsum)[int_genes_down], g_cor[,i][int_genes_down] + x_text_label_factor, 
								   int_genes_down, col = negative_corr_col, cex=cex_text)}
	if (length(int_genes_up)>0) {text(log2(gsum)[int_genes_up],   g_cor[,i][int_genes_up] + y_text_label_factor,
								 int_genes_up, col = positive_corr_col, cex=cex_text)}
	legend(x="bottomright",
		legend=c("corr. trend"), fill=c(trend_col), cex=cex_text,bty="n", text.col=trend_col,
		x.intersp=0.2)
	grid()
	dev.off()
	
	png(paste0(base_dir,i,"_corTrend_example2.png"), width=fig_w, height=fig_h)
	par(mar=c(5,5,1,1))
	plot(log2(gsum), x_gcorn, col="#00000080", pch=19, cex=0.3, 
		ylim = list(range(x_gcorn),ylim_example2)[[ifelse(is.null(ylim_example2),1,2)]],
			cex.lab = cex_text, cex.axis = cex_text,
			ylab = paste0("adjusted correlation to ",i),
			xlab ="total gene UMIs") 
	if (length(int_genes_up)>0) {points(log2(gsum)[int_genes_up], x_gcorn[int_genes_up], col=positive_corr_col, pch=19,cex=1)}
	if (length(int_genes_down)>0) {points(log2(gsum)[int_genes_down], x_gcorn[int_genes_down], col = negative_corr_col, pch=19, cex=1)}
	grid()
	dev.off() 
	}
}

return(list(gcor = g_cor, gcor_n = gcor_n))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



### this function computes adjusted correlation vector for a gene, 
### based on subtraction of its correlation trend (based expression levels) from actual correlation value for each other gene.
#################################################################################################################
cmem_norm_gene = function(x, w, gtot_ord) 
{
n = length(x)
trend = rollmean(x[gtot_ord], w, na.pad=T)
trend[1:w] = trend[w]
trend[(n-w):n] = trend[n-w]
x[gtot_ord] = x[gtot_ord]-trend
return(x)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# This function shuffling expression of a specific gene, then adjusting its corrs by sampling depth of it and all other genes
#################################################################################################################
cmem_permute_gene_adjust_gcor <- function (cells_umi_tab = NULL,
										   min_gtot = 50,
										   min_cell = 1000,
										   max_cell = 40000,
										   rollmean_n_genes = 101,
										   cor_method_spearman = TRUE,
										   symmertrize_adj_corr = FALSE,
										   gene_to_permute = NULL)
{
mat_f = cells_umi_tab[rowSums(cells_umi_tab) >= min_gtot, colSums(cells_umi_tab) <= max_cell & colSums(cells_umi_tab) >= min_cell]
mat_f_perm = mat_f; mat_f_perm[gene_to_permute,] = sample(mat_f[gene_to_permute,], replace = FALSE)

g_cor = tgs_cor(t(as.matrix(mat_f_perm)), spearman = cor_method_spearman)
gsum = rowSums(mat_f_perm)
gtot_ord = order( gsum)

gcor_n = apply(g_cor, 1, cmem_norm_gene, rollmean_n_genes, gtot_ord)
if (symmertrize_adj_corr) { gcor_n = (gcor_n + t(gcor_n)) / 2 }

return(list(perm_gcor = g_cor, perm_gcor_n = gcor_n))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This function plots highest and lowest correlation to a target genes, and barplots over permuted correlations to it,
# also plot distribution of real and shuffled corrs to this gene
#################################################################################################################
cmem_barplot_top_gcor <- function (adj_gcor_mat = NULL,
								   perm_adj_gcor_mat = NULL,
								   gene_name = NULL,
								   n_genes = 10,
								   positive_corr_col = "blue",
								   negative_corr_col = "red",
								   positive_permuted_corr_col = "lightblue",
								   negative_permuted_corr_col = "tomato",
								   positive_corr_margins = c(14,4.5,7,4),
								   negative_corr_margins = c(7,4.5,14,4),
								   ylim_positive_corr = NULL,
								   ylim_negative_corr = NULL,
								   ylim_density = NULL,
								   xlim_ecdf_positive = NULL,
								   ylim_ecdf_positive = NULL,
								   xlim_ecdf_negative = NULL,
								   ylim_ecdf_negative = NULL,
								   text_cex=1.2,
								   barplot_w = 800,
								   barplot_h = 400,
								   density_plot_w = 250,
								   density_plot_h = 250,
								   ecdf_plot_w = 220,
								   ecdf_plot_h = 220,
								   base_dir = "./"
								   )
{
adj_gcor = adj_gcor_mat[,gene_name]
perm_gcor = perm_adj_gcor_mat[,gene_name]

corrs_u= rev(tail(sort(adj_gcor),n_genes+1))[-1]
corrs_d= head(sort(adj_gcor),n_genes)
corrs_u_perm = perm_gcor[names(corrs_u)]
corrs_d_perm = perm_gcor[names(corrs_d)]
corrs_u_tab = rbind(corrs_u,corrs_u_perm)
corrs_d_tab = rbind(rev(corrs_d),rev(corrs_d_perm))

vector_cors_orig = adj_gcor [!names(adj_gcor) %in% gene_name]
vector_cors_perm = perm_gcor [!names(perm_gcor) %in% gene_name]

ecdf_orig=ecdf(vector_cors_orig)
uniq_orig=sort(unique(vector_cors_orig))
ecdf_perm=ecdf(vector_cors_perm)
uniq_perm=sort(unique(vector_cors_perm))

# a) barplots
	png(paste0(base_dir,gene_name,"top_cors_permute.png"),width=barplot_w,height=barplot_h)
	par(mar = positive_corr_margins,mfcol=c(1,2))
	#par(mfcol=c(1,2))
	barplot(corrs_u_tab,
		ylim = list(range(corrs_u_tab),ylim_positive_corr)[[ifelse(is.null(ylim_positive_corr),1,2)]],
		cex.names=1.6,xpd=FALSE,cex.axis=1.8,col=rep(c(positive_corr_col,positive_permuted_corr_col),50),
		space=rep(c(0,1),nrow(corrs_u_tab)-1),cex.lab=2,las=2,beside=TRUE,yaxt="n")
	axis(2,col=positive_corr_col,col.ticks=positive_corr_col,col.axis=positive_corr_col,cex.axis=1.8,las=2)
	par(mar = negative_corr_margins)
	barplot(corrs_d_tab,
		ylim = list(range(corrs_d_tab),ylim_negative_corr)[[ifelse(is.null(ylim_negative_corr),1,2)]],
		cex.names=1.6,xpd=FALSE,cex.axis=1.8,col=rep(c(negative_corr_col,negative_permuted_corr_col),50),
		space=rep(c(0,1),nrow(corrs_d_tab)-1),cex.lab=2,las=2,xlab=NA,
		names=rep("",ncol(corrs_d_tab)),
		beside=TRUE,
		yaxt="n")
	axis(3,at=seq(1.5,(n_genes*3)-0.5,3),labels=colnames(corrs_d_tab),las=2,tick=FALSE,cex.axis=1.6)
	axis(2,col=negative_corr_col,col.ticks=negative_corr_col,col.axis=negative_corr_col,las=2,cex.axis=1.8)
	dev.off()

# b) density
png(paste0(base_dir,gene_name,"_corrs_permut.png"),width=density_plot_w,height=density_plot_h)
par(mar=c(5,5,1,1))
plot(density(vector_cors_orig),col="black",lwd=4,
		ylim = list(range((vector_cors_orig)),ylim_density)[[ifelse(is.null(ylim_density),1,2)]],
xlab=paste0("corr. to ",gene_name),
cex.axis=text_cex,cex.lab=text_cex,main="")
lines(density(vector_cors_perm),col="grey",lwd=4)
grid()
legend(x="topleft",
	   cex=text_cex,
	   bty="n",
	   fill=c("grey","black"),
	   legend=c("permuted","original"),
	   x.inters=0.2,y.intersp=0.8)
dev.off()

# c) ecdf - positive tail
png(sprintf("%s%s_cors_permutation_ecdf_positive.png",base_dir,gene_name),height=ecdf_plot_h,width=ecdf_plot_w)
par(mar=c(5,5,0.1,0.1))
plot(rev(uniq_orig),((1-ecdf_orig(rev(uniq_orig)))*length(vector_cors_orig)),
	xlab="corr. to EpCAM",
	ylab="# genes",
	yaxt="n",
	main="",
	type="l",
	xlim = list(range(rev(uniq_orig)),xlim_ecdf_positive)[[ifelse(is.null(xlim_ecdf_positive),1,2)]],
	ylim = list(range(((1-ecdf_orig(rev(uniq_orig)))*length(vector_cors_orig))),ylim_ecdf_positive)[[ifelse(is.null(ylim_ecdf_positive),1,2)]],
	cex.axis=text_cex,cex.lab=text_cex,
	lwd=4,col=positive_corr_col)	
axis(2,col="black",col.ticks="black",col.axis="black",cex.axis=text_cex,cex.lab=text_cex,cex=text_cex,las=2)
lines(rev(uniq_perm),((1-ecdf_perm(rev(uniq_perm)))*length(vector_cors_orig)),lwd=4,col=positive_permuted_corr_col)
grid(nx=NULL,ny=NULL,col = "lightgrey", lty = "dotted",lwd=par("lwd"))
legend(x="topright",
	   cex=text_cex,
	   bty="n",
	   fill=c(positive_permuted_corr_col,positive_corr_col),
	   legend=c("permuted","original"),
	   x.inters=0.2,y.intersp=0.8)
dev.off()


# d) ecdf - negative tail
png(sprintf("%s%s_cors_permutation_ecdf_negative.png",base_dir,gene_name),height=ecdf_plot_h,width=ecdf_plot_w)
par(mar=c(5,5,0.1,0.1))
plot(uniq_orig,((1-ecdf_orig(rev(uniq_orig)))*length(vector_cors_orig)),
	xlab=paste0("corr. to ",gene_name),
	ylab="# genes",
	yaxt="n",
	main="",
	type="l",
	xlim = list(range(rev(uniq_orig)),xlim_ecdf_negative)[[ifelse(is.null(xlim_ecdf_negative),1,2)]],
	ylim = list(range(((1-ecdf_orig(rev(uniq_orig)))*length(vector_cors_orig))),ylim_ecdf_negative)[[ifelse(is.null(ylim_ecdf_negative),1,2)]],
	cex.axis=text_cex,cex.lab=text_cex,
	lwd=4,col=negative_corr_col)	
axis(2,col="black",col.ticks="black",col.axis="black",cex.axis=text_cex,cex.lab=text_cex,cex=text_cex,las=2)
lines(uniq_perm,((1-ecdf_perm(rev(uniq_perm)))*length(vector_cors_orig)),lwd=4,col=negative_permuted_corr_col)
grid(nx=NULL,ny=NULL,col = "lightgrey", lty = "dotted",lwd=par("lwd"))
legend(x="topleft",
	   cex=text_cex,
	   bty="n",
	   fill=c(negative_permuted_corr_col,negative_corr_col),
	   legend=c("permuted","original"),
	   x.inters=0.2,y.intersp=0.8)
dev.off()

}
#################################################################################################################
#################################################################################################################
#################################################################################################################






# This function downsample a umi table to ds_val
#################################################################################################################
cmem_ds <- function (umi_tab, ds_val)
{
	.downsamp_one=function(v,n)
	{
	hist(sample(rep(1:length(v),times=v),replace=F,size=n),0.5+0:length(v),plot=F)$counts
	}
umi_tab_filt = umi_tab[,colSums(umi_tab)>=ds_val]
umi_tab_ds = apply(umi_tab_filt, 2, .downsamp_one, ds_val)
rownames(umi_tab_ds)=rownames(umi_tab);
return(umi_tab_ds)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This is the MetaCell package downsampling function
##################################################################################################################
cmem_scm_downsamp <- function(umis, n)
{
        umis = umis[,colSums(umis)>= n]
        m = nrow(umis)
        .downsamp_one=function(v,n, replace = F){
                a = tabulate(sample(rep(1:length(v),times=v),replace=replace,size=n),nbins=m)
                return (a)
        }
        max_bin = 16
        doMC::registerDoMC(max_bin)

        max_bin = min(max_bin, ceiling(ncol(umis)/500))

        if(max_bin*10000 < ncol(umis)) {
                max_bin =  round(ncol(umis))/10000
        }
        cell_quant = ceiling(ncol(umis)/max_bin)
        seed = 19
        sub_dsamp = function(x) {
                set.seed(seed)
                i = 1+(x-1)*cell_quant
                j = min(x*cell_quant, ncol(umis))
           ret = Matrix(apply(umis[,i:j], 2, .downsamp_one, n))
           rownames(ret) = rownames(umis)
                return(as(ret,"dgCMatrix"))
        }
        res <- plyr::alply(1:max_bin, 1, sub_dsamp, .parallel=TRUE)
        umis_ds = do.call(cbind, res)
        return(umis_ds)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# This function plots expression sccater of gene or genes over longitudinal tracking of HCT116 cells of 6 single-cell-derived clones
##################################################################################################################
cmem_plot_longitudinal_cells_xy <- function (cells_metadata = NULL,
								 cells_umi_tab_norm = NULL,
								 genes_a = NULL,
								 genes_b = NULL,
								 dot_cex = 0.3,
								 text_cex = 12,
								 lab_a = "lab1",
								 lab_b = "lab2",
								 base_dir = "./",
								 log2_exp_a = FALSE,
								 log2_exp_b = FALSE,
								 ggplot_w = 4.123,
								 ggplot_h = 5.15375,
								 eps = 1)
{
      if(is.null(genes_a) | is.null(genes_b) | (sum(! c(genes_a,genes_b) %in% rownames(cells_umi_tab_norm))>0) ) {
               stop("undefined genes to plot ")
       }

      if(is.null(cells_umi_tab_norm) ) {
               stop("undefined normalized umi-table")
       }
	   
	   if(is.null(cells_metadata) ) {
               stop("undefined longitudinal metadata")
       }	   
	   
if (length(genes_a)==1) { exp_a = cells_umi_tab_norm[genes_a,rownames(cells_metadata)] } else {
						  exp_a = colSums(cells_umi_tab_norm[genes_a,rownames(cells_metadata)]) }
							
if (length(genes_b)==1) { exp_b = cells_umi_tab_norm[genes_b,rownames(cells_metadata)] } else {
						  exp_b = colSums(cells_umi_tab_norm[genes_b,rownames(cells_metadata)]) }

if (log2_exp_a) { exp_a = log2(exp_a + eps) }
if (log2_exp_b) { exp_b = log2(exp_b + eps) }

lt_df = data.frame(clone = cells_metadata[,"clone"],day = cells_metadata[,"day"],color = cells_metadata[,"color"],
				   a = exp_a, b = exp_b)

#browser()
lt_scatter = 
	ggplot(data = lt_df, aes(x = a, y = b, color = color)) +
    geom_point(size = dot_cex) + 
	background_grid() + 
	xlab(lab_a) + 
	ylab(lab_b) + 
    facet_grid(clone ~ day, scales = "fixed") + 
	 scale_colour_manual( values = levels((lt_df$color)) ) + 
	theme(text = element_text(size=text_cex),
		 axis.title.x = element_text(size = text_cex, face = "bold"),
		 axis.title.y = element_text(size = text_cex, face = "bold"),
		 legend.position = "none")
ggsave(filename=paste0(base_dir,lab_a,"_",lab_b,"_longitudinal_XY.png"),plot = lt_scatter,
	   width=ggplot_w, height=ggplot_h)	 		 
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This function plots averaged expression of gene or gene module over 6 time-points 
# during longitudinal tracking of HCT116 cells of 6 single-cell-derived clones
##################################################################################################################
cmem_plot_longitudinal_lines <- function (cells_metadata = cmem_sumi_lt_PAR_metadata,
										  cells_umi_tab = cmem_sumi_lt_PAR,
										  clones_umi_tab1 = cmem_umi_lt_PAR_day10,
										  clones_umi_tab2 = cmem_umi_lt_PAR_day18,
										  genes_to_plot = NULL,
										  line_lwd = 2,
										  lab_y = "lab_y",
										  base_dir = "./",
										  gg_text_sz = 18,
										  clones_line_colors = c("#B79F00","#F564E3","#00BFC4","#F8766D","#619CFF","#00BA38"),
										  error_bar_width = .2,
										  erro_bar_lwd = 1,
										  ggplot_w = 4,
										  ggplot_h = 2.5,
										  gg_legend_position = "none",
										  gg_margins = c(0,0,0,0)
										  )
{
      if( is.null(genes_to_plot) | (sum(!genes_to_plot %in% rownames(cells_umi_tab))>0) ) {
               stop("undefined genes to plot ")
       }
	   
cells_metadata$day = factor(cells_metadata$day,levels = c("d_10","d_18","d_33","d_62","d_98","d_148"))

## get gene UMIs and total UMI counts only from paired clones of day10 and day 18
colnames(clones_umi_tab1)=substring(colnames(clones_umi_tab1),5)
colnames(clones_umi_tab2)=substring(colnames(clones_umi_tab2),5)
if (length(genes_to_plot)==1) { c_exp_day10 = clones_umi_tab1[genes_to_plot,] } else {
							    c_exp_day10 = colSums(clones_umi_tab1[genes_to_plot,]) } 
if (length(genes_to_plot)==1) { c_exp_day18 = clones_umi_tab2[genes_to_plot,] } else {
							    c_exp_day18 = colSums(clones_umi_tab2[genes_to_plot,]) } 
c_target_umis = cbind(c_exp_day10[c("1d12","4e1","7a2","3b3","4b10","7b11")],
						c_exp_day18[c("1d12","4e1","7a2","3b3","4b10","7b11")])					
c_total_umis = cbind(colSums(clones_umi_tab1[,c("1d12","4e1","7a2","3b3","4b10","7b11")]),
					 colSums(clones_umi_tab2[,c("1d12","4e1","7a2","3b3","4b10","7b11")]))
rownames(c_target_umis) = levels(cells_metadata$clone)
rownames(c_total_umis) = levels(cells_metadata$clone)
colnames(c_target_umis) = c("d_10","d_18")
					 
## get gene UMIs and total UMI counts only from longitudinal single cells
if (length(genes_to_plot)==1) { sc_exp = cells_umi_tab[genes_to_plot,rownames(cells_metadata)] } else {
							    sc_exp = colSums(cells_umi_tab[genes_to_plot,rownames(cells_metadata)]) } 
if (length(genes_to_plot)==1) {cells_df = cbind(cells_metadata,t(sc_exp),colSums(cells_umi_tab[,rownames(cells_metadata)])) } else {
							   cells_df = cbind(cells_metadata,sc_exp,colSums(cells_umi_tab[,rownames(cells_metadata)])) }
colnames(cells_df)[ (ncol(cells_df)-1) : ncol(cells_df)] = c("target_umi","total_umi")
sc_target_umis = xtabs( formula = target_umi ~ clone + day, data = cells_df)[,c("d_33","d_62","d_98","d_148")]
sc_total_umis = xtabs( formula = total_umi ~ clone + day, data = cells_df)[,c("d_33","d_62","d_98","d_148")]

## normalize UMI counts of both cells and clones and plot
c_norm_expression = (c_target_umis / c_total_umis) * 1e5
sc_norm_expression = (sc_target_umis / sc_total_umis) * 1e5
combined_norm_expression = cbind(c_norm_expression,sc_norm_expression)
combined_total_umis = cbind(c_total_umis,sc_total_umis)
se = sqrt(( (combined_norm_expression/1e5)* (1-(combined_norm_expression/1e5)) ) / combined_total_umis) * 1e5
se_m = melt(se)[,3]
plot_df = cbind(melt(combined_norm_expression,varnames=c("clone","day"),value.name="target_genes"),se_m)

lt_lines_plot = 
		ggplot(plot_df ,aes(x = day, y = target_genes, colour = clone)) + 
		scale_colour_manual(values = clones_line_colors) +
 		geom_errorbar( aes (ymin = target_genes - se_m, ymax = target_genes + se_m), width=error_bar_width, size=erro_bar_lwd) +
 		ylab(lab_y) + 
 		xlab("") + 
 		theme(legend.position = gg_legend_position,
			  axis.text.x = element_text(size=gg_text_sz,angle = 90, hjust = 1),
			  plot.title = element_text(size=gg_text_sz),
			  legend.text = element_text(size=gg_text_sz),
			  axis.title.x = element_text(size=gg_text_sz),
			  axis.text.y = element_text(size=gg_text_sz),
			  axis.title.y = element_text(size=gg_text_sz),
			  plot.margin = unit(gg_margins, "cm"))+
			  geom_line(aes(group=clone), size=line_lwd) + 
		guides(fill = guide_legend(row=3)) 

 ggsave(lt_lines_plot,file=paste0(base_dir,gsub(" ","_",lab_y),"_longitudinal_lines.png"),w=ggplot_w,h=ggplot_h)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################





# This function gets CpG density table, as well as methylation and coverage table, and plot clonal methylation
##################################################################################################################
cmem_plot_clonalMeth_by_cgCont <- function (meth_table = NULL,
										    cov_table = NULL,
										    cg500_table = NULL,
										    clones_metadata = NULL,
										    dots_cex = 1,
										    text_cex = 1.5,
										    cpg_range_x = NULL,
										    cpg_range_y = NULL,
										    lab_x = "XLAB",
										    lab_y = "YLAB",
										    include_low = FALSE,
										    plot_margins = c(5,5,1,1),
										    plot_h = 280,
										    plot_w = 280,
										    dots_color = "#00000060",
										    fig_nm = NULL,
										    highlight_EMT_clones = FALSE,
											exclude_EMT_clones = FALSE,
										    legend_x = 0.455,
										    legend_y = 0.84,
										    show_EMT_clones_legend = FALSE,
										    EMT_clones_color = "red",
											return_clonal_averages = FALSE
										    )
{

      if( is.null(fig_nm) ) {
               stop("undefined file name to save the plot to")
       }
	   
	   if( is.null(clones_metadata) ) {
               stop("undefined metadata table")
       }
	   if ( (nrow(cg500_table) != nrow(meth_table)) | (nrow(cg500_table) != nrow(cov_table)) ) {
               stop("incompatible CG-content and meth/coverage tables!")	   
	   }
	   
	         if( highlight_EMT_clones & exclude_EMT_clones ) {
               stop("sorry, I can't exclude and highlight the same thing")
       }
	   
cg = cg500_table[,4]
if (!include_low) {
x_avg = colSums(meth_table[cg > cpg_range_x[1] & cg <= cpg_range_x[2],],na.rm=TRUE) /
		colSums(cov_table [cg > cpg_range_x[1] & cg <= cpg_range_x[2],],na.rm=TRUE) ; 
y_avg = colSums(meth_table[cg > cpg_range_y[1] & cg <= cpg_range_y[2],],na.rm=TRUE) /
		colSums(cov_table [cg > cpg_range_y[1] & cg <= cpg_range_y[2],],na.rm=TRUE) ; } else {
x_avg = colSums(meth_table[cg >= cpg_range_x[1] & cg < cpg_range_x[2],],na.rm=TRUE) /
		colSums(cov_table [cg >= cpg_range_x[1] & cg < cpg_range_x[2],],na.rm=TRUE) ; 
y_avg = colSums(meth_table[cg >= cpg_range_y[1] & cg < cpg_range_y[2],],na.rm=TRUE) /
		colSums(cov_table [cg >= cpg_range_y[1] & cg < cpg_range_y[2],],na.rm=TRUE) ; }
				
if (exclude_EMT_clones) {
x_avg = x_avg [!names(x_avg) %in% rownames(clones_metadata)[clones_metadata$VIM_state == "VIM-high"]]
y_avg = y_avg [!names(y_avg) %in% rownames(clones_metadata)[clones_metadata$VIM_state == "VIM-high"]] }

xy_color = rep(dots_color,length(x_avg))
names(xy_color) = names(x_avg)
if (highlight_EMT_clones) { xy_color[rownames(clones_metadata)[clones_metadata$VIM_state == "VIM-high"]] = EMT_clones_color }

png(fig_nm, width=plot_w, height=plot_h)
par(mar = plot_margins)
plot(x_avg, y_avg,
	 pch=19,
	 cex = dots_cex,
	 col = xy_color,
	 ylab = lab_x,
	 xlab = lab_y,
	 cex.lab = text_cex, cex.axis = text_cex); grid();
if (show_EMT_clones_legend){
	legend(x=legend_x,y=legend_y,fill=c(dots_color,EMT_clones_color),
	  legend=c("normal","VIM-high"),cex=text_cex,bty="n",x.intersp=0.2,y.intersp=0.8)}
dev.off()

if (return_clonal_averages){
return(list(x_avg, y_avg))}
}
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This function gets CpG density table, as well as methylation and coverage table, and retrieve both avg and permuted average
##################################################################################################################
cmem_permute_by_cgCont <- function (meth_table = NULL,
									cov_table = NULL,
									cg500_table = NULL,
								    clones_metadata = NULL,
									cpg_range = NULL,
									exclude_EMT_clones = TRUE,
									include_low = FALSE
									)
{
	   if ( (nrow(cg500_table) != nrow(meth_table)) | (nrow(cg500_table) != nrow(cov_table)) ) {
               stop("incompatible CG-content and meth/coverage tables!")	   
	   }

if (exclude_EMT_clones) {
meth_table = meth_table[,!colnames(meth_table) %in% rownames(clones_metadata)[clones_metadata$VIM_state == "VIM-high"]]	   
cov_table = cov_table[,!colnames(cov_table) %in% rownames(clones_metadata)[clones_metadata$VIM_state == "VIM-high"]]	   
}

cg = cg500_table[,4]
if (!include_low) {
meth_table_f = meth_table[which(cg > cpg_range[1] & cg <= cpg_range[2]),]
cov_table_f =  cov_table[which(cg > cpg_range[1] & cg <= cpg_range[2]),] } else {
meth_table_f = meth_table[which(cg >= cpg_range[1] & cg < cpg_range[2]),]
cov_table_f =  cov_table[which(cg >= cpg_range[1] & cg < cpg_range[2]),] } 

cat ("randomly assigning CpGs (in the defined CG-density range) to clones..\n")
n_clones = ncol(meth_table_f)
cgs_all = cbind(meth_table_f, cov_table_f)
cgs_all_perm = t(apply(cgs_all, 1, function(x){tt = sample(1:n_clones); 
						  return( c(x[tt],x[tt+n_clones]))}) )
meth_perm = cgs_all_perm[, 1:n_clones]
cov_perm = cgs_all_perm[, (n_clones+1):ncol(cgs_all)]
colnames(meth_perm) = colnames(meth_table_f)
colnames(cov_perm) = colnames(cov_table_f)

cat ("returning original data and permuted clonal average methylation..\n")
avg_orig = colSums(meth_table_f, na.rm=TRUE) / colSums(cov_table_f, na.rm=TRUE)
avg_perm = colSums(meth_perm, na.rm=TRUE) / colSums(cov_perm, na.rm=TRUE)

return(list( avg_orig = avg_orig, avg_perm = avg_perm))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This function plots variation in average clonal methylation over averages obtained from permuted matrix
##################################################################################################################
cmem_compare_distributions <- function (original = LCG_mat$avg_orig,
							            permuted = LCG_mat$avg_perm,
							            xlim = NULL,
							            original_color = "black",
							            permuted_color = "grey63",
							            plot_margins = c(5,5,1,1),
							            plot_w = NULL,
							            plot_h = NULL,
							            plot_lwd = 4,
							            legend_lwd = 2,
							            plot_legend = FALSE,
							            legend_labels = c("original","perm"),
							            plot_ylim = c(0,280),
							            x_lab = "clonal avg. meth",
							            cex_text = 1.5,
							            fig_nm = NULL)
{
      if( is.null(fig_nm) ) {
               stop("undefined file name to save the plot to")
       }


png(fig_nm, width = plot_w, height = plot_h)
par(mar = plot_margins)
plot(density(original),col="black",lwd = plot_lwd,
	 ylim = plot_ylim,
	 xlab = x_lab,
	 xaxt = "n",
	 cex.axis = cex_text, cex.lab = cex_text, main="")
axis(1, cex.axis=cex_text, cex=cex_text, mgp=c(3, .7, 0))
lines(density(permuted),col=permuted_color, lwd=plot_lwd)
grid()
if (plot_legend) {
legend("topleft", 
 legend=rev(legend_labels),
 col=rev(c(original_color,permuted_color)),
 x.intersp=0.2,
 y.intersp=0.8,
 seg.len=0.9,
 text.col=rev(c(original_color,permuted_color)),
 bty="n", lwd=legend_lwd, cex=cex_text)
}
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################

# this function aggregate CpG methylation stats by given gemonic intervals
#################################################################################################################
cmem_aggregate_CpG_stats_by_intervals <- function (cpg_table = NULL,
												   intervals = NULL,
												   maximal_distance = 100)
{
cg_dist_to_intervs = gintervals.neighbors(cpg_table[,c("chrom","start","end")],intervals)
if (max(cg_dist_to_intervs$dist) > maximal_distance) {
	stop ("at least one CpG is outside the intervals range")}
colnames(cg_dist_to_intervs)[4:6] = c("chrom1","start1","end1")

cgs_to_intervals_convert = cbind(cpg_table[,4:ncol(cpg_table)],cg_dist_to_intervs[,c("chrom1","start1","end1")])
cat (paste0("maximal distance of CpG to interval is ", max(cg_dist_to_intervs$dist)," bases\n"))
cat ("summing CpG stats over intervals..\n ")
cgs_agg_by_intervals = aggregate(.~chrom1+start1+end1,cgs_to_intervals_convert,sum)
intervals_table = cgs_agg_by_intervals[order(cgs_agg_by_intervals$chrom1,cgs_agg_by_intervals$start1),
										   c("chrom1","start1","end1",
										   setdiff(colnames(cgs_agg_by_intervals),c("chrom1","start1","end1")))]
rownames(intervals_table) = paste0(intervals_table$chrom1,"_",intervals_table$start1,"_",intervals_table$end1)

return(intervals_table)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# this function compares and plots average methylation of pooled single cells and clones, filtered by minimal coverage in both
#################################################################################################################
cmem_compare_cells_clones_methylation <- function (min_cov_loci_cells = 50,
												   min_cov_loci_clones = 500,
												   clones_meth_df = NULL,
												   cells_meth_df = NULL,
												   meth_bins = 10,
												   base_dir = "./",
												   boxplot_w = 280,
												   boxplot_h = 280,
												   boxplot_margins = c(5,5,1,1),
												   boxplot_text_cex = 1.5,
												   density_h = 280,
												   density_w = 280,
												   density_margins = c(5,5,1,1),
												   density_text_cex = 1.5,
												   density_lwd = 3,
												   clones_color = "antiquewhite2",
												   cells_color = "black",
												   density_show_legend = TRUE
												   )
{

good_loci = clones_meth_df$cov >= min_cov_loci_clones & cells_meth_df$cov >= min_cov_loci_cells
avg_clones = clones_meth_df$meth[good_loci] / clones_meth_df$cov[good_loci]
avg_cells = cells_meth_df$meth[good_loci] / cells_meth_df$cov[good_loci]

# bin regions by average methylation of single cells, then plot clonal average stratified by it
avg_cells_bins = .bincode(avg_cells,breaks=c(seq(0,1,length.out = meth_bins + 1)),right=TRUE,include.lowest=TRUE)
clones_avg_df = data.frame(avg = avg_clones, singles_bin = as.factor(avg_cells_bins))
png(paste0(base_dir,"clonal_average_over_cells_bins.png"),
	width=boxplot_w,
	height=boxplot_h)
par(mar = boxplot_margins)
boxplot(clones_avg_df$avg ~ clones_avg_df$singles_bin,
		col=clones_color,
		outline=FALSE,
		names = c(1:meth_bins),
		space = rep(2,meth_bins-1),
		yaxt="n",
		cex.axis=boxplot_text_cex,cex.main=boxplot_text_cex,cex.lab=boxplot_text_cex,las=1,
		cex.names=boxplot_text_cex,
		xlab="",
		ylab="")
axis(2, las=2, cex.axis=boxplot_text_cex, cex.lab=boxplot_text_cex, cex=boxplot_text_cex)
title(ylab="clones meth.", line=3.2, cex.lab=boxplot_text_cex, cex.axis=boxplot_text_cex, cex=boxplot_text_cex, outer=FALSE)
title(xlab="cells meth. (bins)", line=2.8, cex.lab=boxplot_text_cex, cex.axis=boxplot_text_cex, cex=boxplot_text_cex, outer=FALSE)
grid(nx=NA,ny=NULL)
dev.off()

# compare distributions of average methylation in covered regions
png(paste0(base_dir,"compare_pools_clones_cells_averageMeth.png"), width = density_w, height = density_h)
par(mar = density_margins)
plot(density(avg_cells),
	 lwd = density_lwd,
	 main = "",
	 cex.lab = density_text_cex, cex.axis = density_text_cex, col= cells_color,
	 xlab= "")
title(xlab = "meth. (capture regions) ", 
	  cex.axis = density_text_cex, cex.lab = density_text_cex, cex = density_text_cex, line=3.5)
lines(density(avg_clones), lwd = density_lwd,
	  col=clones_color)
if (density_show_legend) {
legend(x="topleft",
	   fill=c(clones_color,cells_color),
	   legend=c("clones","cells"), cex = density_text_cex,
	   bty="n",
	   x.intersp=0.2)}
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################


# this function plots distribution of correlations between gene expression and offtarget HCG and LCG methylation in clones
# if return_gcors is TRUE, retrieves correlation to a list object
#################################################################################################################
cmem_compute_methTrend_geneCorr <- function (clones_offTarget_methylation_table = NULL,
											 clones_umi_table = NULL,
											 clones_umi_table_norm = NULL,
											 min_UMI_per_gene = 5,
											 plot_h = 280,
											 plot_w = 280,
											 plot_text_cex = 1.5,
											 plot_lwd = 4,
											 ylim_LCG = c(0,7),
											 ylim_HCG = c(0,6.5),
											 plot_margins = c(5,5,1,1),
											 base_dir = "./",
											 color_original = "black",
											 color_permuted = "grey63",
											 plot_legend = FALSE,
											 legend_lwd = 2,
											 cor_method = "spearman",
											 blacklist_genes = NULL,
											 return_gcors = FALSE
											 )
{
		if (is.null(clones_offTarget_methylation_table)) {
				stop ("undefined table of clonal offtarget methylation HCG and LCG trends") }

		if (is.null(clones_umi_table) | is.null(clones_umi_table_norm)) {
				stop ("undefined clonal UMI tables") }

exp_n = clones_umi_table_norm[rowSums(clones_umi_table[,rownames(clones_offTarget_methylation_table)]) >= min_UMI_per_gene, rownames(clones_offTarget_methylation_table)]
if (!is.null(blacklist_genes)) { exp_n = exp_n[setdiff(rownames(exp_n),blacklist_genes),] }
cat("computing gene correlations to clonal LCG methylation..\n")
LCG_gcor = apply(exp_n, 1, function(x) { cor(x, clones_offTarget_methylation_table$LCG, method = cor_method) } )
LCG_perm_gcor = apply(exp_n, 1, function(x) { cor(x, clones_offTarget_methylation_table$LCG_perm, method = cor_method) } )
cat("computing gene correlations to clonal HCG methylation..\n")
HCG_gcor = apply(exp_n, 1, function(x) { cor(x, clones_offTarget_methylation_table$HCG, method = cor_method) } )
HCG_perm_gcor = apply(exp_n, 1, function(x) { cor(x, clones_offTarget_methylation_table$HCG_perm, method = cor_method) } )

# plot LCG-gene expression correlations distribution over permuted clonal LCG
png(paste0(base_dir,"LCG_gcors_distribution.png"),width=plot_w,height=plot_h)
par(mar = plot_margins)
plot(density(LCG_gcor),col="black",lwd = plot_lwd, ylim = ylim_LCG,
	 xlab="corr. to Low CpG cont.",
	 cex.axis=plot_text_cex, cex.lab=plot_text_cex, main="")
lines(density(LCG_perm_gcor), col=color_permuted, lwd=plot_lwd)
grid()
	if (plot_legend) {
		legend(x="topleft",
		col=rev(c(color_original,color_permuted)),
		cex = plot_text_cex,
		bty = "n",
		lwd = legend_lwd,
		seg.len = 0.9,
		text.col = c(color_permuted,color_original),
		legend=c("permuted","original"),
		x.intersp=0.2, y.intersp=0.8) }
dev.off()

# plot HCG-gene expression correlations distribution over permuted clonal HCG
png(paste0(base_dir,"HCG_gcors_distribution.png"),width=plot_w,height=plot_h)
par(mar = plot_margins)
plot(density(HCG_gcor),col="black",lwd = plot_lwd, ylim = ylim_HCG,
	 xlab="corr. to High CpG cont.",
	 cex.axis = plot_text_cex, cex.lab = plot_text_cex, main="")
lines(density(HCG_perm_gcor), col=color_permuted, lwd=plot_lwd)
grid()
	if (plot_legend) {
		legend(x="topleft",
		cex = plot_text_cex,
		bty = "n",
		lwd = legend_lwd,
		seg.len = 0.9,
		col=rev(c(color_original,color_permuted)),
		text.col=rev(c(color_original,color_permuted)),
		legend=c("permuted","original"),
		x.intersp=0.2, y.intersp=0.8) }
dev.off()

if (return_gcors) {
 return(list(HCG_gcor = HCG_gcor, LCG_gcor = LCG_gcor)) }
}											 
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This function gets clonal methylation and coverage per capture region, and retrieves and average-methylation table,
# only for regions with minimal coverage (also in transcriptome) and minimal methylation variance
#################################################################################################################
cmem_get_regional_average <- function (meth_table = NULL,
									   cov_table = NULL,
									   clones_umi_table = NULL,
									   minimal_regional_coverage_per_clone = 10,
									   minimal_regional_convering_clones = 100,
									   maximal_avg_methylation = 0.98,
									   minimal_avg_methylation = 0.0001,
									   minimal_umi_per_clone = 1e4
									   )
{
# filter by transcriptome and methylome coverage
rna_covered_clones = colnames(clones_umi_table)[colSums(clones_umi_table)>=minimal_umi_per_clone]
rna_and_onTarget_methylation_covered_clones = intersect(rna_covered_clones,colnames(cov_table))
mehtylation_cov_f = apply(cov_table[,rna_and_onTarget_methylation_covered_clones],1,function(x){sum(x >= minimal_regional_coverage_per_clone)})
methylation_covered_regions = mehtylation_cov_f >= minimal_regional_convering_clones
# filter by minimal and maximal average methylation
avg_per_region = rowSums(meth_table[methylation_covered_regions,rna_and_onTarget_methylation_covered_clones],na.rm=T) / 
				 rowSums(cov_table[methylation_covered_regions,rna_and_onTarget_methylation_covered_clones],na.rm=T)
not_boring_regions = names(avg_per_region) [(avg_per_region < maximal_avg_methylation) & (avg_per_region > minimal_avg_methylation)]
return (meth_table[not_boring_regions,rna_and_onTarget_methylation_covered_clones] / 
		cov_table[not_boring_regions,rna_and_onTarget_methylation_covered_clones])
}								 
#################################################################################################################
#################################################################################################################
#################################################################################################################


   
# This function retrieves a boolean matrix, indicating k neareast neighbours to each clone by given features
#################################################################################################################
cmem_get_norm_knn_matrix <- function(feats, k){
    dist_mat <- dist(feats)
    knn_df <- tgs_knn(100-as.matrix(dist_mat), k)
    knn_df <- knn_df %>% mutate(col1 = factor(col1), col2 = factor(col2, levels=levels(col1)))
    
    knn_mat <- sparseMatrix(as.numeric(knn_df$col1), as.numeric(knn_df$col2), x=1)
    rownames(knn_mat) <- levels(knn_df$col1)
    colnames(knn_mat) <- levels(knn_df$col2)
    return(knn_mat)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# This function gets a KNN matrix and average methylation. From each clone, it subtracts the average methylation its K neighbours.
#################################################################################################################
cmem_normalize_knn <- function(raw, knn_mat){
    raw <- raw[, colnames(knn_mat)]
    raw_filt_na <- raw
    raw_filt_na[is.na(raw)] <- 0
    met_exp <- as.matrix(raw_filt_na) %*% as.matrix(t(knn_mat))
    not_na_mat <- !is.na(raw)
    met_exp_n <- as.matrix(not_na_mat) %*% as.matrix(t(knn_mat))
    met_exp_norm <- met_exp / met_exp_n
    met_oe <- raw - met_exp_norm
    return(met_oe)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################



# This function generates a scatter of clonal HCG and LCG for a random clone, highlighting its selected KNN neighbors.
#################################################################################################################
cmem_generate_random_KNN_example <- function (clones_offTarget_methylation_table = NULL,
											  knn_mat = NULL,
											  fig_nm = NULL,
											  plot_w = 280,
											  plot_h = 280,
											  plot_margins = c(5,5,1,1),
											  plot_cex_text = 1.5,
											  selected_clone_dot_cex = 2,
											  plot_legend = FALSE,
											  selected_clone_color = "blue",
											  knn_clones_color = "red",
											  legend_x = 0.67,
											  legend_y = 0.285)
{

      if( is.null(fig_nm) ) {
               stop("undefined file name to save the plot to")
       }
	   
png(fig_nm, width = plot_w, height = plot_h)
par(mar = plot_margins)
plot(clones_offTarget_methylation_table$LCG,
	 clones_offTarget_methylation_table$HCG,
	 cex.lab=plot_cex_text,
	 cex.axis=plot_cex_text,
	 xlab="LCG trend",ylab="HCG trend")
random_clone = sample(rownames(clones_offTarget_methylation_table),size=1)
points(clones_offTarget_methylation_table[colnames(knn_mat)[knn_mat[random_clone,]==1],"LCG"],
	   clones_offTarget_methylation_table[colnames(knn_mat)[knn_mat[random_clone,]==1],"HCG"],pch=19,col=knn_clones_color)
points(clones_offTarget_methylation_table[random_clone,"LCG"],
	   clones_offTarget_methylation_table[random_clone,"HCG"],col=selected_clone_color,cex=selected_clone_dot_cex,pch=19)
if (plot_legend) {
legend(legend = paste0("K = ",sum(knn_mat[random_clone,]==1)),
	   x = legend_x,
	   y = legend_y,
	   x.intersp=0.2,
	   text.col=knn_clones_color,
	   cex=plot_cex_text,
	   bty="n")}
grid()
dev.off()
} 
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This function gets expression and methylation matrices, filter low-variance genes, compute cross correlation, 
# and then filter again genes and regions that not maintain variable or interesting correlations.
#################################################################################################################
cmem_get_ExprMeth_crossCorrelation <- function (exp_mat = NULL,
												exp_mat_norm = NULL,
												exp_ds_value = 5e4,
												meth_mat = NULL,
												min_gene_expression_normVar = 0.15,
												eps = 0.1,
												min_region_cor = 0.195,
												min_gene_cor = 0.28,
												min_region_cc_var = 0.002,
												min_gene_cc_var = 0.003)		
{
	if (sum(!colnames(meth_mat) %in% colnames(exp_mat))>0) {
		stop ("Couldn't find expression for clones defined by the methylation-table")}
		
exp_mat_ds = cmem_ds(exp_mat,exp_ds_value)
exp_v = apply(exp_mat_ds,1,var)
exp_e = rowMeans(exp_mat_ds)
exp_ve = log2( (exp_v+eps)/ (exp_e+eps) )
variable_genes = names(exp_ve)[exp_ve >= min_gene_expression_normVar]

crossCorr_mat = tgs_cor(t(meth_mat) ,
						t(log2(( 7*exp_mat_norm[variable_genes,colnames(meth_mat)])+1)),
						spearman=FALSE,
						pairwise.complete.obs=TRUE)
region_f_cor = apply(abs(crossCorr_mat), 1, function(x){max(x,na.rm=TRUE) >= min_region_cor})
gene_f_cor = apply(abs(crossCorr_mat), 2, function(x){max(x,na.rm=TRUE) >= min_gene_cor})
crossCorr_mat_f = crossCorr_mat[region_f_cor,gene_f_cor]
region_cc_var = apply(crossCorr_mat_f,1,var,na.rm=T)
gene_cc_var = apply(crossCorr_mat_f,2,var,na.rm=T)	
return (crossCorr_mat_f[region_cc_var > min_region_cc_var, gene_cc_var > min_gene_cc_var])
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	


# This function filters (non-Symmetric) correlation matrix, based on minimal objects that maintain interesting negative correlation
#################################################################################################################
#remove boring genes
cmem_filter_cc_mat <- function (crossCorr_mat_to_filter = NULL,
							    max_corr_to_regions = (-0.15),
								min_n_regions = 55,
								max_corr_to_genes = (-0.15),
								min_n_genes = 10)
{								
regions_count = apply(crossCorr_mat_to_filter, 2, function(x){ sum(x < max_corr_to_regions) })
bad_genes = regions_count < min_n_regions
genes_count = apply(crossCorr_mat_to_filter, 1, function(x){ sum(x < max_corr_to_genes) })
bad_regions = genes_count < min_n_genes
return(crossCorr_mat_to_filter[!bad_regions, !bad_genes])
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	


# This function plot a matrix of cross correlations and highlight specific genes
#################################################################################################################
cmem_plot_crossCorr <- function (cc_mat = NULL,
								 genes_to_highlight = NULL,
								 highlight_color = "green1",
								 highlight_title = "module_name",
								 base_dir = "./",
								 fig_nm = "crossCorrelation_matrix.png",
								 plot_w = 1300,
								 plot_h = 1950,
								 plot_cellheight = 20,
								 plot_fontsize = 24,
								 k_genes = 8,
								 k_regions = 5,
								 limit_corr_value = 0.3,
								 retrieve_clusters = FALSE
								 )	
{
      if( is.null(fig_nm) ) {
               stop("undefined file name to save the plot to")
       }

if (!is.null(genes_to_highlight)){
genes_df = data.frame(module = rep(NA,ncol(cc_mat)))
rownames(genes_df) = colnames(cc_mat) 
genes_df$module[rownames(genes_df) %in% genes_to_highlight] = highlight_title
Var1        <- c(highlight_color)
names(Var1) <- c(highlight_title)
anno_colors <- list(module = Var1)
}

cmem_mypheatmap(
			pmin(pmax(t(cc_mat),(-1)*limit_corr_value),limit_corr_value),
			cellheight = plot_cellheight,			
			outdir = base_dir,
		    fn = fig_nm,
		    center0 = TRUE,
		    cluster_rows = TRUE,
		    cluster_cols = TRUE,
		    cutree_cols = k_regions,
		    cutree_rows = k_genes,
		    show_colnames = FALSE,
		    show_rownames = TRUE,
		    treeheight_col = 0,
		    treeheight_row = 0,
		    w = plot_w,
		    h = plot_h,
		    annotation_legend = FALSE,
		    annotation_row = list(NA,genes_df)[[ifelse(is.null(genes_to_highlight),1,2)]],
		    annotation_colors = list(NA,anno_colors)[[ifelse(is.null(genes_to_highlight),1,2)]],
		    fontsize = plot_fontsize)
			
if (retrieve_clusters) {			
ph = pheatmap(pmin(pmax(t(cc_mat),(-1)*limit_corr_value),limit_corr_value),
		      cluster_rows = TRUE,
		      cluster_cols = TRUE,
		      cutree_cols = k_regions,
		      cutree_rows = k_genes,
			  silent = TRUE)
genes_k_vector = cutree(ph$tree_row, k = k_genes)
regions_k_vector = cutree(ph$tree_col, k = k_regions)
return (list (cc_gene_clusters = genes_k_vector,
			  cc_regions_clusters = regions_k_vector)) }
}				 
#################################################################################################################
#################################################################################################################
#################################################################################################################	

# This function plots a scatter of clonal expression for specific genes and clonal avg. methylation of specific regions
#################################################################################################################
cmem_plot_expression_and_methylation <- function (meth_table = NULL,
												  cov_table = NULL,
									              genes_to_plot = NULL,
												  expression_table = NULL,
									              expression_table_norm = NULL,
									              label_expression = "exp.",
									              label_methylation = "meth.",
									              fig_nm = NULL,
									              plot_w = 280,
									              plot_h = 280,
									              plot_text_cex = 1.5,
									              plot_margins = c(5.8,5.8,0.4,0.4),
									              plot_dots_cex = 1,
									              plot_dots_color = "#00000080",
												  min_umi_clone = 5e4,
									              log2_exp = FALSE,
									              reg = 1,
									              corr_method = "pearson",
									              show_corr = FALSE,
									              corr_position = "bottomleft",
									              regions_to_plot = NULL)
{
	  if (sum(!regions_to_plot %in% rownames(meth_table))>0){
		stop ("at least one region is missing from avg. meth table")}
	  
	  if (sum(!genes_to_plot %in% rownames(expression_table))>0){
		stop ("at least one gene is missing from expression table")}
		
      if( is.null(fig_nm) ) {
               stop("undefined file name to save the plot to")
       }
	   
	  if( is.null(regions_to_plot) ) {
               stop("undefined regions to plot")
       }
	   
	  if( is.null(genes_to_plot) ) {
               stop("undefined genes to plot")
       }

exp_n = expression_table_norm[,colSums(expression_table)>=min_umi_clone]	   
good_clones = intersect(colnames(meth_table),colnames(exp_n))
if (length(genes_to_plot)==1) {
clone_exp = exp_n[genes_to_plot,good_clones] } else { clone_exp = colSums(exp_n[genes_to_plot,good_clones]) } 
if (log2_exp) {clone_exp = log2 (clone_exp + reg)}

if (length(regions_to_plot)==1) {
clone_meth = meth_table[regions_to_plot,good_clones];
clone_cov = cov_table[regions_to_plot,good_clones] } else { 
clone_meth = colSums(meth_table[regions_to_plot,good_clones],na.rm=T);
clone_cov = colSums(cov_table[regions_to_plot,good_clones],na.rm=T);
 } 
clone_avg = clone_meth/clone_cov

png(fig_nm, width = plot_w, height = plot_h)
par(mar = plot_margins)
plot(clone_exp,
	 clone_avg,
	 pch = 19,
	 col = plot_dots_color,
	 cex = plot_dots_cex,
	 cex.axis = plot_text_cex,
	 cex.lab = plot_text_cex,
	 xlab = label_expression,
	 ylab = label_methylation)
grid()
if (show_corr) {
legend(legend = paste0("r = ",as.character(round(cor(clone_avg, clone_exp, method = corr_method), digits=2))),
		cex = plot_text_cex,
		bty = "n",
		x.intersp = 0.2,
		x = corr_position)}
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	




#### this function computes the expected methylation variance of given locus based on PERSISTENCY MODEL, based on its average methylation and 6 samplings of it 
#### the model in brief: each founder cell of a clone starts with either 0/0, 0/1 or 1/1 state, based on average methylation of this locus in the population.
#### 					 Once a cell got one of the three combinations, its alleles are faithfully propagated in its clone, just as genetic alleles.
#################################################################################################################
cmem_biallelic_p6 <- function(p)
{
	binom = choose(6,0:6)
	p_onhalf = binom*((1/2)**6)
	devi = (0:6 - 6*p)
	devi2 = sum(p_onhalf * devi * devi)

	v_half = devi2 * 2 * p * (1-p)

	v_onezero = 6*6*2 * (p * (1-p))**2

	return(v_half + v_onezero)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	



#### This function plots the observed methylation-variance, as well as expected variance by the model.
#################################################################################################################
cmem_plot_epipolymorphism <- function(ep_df = NULL,
									  fig_nm = NULL,
									  hot_color = "red",
									  cold_color = "deepskyblue3",
									  normal_color = "#00000050",
									  plot_dots_cex = 0.3,
									  plot_dots_hotCold_cex = 0.8,
									  observed_trend_color = "grey43",
									  persistancy_model_color = "deepskyblue1",
									  mixture_model_color = "coral",
									  plot_ylim = c(0,6.2),
									  plot_w = 380,
									  plot_h = 380,
									  plot_margins = c(5,5,1,1),
									  plot_text_cex=1.7,
									  plot_lwd = 3,
									  plot_legend_models = FALSE,
									  plot_legend_hotCold = FALSE,
									  legend_models_x = -0.085,
									  legend_models_y = 6.65,
									  legend_hotCold_x = 0.6,
									  legend_hotCold_y = 6.65
									  )
{

	  if( is.null(ep_df) ) {
               stop("undefined Epi-polymorphism table")
       }
	   
	  if( is.null(fig_nm) ) {
               stop("undefined file name to save the plot to")
       }
	   
line_mixture_x = seq(0,1,0.001)
line_mixture_y = line_mixture_x*(1-line_mixture_x)*6
line_persistency_x = seq(0,1,0.001)
line_persistency_y = c(); for(p in seq(0,1,0.001)) {line_persistency_y = c(line_persistency_y,cmem_biallelic_p6(p))};
line_observed_x = rollmedian(ep_df$m[order(ep_df$m)],101)
line_observed_y = rollmedian(ep_df$obs[order(ep_df$m)],101)

png(fig_nm, width = plot_w, height = plot_h)
par(mar = plot_margins)
plot(ep_df$m, ep_df$obs,
	cex = plot_dots_cex, pch = 19, col = normal_color,
	cex.main=plot_text_cex, cex.lab=plot_text_cex, cex.axis=plot_text_cex,
	ylim=plot_ylim,
	xlab=bquote(paste("average methylation ",italic("m"))),
	ylab="inter-clonal variance",
	main="")
lines(line_mixture_x, line_mixture_y, col=mixture_model_color,lwd = plot_lwd)
lines(line_persistency_x, line_persistency_y, col = persistancy_model_color, lwd = plot_lwd)
lines(line_observed_x, line_observed_y, col=observed_trend_color, lwd=plot_lwd)
points(ep_df$m [ep_df$regime=="hot"],  ep_df$obs[ep_df$regime=="hot"], cex=plot_dots_hotCold_cex, pch=19, col=hot_color)
points(ep_df$m [ep_df$regime=="cold"], ep_df$obs[ep_df$regime=="cold"], cex=plot_dots_hotCold_cex, pch=19, col=cold_color)
grid()

if (plot_legend_models){
legend(x = legend_models_x,
	   y = legend_models_y,
	   x.intersp = 0.2,
	   y.intersp = 0.8,
	   bty = "n",
	   lty = c(1,1,1),
	   lwd = c(3,3,3),
	   seg.len = c(0.6,0.6,0.6),
	   text.col = c(persistancy_model_color,observed_trend_color,mixture_model_color),
	   col = c(persistancy_model_color,observed_trend_color,mixture_model_color),
	   cex=plot_text_cex,
	   legend=c("persistency","observed","mixture")) }
if (plot_legend_hotCold) {
legend(x=legend_hotCold_x,
	   y=legend_hotCold_y,
	   x.intersp=0.2,
	   y.intersp=0.8,
	   bty="n",
	   pch=c(20,20),
	   text.col = c(cold_color,hot_color),
	   col = c(cold_color,hot_color),
	   cex = plot_text_cex,
	   legend = c("COLD loci","HOT loci"))
}
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	




#### This function retrieves clonal average mehtylation at HOT and HOLD CpGs.
#################################################################################################################
cmem_get_hotCold_clonal_average <- function (min_coverage = 100,
											 meth_tab = NULL,
											 cov_tab = NULL,
											 ep_df = NULL,
											 exclude_EMT_clones = TRUE,
											 HCG_mat = NULL,
											 LCG_mat = NULL,
											 clones_metadata = cmem_clones_metadata
											 ) 
{

	  if( is.null(ep_df) ) {
               stop("undefined Epi-polymorphism table")
       }
	   
hot_cpgs = ep_df$regime=="hot"
cold_cpgs = ep_df$regime=="cold"
hot_meth = colSums(meth_tab[hot_cpgs,4:ncol(meth_tab)],na.rm=T)
hot_cov = colSums(cov_tab[hot_cpgs,4:ncol(cov_tab)],na.rm=T)
cold_meth = colSums(meth_tab[cold_cpgs,4:ncol(meth_tab)],na.rm=T)
cold_cov = colSums(cov_tab[cold_cpgs,4:ncol(cov_tab)],na.rm=T)

good_clones = names(HCG_mat$avg_orig)
hot_cold_df = data.frame (HCG = HCG_mat$avg_orig,
						  LCG = LCG_mat$avg_orig,
						  cold_avg = cold_meth[good_clones]/cold_cov[good_clones],
						  hot_avg = hot_meth[good_clones]/hot_cov[good_clones])
covered_clones = intersect(names(hot_cov)[hot_cov>=min_coverage], names(cold_cov)[cold_cov>=min_coverage])

if (exclude_EMT_clones){
emt_clones = rownames(clones_metadata)[clones_metadata$VIM_state=="VIM-high"]
hot_cold_df = hot_cold_df[!rownames(hot_cold_df) %in% emt_clones,]}

return(hot_cold_df[rownames(hot_cold_df) %in% covered_clones,])
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	


#### This function create boxplots of clonal averages (for example - "cold" and "hot" methylation), stratified by offtarget LCG and HCG clonal trends
#################################################################################################################
cmem_plot_clonal_averages_stratified_by_offTargetTrends <- function (clonal_averages_df = NULL,
																	 plot_w = 180,
																	 plot_h = 300,
																	 plot_margins = c(7,5,1,1),
																	 plot_boxwex = 0.7,
																	 hot_color = "red",
																	 cold_color = "deepskyblue3",
																	 base_dir = "./",
																	 plot_text_cex = 1.5)
{

	  if( is.null(clonal_averages_df) ) {
               stop("undefined clonal averages table")
       }
	   
clonal_averages_df$HCG_group = c("HCG_low","HCG_mid","HCG_high")[.bincode(clonal_averages_df[,"HCG"],breaks=quantile(clonal_averages_df[,"HCG"],c(seq(0,1,length.out = 4))),right=TRUE,include.lowest=TRUE)]
clonal_averages_df$LCG_group = c("LCG_low","LCG_mid","LCG_high")[.bincode(clonal_averages_df[,"LCG"],breaks=quantile(clonal_averages_df[,"LCG"],c(seq(0,1,length.out = 4))),right=TRUE,include.lowest=TRUE)]

png(paste0(base_dir,"LCG_offTar_trends_boxplot.png"), width=plot_w, height=plot_h)
par(mar = plot_margins)
boxplot(clonal_averages_df$cold_avg[clonal_averages_df$LCG_group=="LCG_low"],
		clonal_averages_df$hot_avg[clonal_averages_df$LCG_group=="LCG_low"],
		clonal_averages_df$cold_avg[clonal_averages_df$LCG_group=="LCG_mid"],
		clonal_averages_df$hot_avg[clonal_averages_df$LCG_group== "LCG_mid"],
		clonal_averages_df$cold_avg[clonal_averages_df$LCG_group=="LCG_high"],
		clonal_averages_df$hot_avg[clonal_averages_df$LCG_group== "LCG_high"],
		names = rep("",6),
		xaxt = "n",
		boxwex = plot_boxwex,
		at = c(1:2,3.5:4.5,6:7),
		las = 2, cex.axis=plot_text_cex, cex.lab=plot_text_cex, cex.names=plot_text_cex,
		ylab="clonal average meth.",
		yaxt="n",
		col=rep(c(cold_color,hot_color),3))
axis(2,las=3,cex.axis=plot_text_cex, cex=plot_text_cex)
axis(1,at=c(1.5,4,6.5), cex.lab=plot_text_cex, cex.axis=plot_text_cex,
	 las=2,
	 labels=c("LCG low", "LCG mid", "LCG high"))
dev.off()
		
png(paste0(base_dir,"HCG_offTar_trends_boxplot.png"), width=plot_w, height=plot_h)
par(mar = plot_margins)
boxplot(clonal_averages_df$cold_avg[clonal_averages_df$HCG_group=="HCG_low"],
		clonal_averages_df$hot_avg [clonal_averages_df$HCG_group== "HCG_low"],
		clonal_averages_df$cold_avg[clonal_averages_df$HCG_group=="HCG_mid"],
		clonal_averages_df$hot_avg [clonal_averages_df$HCG_group== "HCG_mid"],
		clonal_averages_df$cold_avg[clonal_averages_df$HCG_group=="HCG_high"],
		clonal_averages_df$hot_avg [clonal_averages_df$HCG_group== "HCG_high"],
		names = rep("",6),
		xaxt = "n",
		boxwex = plot_boxwex,
		at=c(1:2,3.5:4.5,6:7),
		las=2,cex.axis=plot_text_cex, cex.lab=plot_text_cex, cex.names=plot_text_cex,
		ylab="clonal average meth.",
		yaxt="n",
		col=rep(c(cold_color,hot_color),3))
axis(2,las=3,cex.axis=1.5,cex=1.5)
axis(1,at=c(1.5,4,6.5), cex.lab=plot_text_cex, cex.axis=plot_text_cex,
	 las=2,
	 labels=c("HCG low", "HCG mid", "HCG high"))
dev.off()
}			
#################################################################################################################
#################################################################################################################
#################################################################################################################	





######################################################################
#### this function is a modification for pheatmap, where center0 makes 0 values as white
cmem_mypheatmap = function (mat, w=1000, h=1000, fname,outdir, center0 = F,
							my_color = colorRampPalette( rev(brewer.pal(n = 7, name ="RdBu")))(100),
							color_quantile = 0.99, ...) {
  parms = list(... )
  png(sprintf("%s%s", outdir, fname), width=w, height=h)
  if (center0) {
    if (is.null(parms$color)) {
      parms$color =  my_color
    }
    mm = quantile(abs(mat), color_quantile,na.rm=T)
    parms$breaks = sort(unique (c( min(mat), seq(-mm, mm, l=99),max(mat)))) 
  }
  parms$mat = mat
  do.call (pheatmap, parms) 
 dev.off()
}
###############################


# This function will generate a "misha" hg19 database, that will be used throughout the analysis
#################################################################################################################
cmem_generate_hg19_misha_db <- function ()
{
ftp <- "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19" 
gdb.create("hg19",
           paste(ftp, 'chromosomes', paste0('chr', c(1:22, 'X', 'Y'), '.fa.gz'), sep='/'),
           paste(ftp, "database/knownGene.txt.gz", sep = "/"),
           paste(ftp, "database/kgXref.txt.gz", sep = "/"),
           c("kgID", "mRNA", "spID", "spDisplayID", "geneSymbol",
             "refseq", "protAcc", "description", "rfamAcc",
             "tRnaName"))
gdb.init('hg19')
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	




# This function will plot insulation score, as well as the defined TAD borders, over a given genomic region
#################################################################################################################
cmem_plot_insulation_TAD_example <- function (tads_table = NULL,
											  insulation_table = NULL,
											  example_chrom = "chr7",
											  example_start = 1000000,
											  example_end = 8500000,
											  base_dir = "./",
											  plot_w = 937.5,
											  plot_h = 250,
											  insulation_color = "#00000075",
											  TAD_border_color = "red",
											  TAD_border_lwd = 2,
											  coords_to_show = c(2e6,4e6,6e6,8e6),
											  coords_labels = c("2M","4M","6M","8M"),
											  plot_margins = c(5,8,2,2),
											  plot_dots_cex = 1.25,
											  plot_text_cex = 1.5,
											  show_title = FALSE)
{

example_insulation = insulation_table[insulation_table$chrom==example_chrom,]
example_insulation_region = example_insulation[example_insulation$start>example_start & example_insulation$end<example_end,]

example_borders = tads_table[tads_table$chrom==example_chrom,]
example_borders_region = example_borders[example_borders$start>example_start & example_borders$end<example_end,]
browser()
png(paste0(base_dir,example_chrom,"_",example_start,"_",example_end,"_insulation_TADborders_example.png"),
	width = plot_w, height = plot_h)
par(mar = plot_margins)
plot(example_insulation_region$start,
	 -1*(example_insulation_region[,4]),pch=19,
	 cex=plot_dots_cex,col=insulation_color,	 
	 yaxt = "n",
	 xaxt = "n",
	 ylab = "",
	 xlab = "",
	 ylim = c(min(-1*(example_insulation_region[,4]),na.rm=T)-0.3,max(-1*(example_insulation_region[,4]),na.rm=T)),
	 cex.axis=plot_text_cex,cex.lab=plot_text_cex)
apply(example_borders_region,1,
	 function(x){segments(x0 = as.numeric(x[2]), x1 = as.numeric(x[3]), y0 = par("usr")[3]+0.2,y1 = par("usr")[3]+0.2, col = TAD_border_color,lwd = TAD_border_lwd,lty=1);
				 segments(x0 = as.numeric(x[2]), x1 = as.numeric(x[2]), y0 = par("usr")[3]+0.1,y1 = par("usr")[3]+0.3, col = TAD_border_color,lwd = TAD_border_lwd,lty=1);
				 segments(x0 = as.numeric(x[3]), x1 = as.numeric(x[3]), y0 = par("usr")[3]+0.1,y1 = par("usr")[3]+0.3, col = TAD_border_color,lwd = TAD_border_lwd,lty=1);
				 })
abline(h=c(2:5),lwd=1.25)
axis(side=2,cex.axis = plot_text_cex,cex.lab=plot_text_cex, cex = plot_text_cex,las=2)
axis(side=1,cex.axis = plot_text_cex,cex.lab=plot_text_cex, cex = plot_text_cex, at=coords_to_show, label=coords_labels, las=1)
title(ylab="-1 * Insulation",cex.lab=plot_text_cex,cex=plot_text_cex,line=3.5)
title(xlab="Coordinate",cex.lab=plot_text_cex,cex=plot_text_cex,line=2.5)
if (show_title) {
legend(x="topright",
	   x.intersp=0.2,
	   y.intersp=0.8,
	   bty="n",
	   text.col="black",
	   cex=plot_text_cex,
	   legend = example_chrom)
legend(x="topleft",
	   x.intersp=0.2,
	   y.intersp=0.8,
	   bty="n",
	   text.col=TAD_border_color,
	   cex=plot_text_cex,
	   legend="TAD borders")
}
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	


# This function generate stats on gene-to-TAD expression
#################################################################################################################
cmem_generate_gtad_corr_stats <- function (gtad_cor_matrix = NULL,
										   gene_expression_gtad = NULL,
										   tad_to_gene_convertion_table = tad_gene_key
										   )

{
gtad_3rd = rep(-1,nrow(gtad_cor_matrix))
gtad_4th = rep(-1,nrow(gtad_cor_matrix))
gtad_9th = rep(-1,nrow(gtad_cor_matrix))
for (i in 1:nrow(gtad_cor_matrix)){
gtad_3rd[i] = gtad_cor_matrix[rownames(gtad_cor_matrix)[i],as.character(tad_to_gene_convertion_table$TAD[tad_to_gene_convertion_table$gene %in% rownames(gtad_cor_matrix)[i]])]
gtad_4th[i] = (ncol(gtad_cor_matrix)+1) - rank(gtad_cor_matrix[rownames(gtad_cor_matrix)[i],])[as.character(tad_to_gene_convertion_table$TAD[tad_to_gene_convertion_table$gene %in% rownames(gtad_cor_matrix)[i]])]
gtad_9th[i] = sort(unlist(gtad_cor_matrix[rownames(gtad_cor_matrix)[i],!colnames(gtad_cor_matrix) %in% as.character(tad_to_gene_convertion_table$TAD[tad_to_gene_convertion_table$gene %in% rownames(gtad_cor_matrix)[i]])]),decreasing=TRUE)[20]
if (i%%500==0){cat(paste0(i,"\n"))}
}
gtad_df = data.frame(gene_expression_gtad = gene_expression_gtad, gt_cor = gtad_3rd,
					 gt_cor_rank = gtad_4th, 
					 gt_cor_max20_nonself = gtad_9th
					 )
rownames(gtad_df) = rownames(gtad_cor_matrix)
return(gtad_df)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################	


# This function plots the TAD screen ranking histogram od genes to their own TADS
#################################################################################################################
cmem_plot_TAD_screen_rankHist <- function (gtad_stats_table = gtad_stats_hct116,
								  fig_nm = NULL,
								  plot_w = 400,
								  plot_h = 280,
								  plot_margins = c(5,5,1,1),
								  plot_cex_text = 1.5,
								  plot_xlab = "XLAB",
								  plot_breaks = 16)
{
png(fig_nm,width=400,height=280)
par(mar=c(5,5,1,1))
hist(gtad_stats_table$gt_cor_rank,xlab=plot_xlab,cex=1.5,cex.axis=1.5,cex.lab=1.5,
	 breaks=16,
	 main="",
	 #col=")
	 border="grey43",col="grey83")
dev.off()

}							   
#################################################################################################################
#################################################################################################################
#################################################################################################################



# This function compares distributions of surprising gene-TAD expression correlation (when correlation to self-TAD is higher than #20 ranked TAD)
# The comparison is between original values of the screen, and those obtained by the same screen when gene-TAD structure was shuffled
#################################################################################################################
cmem_plot_TAD_screen_over_permutation <- function (fig_nm_ecdf = NULL,
												   fig_nm_scatter = NULL,
												   gtad_stats_orig = NULL,
												   gtad_stats_permute = NULL,
												   FDR5_genes = c(),
												   FDR25_genes = c(),
												   FDR5_color = "deepskyblue3",
												   FDR25_color = "deepskyblue1",
												   plot_ecdf_cex_text = 1.2,
												   ecdf_ylim = c(0,300),
												   ecdf_xlim = c(0.03,0.5),
												   ecdf_lwd = 4,
												   plot_ecdf_w = 280,
												   plot_ecdf_h = 220,
												   plot_ecdf_margins = c(5,5,1,1),
												   show_ecdf_legend = FALSE,
												   plot_xy_w = 220,
												   plot_xy_h = 220,
												   plot_xy_margins = c(5,5,1,1),
												   plot_xy_cex_text = 1.2,
												   plot_xy_xlab = "XLAB",
												   plot_xy_ylab = "YLAB",
												   plot_xy_xlim = c(0,5),
												   eps = 1,
												   plot_xy_dots_cex = 0.8
												   )
{

		if ( is.null(fig_nm_ecdf) | is.null(fig_nm_scatter) ) {
			stop ("undefined file names for figures") }

# ecdf
ecdf_orig_corrTAD=ecdf(gtad_stats_orig$gt_cor-gtad_stats_orig$gt_cor_max20_nonself)
uniq_orig_corrTAD=sort(unique(gtad_stats_orig$gt_cor-gtad_stats_orig$gt_cor_max20_nonself))
ecdf_perm_corrTAD=ecdf(gtad_stats_permute$gt_cor-gtad_stats_permute$gt_cor_max20_nonself)
uniq_perm_corrTAD=sort(unique(gtad_stats_permute$gt_cor-gtad_stats_permute$gt_cor_max20_nonself))
png(fig_nm_ecdf, height=plot_ecdf_h, width=plot_ecdf_w)
par(mar = plot_ecdf_margins)
plot(rev(uniq_orig_corrTAD),((1-ecdf_orig_corrTAD(rev(uniq_orig_corrTAD)))*length(gtad_stats_orig$gt_cor-gtad_stats_orig$gt_cor_max20_nonself)),
	 xlab = "corr. (self TAD - 20th rank)",
	 ylab = "# genes",
	 yaxt = "n",
	 main = "",
	 type = "l",
	 xlim = ecdf_xlim,
	 ylim = ecdf_ylim,	
	 cex.axis = plot_ecdf_cex_text, cex.lab = plot_ecdf_cex_text,
	 lwd = ecdf_lwd, col = "black")	
axis(2,col="black",col.ticks="black",col.axis="black",
	 cex.axis=plot_ecdf_cex_text, cex.lab=plot_ecdf_cex_text, cex=plot_ecdf_cex_text, las=2)
lines(rev(uniq_perm_corrTAD)
	  ,((1-ecdf_perm_corrTAD(rev(uniq_perm_corrTAD)))*length(gtad_stats_permute$gt_cor-gtad_stats_permute$gt_cor_max20_nonself)),
	  lwd=ecdf_lwd,col="grey63")
grid(nx=NULL,ny=NULL,col = "lightgrey", lty = "dotted",lwd=par("lwd"))
if (show_ecdf_legend) {
legend("topright", 
 legend=rev(c("original","randomized")),
 col=rev(c("black","grey53")),
 x.intersp=0.2,
 y.intersp=0.8,
 seg.len=0.3,
 text.col=rev(c("black","grey53")),
 bty="n",lwd=2,cex=plot_ecdf_cex_text) }
 dev.off()
 
# scatter
x_val = log2(gtad_stats_orig$gene_expression_gtad+eps)
y_val = gtad_stats_orig$gt_cor - gtad_stats_orig$gt_cor_max20_nonself
xy_colors = c("#00000030",FDR5_color,FDR25_color)[ifelse(rownames(gtad_stats_orig) %in% FDR5_genes,2,ifelse(
													  rownames(gtad_stats_orig) %in% FDR25_genes,3,1))]
png(fig_nm_scatter, width = plot_xy_w, height = plot_xy_h)
par(mar = plot_xy_margins)
plot(x_val,y_val,pch=19,cex = plot_xy_dots_cex,
	col = xy_colors,
	xlim = plot_xy_xlim,
	cex.lab = plot_xy_cex_text,
	cex.axis = plot_xy_cex_text,
	xlab = plot_xy_xlab,
	ylab = plot_xy_ylab)
grid()
dev.off()
}
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This function just squeezes the generation of gene-TAD correlation stats 
#################################################################################################################
cmem_generate_hic_stats <- function (gene_tad_corr_csv = NULL,
									 gene_expression_csv = NULL,
									 permuted_gene_tad_corr_csv = NULL,
									 permuted_gene_expression_csv = NULL,
									 real_tad_to_gene_convertion_csv = NULL,
									 permuted_tad_to_gene_convertion_csv = NULL)
{
		if (is.null(gene_tad_corr_csv) |
			 is.null(gene_expression_csv) |
		 	  is.null(permuted_gene_tad_corr_csv) |
			   is.null(permuted_gene_expression_csv) |
			    is.null(real_tad_to_gene_convertion_csv) |
			     is.null(permuted_tad_to_gene_convertion_csv) ) {
					stop("Missing .csv path") }
		
# load HCT116 gene-TAD corrs
gtad_cor = read.csv(gene_tad_corr_csv, row.names = 1)
gene_expression_gtad = as.data.frame(read.csv(gene_expression_csv, row.names = 1))[,1]
names(gene_expression_gtad) = rownames(gtad_cor)
#
# load TAD-gene key
tad_gene_key = as.data.frame(apply(read.csv(real_tad_to_gene_convertion_csv, row.names = 1),2,as.character))
rownames(tad_gene_key) = tad_gene_key$gene
#
# generate_stats on gene-TAD correlations
gtad_stats_real = cmem_generate_gtad_corr_stats (gtad_cor_matrix = gtad_cor,
												 gene_expression_gtad = unlist(gene_expression_gtad),
												 tad_to_gene_convertion_table = tad_gene_key[rownames(gtad_cor),])
#
permuted_gtad_cor = read.csv(permuted_gene_tad_corr_csv, row.names = 1)
gene_expression_permuted_gtad = as.data.frame(read.csv(permuted_gene_expression_csv, row.names = 1))[,1]
names(gene_expression_permuted_gtad) = rownames(permuted_gtad_cor)
#
# generate_stats on permuted gene-TAD structure, where TADs have the same number of genes but random ones.
tad_gene_key_perm = as.data.frame(apply(read.csv(permuted_tad_to_gene_convertion_csv, row.names = 1),2,as.character))
rownames(tad_gene_key_perm) = tad_gene_key_perm$gene
gtad_stats_perm = cmem_generate_gtad_corr_stats (gtad_cor_matrix = permuted_gtad_cor,
												 gene_expression_gtad = unlist(gene_expression_permuted_gtad),
												 tad_to_gene_convertion_table = tad_gene_key_perm)
return(list(gtad_stats_orig  = gtad_stats_real, gtad_stats_perm = gtad_stats_perm))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This function plot normalized expression of a given TAD in clones
#################################################################################################################
cmem_tad_exp_plot <- function(tad_id = NULL,
							  expression_table = NULL,
							  expression_table_norm = NULL,
							  min_UMI_per_clone = 1e4,
							  tad_to_gene_convertion = NULL,
							  reg=0.1,
							  heatmap_fontsize = 21,
							  heatmap_cellheight = 18,
							  heatmap_cellwidth = 1.7,
							  base_dir = "./",
							  heatmap_colors = colorRampPalette( rev(brewer.pal(n = 11, name ="BrBG")))(100),
							  annot_columns=NA,
							  annot_columns_colors=NA,
							  show_annotation_legend=TRUE,
							  show_legend=FALSE,
							  lfp_z_lim = 2,
							  min_gene_expression = 0
							  ) 
{
genes_to_plot = as.character(tad_to_gene_convertion[tad_to_gene_convertion$TAD %in% tad_id,"gene"])
genes_to_plot = genes_to_plot[ rowSums(expression_table_norm[genes_to_plot,]) > min_gene_expression]

expression_table_norm = expression_table_norm[, colSums(expression_table) >= min_UMI_per_clone]
mat_fp = t(apply(expression_table_norm[genes_to_plot,],1,function(x) {log2 ( (x+reg) / (median(x)+reg) )}))
mat_lfp = pmin(pmax(mat_fp,(-1)*lfp_z_lim),1*lfp_z_lim)

cmem_mypheatmap(mat_lfp[,order(colSums(expression_table_norm[genes_to_plot,]))],
				my_colors = heatmap_colors,
				legend = show_legend,
				cluster_cols = FALSE,
				center0 = TRUE,
				annotation_col = annot_columns,
				annotation_colors = annot_columns_colors,
				annotation_legend = show_annotation_legend,
				show_colnames=FALSE,
				treeheight_row=0,
				treeheight_col=0,
				cellheight=heatmap_cellheight,
				fontsize=heatmap_fontsize,
				cellwidth=heatmap_cellwidth,
				outdir=base_dir,
				fname=paste0("TAD_expresion_clones_",tad_id,".png"))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################

# This function plot normalized contact landscape of a given TAD and its surrondings (requires SHAMAN PACKAGE!)
#################################################################################################################
cmem_tad_shaman_plot <- function(tad_id = NULL,
								 expand=5e5,
								 TAD_intervals_in_2D_format = NULL,
								 map_point_size = 1,
								 genome_version = "hg19",
								 gene_intervals = NULL,
								 HCT116_observed_contacts_track_name = NULL,
								 base_dir = "./",
								 plot_w = 562.5,
								 plot_h = 562.5,
								 tads_table = tads,
								 plot_annotations_size = 0.7,
								 plot_track_size = 0.1,
								 plot_genes_size = 0.1
								 )
{
	if (is.null(HCT116_observed_contacts_track_name)) {
		stop ("you must download track containing observed contacts, and put it in hg19/ MISHA database directory")}
		
		
expand_kb = expand
rel_interval = TAD_intervals_in_2D_format[tads_table$ID==tad_id,]
rel_interval$start1=rel_interval$start1-min(expand_kb,rel_interval$start1)
rel_interval$start2=rel_interval$start2-min(expand_kb,rel_interval$start2)
rel_interval$end1=rel_interval$end1+expand_kb
rel_interval$end2=rel_interval$end2+expand_kb
point_score = shaman_shuffle_and_score_hic_mat(obs_track_nms = HCT116_observed_contacts_track_name,
											   interval = rel_interval, work_dir = base_dir)

rel_interval_plot = tads_table[tads_table$ID==tad_id,1:3]
rel_interval_plot$start = rel_interval_plot$start-min((rel_interval_plot$start-1),expand_kb)
rel_interval_plot$end = rel_interval_plot$end+expand_kb
								  
shaman_plot_map_score_with_annotations(genome = genome_version, 
									   points_score = point_score$points, 
									   interval_range = rel_interval_plot, 
									   point_size = map_point_size,
									   annotations=list(gene_intervals[,1:3]),
									   a_colors=c("aquamarine3"),
									   add_genes = TRUE,
									   add_ideogram = TRUE,
									   add_axis = TRUE,
									   gene_stacking = "dense",
									   fig_fn = paste0(base_dir,
													  tad_id,"_region.png"),
									   gene_size = plot_genes_size,
									   track_size = plot_track_size,
									   annotation_size = plot_annotations_size,
									   fig_width = plot_w,
									   fig_height = plot_h
										)
										
}		
#################################################################################################################
#################################################################################################################
#################################################################################################################


# These 2 functions generate insulation track bases on a raw HiC contacts track
#################################################################################################################
gtrack.2d.gen_insu_prof = function(track_nm, scale, res, min_diag_d=1000)
{
        #extract - using an iter on

        iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)

        iter_2d = gintervals.2d(chroms1 = iter_1d$chrom,
                        starts1=iter_1d$start,
                        ends1=iter_1d$end,
                        chroms2 = iter_1d$chrom,
                        starts2=iter_1d$start,
                        ends2=iter_1d$end)

        if(length(gvtrack.ls("obs_big")) == 1) {
                gvtrack.rm("obs_big")
        }
        if(length(gvtrack.ls("obs_ins")) == 1) {
                gvtrack.rm("obs_ins")
        }
        gvtrack.create("obs_big", track_nm, "weighted.sum")
        gvtrack.create("obs_ins", track_nm, "weighted.sum")
        gvtrack.iterator.2d("obs_big",
                eshift1=scale, sshift1=-scale,
                eshift2=scale, sshift2=-scale)
        gvtrack.iterator.2d("obs_ins",
                eshift1=0, sshift1=-scale,
                eshift2=scale, sshift2=0)

        message("will iter on ", dim(iter_2d)[1])
        ins = gextract("obs_big", "obs_ins",
                        gintervals.2d.all(),
                         iterator=iter_2d, band=c(-scale*2,0))
        ins_diag = gextract("obs_big", "obs_ins",
                        gintervals.2d.all(),
                         iterator=iter_2d, band=c(-min_diag_d,0))
        ins$obs_big = ins$obs_big - ins_diag$obs_big
        ins$obs_ins = ins$obs_ins - ins_diag$obs_ins

        message("will retrun ins with ", dim(ins)[1], " rows")

        return(ins)
}
gtrack.2d.gen_insu_track = function(track_nm, scale, res, min_diag_d=1000, new_track, description="")
{
        k_reg = 10
        prof = gtrack.2d.gen_insu_prof(track_nm, scale, res, min_diag_d)
        message("names ", paste(names(prof), collapse=","))
        names(prof)[1] = "chrom"
        names(prof)[2] = "start"
        names(prof)[3] = "end"

        gtrack.create_sparse(track=new_track, prof[,c(1,2,3)], log2(prof$obs_ins/(prof$obs_big+k_reg)), description=description)
}
#################################################################################################################
#################################################################################################################
#################################################################################################################




# This function generates supporting bar scales for the regional HiC + clonal-expression plots
#################################################################################################################
cmem_plot_colorBars_HiC_examples <- function (hic_vals = c(0:200),
											  hic_colors = shaman_hic_colors,
											  hic_values_to_show = c(1,51,101,151,201),
											  hic_values_Xcoord = 15,
											  hic_bar_title = NULL,
											  hic_fig_nm = NULL,
											  expression_vals = round(seq(-2,2,length.out=100),digits=1),
											  expression_colors = expression_colors,
											  expression_values_to_show = c(1,50,100),
											  expression_values_Xcoord = 15,
											  expression_bar_title = NULL,
											  expression_fig_nm = NULL
											  )
{
# contacts
cmem_plot_color_bar (vals = hic_vals, cols = hic_colors, fig_fn = hic_fig_nm,
					 show_vals_ind = hic_values_to_show,
					 values_xcoord = hic_values_Xcoord,
					 title = hic_bar_title)
					 
# exprssion					 
cmem_plot_color_bar(vals = expression_vals, cols = expression_colors,fig_fn = expression_fig_nm,
					show_vals_ind = expression_values_to_show,
					values_xcoord = expression_values_Xcoord,
					title =expression_bar_title )					 
}
#################################################################################################################
#################################################################################################################
#################################################################################################################


# This is a general function that draws color-bars
#################################################################################################################
cmem_plot_color_bar <- function (vals, cols, fig_fn=NULL, title="",
								 show_vals_ind=NULL,
								 plot_w=400,
								 plot_h=400,
								 values_angle=0,
								 values_xcoord=19,
								 pos_bar_in_plot = c(0,100),
								 cex_values_key=2,cex_ylab=2,ylab_xcoord=1
								 )
{
  if (!is.null(fig_fn)) {
   png(fig_fn, width = plot_w, height = plot_h)
		}
		plot.new()
		plot.window(xlim=pos_bar_in_plot, 
					ylim=c(0, length(cols) + 3))
		rect(9, 1:length(cols), 15, 1:length(cols) + 1, border=NA, col=cols)
		rect(9, 1, 15, length(cols)+1, col=NA, border = 'black')
		if (is.null(show_vals_ind)) {
			show_vals_ind = rep(T, length(cols))
			}
		text(values_xcoord, (1:length(cols))[show_vals_ind] + 0.5, labels = vals[show_vals_ind], 
		pos = 4, cex = cex_values_key, srt = values_angle)
		text(ylab_xcoord, length(cols)/2 + 1, labels = title, srt = 90, cex = cex_ylab)
		
		if (!is.null(fig_fn)) {
    dev.off()
  }
}
#################################################################################################################
#################################################################################################################
#################################################################################################################

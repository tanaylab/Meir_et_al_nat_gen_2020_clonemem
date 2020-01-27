##
##
## HiC.r
## Date: 26.01.2020
## This file summarizes main steps used for screen for TAD-memory in HCT116 and other cells lines, as part of the Luria & Delbruck memory screen comparing single cells and single-cell-derived clonal populations (Nat. gen. Meir et al. 2020)
##
##


# load dependency packages and project functions
source("Meir_et_al_2020_nat_gen_functions.r")

# download and initialize a hg19 environment for analysis with "misha" packages
#  NOTE: this will download the assembly and other files into /hg19 directory, that will be generated in current directory.
#  NOTE: this will take ~5 minutes.
#  NO WORRIES: If you already have an hg19 database, it wouldn't download it again.
cmem_generate_hg19_misha_db ()
# now download 2 tracks (a smaller one with insulation score per region, and a bigger one with raw contact counts)
name_of_insulation_data_track = "HCT116_hic_notMerged_insulation_1e3"
name_of_raw_data_track = "HCT116_hic_notMerged"
download.file(url = "http://www.wisdom.weizmann.ac.il/~zoharme/Meir_et_al_NatGen_2020/misha_tracks/HCT116_hic_notMerged_insulation_1e3.track.zip", 
			  destfile = paste0("hg19/",name_of_insulation_data_track,".track.zip"))
download.file(url = "http://www.wisdom.weizmann.ac.il/~zoharme/Meir_et_al_NatGen_2020/misha_tracks/HCT116_hic_notMerged.track.zip", 
			  destfile = paste0("hg19/",name_of_raw_data_track,".track.zip"))		  
unzip(zipfile = paste0("hg19/",name_of_insulation_data_track,".track.zip"),
	  file = paste0("hg19/",name_of_insulation_data_track,".track"))
unzip(zipfile = paste0("hg19/",name_of_raw_data_track,".track.zip"),
	  file = paste0("hg19/",name_of_raw_data_track,".track"))
gsetroot("hg19/")
gdb.reload()

# load TAD borders and insulation as obtained eventually from module 3.1 and alaysis of tracks
tads_table = read.csv ("HiC_data/HCT116_TAD_borders.csv", row.names = 1)
insulation_table = read.csv ("HiC_data/HCT116_insulation_table.csv", row.names = 1)



################################################################################################
######
###### part 3: Analysis of clonal TAD-memory
######
######		module 3.1 - define TAD-borders based on HCT116 HiC Data (raw contacts downloaded from Rao et al 2017)
###### 		module 3.2 - screen for clonal TAD-memory over shuffled controls in HCT116, A549 and WI38 cells.
################################################################################################



################################################################################################
######
###### module 3.1 - define TAD-borders based on HCT116 HiC Data (raw contacts downloaded from Rao et al 2017)
######
######	step 3.1.1 - generate insulation track from contacts.
######	step 3.1.2 - Use insulation track to define TAD borders.
######  step 3.1.3 - plot examples for insulation and defined TADs over genomic regions.
######
################################################################################################

#
######	step 3.1.1 - generate insulation track from contacts.
#
# IMPORTANT: this step will only work if you download MISHA and SHAMAN packages, as well as the raw contacts track
gtrack.2d.gen_insu_track(track_nm = name_of_raw_data_track,
						 scale=2e5,
						 new_track = "HCT116_hic_notMerged_insulation_1e3",
						 description = "Insulation track of HCT116 Hic taken from Rao et al. 2017",
						 res=1e3)

#
######	step 3.1.2 - Use insulation track to define TAD borders.
#
gtrack.2d.get_insu_doms <- function(insu_track, thresh, iterator=500){
    doms <- gscreen(sprintf("is.na(%s) | %s > %f", insu_track, insu_track, thresh), iterator = iterator)
    return(doms)
}
#
cmem_hic_domains <- function(insulation_track, thresh, min_domain_size = 50000) {
    domains <- gtrack.2d.get_insu_doms(insulation_track, thresh)
    domains <- domains %>% filter(end - start > min_domain_size)    
    message(glue('relying on {dim(domains)[1]}domains'))
    return(domains)
}
#
cmem_divide_genome_into_doms <- function(insu_track, min_domain_size=50000, ins_dom_thresh=0.1){
    domains <- cmem_hic_domains(insu_track, thresh=gquantiles(insu_track, percentiles=c(ins_dom_thresh)), min_domain_size=min_domain_size)
    return(domains)    
}
#
insulation_track_name = "HCT116_hic_notMerged_insulation_1e3"
tads = cmem_divide_genome_into_doms(insu_track = insulation_track_name,
									 min_domain_size = 50000, 
									 ins_dom_thresh = 0.1)

									 
#	
######  step 3.1.3 - plot examples for insulation and defined TADs over genomic regions.
#
hct_ins_1e3 = gextract(insulation_track_name, intervals = gintervals.all())
#
cmem_plot_insulation_TAD_example (tads_table = tads,
								  insulation_table = hct_ins_1e3,
								  example_chrom = "chr7",
								  example_start = 1000000,
								  example_end = 8500000,
								  insulation_color = "#00000075",
								  TAD_border_color = "red",
								  coords_to_show = c(2e6,4e6,6e6,8e6),
								  coords_labels = c("2M","4M","6M","8M"),
								  show_title = TRUE)
#
cmem_plot_insulation_TAD_example (tads_table = tads,
								  insulation_table = hct_ins_1e3,
								  example_chrom = "chr12",
								  example_start = 1000000,
								  example_end = 10000000,
								  insulation_color = "#00000075",
								  TAD_border_color = "red",
								  coords_to_show = c(2e6,4e6,6e6,8e6),
								  coords_labels = c("2M","4M","6M","8M"),
								  show_title = TRUE)



################################################################################################
######
###### module 3.2 - screen for clonal TAD-memory over shuffled controls in HCT116, H1299, A549 and WI38 cells.
######
######		in each of the cell lines (based on the HCT116 TAD-borders):
######	step 3.2.1 - filter TADs with low expression, variance or few genes, and compute correlation between gene expression and their total TADs output
######              (after subtracting each gene's expression from its TAD when computing correlation between them).
######	step 3.2.2 - Compare distribution of gene-to-TAD expression correlation, over shuffled control where the TAD-gene structure is randomized.
######  step 3.2.3 - Plot interesting examples with Shaman package.
######
################################################################################################



#
######	step 3.2.1 - filter TADs with low expression, variance or few genes, and compute correlation between gene expression and their total TADs output
#              (after subtracting each gene's expression from its TAD when computing correlation between them).
#
# load HCT116 gene-TAD corrs
gtad_cor_hct116 = read.csv("HiC_data/HCT116_clones_gene_to_TAD_correlation.csv", row.names = 1)
gene_expression_gtad_hct116 = as.data.frame(read.csv("HiC_data/HCT116_clones_TAD_gene_expression.csv", row.names = 1))[,1]
names(gene_expression_gtad_hct116) = rownames(gtad_cor_hct116)
#
# load TAD-gene key
tad_gene_key = as.data.frame(apply(read.csv("HiC_data/gene_to_TAD_key.csv", row.names = 1),2,as.character))
rownames(tad_gene_key) = tad_gene_key$gene
#
# generate_stats on gene-TAD correlations
gtad_stats_hct116 = cmem_generate_gtad_corr_stats (gtad_cor_matrix = gtad_cor_hct116,
							   gene_expression_gtad = unlist(gene_expression_gtad_hct116),
							   tad_to_gene_convertion_table = tad_gene_key[rownames(gtad_cor_hct116),])
#
# plot correlation ranks distribution (of genes to their self TADs vs all other TADs)
cmem_plot_TAD_screen_rankHist (gtad_stats_table = gtad_stats_hct116,
							   fig_nm = "gene_to_TAD_correlationRanks_HCT116.png",
							   plot_cex_text = 1.5,
							   plot_xlab = "rank of corr. to self-TAD (HCT116)",
							   plot_breaks = 16)

							   
							   
							   
######	step 3.2.2 - Compare distribution of gene-to-TAD expression correlation, over shuffled control where the TAD-gene structure is randomized.
#
# generate stats on gene-to TAD correlation, when TAD structure is randomized
# load HCT116 gene-TAD corrs
permuted_gtad_cor_hct116 = read.csv("HiC_data/HCT116_clones_gene_to_permuteTAD_correlation.csv", row.names = 1)
gene_expression_permuted_gtad_hct116 = as.data.frame(read.csv("HiC_data/HCT116_clones_permuteTAD_gene_expression.csv", row.names = 1))[,1]
names(gene_expression_permuted_gtad_hct116) = rownames(permuted_gtad_cor_hct116)
#
# generate_stats on permuted gene-TAD structure, where TADs have the same number of genes but random ones.
tad_gene_key_hct_perm = as.data.frame(apply(read.csv("HiC_data/gene_to_TAD_key_HCT116permute.csv", row.names = 1),2,as.character))
rownames(tad_gene_key_hct_perm) = tad_gene_key_hct_perm$gene
gtad_stats_hct116_permute = cmem_generate_gtad_corr_stats (gtad_cor_matrix = permuted_gtad_cor_hct116,
							   gene_expression_gtad = unlist(gene_expression_permuted_gtad_hct116),
							   tad_to_gene_convertion_table = tad_gene_key_hct_perm)							   
# plot permuted correlation ranks distribution (of genes to their self TADs vs all other TADs)
cmem_plot_TAD_screen_rankHist (gtad_stats_table = gtad_stats_hct116_permute,
							   fig_nm = "gene_to_TAD_correlationRanks_HCT116_permute.png",
							   plot_cex_text = 1.5,
							   plot_xlab = "rank of corr. to self-TAD (HCT116 permuted)",
							   plot_breaks = 16)
# plot interesting gene-to-TAD correlations over permuted screen
TADS_screen_metadata = read.csv("HiC_data/Supplementary_Table9_TADS_metadata.csv",row.names=1)
cmem_plot_TAD_screen_over_permutation (fig_nm_ecdf = "HCT116_TAD_screen_over_permutation_ecdf.png",
									   fig_nm_scatter = "HCT116_TAD_screen_intersting_corrs.png",
									   gtad_stats_orig = gtad_stats_hct116,
									   gtad_stats_permute = gtad_stats_hct116_permute,
									   FDR5_genes = rownames(TADS_screen_metadata)[TADS_screen_metadata$HCT116_FDR_0_05],
									   FDR25_genes = rownames(TADS_screen_metadata)[TADS_screen_metadata$HCT116_FDR_0_25],
									   FDR5_color = "deepskyblue3",
									   FDR25_color = "deepskyblue1",
									   plot_ecdf_cex_text = 1.2,
									   ecdf_ylim = c(0,300),
									   ecdf_xlim = c(0.03,0.5),
									   ecdf_lwd = 4,
									   plot_ecdf_w = 280,
									   plot_ecdf_h = 220,
									   plot_ecdf_margins = c(5,5,1,1),
									   show_ecdf_legend = TRUE,
									   plot_xy_w = 220,
									   plot_xy_h = 220,
									   plot_xy_margins = c(5,5,1,1),
									   plot_xy_cex_text = 1.2,
									   plot_xy_xlim=c(0,5),
									   plot_xy_xlab = "gene expression",
									   plot_xy_ylab = "gene-TAD corr.",
									   plot_xy_dots_cex = 0.3
									   )

# repeat TAD-memory screen for A549 clones
gtad_stats_a549 = cmem_generate_hic_stats (gene_tad_corr_csv = "HiC_data/A549_clones_gene_to_TAD_correlation.csv",
										   gene_expression_csv = "HiC_data/A549_clones_TAD_gene_expression.csv",
										   permuted_gene_tad_corr_csv = "HiC_data/A549_clones_gene_to_permuteTAD_correlation.csv",
										   permuted_gene_expression_csv = "HiC_data/A549_clones_permuteTAD_gene_expression.csv",
										   real_tad_to_gene_convertion_csv = "HiC_data/gene_to_TAD_key.csv",
										   permuted_tad_to_gene_convertion_csv = "HiC_data/gene_to_TAD_key_A549permute.csv")
#
# plot correlation ranks distribution (of genes to their self TADs vs all other TADs)
cmem_plot_TAD_screen_rankHist (gtad_stats_table = gtad_stats_a549$gtad_stats_orig,
							   fig_nm = "gene_to_TAD_correlationRanks_A549.png",
							   plot_cex_text = 1.5,
							   plot_xlab = "rank of corr. to self-TAD (A549)",
							   plot_breaks = 16)
#
cmem_plot_TAD_screen_rankHist (gtad_stats_table = gtad_stats_a549$gtad_stats_perm,
							   fig_nm = "gene_to_TAD_correlationRanks_permuted_A549.png",
							   plot_cex_text = 1.5,
							   plot_xlab = "rank of corr. to self-TAD (A549 permuted)",
							   plot_breaks = 16)							   
#							   
cmem_plot_TAD_screen_over_permutation (fig_nm_ecdf = "A549_TAD_screen_over_permutation_ecdf.png",
									   fig_nm_scatter = "A549_TAD_screen_intersting_corrs.png",
									   gtad_stats_orig = gtad_stats_a549$gtad_stats_orig,
									   gtad_stats_permute = gtad_stats_a549$gtad_stats_perm,
									   FDR5_genes = rownames(TADS_screen_metadata)[TADS_screen_metadata$A549_FDR_0_05],
									   FDR25_genes = rownames(TADS_screen_metadata)[TADS_screen_metadata$A549_FDR_0_25],
									   FDR5_color = "deepskyblue3",
									   FDR25_color = "deepskyblue1",
									   plot_ecdf_cex_text = 1.2,
									   ecdf_ylim = c(0,300),
									   ecdf_xlim = c(0.03,0.5),
									   ecdf_lwd = 4,
									   plot_ecdf_w = 280,
									   plot_ecdf_h = 220,
									   plot_ecdf_margins = c(5,5,1,1),
									   show_ecdf_legend = TRUE,
									   plot_xy_w = 220,
									   plot_xy_h = 220,
									   plot_xy_margins = c(5,5,1,1),
									   plot_xy_cex_text = 1.2,
									   plot_xy_xlim=c(0,5),
									   plot_xy_xlab = "gene expression",
									   plot_xy_ylab = "gene-TAD corr.",
									   plot_xy_dots_cex = 0.3
									   )

# repeat TAD-memory screen for WI38 clones
gtad_stats_wi38 = cmem_generate_hic_stats (gene_tad_corr_csv = "HiC_data/WI38_clones_gene_to_TAD_correlation.csv",
										   gene_expression_csv = "HiC_data/WI38_clones_TAD_gene_expression.csv",
										   permuted_gene_tad_corr_csv = "HiC_data/WI38_clones_gene_to_permuteTAD_correlation.csv",
										   permuted_gene_expression_csv = "HiC_data/WI38_clones_permuteTAD_gene_expression.csv",
										   real_tad_to_gene_convertion_csv = "HiC_data/gene_to_TAD_key.csv",
										   permuted_tad_to_gene_convertion_csv = "HiC_data/gene_to_TAD_key_WI38permute.csv")
#
# plot correlation ranks distribution (of genes to their self TADs vs all other TADs)
cmem_plot_TAD_screen_rankHist (gtad_stats_table = gtad_stats_wi38$gtad_stats_orig,
							   fig_nm = "gene_to_TAD_correlationRanks_WI38.png",
							   plot_cex_text = 1.5,
							   plot_xlab = "rank of corr. to self-TAD (WI38)",
							   plot_breaks = 16)
#
cmem_plot_TAD_screen_rankHist (gtad_stats_table = gtad_stats_wi38$gtad_stats_perm,
							   fig_nm = "gene_to_TAD_correlationRanks_permuted_WI38.png",
							   plot_cex_text = 1.5,
							   plot_xlab = "rank of corr. to self-TAD (WI38 permuted)",
							   plot_breaks = 16)							   
#							   
cmem_plot_TAD_screen_over_permutation (fig_nm_ecdf = "WI38_TAD_screen_over_permutation_ecdf.png",
									   fig_nm_scatter = "WI38_TAD_screen_intersting_corrs.png",
									   gtad_stats_orig = gtad_stats_wi38$gtad_stats_orig,
									   gtad_stats_permute = gtad_stats_wi38$gtad_stats_perm,
									   FDR5_genes = rownames(TADS_screen_metadata)[TADS_screen_metadata$WI38_FDR_0_05],
									   FDR25_genes = rownames(TADS_screen_metadata)[TADS_screen_metadata$WI38_FDR_0_25],
									   FDR5_color = "deepskyblue3",
									   FDR25_color = "deepskyblue1",
									   plot_ecdf_cex_text = 1.2,
									   ecdf_ylim = c(0,200),
									   ecdf_xlim = c(0.03,0.35),
									   ecdf_lwd = 4,
									   plot_ecdf_w = 280,
									   plot_ecdf_h = 220,
									   plot_ecdf_margins = c(5,5,1,1),
									   show_ecdf_legend = TRUE,
									   plot_xy_w = 220,
									   plot_xy_h = 220,
									   plot_xy_margins = c(5,5,1,1),
									   plot_xy_cex_text = 1.2,
									   plot_xy_xlim=c(0,5),
									   plot_xy_xlab = "gene expression",
									   plot_xy_ylab = "gene-TAD corr.",
									   plot_xy_dots_cex = 0.3
									   )

								   
#
######  step 3.2.3 - Plot interesting examples with Shaman package.
#
# showing one example - the B-globins region
tad_to_plot = "TAD_523"
# plot expression in clones
cmem_tad_exp_plot (tad_id = tad_to_plot,
				   expression_table = cmem_umi_PAR,
				   expression_table_norm = cmem_umi_PAR_norm,
				   min_UMI_per_clone = 1e4,
				   tad_to_gene_convertion = tad_gene_key,
				   reg = 0.1,
				   heatmap_fontsize = 21,
				   heatmap_cellheight = 18,
				   heatmap_cellwidth = 1.7,
				   heatmap_colors = colorRampPalette( rev(brewer.pal(n = 11, name ="BrBG")))(100),
				   base_dir = "./",
				   lfp_z_lim = 2,
				   min_gene_expression = 0
				   ) 
# plot normalized contacts-map of the region (over shuffling, considering the local cis-decay)
# IMPROTANT: you must install "shaman" package in order to run this part of analysis, and also download the observed contacts track!
#
intervs_2d_format_tads = read.csv ("HiC_data/TAD_inetrvals_2D_format.csv",row.names = 1)
hg19_gene_intervs = read.csv ("HiC_data/hg19_gene_intervals.csv", row.names = 1)
raw_HCT116_observed_contacts_track_name = "HCT116_hic_notMerged"
cmem_tad_shaman_plot	(tad_id = tad_to_plot,
						 expand=5e5,
						 TAD_intervals_in_2D_format = intervs_2d_format_tads,
						 map_point_size = 1,
						 genome_version = "hg19",
						 gene_intervals = hg19_gene_intervs,
						 HCT116_observed_contacts_track_name = raw_HCT116_observed_contacts_track_name,
						 base_dir = "./",
						 tads_table = tads_table
						 )
#
## plot color-keys scale-bar for examples of HiC SHAMAN-maps over expression of specific regions
expression_lfp_z_lim = 2
expression_colors = colorRampPalette( rev(brewer.pal(n = 11, name ="BrBG")))(100)
shaman_hic_colors = shaman_score_pal()
cmem_plot_colorBars_HiC_examples (hic_vals = c(0:200),
								  hic_colors = shaman_hic_colors,
								  hic_values_to_show = c(1,51,101,151,201),
								  hic_values_Xcoord = 15,
								  hic_bar_title = "norm. contacts score",
								  hic_fig_nm = "local_region_shaman_colorkey.png",
								  expression_vals = round(seq( (-1)*expression_lfp_z_lim,
																1*expression_lfp_z_lim,length.out=100),digits=1),
								  expression_colors = expression_colors,
								  expression_values_to_show = c(1,50,100),
								  expression_values_Xcoord = 15,
								  expression_bar_title = "norm. expression",
								  expression_fig_nm ="local_region_expression_colorkey.png")
								 
								  
								 



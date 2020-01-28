##
##
## methylation.r
## Date: 22.01.2020
## This file summarizes main steps use for DNA-methylation analysis of HCT116 single-cells (wgPBAT) and single-cell-derived clonal populations (capture-PBAT) in Nat. gen. Meir et al. 2020
##
##

#load dependency packages and project functions
source("Meir_et_al_2020_nat_gen_functions.r")

# Download tables used for this analysis from our repositories
cmem_download_meth_tables()

#load whole genome DNA-methylation of clones
cmem_clones_meth = read.table("methylation_data/clonesMethylome_captPBAT_HCT116_WT_methCalls.txt", sep="\t", header = TRUE)
cmem_clones_cov = read.table("methylation_data/clonesMethylome_captPBAT_HCT116_WT_coverage.txt", sep="\t", header = TRUE)

#load CRC_V3 capture intervals
cmem_crcV3_intervals = read.csv("/net/mraid14/export/data/users/zoharme/cloneMem/scripts/gitHub_scripts/methylation_data/CRC_v3_probeSet_intervals_hg19.csv", row.names = 1)

#load averaged genomic CG content per 500 bases
cg500 = read.csv("methylation_data/CGcontent500_per_CpG_hg19.csv", row.names = 1)

#load wg-sc methylation marginals per cell
cmem_singles_wgPBAT = read.csv("methylation_data/scMethylome_wgPBAT_HCT116.csv")

#load clones bulk MARS-seq umi_tables
cmem_load_clones()

#load HCT116 clones metadata
cmem_clones_metadata = read.csv("methylation_data/Supplementary_Table5_clonesMetadata.csv",row.names = 1)

#load cc independent gene modules
cc_independent_geneModules = read.csv("expression_data/Supplementary_Table3_ccIndependet_geneModules.csv",row.names="gene")


#download and initialize a hg19 environment for analysis with "misha" packages
# NOTE: this will download the assembly and other files into /hg19 directory, that will be generated in current directory.
# NOTE: this will take ~5 minutes.
# NO WORRIES: If you already have an hg19 database, it wouldn't download it again.
cmem_generate_hg19_misha_db ()

################################################################################################
######
###### part 2: Analysis of single-cell-derived clones DNA methylation
######
######		module 2.1 - analyze clonal average methylation 
###### 		module 2.2 - analyze methylation dynamics of specific loci
################################################################################################




################################################################################################
######
###### module 2.1: analyze clonal average methylation 
######
######	step 2.1.1 - Analyze off-target methylation (output of loci not covered by the CRC-oriented probe set), defining HCG and LCG methylation
######  step 2.1.2 - Analyze on-target methylation (regions covered by multiple clones)
######	step 2.1.3 - compute expression correlation to clonal LCG and HCG, KNN normalization by HCG trend
######  step 2.1.4 - Cluster cross-correlation of clonal methylation and transcription.
######
################################################################################################

#
######	step 2.1.1 - Analyze clonal off-target methylation (output of loci not covered by the CRC-oriented probe set)
#
# get high-CpG (HCG) and low-CpG (LCG) off-target clonal methylation
minimal_dist_from_probe = 100
minimal_offTarget_coverage = 5e4;
distances_to_pbat_intervs = gintervals.neighbors (cmem_clones_cov[, c("chrom","start","end")], cmem_crcV3_intervals)
cmem_clones_cov_offtar = cmem_clones_cov [distances_to_pbat_intervs$dist > minimal_dist_from_probe, ]
cmem_clones_meth_offtar = cmem_clones_meth [distances_to_pbat_intervs$dist > minimal_dist_from_probe, ]
cg500_offtar = cg500[distances_to_pbat_intervs$dist > minimal_dist_from_probe,]
cmem_clones_meth_offtar_f = cbind(cmem_clones_meth_offtar[, 1:4],
								  cmem_clones_meth_offtar[, 5:ncol(cmem_clones_meth_offtar)] [,colSums(cmem_clones_cov_offtar[, 5:ncol(cmem_clones_cov_offtar)],na.rm=T) >= minimal_offTarget_coverage])
cmem_clones_cov_offtar_f = cbind(cmem_clones_cov_offtar[, 1:4],
								  cmem_clones_cov_offtar[, 5:ncol(cmem_clones_cov_offtar)] [,colSums(cmem_clones_cov_offtar[, 5:ncol(cmem_clones_cov_offtar)],na.rm=T) >= minimal_offTarget_coverage])
offtarget_meth_tab = cmem_clones_meth_offtar_f[,5:ncol(cmem_clones_meth_offtar_f)]	
offtarget_cov_tab = cmem_clones_cov_offtar_f[,5:ncol(cmem_clones_cov_offtar_f)]		
#			
#
# plot LCG clonal methylation (highlight VIM-high clones)
#
cmem_plot_clonalMeth_by_cgCont (meth_table = offtarget_meth_tab,
								cov_table = offtarget_cov_tab,
								cg500_table = cg500_offtar,
								clones_metadata = cmem_clones_metadata,
								dots_cex = 1,
								text_cex = 1.5,
								cpg_range_x = c(0,0.015),
								cpg_range_y = c(0.015,0.03),
								lab_x = "0-1.5% (CG cont.)",
								lab_y = "1.5-3% (CG cont.)",
								include_low = FALSE,
								plot_margins = c(5,5,1,1),
								plot_h = 280,
								plot_w = 280,
								dots_color = "#00000060",
								fig_nm = "./LCG_methylation.png",
								highlight_EMT_clones = TRUE,
								legend_x = 0.455,
								legend_y = 0.84,
								show_EMT_clones_legend = TRUE,
								EMT_clones_color = "red",
								return_clonal_averages = FALSE
								)						
#
# plot HCG clonal methylation (highlight VIM-high clones)
#
cmem_plot_clonalMeth_by_cgCont (meth_table = offtarget_meth_tab,
								cov_table = offtarget_cov_tab,
								cg500_table = cg500_offtar,
								clones_metadata = cmem_clones_metadata,
								dots_cex = 1,
								text_cex = 1.5,
								cpg_range_x = c(0.07,0.09),
								cpg_range_y = c(0.09,0.15),
								lab_x = "7-9% (CG cont.)",
								lab_y = "9-15% (CG cont.)",
								include_low = FALSE,
								plot_margins = c(5,5,1,1),
								plot_h = 280,
								plot_w = 280,
								dots_color = "#00000060",
								fig_nm = "./HCG_methylation.png",
								highlight_EMT_clones = TRUE,
								legend_x = 0.455,
								legend_y = 0.84,
								show_EMT_clones_legend = TRUE,
								EMT_clones_color = "red",
								return_clonal_averages = FALSE
								)	
#
# plot lack of corr. between clonal LCG and HCG
#
cmem_offtar_trends = cmem_plot_clonalMeth_by_cgCont (meth_table = offtarget_meth_tab,
													 cov_table = offtarget_cov_tab,
													 cg500_table = cg500_offtar,
													 clones_metadata = cmem_clones_metadata,
													 dots_cex = 1,
													 text_cex = 1.5,
													 cpg_range_x = c(0,0.03),
													 cpg_range_y = c(0.07,0.15),
													 lab_x = "0-3% (LCG)",
													 lab_y = "7-15% (HCG)",
													 include_low = FALSE,
													 plot_margins = c(5,5,1,1),
													 plot_h = 280,
													 plot_w = 280,
													 dots_color = "#00000060",
													 fig_nm = "./HCG_vs_LCG_clonalMethylation.png",
													 highlight_EMT_clones = FALSE,
													 exclude_EMT_clones = TRUE,
													 legend_x = 0.455,
													 legend_y = 0.84,
													 show_EMT_clones_legend = FALSE,
													 return_clonal_averages = TRUE
													 )	
#
# test clonal variation in HCG and LCG over shuffled matrix (meth calls are randomly assigned to clones)
#
HCG_mat = cmem_permute_by_cgCont (meth_table = offtarget_meth_tab,
								  cov_table = offtarget_cov_tab,
								  cg500_table = cg500_offtar,
								  clones_metadata = cmem_clones_metadata,
								  cpg_range = c(0.07,0.15),
								  exclude_EMT_clones = TRUE,
								  include_low = FALSE
								  )
#
LCG_mat = cmem_permute_by_cgCont (meth_table = offtarget_meth_tab,
								  cov_table = offtarget_cov_tab,
								  cg500_table = cg500_offtar,
								  clones_metadata = cmem_clones_metadata,
								  cpg_range = c(0,0.03),
								  exclude_EMT_clones = TRUE,
								  include_low = FALSE
								  )
#
cmem_compare_distributions (original = LCG_mat$avg_orig,
							permuted = LCG_mat$avg_perm,
							original_color = "black",
							permuted_color = "grey63",
							plot_margins = c(5,5,1,1),
							plot_w = 280,
							plot_h = 280,
							plot_lwd = 4,
							legend_lwd = 2,
							plot_legend = TRUE,
							legend_labels = c("Low CpG","Low CpG (perm.)"),
							plot_ylim = c(0,300),
							cex_text = 1.5,
							fig_nm = "./LCG_clonal_variation.png")
#							
cmem_compare_distributions (original = HCG_mat$avg_orig,
							permuted = HCG_mat$avg_perm,
							original_color = "black",
							permuted_color = "grey63",
							plot_margins = c(5,5,1,1),
							plot_w = 280,
							plot_h = 280,
							plot_lwd = 4,
							legend_lwd = 2,
							plot_legend = TRUE,
							legend_labels = c("High CpG","High CpG (perm.)"),
							plot_ylim = c(0,220),
							cex_text = 1.5,
							fig_nm = "./HCG_clonal_variation.png")
							

							
#
######  step 2.1.2 - Analyze on-target methylation (regions covered by multiple clones)
#
# load on-target methylation of clones and single-cells (pooled)
cmem_clones_cov_ontar = read.csv("methylation_data/clonesMethylome_captPBAT_HCT116_WT_onTarget_coverage.csv", row.names=1)
cmem_clones_meth_ontar = read.csv("methylation_data/clonesMethylome_captPBAT_HCT116_WT_onTarget_methCalls.csv", row.names=1)
sc_meth_pool = read.csv("methylation_data/cellsMethylome_wgPBAT_HCT116_WT_onTarget_marginals.csv", row.names = 1)
#
# compare average methylation of HCT116 cells and clones over capture regions where both are covered
#
min_ontar_coverage = 5e4
covered_ontar_clones = colnames(cmem_clones_cov_ontar[4:ncol(cmem_clones_cov_ontar)])[colSums(cmem_clones_cov_ontar[,4:ncol(cmem_clones_cov_ontar)]) >= min_ontar_coverage]
cmem_clones_cov_ontar_f = cbind(cmem_clones_cov_ontar[,1:3],cmem_clones_cov_ontar[,covered_ontar_clones])
cmem_clones_meth_ontar_f = cbind(cmem_clones_meth_ontar[,1:3],cmem_clones_meth_ontar[,covered_ontar_clones])
#
clones_pooled_region_meth = cmem_aggregate_CpG_stats_by_intervals(cpg_table = cmem_clones_meth_ontar_f, intervals = crc_ints)
clones_pooled_region_cov = cmem_aggregate_CpG_stats_by_intervals(cpg_table = cmem_clones_cov_ontar_f, intervals = crc_ints) 																  
covered_regions = intersect(paste0(sc_meth_pool$chrom,"_",sc_meth_pool$start,"_",sc_meth_pool$end),rownames(clones_pooled_region_cov))
clones_meth_regions = data.frame(meth = rowSums(clones_pooled_region_meth[covered_regions,setdiff(colnames(clones_pooled_region_meth),c("chrom1","start1","end1"))]),
								 cov = rowSums(clones_pooled_region_cov[covered_regions,setdiff(colnames(clones_pooled_region_cov),c("chrom1","start1","end1"))]))
cells_meth_regions = sc_meth_pool [match(paste0(sc_meth_pool$chrom,"_",sc_meth_pool$start,"_",sc_meth_pool$end),covered_regions),]
cmem_compare_cells_clones_methylation (min_cov_loci_cells = 50,
									   min_cov_loci_clones = 500,
									   clones_meth_df = clones_meth_regions,
									   cells_meth_df = cells_meth_regions,
									   meth_bins = 10,
									   base_dir = "./",
									   clones_color = "antiquewhite2",
									   cells_color = "black",
									   density_show_legend = TRUE
									   )


									   
#									   
######	step 2.1.3 - compute expression correlation to clonal LCG and HCG, KNN normalization by HCG trend
#
# get clones with off-target coverage and expression coverage (excluding the distinctly hypomethylated VIM-high clones)
clones_covered_by_expression_and_offTargetMethylation = intersect(names(HCG_mat$avg_orig),colnames(cmem_umi_PAR))
expression_offtarMethylation_df = data.frame(clone = clones_covered_by_expression_and_offTargetMethylation,
											 HCG = HCG_mat$avg_orig[clones_covered_by_expression_and_offTargetMethylation],
											 HCG_perm = HCG_mat$avg_perm[clones_covered_by_expression_and_offTargetMethylation],
											 LCG = LCG_mat$avg_orig[clones_covered_by_expression_and_offTargetMethylation],
											 LCG_perm = LCG_mat$avg_perm[clones_covered_by_expression_and_offTargetMethylation])
cmem_compute_methTrend_geneCorr (clones_offTarget_methylation_table = expression_offtarMethylation_df,
								 clones_umi_table = cmem_umi_PAR,
								 clones_umi_table_norm = cmem_umi_PAR_norm,
								 min_UMI_per_gene = 5,
								 base_dir = "./",
								 color_original = "black",
								 color_permuted = "grey63",
								 plot_legend = TRUE,
								 legend_lwd = 2,
								 cor_method = "spearman",
								 blacklist_genes = grepl("^MTRNR2|MALAT1|NEAT1|BC010924",rownames(cmem_umi_PAR)),
								 return_gcors = FALSE
								 )
#
# generate per-region stats on each clone, keeping only regions with minimal coverage and variance
#
avg_agg = cmem_get_regional_average (meth_table = clones_pooled_region_meth[4:ncol(clones_pooled_region_meth)],
									 cov_table = clones_pooled_region_cov[4:ncol(clones_pooled_region_cov)],
									 clones_umi_table = cmem_umi_PAR,
									 minimal_regional_coverage_per_clone = 10,
									 minimal_regional_convering_clones = 80,
									 maximal_avg_methylation = 0.98,
									 minimal_avg_methylation = 0.0001,
									 #minimal_avg_methylation = 0.1,
									 minimal_umi_per_clone = 1e4
									 )	

#
# normalize offtarget effect by subtracting from each clone the average methylation of its K=25 nearest neighbors in LCG and HCG space
#
clones_covered_by_expression_and_offTargetMethylation = intersect( names(HCG_mat$avg_orig), colnames(avg_agg))
offTarget_features = data.frame(HCG = HCG_mat$avg_orig[clones_covered_by_expression_and_offTargetMethylation],
					            LCG = LCG_mat$avg_orig[clones_covered_by_expression_and_offTargetMethylation]) %>% t() %>% t() %>% scale(center=TRUE, scale=TRUE)
offTarget_knn_mat = cmem_get_norm_knn_matrix(offTarget_features, k=25)
#plot 3 examples for KNN logic
for (i in 1:3){
cmem_generate_random_KNN_example (clones_offTarget_methylation_table = expression_offtarMethylation_df [clones_covered_by_expression_and_offTargetMethylation,],
								  knn_mat = offTarget_knn_mat,
								  fig_nm = paste0("./KNN_feature_Selection_example",as.character(i),".png"),
								  plot_legend = TRUE,
								  selected_clone_color = "blue", 
								  knn_clones_color = "red")}
avg_agg_f = avg_agg[, intersect(colnames(avg_agg),rownames(offTarget_features))]
avg_agg_norm = cmem_normalize_knn(avg_agg_f, offTarget_knn_mat)



#
######  step 2.1.4 - Cluster cross-correlation of clonal methylation and transcription, load cross-correlation use for plot.
#
# compute cross correlation between methylation and expression, and filter boring regions or genes
cc_mat = cmem_get_ExprMeth_crossCorrelation (exp_mat = cmem_umi_PAR[,colnames(avg_agg_norm)],
											 exp_mat_norm = cmem_umi_PAR_norm[,colnames(avg_agg_norm)],
											 exp_ds_value = 1e5,
											 meth_mat = avg_agg_norm,
											 min_gene_cor = 0.28,
											 min_gene_expression_normVar = 0.15,
											 eps = 0.1)								
# filter out boring clusters of gene and regions after initial clustering
cc_mat_f = as.matrix(read.csv ("/net/mraid14/export/data/users/zoharme/cloneMem/scripts/gitHub_scripts/methylation_data/HCT116clones_methylation_expression_crossCorrelation.csv", row.names = 1, check.names = FALSE))
# filter out genes and regions that do not maintain numerous interesting corrs
cc_mat_f_plot = cmem_filter_cc_mat (crossCorr_mat_to_filter = cc_mat_f,
							        max_corr_to_regions = (-0.15),
							        min_n_regions = 55,
							        max_corr_to_genes = (-0.15),
							        min_n_genes = 10)
# plot cross-correlation matrix
cc_clusters = cmem_plot_crossCorr (cc_mat = cc_mat_f_plot,
								   genes_to_highlight = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$HCT116_WT=="Epithelial"]),"[[",1)),
								   highlight_color = "green1",
								   highlight_title = "Epithelial Module",
								   base_dir = "./",
								   fig_nm = "crossCorrelation_matrix.png",
								   plot_w = 1300,
								   plot_h = 1950,
								   plot_cellheight = 20,
								   plot_fontsize = 24,
								   k_genes = 8,
								   k_regions = 5,
								   limit_corr_value = 0.3,
								   retrieve_clusters = TRUE
								   )	
# plot example of actual DNA methylation (without KNN-norm) at epithelial-related regions, over epithelial expression 
ep_regions = names(cc_clusters$cc_regions_clusters)[cc_clusters$cc_regions_clusters==4]
cmem_plot_expression_and_methylation (meth_table = clones_pooled_region_meth[,colnames(avg_agg_f)],
									  cov_table = clones_pooled_region_cov[,colnames(avg_agg_f)],
									  genes_to_plot = unlist(lapply(na.omit (rownames(cc_independent_geneModules) [cc_independent_geneModules$HCT116_WT=="Epithelial"]),"[[",1)),,
									  expression_table = cmem_umi_PAR,
									  expression_table_norm = cmem_umi_PAR_norm,
									  label_expression = "Ep. module expression",
									  label_methylation = "Ep. regions methylation",
									  fig_nm = "epithelial_regions_methylation_over_expression.png",
									  reg = 1,
									  show_corr = TRUE,
									  corr_position = "bottomleft",
									  regions_to_plot = ep_regions)
									  
									
									
									





						
################################################################################################
######
###### module 2.2: analyze methylation dynamics of specific loci
######
######	step 2.2.1 - For most covered CpGs, compute epi-polymorphism and compare to mixture and persistency models.
######  step 2.2.2 - Define "hot" and "cold" regions.
######	step 2.2.3 - Stratify "hot" and "cold" regimes by LCG and HCG off-target trends in clones.
######
################################################################################################


#
######	step 2.2.1 - For most covered CpGs, compute epi-polymorphism and compare to mixture and persistency models.
#
# generate a matrix of most covered CpGs in clones (sampled at least 6 times in at least 60 clones, downsampled to 6 samples in 60 random clones)
meth6_60 = read.csv("methylation_data/clonesMethylome_captPBAT_HCT116_WT_downsampled_6calls_60clones.csv", row.names = 1)
meth6_60_clean = apply(cmem_meth_ds6_60[,4:ncol(cmem_meth_ds6_60)], 2, as.numeric)
meth6_60_metadata = read.csv("methylation_data/clonesMethylome_captPBAT_HCT116_WT_downsampled_6calls_60clones_metadata.csv", row.names = 1)
# according to average_methylation of each locus, compute OBSERVED variance, as well as expected variance according to MIXTURE (no-clonality) and PERSISTENCY (Hardy-Weinberg mode of epialleles heritablility) models.
average_m = rowMeans(meth6_60_clean)/6
var_observed = apply(meth6_60_clean,1,var)
var_persistency = sapply(average_m,cmem_biallelic_p6)
var_mixture = average_m*(1-average_m)*6


#
######  step 2.2.2 - Define "hot" (non-persistent) and "cold" (persistent) CpGs and plot epi-polymorphism.
#
deviation_from_persistency = var_persistency - var_observed
epipolymorphism_df = data.frame(m = average_m, 
						 obs = var_observed,
						 mix = var_mixture,
						 per = var_persistency,
						 dev_per = deviation_from_persistency)
partially_methylated_cpgs =  (meth6_60_metadata$pool_e > 0.3) & (meth6_60_metadata$pool_e < 0.7)
hot_cpgs =  deviation_from_persistency > 3   & partially_methylated_cpgs
cold_cpgs = (abs(deviation_from_persistency) < 0.3) & partially_methylated_cpgs
epipolymorphism_df$regime = c("normal","hot","cold")[ifelse(hot_cpgs,2,ifelse(cold_cpgs,3,1))]
cmem_plot_epipolymorphism (ep_df = epipolymorphism_df,
						   fig_nm = "./comparison_of_DNAmethylation_variance_models.png",
						   hot_color = "red",
						   cold_color = "deepskyblue3",
						   normal_color = "#00000050",
						   observed_trend_color = "grey43",
						   persistancy_model_color = "deepskyblue1",
						   mixture_model_color = "coral",
						   plot_text_cex=1.7,
						   plot_lwd = 3,
						   plot_legend_models = TRUE,
						   plot_legend_hotCold = TRUE
						   )


#
######	step 2.2.3 - Stratify "hot" and "cold" regimes by LCG and HCG off-target trends in clones.
#
# get all methylation calls of hot and cold loci
meth6_60_all_meth = read.csv("methylation_data/clonesMethylome_captPBAT_HCT116_WT_downsampled_6calls_60clones_allMethCalls.csv", row.names = 1)
meth6_60_all_cov = read.csv("methylation_data/clonesMethylome_captPBAT_HCT116_WT_downsampled_6calls_60clones_allCoverage.csv", row.names = 1)
hot_cold_df = cmem_get_hotCold_clonal_average (min_coverage = 100,
											   meth_tab = meth6_60_all_meth,
											   cov_tab = meth6_60_all_cov,
											   ep_df = epipolymorphism_df,
											   exclude_EMT_clones = TRUE,
											   HCG_mat = HCG_mat,
											   LCG_mat = LCG_mat,
											   clones_metadata = cmem_clones_metadata) 
#											 
cmem_plot_clonal_averages_stratified_by_offTargetTrends (clonal_averages_df = hot_cold_df,
														 plot_boxwex = 0.7,
														 hot_color = "red",
														 cold_color = "deepskyblue3",
														 base_dir = "./",
														 plot_text_cex = 1.5)



library(purrr)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)


# prerequisites:
# pr3c00353_si_002.xlsx


if(!dir.exists("figures_tables"))
   dir.create("figures_tables")

files_to_execute <- c(
   "Fig1_basic_stats.R",
   "Fig2_TD_BUP_correlation.R",
   "Fig3_punchcard_heatmaps.R",
   "Fig4_abeta_species.R",
   "Fig5_4abc_APP-7_test_results_spectral_counts_v2.R",
   "Fig6_VGF_cleavage_sites.R",
   "Fig7_WGCNA.R",
   "FigS1_proteome_coverage.R",
   "FigS2_patrie_sol_insol.R",
   "FigS5_haptoglobin.R",
   "FigS6_densityplots_pLDDT_alphafold.R",
   "FigS8_SRM_TMT_corr.R",
   "FigS9_WGCNA_APP_VGF_cluster_abundances.R",
   "TableS1S2_main_supplementary_tables.R")


files_to_execute %>% 
   map(\(x) file.path(script_dir, x)) %>%
   walk(source, echo = TRUE, .progress = TRUE)


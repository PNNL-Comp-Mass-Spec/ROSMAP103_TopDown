
library(purrr)

script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)


# check if the prerequisite files exist in the ./source_data
# toppic_output.RData
# TopPIC_Dynamic_Mods.txt
# ID_008032_8627C6BD.fasta
# final_r1r2r3_peptide_srm_data.RData
# processed_protein_level_msnset.RData


if(!dir.exists("output_data"))
   dir.create("output_data")

files_to_execute <- c(
   "0-common_first_steps.R",
   "1a-final_intensity_pipeline.R",
   "1b-final_spectral_counting_pipeline.R",
   "1c-remove_tailing_APPs.R",
   "2-annotate_with_mods.R",
   "3-add_counts.R",
   "4a-add_test_results_int.R",
   "4b-add_test_results_spectral_counts.R",
   "6a-blanket_test_with_intensities.R",
   "6b-blanket_test_with_spectral_counts.R")

files_to_execute %>% 
   map(\(x) file.path(script_dir, x)) %>%
   walk(source, echo = TRUE, .progress = TRUE)


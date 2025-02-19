
# R code repository for the manuscript on Top-Down LC-MS analysis of the 103 subjects from ROS/MAP cohorts


* The raw mass spectrometry data can be downloaded at [MassIVE](https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp) repository by ***MSV000093728*** ID.

* The clinical and pathological metadata is protected by Data Use Agreement (DUA) of Rush Alzheimer's Disease Center ([RADC](https://www.radc.rush.edu/))

* The current repository provides the code for the trasparency of the data analysis.


## The intended use:

1. `run_processing_pipeline.R` takes the results of [TopPIC](https://www.toppic.org/software/toppic/index.html) MS/MS search and feature finding, FASTA file and number of additional data files and produces two Bioconductor's [`MSnSet`](https://lgatto.github.io/MSnbase/reference/MSnSet-class.html) objects with proteoforms quantified using 1) MS1 intensities and 2) spectral counts.

2. `run_to_make_figures_and_tables.R` recreates most of the figures (both in the main text and the supplementary) and two supplementary tables with the quantitative data on the proteoforms.



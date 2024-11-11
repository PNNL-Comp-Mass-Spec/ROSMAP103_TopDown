library(MSnSet.utils)
library(pcaMethods)
library(tidyverse)
library(glue)
library(magrittr)
library(WGCNA)
options(stringsAsFactors = FALSE)
disableWGCNAThreads() # on RStudio for some reason multithreading does not work


MEDissThres <- 0.15
softPower <- 20
minModuleSize <- 5 # default is 20
deepSplit <- 1 # default is 1
networkType <- "signed" # default
TOMType <- "unsigned" # default
adj_type <- "signed" # default





# Import ------------------------------------------------------------------


(load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData"))




# Prep data ---------------------------------------------------------------

m <- correct_batch_effect_empiricalBayesLM(m,
                                           removed_cov_name = "pmi")

m1 <- m

# Impute
exprs(m1) <- t(completeObs(pca(as(m1,"ExpressionSet"), method="svdImpute",
                               nPcs=min(dim(m1)), center = TRUE)))

# Convert to z-score
exprs(m1) <- sweep(exprs(m1),
                   MARGIN = 1,
                   STATS = apply(exprs(m1), 1, mean, na.rm = TRUE),
                   FUN = "-")

exprs(m1) <- sweep(exprs(m1),
                   MARGIN = 1,
                   STATS = apply(exprs(m1), 1, sd, na.rm = TRUE),
                   FUN = "/")



# WGCNA -------------------------------------------------------------------

## Select transform power ----
powers = 1:30
sft <- pickSoftThreshold(t(exprs(m1)),
                         powerVector = powers,
                         networkType = networkType,
                         verbose = 0)


## Clustering (first round) ----
adjacency.mat <- adjacency(t(exprs(m1)), power = softPower, type = adj_type)

# convert to topological overlap
TOM <- TOMsimilarity(adjacency.mat, TOMType = TOMType)

# convert to distance
dissTOM = 1 - TOM

# hclust
geneTree <- hclust(as.dist(dissTOM), method = "average")



## Dynamic tree cutting ----
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = deepSplit,
                             minClusterSize = minModuleSize)


## Color the modules ----
# Convert numeric lables into colors
library(RColorBrewer)
# dynamicColors <- labels2colors(dynamicMods, colorSeq = brewer.pal(11,'Spectral'))
# this `standardColors()` is a bit nicer as the colors come with names
dynamicColors <- labels2colors(dynamicMods, colorSeq = standardColors())
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


## Decision on Merging modules ----

# Calculate eigengenes
MEList <- moduleEigengenes(t(exprs(m1)), colors = dynamicColors)
MEs <- MEList$eigengenes

# calculate module distances
MEDiss <- 1 - cor(MEs)

# cluster eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# plot eigengene clusters
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")


## Merge modules ----

# merge
merge <- mergeCloseModules(t(exprs(m1)),
                           dynamicColors, 
                           cutHeight = MEDissThres, 
                           verbose = 3)
moduleColors <- merge$colors
fData(m1)$cluster <- moduleColors 


###########################################################################################

mc2 <- MSnSet(as.matrix(t(merge$newMEs)), pData = pData(m1))
featureNames(mc2) <- sub("ME", "", featureNames(mc2))
fData(mc2) <- data.frame(Cluster = featureNames(mc2), row.names = featureNames(mc2))


# Run limma - clusters as single features
phenotype <- c("anye4", "caa_4gp", "ci_num2_gct", "cogn_global_lv",
               "cogng_demog_slope", "cogng_path_slope", "hip_scl_3reg_yn",
               "lb_neo", "sqrt_amyloid", "sqrt_tangles")


phenotype <- c("sqrt_tangles",
         "sqrt_amyloid",
         "lb_neo",
         "dxpark",
         "hip_scl_3reg_yn",
         "cvda_4gp2",
         "cogn_global_lv",
         "cogng_demog_slope",
         # "cogng_path_slope",
         # "cogdx_3gp",
         "ci_num2_gct",
         "ci_num2_mct",
         "caa_4gp",
         "anye2",
         "anye4",
         "age_death",
         "msex",
         "hspath_3reg",
         "tdp_st4",
         "arteriol_scler",
         "parksc_lv")

# same thing, but just adds arteriol_scler


res_list <- vector("list", length(phenotype))
names(res_list) <- phenotype

for (i in phenotype) {
   res <- limma_gen(mc2, model.str = glue("~ {i}"), coef.str = i)
   res_list[[i]] <- res %>%
      rownames_to_column("Cluster") %>%
      mutate(phenotype = i, .before = Cluster)
}

# Limma results
res_list_df <- bind_rows(res_list) %>%
   filter(adj.P.Val < 0.05) %>%
   group_by(Cluster) %>%
   arrange(Cluster, adj.P.Val, .by_group = F)



clstrs <- distinct(res_list_df, Cluster) %>%
   dplyr::rename(cluster = Cluster) %>%
   inner_join(distinct(select(fData(m1), proteoform_id, cluster, firstAA, lastAA)))


length(unique(res_list_df$Cluster))
length(unique(fData(m1)$cluster))






sig_modules <- paste0("ME", unique(res_list_df$Cluster))

# orderings
sig_modules <- c("MEcoral","MEdarkolivegreen4","MEmediumpurple1",
                  "MElightblue4", "MEplum4",
                  "MEyellowgreen","MElightcoral", 
                  "MEdarkmagenta",  "MEdarkseagreen3",
                  "MEskyblue2")

# interpretation
sig_modules_reps <- c("APP-1","APP-2","APP-3",
                 "MAP2", "STMN1",
                 "VGF-1","VGF-2", 
                 "TMSB4X",  "CALM1", "GFAP")

row_lbls <-  paste(sig_modules,
                     "\n(",
                     sig_modules_reps,
                     ")", sep = "")  
   
###########################################################################################


## Module-trait Association ----

# prep traits data
ptypes <- c("sqrt_amyloid", "sqrt_tangles", "anye4", "caa_4gp", "ci_num2_gct", 
            "hip_scl_3reg_yn",
            "cogn_global_lv",
               "cogng_demog_slope", 
            # "cogng_path_slope",
            "arteriol_scler",
               "lb_neo")

datTraits <- pData(m1)[,ptypes]


# calculate correlation p-values
# moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor = cor(merge$newMEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, ncol(m1));

moduleTraitPvalue <- moduleTraitPvalue %>%
   as.data.frame() %>%
   mutate(across(everything(), .fns = \(x){p.adjust(x, n = nrow(moduleTraitPvalue), method = "BH")})) %>%
   as.matrix()

moduleTraitCor <- moduleTraitCor[sig_modules,]
moduleTraitPvalue <- moduleTraitPvalue[sig_modules,]

# trait association heatmap
textMatrix <-  paste(signif(moduleTraitCor, 2),
                     "\n(",
                     signif(moduleTraitPvalue, 1),
                     ")", sep = "");
textMatrix <- matrix(textMatrix, ncol = ncol(moduleTraitCor))



save(clstrs, file = "./output_data/wgcna_clusters.RData")



labelleur <- c("sqrt_amyloid" = "amyloid",
               "sqrt_tangles" = "tangles",
               "anye4" = expression(italic("APOE")~ε4),
               "caa_4gp" = "cerebral amyloid angiopathy",
               "ci_num2_gct" = "macroscopic infarcts",
               "hip_scl_3reg_yn" = "hippocampal sclerosis",
               "cogn_global_lv" = "global cognition",
               "cogng_demog_slope" = "slope of cognitive decline",
               # "cogng_path_slope" = "slope of cognitive decline\nadjusted by pathologies",
               "arteriol_scler" = "arteriolosclerosis",
               "lb_neo" = "neocortical Lewy bodies")
               
png("./figures_tables/Fig7_WGCNA.png", width = 1100, height = 800, pointsize = 20)
par(mar = c(10, 10, 3, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = labelleur[colnames(datTraits)],
               yLabels = sig_modules,
               ySymbols = row_lbls,
               colorLabels = FALSE,
               # colors = greenWhiteRed(50),
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, 
               setStdMargins = FALSE,
               cex.text = 0.75,
               zlim = c(-1,1),
               main = paste("Module-Trait Association (adjusted p-value)"))
dev.off()



# Anye4: APOE E4
# Caa_4gp: cerebral amyloid angiopathy
# Ci_num2_gct: macroscopic infarcts
# Hip_scl..: hippocampal sclerosis
# Cogn_global_lv: Global cognition
# Cogng_demog…: slope of cognitive decline adjusted by demographics
# Cogng_path…: slope of cognitive decline adjusted by pathologies (I suggest removing this as it confuses more than it helps).
# Arteriol_scler: Arteriolosclerosis
# Lb_neo: neocortical Lewy bodies


























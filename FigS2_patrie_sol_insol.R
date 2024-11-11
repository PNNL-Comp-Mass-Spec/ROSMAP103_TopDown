
library(readxl)
library(dplyr)
library(tidyverse)
library(ggbreak)
library(tidyverse)
library(TopPICR) # remotes::install_github("evanamartin/TopPICR", ref="0.0.3")
library(memoise)
library(MSnbase)
library(MSnSet.utils)
library(PNNL.DMS.utils)

############

####Intensity based check 

#############


si <- read_excel("./source_data/pr3c00353_si_002.xlsx", sheet = 2)

(load("./output_data/shotgun_topdown_int_20240730_modann_cnt.RData"))
mi <- m

(load("./output_data/shotgun_topdown_int_20240730_beforecorrections_mid1a.RData"))

m0 <- m
#for raw intensity needs to be division not subtraction 
#for log transformed needs to be substraction

#only need to do this once, saved rollup as .Rdata below 
if(!file.exists("./output_data/rollup_forpatrie.RData")){
   m <- rrollup(m, "proteoform_id", rollFun = "/", verbose = FALSE)
   save(m, file= "./output_data/rollup_forpatrie.RData")
}else{
   load("./output_data/rollup_forpatrie.RData")
}


# let's filter by count to robustly observed species
mi <- mi[fData(mi)$count > 20,]
m <- m[featureNames(m) %in% featureNames(mi),]






selected_features <- m %>%
   exprs() %>%
   as.data.frame() %>%
   rownames_to_column("feature_name") %>%
   pivot_longer(cols = -feature_name, names_to = "sample_name", values_to = "intensity") %>%
   inner_join(pData(m)) %>%
   select(feature_name, sample_name, intensity, batch) %>%
   filter(intensity > 0) %>%
   group_by(feature_name, batch) %>%
   tally() %>%
   filter(n >= 2) %>%
   group_by(feature_name) %>%
   tally() %>%
   filter(n >= 3) %>%
   pull(feature_name) 

m <- m[selected_features,]
m$batch <- as.factor(m$batch)

#get median intensity for each pfr and add back to fData 
medianint <- apply(exprs(m), 1, median, na.rm=T) %>%
   data.frame() %>%
   rownames_to_column(var="proteoform_id") 

extract_mods <- function(pform){
   mods <- str_extract_all(pform, "(?<=\\[)[^\\]\\[]*(?=\\])", simplify = TRUE)
   posi <- str_locate_all(pform, "(?<=\\[)[^\\]\\[]*(?=\\])")[[1]][,1]
   posi_adj <- (seq_along(mods) * 4 + 2) + c(0, nchar(mods)[-length(mods)])
   posi <- posi - posi_adj
   mod_str <- paste(map2_chr(mods, posi, paste, sep="@"), collapse=", ")
   return(mod_str)
}
extract_mods <- Vectorize(extract_mods)




# what is the meaning of this pipe? and why count (former spectralCount) filter?
# retaining only unmodified APP species plus some housekeeping recalculations
appdf <- fData(m0) %>%
   filter(Gene == "APP") %>%
   # filter(count > 20) %>%
   mutate(firstAA = firstAA - 671,
          lastAA = lastAA - 671) %>%
   mutate(mod_str = extract_mods(Proteoform)) %>%
   as_tibble() %>%
   #dplyr::select(-mods) %>%
   filter(mod_str == "") %>%
   mutate(aa = paste(firstAA, lastAA)) %>%
   select(proteoform_id, aa, firstAA, lastAA) %>%
   distinct() %>%
   left_join(medianint, by ="proteoform_id") %>%
   dplyr::select("firstAA" , "lastAA", "aa", ".") 



##for loop that makes correlation data for our TD data versus Steve's mouse data, each loop is using a different time point in Steve's data 

plottindf <- c()

for(x in unique(si$`Treatment 2`)){
   
   filter <- x
summary <- si %>%
   filter(`Treatment 2` == filter) %>%
   group_by(PFR,Abetaproteoform, Treatment) %>%
   #summarize(sum= sum(Raw), sd= sd(Raw)) %>%
   summarize(sum= median(Raw), sd= sd(Raw)) %>%
   #summarize(sum= mean(Raw), sd= sd(Raw)) %>%
   separate(Abetaproteoform, c("firstAA", "lastAA"), sep="-") %>%
   mutate(lastAA = as.numeric(str_extract(lastAA, "[[:digit:]]+"))) %>%
   mutate(firstAA = as.numeric(str_extract(firstAA, "[[:digit:]]+")))%>%
   mutate(aa = paste(firstAA, lastAA))

merger <- merge(summary, appdf, by ="aa") %>%
   mutate(sum = log10(sum),
          sd = log10(sd),
          logmed = log10(.))

insolcor <- cor(((merger %>% filter(Treatment == "Insol"))$sum), 
                ((merger %>% filter(Treatment == "Insol"))$logmed),
                use = "pairwise.complete.obs")


solcor <- cor(((merger %>% filter(Treatment == "Sol"))$sum), 
              ((merger %>% filter(Treatment == "Sol"))$logmed),
              use = "pairwise.complete.obs")

temp <- merger %>%
   mutate(insolcor =insolcor,
          solcor= solcor,
          timelabel = paste("Mouse age: ", as.numeric(str_extract(x, "[[:digit:]]+")), " months", sep=""),
          time = as.numeric(str_extract(x, "[[:digit:]]+"))
)

plottindf <- rbind(plottindf, temp)

}

library(ggpubr)
library(ggrepel)

plottindf <- plottindf %>%
   mutate(Fraction = case_when(Treatment == "Insol" ~ "Insoluble", Treatment == "Sol" ~ "Soluble")) 

ggscatter(plottindf, x="sum", 
          y="logmed", 
          add="reg.line", 
          conf.int=T, 
          fill="Fraction", 
          cor.method ="spearman",
          color="Fraction")+
   labs(x="log10(Median intensity) - Mouse", y="log10(Median intensity) - Human")+
   #ggtitle(paste(filter, "; insoluble cor =", round(insolcor, 2), "; soluble cor =", round(solcor, 2)))+ 
   #theme(plot.title = element_text(size = 20, face = "bold"))+
   #geom_text_repel(aes(label = aa))+
   geom_text(size=7, x = 10, y = 8.5, aes(label = paste(format(round(insolcor, 2), nsmall = 2))), 
             data = plottindf, 
             color="#F8766D",
             check_overlap = TRUE)+
   geom_text(size=7, x = 4.5, y = 8.5, aes(label = paste(format(round(solcor, 2), nsmall = 2))), 
             data = plottindf, 
             color="#00BFC4",
             check_overlap = TRUE)+
   xlim(4,10.5)+
   ylim(6, 9)+
   theme_classic(base_size = 18)+
   theme(legend.position="bottom")+
   #facet_wrap(~timelabel, ncol = 1)+
   facet_wrap(~factor(timelabel, levels=c('Mouse age: 12 months', 
                                          'Mouse age: 8 months', 
                                          'Mouse age: 5 months', 
                                          'Mouse age: 2 months')), ncol=1)
   # guides(color = guide_legend(title = "Fractions *Kandi et al."),
   #        fill = guide_legend(title = "Fractions *Kandi et al.")) 
   # 
   
ggsave(plot = last_plot(), 
       path = paste(getwd(), "/figures_tables", sep=""),
       filename = "FigS2_intensity_correlation_solvinsol.png", 
       width = 6,
       height = 13,
       units = c("in"))


   

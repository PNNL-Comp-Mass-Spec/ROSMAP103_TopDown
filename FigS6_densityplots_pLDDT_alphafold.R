library(TopPICR)
library(tidyverse)
library(MSnSet.utils)
library(protti)
library(patchwork)
library(drawProteins)
library(ggplot2)
library(ggnewscale)

load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")

#load or retrieve alphafold predications for detected accessions 
if("./source_data/pLDDT_ref_human.RData" %in% list.files()){
   load("./source_data/pLDDT_ref_human.RData")
}else{
   # ### Retrive UNIPROT accessions
   UNIPROT <- unique(fData(m)$UniProtAcc)

   ### Fetch pLDDT scores for each UNIPROT accession
   pLDDT_ref_human <- fetch_alphafold_prediction(
      uniprot_ids = UNIPROT,
      #  organism_name = "Homo sapiens",
      version = "v4",
      timeout = 3600,
      return_data_frame = T,
      show_progress = TRUE
   ) %>%
      distinct(uniprot_id, prediction_score, label_seq_id,auth_comp_id) %>%
      mutate(disc_pLDDT = case_when(
         prediction_score < 50 ~ "very low",
         prediction_score >= 50 & prediction_score < 70 ~ "low",
         prediction_score >= 70 & prediction_score < 90 ~ "high",
         prediction_score >= 90 ~ "very high"
      )) %>%
      mutate(disc_pLDDT = factor(disc_pLDDT, levels = c("very low", "low", "high", "very high")))

   save(pLDDT_ref_human, file="./source_data/pLDDT_ref_human.RData")

}

#####################
#mapping cleavage sites to local prediction scores stored in pLDDT_ref_human.RData

# finds all proteforms that are cleaved at least 3-4 aa away from termini
# so they genuinely be considered cleavages
subset <- fData(m) %>%
   filter(firstAA > 3) %>%
   filter(lastAA < protLength-4)%>%
   distinct(proteoform_id, .keep_all = T)

#range over which pLDDT is averaged 
AArange <- 3

summary <- c()

for(x in 1:nrow(subset)){
   
   ref <- filter(pLDDT_ref_human, uniprot_id == subset[x,]$UniProtAcc)
   
   subset[x,]$firstAA
   
   subset[x,]$lastAA
   
   Navg <- ref %>%
      filter(between(label_seq_id, subset[x,]$firstAA - AArange, subset[x,]$firstAA + AArange-1)) %>%
      summarise(mean = mean(prediction_score), sd=sd(prediction_score))
   #-1 is to make it symmetric, because cut is Nterm to whatever that residue is 
   
   Cavg <- ref %>%
      filter(between(label_seq_id, subset[x,]$lastAA - AArange+1, subset[x,]$lastAA + AArange))%>%
      summarise(mean = mean(prediction_score), sd=sd(prediction_score))
   
   temp <- data.frame(UniProtAcc = subset[x,]$UniProtAcc,
                      Nmean =  Navg$mean,
                      Nsd = Navg$sd,
                      Nloc = 100*subset[x,]$firstAA/subset[x,]$protLength, 
                      Cmean =  Cavg$mean,
                      Csd = Cavg$sd,
                      Cloc =100*subset[x,]$lastAA/subset[x,]$protLength,
                      proteoform_id = subset[x,]$proteoform_id)
   
   summary <- rbind(summary, temp)
   
}  

#making random distribution

Ncuts <- subset %>%
   group_by(UniProtAcc) %>%
   distinct(firstAA) %>%
   tally()

Ccuts <- subset %>%
   group_by(UniProtAcc) %>%
   distinct(lastAA) %>%
   tally()

merge(Ncuts, Ccuts, by ="UniProtAcc") %>%
   mutate(total = n.x+ n.y)


##############
#plotting 

Nn <- summary %>%
   dplyr::select(Nmean, Nsd) %>%
   dplyr::rename(mean = Nmean, sd = Nsd) %>%
   #mutate(data = "N-terminal cleavages")
   mutate(data = "Cleavage")


Cc <- summary %>%
   dplyr::select(Cmean, Csd) %>%
   dplyr::rename(mean = Cmean, sd = Csd) %>%
   #mutate(data = "C-terminal cleavages")
   mutate(data = "Cleavage")

forplotting <- rbind(Nn, Cc)

cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

############plot for manuscript 
library(grid)

forplotting2 <- rbind(forplotting, pLDDT_ref_human %>%
                         filter(uniprot_id %in% unique(fData(m)$UniProtAcc)) %>%
                         mutate(mean = prediction_score, sd = 0, data = "Reference") %>%
                         select(mean, sd, data))

##########
##Kolmogorov–Smirnov test 

#http://www.matf.bg.ac.rs/p/files/69-[Michael_J_Panik]_Advanced_Statistics_from_an_Elem(b-ok.xyz)-776-779.pdf
#You can compare the statistic D to critical values of the D distribution, which appear in tables. 
#If the statistic is greater than the critical value, you reject the null hypothesis that the sample came from the reference distribution.
#if D > Dtable then sample != reference 

kst <- ks.test(
   (forplotting2 %>% filter(data == "Reference"))$mean, c(summary$Nmean, summary$Cmean))

densityplot <- ggplot(forplotting2 %>% filter(data != "Random sampling"), aes(mean, fill = data)) + 
   annotate(geom = "rect",
            xmin = 0, xmax = 50,
            ymin = -Inf,   ymax = Inf,
            fill = "dodgerblue4", alpha = 0.05) +
   annotate(geom = "rect",
            xmin = 50, xmax = 70,
            ymin = -Inf,   ymax = Inf,
            fill = "dodgerblue4", alpha = 0.1) +
   annotate(geom = "rect",
            xmin = 70, xmax = 90,
            ymin = -Inf,   ymax = Inf,
            fill = "dodgerblue4", alpha = 0.2) +
   annotate(geom = "rect",
            xmin = 90, xmax = 100,
            ymin = -Inf,   ymax = Inf,
            fill = "dodgerblue4", alpha = 0.3) +
   geom_density(alpha = 0.5)+
   theme_classic(base_size = 24)+
   theme(legend.position="bottom",  plot.margin = margin(t = 20,  # Top margin
                                                         r = 50,  # Right margin
                                                         b = 40,  # Bottom margin
                                                         l = 10))+ # Left margin)+
   labs(x="pLDDT", y="Density", fill="")+
   scale_x_continuous(expand=c(0, 0), limits=c(0, 100)) +
   scale_y_continuous(expand=c(0, 0))+
   scale_fill_manual(values = c("#E69F00", "#999999"))+
   #geom_vline(xintercept = 50, colour="black", linetype = "longdash")+
   annotate("text", x=65, y=0.040, label= paste(kst$method, "\n", "D =", round(kst$statistic, digits=4), "\n", "p-value < 2.2e-16"), size = 6)+
   annotate("text", x=95, y=0.058, label= paste("Very high \n (pLDDT > 90)"), color = "dodgerblue4")+
   annotate("text", x=80, y=0.058, label= paste("High \n (90 > pLDDT > 70)"), color = "dodgerblue3")+
   annotate("text", x=60, y=0.058, label= paste("Low \n (70 > pLDDT > 50)"), color = "dodgerblue2")+
   annotate("text", x=25, y=0.058, label= paste("Very low \n (pLDDT < 50)"), color = "dodgerblue")
densityplot  

ggsave(densityplot, filename = "./figures_tables/FigS6_densityplots_pLDDT_alphafold.png",
       width = 14, height = 6,
       bg = "white")



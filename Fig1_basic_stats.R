library(tidyverse)
library(ggbreak)
library(extrafont)
library(TopPICR) # remotes::install_github("evanamartin/TopPICR", ref="0.0.3")
library(memoise)
library(MSnbase)
library(MSnSet.utils)
library(PNNL.DMS.utils)
library(plotly)
library(RColorBrewer)
library(ggsci)
library(grid)
library(ggpubr)
library(MSnSet.utils)
library(tidyverse)
library(extrafont)
loadfonts()


# just in case
if(!dir.exists("figures_tables"))
   dir.create("figures_tables")



load("./output_data/shotgun_topdown_int_20240730_modann.RData")


x <- exprs(m) %>%
   as.data.frame()

meta <- m@featureData@data

meta_short <- meta %>%
   dplyr::select(Gene, proteoform_id, UniProtAcc, mass)

x_long <- x %>%
   rownames_to_column(var = "proteoform_id") %>%
   remove_rownames() %>%
   pivot_longer(-proteoform_id, names_to = "SubjectID",
                values_to = "Intensity") %>%
   inner_join(.,meta_short)

#these are averages 
pf_count <- x_long %>%
   filter(!is.na(Intensity)) %>%
   group_by(SubjectID) %>%
   tally() %>%
   summarize(mean = mean(n),
             SD = sd(n)) %>%
   mutate(Type = "Proteoforms")

#these are averages
gene_count <- x_long %>%
   filter(!is.na(Intensity)) %>%
   distinct(Gene, SubjectID) %>%
   group_by(SubjectID) %>%
   tally() %>%
   summarize(mean = mean(n),
             SD = sd(n)) %>%
   mutate(Type = "Genes")

combo <- full_join(pf_count, gene_count)

panela <- combo %>%
   ggplot()+
   aes(x = Type, y = mean, fill = Type)+
   geom_bar(stat = "identity")+
   theme_bw(base_size = 22) +
   theme(legend.position = "none",
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 22),
         axis.text.x=element_text(angle = 0,
                                  vjust = 0.5, 
                                  hjust= 0.5, 
                                  color='black',
                                  size = 22),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         text=element_text(family="Helvetica"),
         axis.line = element_line())+
   xlab("")+
   scale_fill_manual(values = c("#88419d", "#8c96c6"))+
   geom_errorbar( aes(x=Type,
                      ymin=mean-SD, ymax=mean+SD),
                  width=0.4, colour="black", alpha=1, size=1.2)+
   scale_y_break(c(1000,5000), space = 0.6)+
   ylab("Mean Observations (n)")

unique_genes <- x_long %>%
   filter(!is.na(Intensity)) %>%
   distinct(Gene) %>% 
   tally() 

unique_pfr <- x_long %>%
   filter(!is.na(Intensity)) %>%
   distinct(proteoform_id) %>% 
   tally() 

unique_acc <- x_long %>%
   filter(!is.na(Intensity)) %>%
   distinct(UniProtAcc) %>% 
   tally() 

panelb <- combo %>%
   mutate(mean = c(as.numeric(unique_pfr),as.numeric(unique_genes))) %>%
   ggplot()+
   aes(x = Type, y = mean, fill = Type)+
   geom_bar(stat = "identity")+
   theme_bw(base_size = 22) +
   theme(legend.position = "none",
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 22),
         axis.text.x=element_text(angle = 0,
                                  vjust = 0.5, 
                                  hjust= 0.5, 
                                  color='black',
                                  size = 22),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         text=element_text(family="Helvetica"),
         axis.line = element_line())+
   xlab("")+
   scale_fill_manual(values = c("#88419d", "#8c96c6"))+
   scale_y_break(c(2000,8000),  space = 0.6)+
   ylab("Total Observations (n)")

numberofsamples <- 103 

panelc <- x_long %>%
   filter(!is.na(Intensity)) %>%
   group_by(proteoform_id) %>%
   tally() %>% 
   mutate(n = n/numberofsamples) %>% 
   mutate(Group = case_when(between(n, 0, 0.25) ~ "25-0%",
                            between(n, 0.25, 0.5) ~ "50-25%",
                            between(n, 0.5, 0.75) ~ "75-50%",
                            between(n, 0.75, 1) ~ "100-75%")) %>%
   group_by(Group) %>%
   tally() %>%
   ggplot()+
   aes(x= fct_relevel(Group,"100-75%",
                      "75-50%", "50-25%", "25-0%" ),
       y = n, fill = fct_relevel(Group,"100-75%",
                                 "75-50%", "50-25%", "25-0%"))+
   geom_bar(stat= "identity")+
   theme_bw(base_size = 22) +
   theme(legend.position = "none",
         panel.background = element_rect(fill= 'white'),
         axis.text.y=element_text(color = 'black', size = 22),
         axis.text.x=element_text(angle = 0,
                                  vjust = 0.5, 
                                  hjust= 0.5, 
                                  color='black',
                                  size = 22),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         text=element_text(family="Helvetica"),
         axis.line = element_line())+
   scale_fill_brewer(palette = "Blues", direction = -1)+
   xlab("Data Completeness")+
   ylab("Total (n)")



x <- fData(m)

#' # getting the stats for all accessions
mod_stats <- map(unique(x$UniProtAcc), ~ get_mods_counts(x,.x))
mod_vec <- map(mod_stats, names) %>% unlist() %>% unique()
mod_vec <- vector("numeric", length(mod_vec)) %>% setNames(mod_vec)
for(i in mod_stats){
   mod_vec[names(i)] <- mod_vec[names(i)] + i
}
mod_vec <- sort(mod_vec, decreasing = TRUE)

#' barplot of top 50
mod_vec_df <- data.frame(mod_name = names(mod_vec), mod_count = mod_vec) %>%
   mutate(mod_name = ordered(mod_name, levels=rev(unique(mod_name))))

ggplot(head(mod_vec_df, 50)) +
   geom_bar(aes(y=mod_name, x=mod_count), stat="identity", orientation = "y") +
   xlab("count") +
   ylab("modification")


#' # removing some artifact mods and redoing the plot and N-term acetylation
x2 <- x

x2 <- remove_a_mod(x2, "+1C13")
x2 <- remove_a_mod(x2, "-1C13")
x2 <- remove_a_mod(x2, "IronAdduct")
x2 <- remove_a_mod(x2, "Cation:K")
x2 <- remove_a_mod(x2, "Carbamyl")
x2 <- remove_a_mod(x2, "N-Acetyl")


mod_stats <- map(unique(x2$UniProtAcc), ~ get_mods_counts(x2,.x))
mod_vec <- map(mod_stats, names) %>% unlist() %>% unique()
mod_vec <- vector("numeric", length(mod_vec)) %>% setNames(mod_vec)
for(i in mod_stats){
   mod_vec[names(i)] <- mod_vec[names(i)] + i
}
mod_vec <- sort(mod_vec, decreasing = TRUE)

#' barplot of top 50
mod_vec_df <- data.frame(mod_name = names(mod_vec), mod_count = mod_vec) %>%
   mutate(mod_name = case_when(mod_name == "Plus1Oxy" ~ "Oxidation",
                               #mod_name == "N-Acetyl" ~ "N-Acetylation", # removed
                               mod_name == "Acetyl" ~ "Acetylation",
                               mod_name == "Phosph" ~ "Phosphorylation",
                               mod_name == "Deamide" ~ "Deamidation",
                               mod_name == "Dimethyl" ~ "Dimethylation",
                               mod_name == "Plus2Oxy" ~ "Dioxidation",
                               mod_name == "Glu->pyro-Glu" ~ "pyro-Glu",
                               mod_name == "Methyl" ~ "Methylation",
                               TRUE ~ mod_name))%>%
   mutate(mod_count2 = mod_count/11507 *100) %>%
   mutate(mod_name = ordered(mod_name, levels=rev(unique(mod_name)))) 

paneld <- ggplot(head(mod_vec_df, 10)) +
   geom_bar(aes(y=mod_name, x=mod_count2, fill = ""), stat="identity", orientation = "y",
            show.legend = FALSE, alpha = 0.5) +
   xlab("Percentage of Proteoforms (%)") +
   ylab("Modifications")+
   scale_fill_manual(values = "#238b45") +
   #scale_x_continuous(breaks=seq(0,25, by =5))+
   xlim(0,5)+
   theme_minimal(base_size = 20)+
   theme(
      #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      text=element_text(family="Helvetica"))

library(patchwork)
#need to manually add tags A-D, for some reason sizing in non-linear 
panela + panelb + panelc + paneld
plot_layout(widths = c(0.35, 0.4, 1),
            guides = "collect")

ggsave(filename = "./figures_tables/Fig1_basic_stats.png", width = 16, height = 14)



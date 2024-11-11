
library(RColorBrewer)
library(tidyverse)
library(ggsci)
library(MSnSet.utils)
library(grid)
library(TopPICR)
library(tidyr)
library(dplyr)


load("./output_data/shotgun_topdown_sc_20240730_modann_wres.RData")
m_spc <- m
load("./output_data/spectral_counting_test_results.RData")
res_spc <- results

app_c <- grep("^APP_", featureNames(m_spc), value=TRUE) 


# ID, first, last, mod, counts, sig_sqrt_amyloid, caa_4gp, sig_sqrt_amyloid, caa_4gp
##USES SC data and not intensity 
ab <- fData(m_spc[app_c,]) %>%
  select(-protLength, -mass, -rt, -Gene, -pcGroup,
         -proteoform_id, -UniProtAcc, -collision)

extract_mods <- function(pform){
  mods <- str_extract_all(pform, "(?<=\\[)[^\\]\\[]*(?=\\])", simplify = TRUE)
  posi <- str_locate_all(pform, "(?<=\\[)[^\\]\\[]*(?=\\])")[[1]][,1]
  posi_adj <- (seq_along(mods) * 4 + 2) + c(0, nchar(mods)[-length(mods)])
  posi <- posi - posi_adj
  mod_str <- paste(map2_chr(mods, posi, paste, sep="@"), collapse=", ")
  return(mod_str)
}
extract_mods <- Vectorize(extract_mods)


ab <- ab %>%
  mutate(mod_str = extract_mods(Proteoform)) %>%
  mutate(firstAA = firstAA - 671,
         lastAA = lastAA - 671) %>%
  # arrange(-count) %>%

  # arrange(firstAA) %>%
  as_tibble() %>%
  dplyr::select(-mods)

#####

#kind of a table of proteoforms, need this for SI

######
#write.csv(ab, "manual_adjustment_unedited_vladversion.csv")

#before running next lines, you need to open and manually edit this .csv to change mod names to something more palatable, eg -18@## = "WaterLoss@##"
#Notes
#Carbamyl is actually Acetyl should be at K16
#16.9935@24 for Ab 0-42 is N-term oxidation?
#all pyro glus must be N-terminal!!!!!Make sure residue with pyroglu = firstAA
#-17 is probably a gln -> pyroglu@15, looks like ammonia loss 
#-18 near Ntermini is porbably glu-> pyroglu@3, looks likem water loss
#111.0376@4 for AB4-43 is actually  AB3-43 with pyroglu@3 
#1585 is 1-42 ACetyl@K16, it thinks its a truncation with huge mod 
#-89@13 for 2-42 is actually 3-42 with pyroglu, pyroglu+alanine = 89 DA 
# 
#  ab <-  read.csv("C:/Users/ives435/OneDrive - PNNL/Desktop/Alz TD PTMS/Evan data/data/manual_adjustment_APP_JMFcopy.csv") %>%
#     arrange(lastAA, firstAA)



# let's block this for now
################################################################################
# ab <-  read.csv("manual_adjustment_edited_vladversion.csv") %>%
#    arrange(lastAA, firstAA)
################################################################################


# write_csv(ab, file = "current_abetas.csv")

# scripting the manual interpretation of the Abeta species
ab <- as.data.frame(ab)
rownames(ab) <- ab$feature_name

ab["APP_106","mod_str"] <- "Acetyl@16"
ab["APP_18","mod_str"] <- "WaterLoss@1-6"
ab["APP_29","mod_str"] <- "PyroGlu@3"
ab["APP_37","mod_str"] <- "Oxidation@35"
ab["APP_45","mod_str"] <- "Oxidation@1'"
ab["APP_46","mod_str"] <- "PyroGlu@3;Deamidation@27"
ab["APP_50","mod_str"] <- "Oxidation@35"
ab["APP_51","mod_str"] <- "Acetyl@16"
ab["APP_57","mod_str"] <- "WaterLoss@40"
ab["APP_61","mod_str"] <- "PyroGlu@3"; ab["APP_61","firstAA"] <- 3
ab["APP_67","mod_str"] <- "PyroGlu@3"; ab["APP_67","lastAA"] <- 43
ab["APP_7","mod_str"] <- "PyroGlu@3"; ab["APP_7","firstAA"] <- 3




ab2 <- res_spc %>%
  dplyr::rename(feature_name = proteoform) %>%
  semi_join(ab) %>%
  filter(variable %in% c("sqrt_amyloid", "caa_4gp", "anye4", "cogng_demog_slope", "cogn_global_lv")) %>%
  select(feature_name, adj.P.Val, variable) %>%
  pivot_wider(names_from = variable, values_from = adj.P.Val) %>%
  left_join(ab,.)

ab2 <- ab2 %>%
  select(-Proteoform) %>%
  filter(count >= 21) %>%
   #watch out for this step!!! turning data remakes some pfrs get lost if counts change and drop below 30 (old threshold)
  mutate(feature_name = ordered(feature_name, levels=rev(feature_name))) %>%
  mutate(is_modified = mod_str != "") %>%
  mutate(mod_str2 = case_when(mod_str == "" ~ "None", TRUE ~ mod_str)) %>%
  mutate(modification = ordered(.$mod_str2, levels = unique(.$mod_str2)))

# ordering proteoforms according to N- and C- termini
ab2 <- ab2 %>%
  arrange(desc(lastAA), desc(firstAA)) %>%
  mutate(feature_name = ordered(feature_name, levels= unique(feature_name)))

#checking if p values changed a lot after 4 kda search for cogn_global of trunco forms 

temp <- ab2 %>%
  pivot_longer(cols = c(caa_4gp, sqrt_amyloid, anye4, cogng_demog_slope, cogn_global_lv),
               names_to = "variable", values_to = "adj.P.Val") %>%
  mutate(adj.P.Val = case_when(adj.P.Val > 0.05 ~ NA_real_,
                               TRUE ~ adj.P.Val)) %>%
   #this is filtering for P.val > 0.05
  mutate(specie_lbl = paste(firstAA, lastAA))

ABlabels <- temp %>%
  distinct(feature_name ,.keep_all = TRUE) %>%
  # mutate(feature_name2 = paste0("A\U03B2","(",firstAA,"-", lastAA,")", if_else(is_modified  == TRUE,"*",""))) %>% #adds afsterik to mod forms
   mutate(feature_name2 = paste0("A\U03B2","(",firstAA,"-", lastAA,")", if_else(is_modified  == TRUE,"",""))) %>%
  mutate(feature_name3 = fct_rev(feature_name)) %>%
  select(feature_name, feature_name2) %>%
  arrange((feature_name))

temp$count_variable <- "Abundance"

#spectral counts sideways 

p0 <- ggplot(temp) +
  aes(y=feature_name, x=count) +
   #before position dodge it was stacking counts multiple times because of different tests, so values looked way higher 
  geom_bar(stat="identity", color = NA, fill = "black", position = "dodge") +
#scale_x_continuous(limits = c(0, 8e3), )+
  scale_x_reverse(labels = scales::scientific, breaks = c(0, 2e3, 4e3, 6e3, 8e3, 1e4)) +
  scale_fill_gradientn(#colors = MSnSet.utils::hot2.colors(50),
    colors = terrain.colors(50),
    name = expression(paste("log"[10],"(count)"))) +
  guides(size = guide_legend(title = expression(paste("-log"[10],"(adj p-value)")))) +
  facet_grid(~ count_variable, scales = "free_x", drop = T) +
  theme_bw(base_size=12) +
  theme(axis.text.y.right = element_text(hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "left",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        # panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background =element_rect(fill="white"),
        # strip.text = element_text(face="bold"),
        strip.text = element_text(face="bold", size = 12))+
  scale_y_discrete(name = NULL, expand = expansion(0.02),
                   labels= ABlabels$feature_name)
                   #labels= ABlabels$feature_name2)
p0

#significance plot, dots on grid
p1 <- temp %>%
  mutate(facet_var = "Significance") %>%
  ggplot() +
  aes(y=feature_name, x=variable, size=-log10(adj.P.Val),
      fill = -log10(adj.P.Val)) +
  geom_point(shape=21, stroke = 0.5#, fill="#B0BEC5"
             ) +
   scale_fill_viridis_c(guide = "legend", breaks = c(4,8,12,16), direction =-1) +
   scale_size_continuous(breaks = c(4,8,12,16))+
  # scale_fill_gradientn(colors = viridis::viridis(4, direction = -1), breaks = c(4,8,12,16)) +
  guides(size = guide_legend(title = expression(paste("-log"[10],"(adj p-value)")),
                             title.position = "top"),
  fill = guide_legend(title = expression(paste("-log"[10],"(adj p-value)")),
                       title.position = "top")) +
  facet_grid(cols = vars(facet_var)) +
  # facet_grid(~ variable, scales = "free_x", drop = T) +
  theme_bw(base_size=12) +
  theme(axis.text.y.right = element_text(hjust = 1),
        axis.title.x = element_blank(),
        #
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # panel.grid.major.x = element_blank(),
        #
        # legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        strip.background =element_rect(fill="white"),
        # strip.text = element_text(face="bold"),
        strip.text = element_text(face="bold", size = 12)) +
   
  scale_y_discrete(name = NULL, expand = expansion(0.02))+
   scale_x_discrete(labels=c(#"anye4" = expression(italic("APOE")~e4), 
                             "anye4" = expression("" *italic("APOE")* "\u03B54"),
                             "caa_4gp" = "cerebral amyloid angiopathy",
                             "cogn_global_lv" = "global cognition", 
                             "cogng_demog_slope" = "slope of cognitive decline", 
                             "sqrt_amyloid" = "amyloid"))
p1

# plot Abeta fragments

temp <- ab2 %>%
  select(feature_name, firstAA, lastAA, is_modified, modification) %>%
  mutate(firstAA = firstAA,
         lastAA = lastAA) %>%
  pivot_longer(cols = c(firstAA, lastAA), names_to = "terminus", values_to = "deviation")

temp <- temp %>%
  mutate(refAA = case_when(terminus == "firstAA" ~ 1, TRUE ~ 40))

temp2 <- temp %>%
  select(-refAA) %>%
  pivot_wider(id_cols = c(feature_name, is_modified, modification),
              names_from = "terminus", values_from = deviation) %>%
  mutate(ymin = as.numeric(feature_name) - 0.2,
         ymax = as.numeric(feature_name) + 0.2)

temp2$modification <- ordered(temp2$modification, levels = c("None", 
                                                             "Acetyl@16", 
                                                             "IronAdduct@23",
                                                             "PyroGlu@3", 
                                                             "PyroGlu@3;Deamidation@27", 
                                                             "Oxidation@1'", 
                                                             "Oxidation@35", 
                                                             "WaterLoss@1-6",
                                                             "WaterLoss@40"))

mod_cols <- c(
   "#BDBDBD",
   "#FDD835",
   "#F4511E",
   "#8E24AA",
   "#5E35B1",
   "#00897B",
   "#7CB342",
   "#00ACC1",
   "#1E88E5")

#C-term + body of abeta frags 
p3 <- temp2 %>%
  mutate(facet_var = "Proteoform") %>%
  ggplot() +
  aes(ymin = ymin, ymax = ymax,
      xmin = firstAA, xmax = lastAA,
      fill=modification) +
  geom_rect(stat = "identity", color="black", size = 0.1) +
  facet_grid(cols = vars(facet_var)) +
  geom_vline(xintercept = c(1, 40)) +
  scale_fill_manual(breaks = levels(temp2$modification),
                    values = mod_cols,
                    name = "Modification") +
  scale_y_continuous(breaks = seq_along(temp2$feature_name),
                     labels = ABlabels$feature_name2,
                     expand = expansion(0.015), position = "right") +
  # ggbreak::scale_x_break(breaks = c(12,36),
  #                        expand = expansion(0),
  #                        space = 0.5) +
  scale_x_continuous(breaks = c(0, seq(4,44,4)), labels = c("1'", as.character(seq(4,44,4)))) +
  # annotation_custom(grob = textGrob(label = "Residue",
  #                                   gp = gpar(fontsize = 11),
  #                                   y = -20, default.units = "pt")) +
   annotation_custom(grob = textGrob(label = "Amino Acid Residue",
                                     gp = gpar(fontsize = 10),
                                     y = -30, default.units = "pt")) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "lightgrey",
                                          size = 0.4,
                                          linetype = 1),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size = 12),
        # strip.text = element_text(size = 14),
        legend.position = "right",
        axis.title.y = element_blank())
p3


library(patchwork)
p0 + p1 + p3 +
  plot_layout(widths = c(0.35, 0.25, 1),
              guides = "collect")

ggsave(path = "figures_tables", filename = "Fig4_abeta_species.png" , width = 10, height = 6, scale = 1.2)





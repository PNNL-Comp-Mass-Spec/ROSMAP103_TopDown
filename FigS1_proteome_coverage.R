
library(tidyverse)
library(Biostrings)
library(MSnSet.utils)

#FASTA @ Documents/R/Large Sample Practice/Evan Files/_Brain_TD_103_manuscript/R Data
#need to change this from local 
# x <- Biostrings::readAAStringSet(
#    file.path(r"(\\gigasax\DMS_FASTA_File_Archive\Dynamic\Forward\)",
#              "ID_008032_8627C6BD.fasta"))

x <- Biostrings::readAAStringSet("./source_data/ID_008032_8627C6BD.fasta")



sp <- data.frame(protLen = width(x), protNames = names(x))

sp <- sp %>%
   filter(!grepl("^DECOY", protNames)) %>%
   filter(grepl("_HUMAN", protNames)) %>%
   mutate(protNames = sub("(^[^ ]*) .*","\\1",protNames)) %>%
   mutate(annType = case_when(grepl("-", protNames) ~ "isoforms",
                              !grepl("-", protNames) & grepl("^sp", protNames) ~ "Swiss-Prot Database",
                              grepl("^tr", protNames) ~ "tentative",
                              TRUE ~ "stuff"))

load("./output_data/shotgun_topdown_sc_20240730_modann.RData")

covered_proteins <- fData(m) %>%
   distinct(UniProtAcc, protLength) %>%
   dplyr::rename(protLen = protLength) %>%
   mutate(annType = "Proteins \n (This Study)")

all_fragments <- fData(m) %>%
   mutate(fragLength = lastAA - firstAA + 1) %>%
   select(UniProtAcc, fragLength) %>%
   dplyr::rename(protLen = fragLength) %>%
   mutate(annType = "Proteoforms \n (This Study)")

max_fragments <- all_fragments %>%
   group_by(UniProtAcc) %>%
   summarise(protLen = max(protLen)) %>%
   mutate(annType = "Proteoforms \n (Maximum)")


# merging it all together
sp <- select(sp, annType, protLen) %>% as_tibble()
covered_proteins <- select(covered_proteins, annType, protLen) %>% as_tibble()
all_fragments <- select(all_fragments, annType, protLen) %>% as_tibble()
max_fragments <- select(max_fragments, annType, protLen) %>% as_tibble()

x <- bind_rows(sp, covered_proteins, all_fragments, max_fragments) %>%
   mutate(kDa = protLen * 0.11) %>%
   mutate(annType = ordered(annType, levels = c("Swiss-Prot Database",
                                                "isoforms","tentative",
                                                "Proteins \n (This Study)",
                                                "max fragments",
                                                "Proteoforms \n (This Study)")))

xx <- x %>%
   filter(annType == "Swiss-Prot Database" | 
             annType == "Proteins \n (This Study)" | 
             annType == "Proteoforms \n (This Study)") %>%
   mutate(annType = ordered(annType, levels = c("Proteoforms \n (This Study)", "Proteins \n (This Study)", "Swiss-Prot Database")))

# 110 Da
sumplot <- ggplot(xx) +
   aes(x = kDa) +
   geom_histogram(binwidth = 2.5, show.legend = FALSE, 
                  color = "#212121",
                  fill = "#757575") +
   scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0,0), add = 0)) +
   scale_y_continuous(expand = expansion(mult = c(0,0.1), add = 0)) +
   facet_wrap(~annType, scales = "free_y") +
   theme_bw(base_size = 16) +
   theme(strip.background = element_rect(fill="#FAFAFA"),
         plot.margin = margin(t = 2, r = 16, b = 1, l = 1, unit = "pt"))

ggsave(plot= sumplot, filename = "./figures_tables/FigS1_proteome_coverage.png", width = 11, height = 6, scale = 0.7)






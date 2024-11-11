


library(plotly)
library(tidyverse)
library(ggplot2)

load("./output_data/pre_identifications.RData")

#loads in x_meta and x_grp 

w_cluster <- x_grp

w_prot <- w_cluster %>%
   filter(Gene == "HP") %>%
   mutate

#Makes Figure S4 
plot_ly(data = w_prot,
        x = ~RTalign,
        y = ~RecalMass,
        type = "scatter",
        mode = "markers",
        color = ifelse(w_prot$RTalign >= 3750, "#ef8a62", "#67a9cf"), #colors are green and purple? 
        size = ifelse(w_prot$pcGroup == 0, 0, 1)) %>%
   layout(xaxis = list(range=c(3400,4100), title = "Retention time (seconds)"),
          title="Haptoglobin (HP)", 
          yaxis = list(title = "Mass (Da)",range=c(8500, 9500)))

#remake in ggplot
hpplot <- w_prot %>%
   mutate(range = case_when(RTalign <= 3750 ~ "HP*1F", RTalign > 3750 ~ "HP*1S")) %>%
   ggplot() +
   aes(x=RTalign, y= RecalMass, color=range) +
   geom_point()+
   ylim(8800, 9400)+
   xlim(3400, 4100) +
   labs(x="Retention time (sec)", y="Mass (Da)", title= "Haptoglobin (HP)", color="Allele")+
   theme_bw(base_size = 16)+
   theme(text=element_text(family="Helvetica"),
         panel.border = element_blank())
hpplot

ggsave(plot = hpplot, 
       path = paste(getwd(), "/figures_tables", sep=""),
       filename = "FigS5_haptoglobin.png", 
       width = 8,
       height = 6,
       units = c("in"))


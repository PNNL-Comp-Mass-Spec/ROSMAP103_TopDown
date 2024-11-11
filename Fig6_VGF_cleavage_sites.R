

library(tidyverse)
library(MSnSet.utils)
library(TopPICR)
library(viridis)
library(stringr)



vgf_seq <- "
MKALRLSASALFCLLLINGLGAAPPGRPEAQPPPLSSEHKEPVAGDAVPGPKDGSAPEVR
GARNSEPQDEGELFQGVDPRALAAVLLQALDRPASPPAPSGSQQGPEEEAAEALLTETVR
SQTHSLPAPESPEPAAPPRPQTPENGPEASDPSEELEALASLLQELRDFSPSSAKRQQET
AAAETETRTHTLTRVNLESPGPERVWRASWGEFQARVPERAPLPPPAPSQFQARMPDSGP
LPETHKFGEGVSSPKTHLGEALAPLSKAYQGVAAPFPKARRPESALLGGSEAGERLLQQG
LAQVEAGRRQAEATRQAAAQEERLADLASDLLLQYLLQGGARQRGLGGRGLQEAAEERES
AREEEEAEQERRGGEERVGEEDEEAAEAEAEAEEAERARQNALLFAEEEDGEAGAEDKRS
QEETPGHRRKEAEGTEEGGEEEDDEEMDPQTIDSLIELSTKLHLPADDVVSIIEEVEEKR
KRKKNAPPEPVPPPRAAPAPTHVRSPQPPPPAPAPARDELPDWNEVLPPWDREEDEVYPP
GPYHPFPNYIRPRTLQPPSALRRRHYHHALPPSRHYPGREAQARRAQEEAEAEERRLQEQ
EELENYIEHVLLRRP"
vgf_seq <- gsub("\\n", "", vgf_seq)
basic_sites <- str_locate_all(vgf_seq, "[KR]{2,}")[[1]]



load("./output_data/shotgun_topdown_int_20240730_modann_cnt_wres.RData")

x <- fData(m) %>%
   mutate(ProtLen = protLength) %>%
   mutate(First_AA = firstAA,
          Last_AA = lastAA) %>%
   mutate(feature_name = proteoform_id)


xg <- dplyr::filter(x, Gene == "VGF") %>%
   mutate(SigLvl = case_when(`cogng_demog_slope.apv` <= 0.05 ~ "adjusted",
                            `cogng_demog_slope.pv` <= 0.05 ~ "nominal",
                            TRUE ~ "not significant"))




library(patchwork)

anno <-  data.frame(entry = c("1_22_Signal_1", 
"177_206_NERP-3_1",
"281_306_NERP-1_1.5",
"310_347_NERP-2_1",
"419_427_TPGH_2_1.5",
"485_615_NAPP129_1.5",
"485_503_NERP-4_1",
"554_615_TLQP-62_2",
"554_577_AMP_1",
"586_615_AQEE-30_2.5",
"554_574_TLQP-21_3",
"597_615_LQEQ-19_1")) %>% 
   separate(entry, into=c("firstAA", "lastAA", "entry", "rowname"), sep="_") %>%
   mutate(firstAA = as.numeric(str_extract(firstAA, "[[:digit:]]+"))) %>%
   mutate(lastAA = as.numeric(str_extract(lastAA, "[[:digit:]]+"))) 

vgf_known <-
   tribble(
      ~entry,   ~firstAA, ~lastAA, ~ymin, ~ymax,
      "Signal",              1,           22, -0.22,  -0.26, 
      "NERP-3",            177,          206, -0.22,  -0.26, 
      "NERP-1",            281,          306, -0.22,  -0.26,
      "NERP-2",            310,          347, -0.32,  -0.36, 
      "TPGH",              419,          427, -0.22, -0.26, 
      "NERP-4",            485,          503, -0.32,  -0.36, 
      "NAPP129",           485,          615, -0.22, -0.26, 
      "TLQP-21",           554,          574, -0.52,  -0.56, 
      "AMP",               554,          577, -0.42,  -0.46, 
      "TLQP-62",           554,          615, -0.32,  -0.36, 
      "AQEE-30",           586,          615, -0.42,  -0.46, 
      "LQEQ-19",           597,          615, -0.52,  -0.56) 



# version from Becky Carlyle and Steve Arnold review
vgf_known <-
   tribble(
      ~entry,   ~firstAA, ~lastAA, ~ymin, ~ymax,
      "Signal",              1,           22, -0.22,  -0.26, 
      "NERP-3",            177,          206, -0.22,  -0.26, 
      "NERP-1",            281,          306, -0.22,  -0.26,
      "NERP-2",            310,          347, -0.32,  -0.36, 
      "PGH",              419,          427, -0.22, -0.26, 
      "NERP-4",            485,          507, -0.32,  -0.36, 
      "NAPP-129",           485,          615, -0.22, -0.26, 
      "TLQP-62",           554,          615, -0.32,  -0.36, 

      "NAPP-19",           485,          503, -0.42,  -0.46, 
      "NYPG-41",           575,          615, -0.42,  -0.46,
      
      "AMP",               554,          577, -0.52,  -0.56, 
      "AQEE-30",           586,          615, -0.52,  -0.56, 
      
      "TLQP-21",           554,          574, -0.62,  -0.66, 
      "LQEQ-19",           597,          615, -0.62,  -0.66
      
      ) 




vgf_known$ymin <- vgf_known$ymin + 0.1
vgf_known$ymax <- vgf_known$ymax + 0.1

cleavages <- c(
   "23,s",
   "177,k",
   "206,k",
   "310,k",
   "347,k",
   "419,k",
   "430,k",
   #"481,k",
   "485,k",
   "505,k",
   "554,k",
   "576,k",
   "586,k",
   "597,k",
   "64,u",
   "121,u",
   "373,u",
   "564,u",
   "281,k"
)

cleavagesdf <- data.frame(sites = cleavages) %>%
   separate(sites, into=c("site", "status"), sep=",") %>%
   mutate(site = as.numeric(str_extract(site, "[[:digit:]]+"))) 
   

breaksdf <- data.frame(breaks= c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600),
                       labels = c("-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"))

library(simplecolors)
#remotes::install_github("hrbrmstr/ggchicklet")
library(ggchicklet)

vgfplot <- TopPICR::plot_accession_ptm(xg, "O15240", fill_by = "SigLvl", mods_to_name = "top10") +
   scale_fill_manual(values = alpha(c("#EF5350", "#66BB6A", "#BDBDBD"), 0.8)) +
   scale_shape_manual(name = "Modification:", 
                        labels = c("Water loss", "Iron Adduct", "Met -> Asn", "Gln -> pyro-Glu", "Phospho"),
               values=c(49, 50, 51, 52, 53)) +
   labs(fill="") +
     geom_vline(xintercept = (filter(cleavagesdf, status == "u"))$site, color="red", alpha=0.2, linetype = "longdash", size=0.7)+
     geom_vline(xintercept = (filter(cleavagesdf, status == "k"))$site, color="blue", alpha=0.2, linetype = "longdash", size=0.7)+
     geom_vline(xintercept = (filter(cleavagesdf, status == "s"))$site, color="forestgreen", alpha=0.2, linetype = "longdash", size=0.7)+
     # ylim(-0.50, 0.4)+
     ylim(-0.60, 0.4)+
     xlim(0, 620) +
     scale_x_continuous(breaks= c(1, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600)) +
   #geom_text(aes(x=firstAA, y=ymin, label=ymin))+
   geom_rect(fill = NA, data=vgf_known, inherit.aes = FALSE,
            aes(xmin=firstAA, xmax=lastAA, ymin=(ymin+ymax)/2, ymax=ymax),
            color="black", size=0.5, alpha=0.1)+
        geom_text(data=vgf_known, aes(label=entry, x=(firstAA + lastAA)/2, y=ymax-0.02), size=4) +
   
   geom_rect(fill = "yellow", data = data.frame(basic_sites), inherit.aes = FALSE,
             aes(xmin=start-1, xmax=end+1), ymin= 0, ymax= +0.025,
             color="black", size=0.2) +
   
   geom_text(data=vgf_known, aes(label=firstAA, x=firstAA-0, y=ymin, angle=90), size=4) +
   geom_text(data=vgf_known, aes(label=lastAA, x=lastAA+0, y=ymin , angle=90), size=4) +
     #brute force way to center x axis marks 
     geom_text(data=breaksdf, aes(label=breaks, x=breaks+0.5, y=-0.022 , angle=90), size=4)+
     geom_text(data=breaksdf, aes(label=labels, x=breaks, y=-0.00 , angle=90), size=4)+
     annotate("text", x = 320, y = -0.05, label = "residue")+
       theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank())
  
# vgfplot


ggsave(plot = vgfplot, path = "figures_tables", filename = "Fig6_VGF_cleavage_sites.png", width = 12.5, height = 8)





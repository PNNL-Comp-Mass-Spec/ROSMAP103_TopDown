

library(tidyverse)
library(MSnSet.utils)

# spectral counting data
(load("./output_data/shotgun_topdown_sc_20240730_modann_wres.RData"))


plot(m$sqrt_amyloid, as.numeric(exprs(m)["APP_6",])) # reference 1-42 
plot(m$sqrt_amyloid, as.numeric(exprs(m)["APP_7",])) # 3-42 pyro-Glu


# with 1-40
apps <- exprs(m)[c("APP_1", "APP_6","APP_7", "APP_10"),] %>%
   t() %>%
   as.data.frame() %>%
   rownames_to_column("ProjID") %>%
   inner_join(select(pData(m), ProjID, sqrt_amyloid, cogn_global_lv)) %>%
   pivot_longer(cols = c(APP_1, APP_6, APP_7, APP_10))

appsl <- apps %>%
   pivot_longer(cols = c(sqrt_amyloid, cogn_global_lv), names_to = "variable", values_to = "clin_val") %>%
   mutate(name = ordered(name, levels = c("APP_1", "APP_6", "APP_7", "APP_10"))) %>%
   mutate(variable = ordered(variable, levels = c("cogn_global_lv", "sqrt_amyloid")))




appsl <- appsl %>%
   mutate(variable0 = variable,
          name0 = name) %>%
   mutate(variable = ordered(variable, 
                             levels = c("cogn_global_lv", "sqrt_amyloid"), 
                             labels = c("Global~Cognition", "Amyloid")),
          name = ordered(name, 
                         levels = c("APP_1", "APP_6", "APP_7", "APP_10"),
                         labels = c(expression(atop(paste("A",beta[1-40]), scriptstyle("APP_1"))),
                                    expression(atop(paste("A",beta[1-42]), scriptstyle("APP_6"))),
                                    expression(atop(paste("A",beta[pE3-42]), scriptstyle("APP_7"))),
                                    expression(atop(paste("A",beta[9-42]), scriptstyle("APP_10"))))))

load("./output_data/spectral_counting_test_results.RData") # results
results <- results %>%
   filter(proteoform %in% c("APP_1", "APP_6", "APP_7", "APP_10"),
          variable %in% c("sqrt_amyloid", "cogn_global_lv")) %>%
   dplyr::rename(variable0 = variable,
                 name0 = proteoform)

# for Y coord
y_pos_tbl <- appsl %>% 
   group_by(name0) %>% 
   summarise(y_pos = 0.9 * max(value)) %>%
   as.data.frame()
rownames(y_pos_tbl) <- y_pos_tbl$name0

y_pos_tbl["APP_1","y_pos"] <- 128
y_pos_tbl["APP_6","y_pos"] <- 150
y_pos_tbl["APP_7","y_pos"] <- 24
y_pos_tbl["APP_10","y_pos"] <- 7.3


# for X coord
x_pos_tbl <- appsl %>% 
   group_by(variable0) %>% 
   summarise(maxVal = max(clin_val), 
             minVal = min(clin_val)) %>%
   mutate(x_pos = case_when(variable0 == "cogn_global_lv" ~ maxVal - 0.33 * (maxVal - minVal) + 0.3,
                            variable0 == "sqrt_amyloid" ~ minVal + 0.33 * (maxVal - minVal)) - 1.2) %>%
   select(variable0, x_pos)



results <- results %>%
   inner_join(x_pos_tbl) %>%
   inner_join(y_pos_tbl)

results <- inner_join(results, distinct(appsl, variable0, name0, variable, name))
results <- results %>%
   mutate(is_sig = adj.P.Val < 0.05)
sigs <- cut(results$adj.P.Val, breaks = c(0, 0.001, 0.01, 0.05, 1))
levels(sigs) <- c("***", "**", "*", "ns")
results$sigs <- sigs
results <- results %>%
   mutate(adj.P.Val = sprintf('%.01e', adj.P.Val)) %>%
   mutate(adj.P.Val = sprintf("%s (%s)", adj.P.Val, sigs))



ggplot(appsl) +
   aes(x = clin_val, y = value) +
   geom_point(shape=21, size=2.5, fill='#616161', color = "black") +
   geom_smooth(method = "glm", method.args = list(family = "quasipoisson"), se = T) +
   geom_label(data = results, aes(x = x_pos, y = y_pos, label = `adj.P.Val`), size = 5, fill = "#FAFAFA", hjust = "left") +
   facet_grid(name ~ variable, scales = "free",
              labeller = label_parsed) +
   theme_bw(base_size = 16) +
   theme(strip.text = element_text(size = 20,
                                   face = "plain"),
         strip.text.y.right = element_text(angle = -90),
         strip.background = element_blank(),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank()) +
   ylab("Spectral Counts") +
   xlab("Arbitrary Units")


ggsave(path = "figures_tables", filename = "Fig5_abeta_scatterplots.pdf", width = 7.6, height = 12)
ggsave(path = "figures_tables", filename = "Fig5_abeta_scatterplots.png", width = 7.6, height = 12)











# geom_text(data=test_results, aes(x=-0.5, y=-1.3, label=est_lbl),
#           colour="#606060",
#           fontface=4,
#           size=5, inherit.aes=FALSE, parse=T, adj=0) +
#    geom_text(data=test_results, aes(x=-0.5, y=-1.7, label=pval_lbl),
#              colour="#606060",
#              fontface=4,
#              size=5, inherit.aes=FALSE, parse=T, adj=0) +





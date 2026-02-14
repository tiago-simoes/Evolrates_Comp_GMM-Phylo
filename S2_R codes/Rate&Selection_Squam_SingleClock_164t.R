library(EvoPhylo)
library(ggtree)
library(ggplot2)
library(deeptime)
library(ape)

setwd("F:/Programas/Rcodes/PhyloParameters/Squamata_2024/1Mol1Mor_164t")

## Accessory function
#Ln and Min-Max Transf; range = [-1,+1]
normalize_ln_minmax <- function(x) {
  x_log<- log(x+ 1e-9)  # Adding a small constant to avoid log(0)
  x_normalized <-2 * ((x_log - min(x_log)) / (max(x_log) - min(x_log))) - 1
  return(x_normalized)
}

## Evolutionary Rates Statistics and Plots

### 1. Get rates from the clock tree and create a rate table

## Import summary tree produced by Mr. Bayes
tree<-treeio::read.mrbayes("Lepido_COMB_BayesCal_164t_SFBD_noSA_ILN_MCT_edited.t.con.tre") #edited as treeio misreads trees with multiple partitions within a single clock

#Convert MrB parameters class
tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)
tree@data <- tree@data %>% dplyr::mutate_if(is.character,as.numeric)
tidytree::get.fields(tree)

#Normalize tree rate means
#ln and min-max normalize rates
tree@data <- tree@data %>% mutate(across(where(is.numeric) & contains("rate"), normalize_ln_minmax))


### 2. Get and export the rate table
## Get table of clock rates with summary stats for each node in the tree for each relaxed clock partition
RateTable_Medians_no_clades <- get_clockrate_table_MrBayes(tree, summary = "median")
RateTable_Means_no_clades <- get_clockrate_table_MrBayes(tree, summary = "mean")

## Export the rate tables
write.csv(RateTable_Medians_no_clades, file="RateTable_Medians_Norm.csv", row.names = FALSE)
write.csv(RateTable_Means_no_clades, file="RateTable_Means_Norm.csv", row.names = FALSE)


### 3. Plot tree node labels

## Plot tree node labels
tree_nodes<-ggtree::ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
                 branch.length="none", size = 0.05)+
  geom_tiplab(size=2, linesize = 0.01, color="black",  offset = 0.5)+
  geom_label(aes(label=node), size=2, color="purple", position = "dodge")
tree_nodes

## Save your plot to your working directory as a PDF
ggplot2::ggsave	("Tree_nodes.pdf", width=8, height=20)



### 4. Get summary statistics table and plots

## Import rate table with clade membership (after new "clade" column added manually)
RateTable_Means<- read.csv("RateTable_Means_Clades1.csv", header = TRUE)
RateTable_MeansNorm<- read.csv("RateTable_Means_Norm_Clades.csv", header = TRUE)

## Get summary statistics table for each clade by clock
clockrate_summary(RateTable_Means, "Sum_RateTable_Means_Clades.csv", digits=2)
clockrate_summary(RateTable_MeansNorm, "Sum_RateTable_MeansNorm_Clades.csv", digits=2)

### 5. Plot rates by clock partition and clade

## Norm rates - Stacked plots with viridis color scale
clockrate_dens_plot(RateTable_MeansNorm, stack = TRUE, clock = 2, nrow = 2, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")

## Save your plot to your working directory as a PDF
ggplot2::ggsave	("RatesByClade_MeanNorm_Morpho_Stack.pdf", width=6, height=4.5)




##################  Selection mode inference ##########

### 1. Import and transform table

## Import rate table with clade membership
RateTable_Means<- read.csv("RateTable_Means_Clades1.csv", header = TRUE)

## Transform table from wide to long format
RatesByClade <- clock_reshape(RateTable_Means)


### 2. Import combined log file from all runs.

## Import all log (.p) files from all runs and combine them, with burn-in = 25% and downsampling to 2.5k trees in each log file
Comb_posterior <- combine_log("E:/Programas/Rcodes/PhyloParameters/Amniotes/TK02(Exp1)_6l_2p_Fix_(EvoPhylo)/P_files"
                              , burnin = 0.25, downsample = 2500)

### 3. Pairwise t-tests of Rate values

## Get table of pairwise t-tests for difference between the posterior mean and the rate for each tree node
RateSign_tests<- get_pwt_rates_MrBayes(RateTable_Means, Comb_posterior)

## Export the table
write.csv(RateSign_tests, file="RateSign_tests.csv")


### 5. Plot selection gradient on the summary tree (width=15, height=20)

## Plot tree using various thresholds for each clock partition
S1<- plot_treerates_sgn(type = "MrBayes",
                   tree, Comb_posterior, trans = "none",
                   clock = 1,               #Show rates for clock partition 1
                   summary = "mean",        #sets summary stats to get from summary tree nodes
                   branch_size = 1, tip_size = 3,                      #sets size for tree elements
                   xlim=c(-350,60), nbreaks = 10, geo_size=list(3, 3),   #sets limits and breaks for geoscale
                   geo_skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
                   threshold = c("1 SD", "3 SD"))+                       #sets threshold for selection mode
    labs(title = "Selection for cranial partition")

S1

S2<- plot_treerates_sgn(type = "MrBayes",
                   tree, Comb_posterior, trans = "none",
                   clock = 2,               #Show rates for clock partition 1
                   summary = "mean",        #sets summary stats to get from summary tree nodes
                   branch_size = 1, tip_size = 3,                      #sets size for tree elements
                   xlim=c(-350,60), nbreaks = 10, geo_size=list(3, 3),   #sets limits and breaks for geoscale
                   geo_skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
                   threshold = c("1 SD", "3 SD"))                 #sets threshold for selection mode
S2


S1|S2


## Save your plot to your working directory as a PDF
ggsave("125t_Rates_2P_TSign.pdf", width=30, height=20)


######################### SELECTION MODE #######################################
### Full TREE (mirrored) (5 x 8)
#Large mirrored multiclock tree (custom ggtree plot script)
#Using threshold for selection mode from outputted from S1 and S3 above

tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)
tree@data <- tree@data %>% dplyr::mutate_if(is.character,as.numeric)




#Postcranial rates
Pl2b<-ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
             aes(color=`rateTK02Brlens{2}_median`),
             branch.length="none",
             size=1)%>% flip(188, 151)%>% flip(189, 200)%>% flip(156,172)%>% flip(157,163)+
  #geom_tiplab(size=2, linesize = 0.01, color="black", fontface = "italic", offset=0, align=F)+
  scale_colour_steps2("Rate threshold",
                      low ="blue", mid= "gray", high ="red",
                      breaks=c(0.5550208, 0.8516736, 1.1483264, 1.4449792),
                      labels =c("-3 SD", "-1 SD", "+1 SD", "+3 SD"),
                      limits = c(0,2),
                      midpoint = 1)+
    labs(title = "Postcranial") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position="none")+
    scale_x_reverse()
Pl2b

#Combine plots (10x10)
par(mar = c(7, 7, 7, 7))
A1<-cowplot::plot_grid(Pl2a, Pl2b, ncol = 2)
A1

ggsave("125t_Rates_2P_TSign.pdf", width=10, height=10)





##############  EVOL RATE PLOTS


setwd("E:/Programas/Rcodes/PhyloParameters/Squamata_2024/1Mol1Mor_164t")

## Import summary tree produced by Mr. Bayes (1Mol1Mor)
tree<-treeio::read.mrbayes("Lepido_COMB_BayesCal_164t_SFBD_noSA_ILN_MCT_edited.t.con.tre") #edited as treeio misreads trees with multiple partitions within a single clock

#Normalize tree rate means
#ln and min-max normalize rates
tree@data <- tree@data %>% mutate(across(where(is.numeric) & contains("rate"), normalize_ln_minmax))

tidytree::get.fields(tree)

#Rect Node labels
ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE, 
       branch.length="none", size = 0.05)+
  geom_tiplab(size=2.5, linesize = 0.01, color="black",  offset = 0.2)+
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIgrBrlens_Morph_median, 0.5)), vjust=-.5, size=1)+
  geom_label(aes(label=node), size=3, color="purple", position = "dodge")

ggsave("MCT_NodeNo.pdf", width=20, height=22)

#Posterior only (dropped tips)

par(mar = c(7, 7, 7, 7))
Pl1 <- ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
              aes(size=I(prob))) +
  geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic") +
  geom_nodelab(aes(x=branch, label=round(prob, digits=2)), vjust=-.4, size=3) +
  coord_geo(xlim=c(-330,60), ylim=c(0,Ntip(tree)+1), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE), 
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"), alpha = 1, height = unit(1, "line"),
            rot = 0, size = list(3,4), neg = TRUE) +
  scale_x_continuous(breaks=seq(-330,60,20), labels=abs(seq(-330,60,20))) +
  theme_tree() +
  xlab("Time")
Pl1 <- revts(Pl1)
Pl1

ggsave("SFBD_noSA_ILN_1Mol1Mor_Probs.pdf", width=20, height=25)

######## Ages/ Age+Probs / Rates

#Age + Probs (FULL Tree)
Pl1a <- ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE, size=0.1) +
  geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic", offset=0.5) +
  geom_range(range='age_0.95HPD', color='purple', alpha=.6, size=2) +
  geom_nodelab(aes(label=round(age_median, digits = 1)), 
               vjust=+1.4, size=3) +
  geom_nodelab(aes(label=round(prob, digits=2)), color="black", fontface = "bold", 
               vjust=-0.5, size=3) +
  coord_geo(xlim=c(-330,60), ylim=c(0,Ntip(tree)+1), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE), 
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"), alpha = 1, height = unit(1.5, "line"),
            rot = 0, size = list(4,5), neg = TRUE) +
  scale_x_continuous(breaks=seq(-330,0,20), labels=abs(seq(-330,0,20))) +
  theme_tree() +
  geom_cladelabel(node=167, label="Early Archosauriformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA, offset=-200, offset.text=2) +
  geom_cladelabel(node=172, label="Stem Lepidosauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA, offset=-200, offset.text=2) +
  geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=244, label="Amph.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA, offset = -55, offset.text=4) +
  geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=311, label="'Haenophidia'", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +


Pl1a <- revts(Pl1a)
Pl1a

ggsave("SFBD_noSA_ILN_1Mol1Mor_Age&Probs_Clades.pdf", width=20, height=25)



#Single clock
#Ladder Tree - Rate log


Pl1b<-ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
             aes(color=rateIlnBrlensMorpho_mean), size=2)+
  geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic", offset=1)+
  geom_tiplab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50", 
                         name = "Ln(rel rate mean)\n Overall Morphology")+
 #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  geom_nodelab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  coord_geo(xlim=c(-330,60), ylim=c(0,Ntip(tree)+1), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE), 
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1.5, "line"),
            rot = 0, size = list(4,5), neg = TRUE) +
  scale_x_continuous(position = "bottom", breaks=seq(-360,30,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  theme(legend.position = "inside",
        legend.position.inside=c(.15, y=.1),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=14), 
        legend.key.height=unit(1.3, "cm"),
        legend.key.width=unit(1.3, "cm"))

Pl1b <- revts(Pl1b)
Pl1b

ggsave("SFBD_noSA_ILN_1Mol1Mor_RateMeanNorm(Rect).pdf", width=20, height=25)



#Full TREE (Circular)- Rate Log - no BL

Pl2a<-ggtree(tree,layout="fan", open.angle=10, ladderize=TRUE, right=TRUE,
             aes(color=rateIlnBrlensMorpho_mean), branch.length="none", size=2)+
  geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic", offset=0.2)+
  geom_tiplab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50", 
                         name = "Ln (rel rate mean)\n Overall Morphology")+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_x_continuous(position = "bottom", breaks=seq(-360,30,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  theme(legend.position = "inside",
        legend.position.inside=c(.15, y=.1),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=8), 
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))
Pl2a
ggsave("SFBD_noSA_ILN_1Mol1Mor_RateMeanNorm_fan_noBL.pdf", width=16, height=16)



#Full TREE (Circular)- Rate Log - BL (species names)

Pl2a<-ggtree(tree,layout="fan", open.angle=0, ladderize=TRUE, right=TRUE,
             aes(color=rateIlnBrlensMorpho_mean), size=3)+
  geom_tiplab(size=4,  color="black", fontface = "italic", offset=2,
              align=TRUE, linetype = "dotted", linesize = 0.001, alpha = 0.5)+
  #geom_tiplab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50", 
                         name = "Ln (rel rate mean)\n Overall Morphology")+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_x_continuous(position = "bottom", breaks=seq(-360,30,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.2),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=8), 
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))+
    geom_cladelabel(node=167, label="Archos.", align=TRUE, 
                    barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=172, label="E.Lepid.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=244, label="Amphisb.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
                    geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 22, taxa2= 25, label="E. Squam.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
               geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 142, taxa2= 149, label="'Haenophidia'", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 39, horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 132, taxa2= 136, label="E. Snakes", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
               geom = "text", color = "black", fill = NA)

Pl2a
ggsave("SFBD_noSA_ILN_1Mol1Mor_RateMeanNorm_fan_SppClades_BL.pdf", width=16, height=16)


#Full TREE (Circular) - Clade Names Only

Pl2b<-ggtree(tree, layout = "circular", ladderize=TRUE, right=TRUE,
             aes(color=rateIlnBrlensMorpho_mean), branch.length="none", size=3)+
  #geom_tiplab(size=4,  color="black", fontface = "italic", offset=2, align=TRUE, linetype = "dotted", linesize = 0.001, alpha = 0.5)+
  #geom_tiplab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50", 
                         name = "Ln(rate mean)\n Overall Morphology")+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIlnBrlensMorpho_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_x_continuous(position = "bottom", breaks=seq(-360,30,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  theme(legend.position = "inside",
        legend.position.inside=c(.05, y=.2),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=8), 
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))+
  geom_cladelabel(node=167, label="Archos.", align=TRUE, 
                  barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=172, label="E.Lep.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=244, label="Amphisb.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                  geom = "text", color = "black", fill = NA) +
  geom_strip(taxa1= 22, taxa2= 25, label="E. Squam.", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset.text =1.5,
             geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 142, taxa2= 149, label="'Haenophidia'", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 39, horizontal = FALSE, hjust='center', offset.text=1,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 132, taxa2= 136, label="E. Serp.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset.text=1.5,
               geom = "text", color = "black", fill = NA)

Pl2b

ggsave("SFBD_noSA_ILN_1Mol1Mor_RateMeanNorm_fan_Clades_noBL.pdf", width=14, height=14)





#Full TREE (mirrored) (5 x 8)

Pl2a<-ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
             aes(color=rateTK02Brlens1_median), 
             branch.length="none",
             size=1)%>% flip(188, 151)%>% flip(189, 200)%>% flip(156,172)%>% flip(157,163)+
  #geom_tiplab(size=2, linesize = 0.01, color="black", fontface = "italic", offset=0, align=F)+
  scale_x_continuous(breaks=seq(-360,0,50), labels=abs(seq(-360,0,50)))+
  scale_colour_gradient2(low ="blue", mid = "orange",high ="red",
                         midpoint = 1.5, space = "Lab", na.value = "grey50", 
                         guide = "colourbar",aesthetics = "colour")+
  theme(legend.position="none")+
  geom_cladelabel(node=128, label="Diadectidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA, hjust=0)+
  geom_cladelabel(node=246, label="Captorhinidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=249, label="Araeoscelidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  
  geom_cladelabel(node=211, label="Caseasauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=220, label="Varanopidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=229, label="Ophiacodontidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=241, label="Edaphosauridae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=232, label="Sphenacodontidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  
  geom_cladelabel(node=133, label="Acleistorhinidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=207, label="Milleretidae", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=138, label="Ankyromorpha", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=144, label="Younginiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=152, label="Testudines", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=169, label="Archosauriformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=164, label="Rhyncho+Allok", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=157, label="Protorosauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=173, label="Ichthyosauromorpha", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=179, label="Thalattosauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=181, label="Sauropterygia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=189, label="Sphenodontia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)+
  geom_cladelabel(node=200, label="Squamata", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 0,
                  geom = "text", color = "black", fill = NA)
Pl2a


Pl2b<-ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
             aes(color=rateTK02Brlens2_median), 
             branch.length="none",
             size=1)%>% flip(188, 151)%>% flip(189, 200)%>% flip(156,172)%>% flip(157,163)+
  #geom_tiplab(size=2, linesize = 0.01, color="black", fontface = "italic", offset=0, align=F,aes(angle=180))+
  scale_colour_gradient2(low ="blue", mid = "orange",high ="red",
                         midpoint = 1.5, space = "Lab", na.value = "grey50", 
                         guide = "colourbar",aesthetics = "colour")+
  theme(legend.position="none")+
  scale_x_reverse()
Pl2b

#Combine plots (10x10)
par(mar = c(7, 7, 7, 7))
A1<-plot_grid(Pl2a, Pl2b, ncol = 2)  
A1


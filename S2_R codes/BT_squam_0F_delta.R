library(ape)
library(phytools)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(deeptime)
library(cowplot)


## Accessory function
#Ln and Min-Max Transf; range = [-1,+1]
normalize_ln_minmax <- function(x) {
  x_log<- log(x+ 1e-9)  # Adding a small constant to avoid log(0)
  x_normalized <-2 * ((x_log - min(x_log)) / (max(x_log) - min(x_log))) - 1
  return(x_normalized)
}


##This code is used to import the results of BayesTraits RJ-MCMC analysaes

setwd("D:/Programas/BayesTraits/Datasets/Squam(122t)/PCs_95var/delta/Squam_0F")

tree_gmm<- read.nexus("pruned_tree0.nex")

Squam0F_rjpp <- rjpp(rjlog = "squam_0F_PCscores_run4.txt.VarRates.txt",
                     rjtrees =  "squam_0F_PCscores_run4.txt.Output.trees",
                     tree=tree_gmm,
                     burnin = 0, thinning = 10)

gmm_data0F<- Squam0F_rjpp$data

write.csv(gmm_data0F, file = "./results_plots/Squam0F_BT_phyloPC_stats_raw.csv", row.names = FALSE)




############# calculate additional variables for data table

# Corrected sum stats table (if necessary)

setwd("D:/Programas/BayesTraits/Datasets/Squam(122t)/PCs_95var/delta/Squam_0F/results_plots")

gmm_data<-read.csv("./Squam0F_BT_phyloPC_stats.csv")

#adjust gmm values
gmm_data<-gmm_data[-c(1),] #remove root stats (all 0 and not in tree_gmm)
gmm_data$node<-gmm_data$descNode

#Recalculate re-scaled mean branch rates (default rate output contains many 0s for large datasets)
#drop problematic estimates and recalculate rates
gmm_data <- subset(gmm_data, select=-c(meanRate,medianRate,modeRate))
gmm_data$meanRate <- (gmm_data$meanBL/gmm_data$orgBL)*gmm_data$meanDelta

# provide relative gmm rates (for normalization and comparability to clock rates)
gmm_data$RelmeanRate <- gmm_data$meanRate/mean(gmm_data$meanRate)

#Normalize tree rate means
#min-max normalize rates [-1,+1]
gmm_data$RelmeanRateNorm <- normalize_ln_minmax(gmm_data$RelmeanRate)

#ln transf rates
gmm_data$RelmeanRateLn <- log(gmm_data$RelmeanRate)

#Write post summary stats
write.csv(gmm_data, file = "./Squam0F_BT_phyloPC_stats_Final.csv", row.names = FALSE)



######### Get clade Membership for each node in the phylogenetic tree
#(necessary for downstream analyses)

#gmm_data<-read.csv("./Squam0F_BT_phyloPC_stats_Final.csv")
#tree_gmm<- read.nexus("pruned_tree0.nex")

#Rect Node labels
ggtree(tree_gmm, layout = "rectangular", ladderize=TRUE, right=TRUE,
       branch.length="none", size = 0.05)+
  geom_tiplab(size=2.5, linesize = 0.01, color="black",  offset = 0.2)+
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIgrBrlens_Morph_median, 0.5)), vjust=-.5, size=1)+
  geom_label(aes(label=node), size=3, color="purple", position = "dodge")

ggsave("Tree_92t_NodeNo.pdf", width=20, height=22)


# 1.Specify the ancestral nodes and their string names
ancestral_nodes <- list(
  Other_Lepidosauria = 76,
  Gekkota = 171,
  Scincoidea = 160,
  Lacertoidea = "148 - 156",
  Amphisbaenia = 156,
  Anguiformes = 114,
  Iguania = 100,
  Caenophidia = 127,
  Early_Serpentes = "123 - 127"
)


# Generate the data frame
Nodes_Clade_Table <- clade_membership(
  tree = "./pruned_tree0.nex",
  ancestral_nodes = ancestral_nodes,
  other_nodes_label = "Others")

## Merge data and clade table

gmm_data_clades<- merge(gmm_data, Nodes_Clade_Table)

#Write post summary stats
write.csv(gmm_data_clades, file = "Squam0F_stats_BT_phyloPC_Final_Clades.csv", row.names = FALSE)



################## Evolutionary Rates Statistics and Plots ##############


############## Evolutionary Rates Statistics

#gmm_data_clades<-read.csv("./Squam0F_BT_phyloPC_stats_Final_Clades.csv")

#Check rate dist
pl_rate<-  ggplot(gmm_data_clades, aes(x=RelmeanRate)) +
  geom_histogram(bins = 50, colour="black", fill="white") +
  geom_density(alpha=.2, fill="cyan") +
  geom_vline(aes(xintercept=mean(RelmeanRate)), color="red", linetype="dashed", linewidth=1) +
  labs(x = "Rates", y = "Density") +
  #xlim(c(0,quantile(gmm_data$RelmeanRate, 0.95))) +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold", hjust = 0.5))
pl_rate

ggplot2::ggsave	("rates_0F_RelmeanRate_dist.pdf", width=10, height=5)


#Check normalized rate dist
pl_rate_norm<-  ggplot(gmm_data_clades, aes(x=RelmeanRateNorm)) +
  geom_histogram(bins = 50, colour="black", fill="white") +
  geom_density(alpha=.2, fill="cyan") +
  geom_vline(aes(xintercept=mean(RelmeanRateNorm)), color="red", linetype="dashed", linewidth=1) +
  labs(x = "Rates", y = "Density") +
  #xlim(c(-quantile(gmm_data_norm$RelmeanRateNorm, 0.95),quantile(RateTable_Means_noCl_norm$RelmeanRateNorm, 0.95))) +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold", hjust = 0.5))
pl_rate_norm

ggplot2::ggsave	("rates_0F_RelmeanRateNorm_dist.pdf", width=10, height=5)


### Get summary statistics table and plots

#convert gmm_data to EvoPhylo stly tables
RateTable_Means_Clades<-data.frame(nodes = gmm_data_clades$node,
                            age_median = gmm_data_clades$end,
                            rates = gmm_data_clades$RelmeanRate,
                            clade = gmm_data_clades$clade)
RateTable_MeansNorm_Clades<-data.frame(nodes = gmm_data_clades$node,
                            age_median = gmm_data_clades$end,
                            rates = gmm_data_clades$RelmeanRateNorm,
                            clade = gmm_data_clades$clade)

write.csv(RateTable_Means_Clades, file="RateTable_Means_Clades.csv", row.names = FALSE)
write.csv(RateTable_MeansNorm_Clades, file="RateTable_MeansNorm_Clades.csv", row.names = FALSE)

## Get summary statistics table for each clade by clock
clockrate_summary(RateTable_Means_Clades, "Sum_RateTable_Means_Clades.csv", digits=2)
clockrate_summary(RateTable_MeansNorm_Clades, "Sum_RateTable_MeansNorm_Clades.csv", digits=2)

#
#RateTable_Means_Clades<- read.csv("./RateTable_Means_Clades.csv", header = TRUE)
#RateTable_MeansNorm_Clades<- read.csv("./RateTable_MeansNorm_Clades.csv", header = TRUE)

### Plot rates by clock partition and clade

#untransf rates
## Stacked plots with viridis color scale - Log10 Scale
BT_92t_Rel<- clockrate_dens_plot(RateTable_Means_Clades, stack = TRUE, nrow = 2, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")+
  ggplot2::scale_x_log10()+
  labs(x = "Log10 (relative rate mean)")
BT_92t_Rel
#ggplot2::ggsave	("BT_phyloPC_92t_RatesByClade_SkullSurface.pdf", width=8, height=6)


#Norm rates
## Stacked plots with viridis color scale
BT_92t_N<- clockrate_dens_plot(RateTable_MeansNorm_Clades, stack = TRUE, nrow = 4, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")+
  labs(x = " Relative rate mean (normalized)")
BT_92t_N
ggplot2::ggsave	("BT_phyloPC_92t_Rates_MeanNorm_ByClade_SkullSurface.pdf", width=8, height=6)



### Merge all Plots for BT 122-106-92t

cowplot::plot_grid(BT_122t_N, BT_106t_N, BT_92t_N,
          labels = c('A', 'B', 'C'),
          label_size = 8, label_fontface = "bold",
          align = "hv", ncol = 1)

ggplot2::ggsave	("D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons/BT_phyloPC_All_Rates_MeanNorm_ByClade_SkullSurface.pdf", width=4, height=12)


############ Tree rate plots

#merge time tree and gmm rate & data stats
tree_gmm_data <- full_join(tree_gmm, gmm_data, by = 'node')


#rel rate means - BL & tip labels
p1 <- ggtree(tree_gmm_data, aes(color=RelmeanRate), layout = 'circular',
        ladderize = FALSE, size=2) +
    scale_color_gradientn(colours=c('blue','cyan3','orange','red'),
                          transform = "log", na.value = "grey50",
                          name = "Rel rate mean)") +
    geom_tiplab(size=3, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
    theme_tree() +
      theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
          legend.title = element_text(size = 12, face = "bold", hjust=0.5),
          legend.text=element_text(size=12),
          legend.key.height=unit(1, "cm"),
          legend.key.width=unit(1, "cm"))
p1
#ggsave("BT_phyloPC_Tree92t_GMM_BL_Tips_RateMean_RelAbsolute(Circ).pdf", width=16, height=16)

#(circular) ln (rel rate means) - BL & tip labels
p2 <- ggtree(tree_gmm_data, aes(color=RelmeanRateLn), layout = 'circular',
        ladderize = FALSE, size=3) +
        geom_tiplab(size=3.5, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
       scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50",
                         name = "Relative rate mean \n (normalized)")+
      theme_tree() +
      theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
          legend.title = element_text(size = 12, face = "bold", hjust=0.5),
          legend.text=element_text(size=12),
          legend.key.height=unit(1, "cm"),
          legend.key.width=unit(1, "cm"))
p2
ggsave("BT_phyloPC_Tree92t_GMM_BL_Tips_RateMean_RelLn(Circ).pdf", width=16, height=16)

#(circular) rel rate means normalized [-1,+1] - BL & tip labels
p3 <- ggtree(tree_gmm_data, aes(color=RelmeanRateNorm), layout = 'circular',
        ladderize = FALSE, size=3) +
        geom_tiplab(size=3.5, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
       scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50",
                         name = "Relative rate mean \n (normalized)")+
      theme_tree() +
      theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
          legend.title = element_text(size = 12, face = "bold", hjust=0.5),
          legend.text=element_text(size=12),
          legend.key.height=unit(1, "cm"),
          legend.key.width=unit(1, "cm"))
p3
ggsave("BT_phyloPC_Tree92t_GMM_BL_Tips_RateMean_RelNorm(Circ).pdf", width=16, height=16)

#(Ladder Tree) rel rate means normalized [-1,+1] - BL & tip labels
Pl1a<-ggtree(tree_gmm_data, layout = "rectangular", ladderize=TRUE, right=TRUE,
                 aes(color=RelmeanRateLn), size=2)+
   geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
   geom_tiplab(aes(x=branch, label=round(RelmeanRateLn, 0.5)), vjust=-.5, size=3, color="black")+
   scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50",
                         name = "ln(relative rate mean)")+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  geom_nodelab(aes(x=branch, label=round(RelmeanRateLn, 0.5)), vjust=-.5, size=3, color="black")+
  coord_geo(xlim=c(-330,60), ylim=c(0,Ntip(tree_gmm_data)+1), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE),
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1.5, "line"),
            rot = 0, size = list(4,5), neg = TRUE) +
  scale_x_continuous(position = "bottom", breaks=seq(-330,60,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  #ggtitle("Evol Rates") +
  #theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust=-2))+
  theme(legend.position = "inside",
        legend.position.inside=c(.15, y=.1),
        legend.title = element_text(size =16, face = "bold", hjust=0.5),
        legend.text=element_text(size=14),
        legend.key.height=unit(1.3, "cm"),
        legend.key.width=unit(1.3, "cm"))+
  geom_cladelabel(node=171, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=160, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=149, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=156, label="Amph.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=114, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=112, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=101, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=146, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=138, label="'Haenophidia'", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=127, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4)


Pl1a <- revts(Pl1a)
Pl1a
ggsave("BT_phyloPC_Tree92t_GMM_BL_Tips_RateMean_RelLn_Tip+Clade(Rect).pdf", width=25, height=25)




# Tree circular plot with mean rates - BL
pl1c <- ggtree(tree_gmm_data, layout = "circular", ladderize=TRUE, right=TRUE,
                 aes(color=RelmeanRateNorm), size=3) +
    scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50",
                         name = "Relative rate mean \n (normalized)")+
    theme_tree() +
    theme(legend.position=c(.07, y=.2),
          legend.title = element_text(size = 14, face = "bold", hjust=0.5),
          legend.text=element_text(size=12),
          legend.key.height=unit(1, "cm"),
          legend.key.width=unit(1, "cm")) +
    ggtitle(paste0("Skull Surface GMM\n (100% fossils) ")) +
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust=-10))+
 geom_cladelabel(node=171, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=160, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=149, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=156, label="Amphisb.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=114, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=112, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=101, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=146, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=127, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 7,
                  geom = "text", color = "black", fill = NA) +
  geom_strip(taxa1= 85, taxa2= 20, label="'Haenophidia'", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = -17, horizontal = FALSE, hjust='center', offset.text = 8,
             geom = "text", color = "black", fill = NA)
pl1c

ggsave("BT_phyloPC_Tree92t_GMM_Tips_RateMean_RelLn_Clade(Circ_BL).pdf", width=14, height=14)






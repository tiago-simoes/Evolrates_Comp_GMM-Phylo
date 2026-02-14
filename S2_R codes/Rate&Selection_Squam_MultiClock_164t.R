library(EvoPhylo)
library(ggtree)
library(treeio)
library(dplyr)
library(ggplot2)
library(deeptime)
library(cowplot)
library(ape)

setwd("D:/Programas/Rcodes/PhyloParameters/Squamata_2024/4Mol4Mor_164t")

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
tree<-treeio::read.mrbayes("Comb_164t_SFBD_noSA_ILN_4Mol4Mor_FixTree_MCT.t.con.tre") #edited as treeio misreads trees with multiple partitions within a single clock

#### Adjust MrB parameters names and class for data analysis
tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)
tree@data <- tree@data %>% dplyr::mutate_if(is.character,as.numeric)
names(tree@data)  <- gsub("\\{|\\}", "", names(tree@data)) 
tidytree::get.fields(tree)


### 2. Get and export the rate table
## Get table of clock rates with summary stats for each node in the tree for each relaxed clock partition
RateTable_Medians <- get_clockrate_table_MrBayes(tree, summary = "median")
RateTable_Means <- get_clockrate_table_MrBayes(tree, summary = "mean")

## Export the rate tables
write.csv(RateTable_Medians, file="RateTable_Medians.csv", row.names = FALSE)
write.csv(RateTable_Means, file="RateTable_Means.csv", row.names = FALSE)


#Check rate dist
pl_rate<-  ggplot(RateTable_Means, aes(x=rates5)) +
  geom_histogram(bins = 50, colour="black", fill="white") +
  geom_density(alpha=.2, fill="cyan") +
  geom_vline(aes(xintercept=mean(rates5)), color="red", linetype="dashed", linewidth=1) +
  labs(x = "Rates", y = "Density") +
  #xlim(c(0,quantile(RateTable_Means$rates5, 0.95))) +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold", hjust = 0.5))
pl_rate

ggplot2::ggsave	("rates5.pdf", width=10, height=5)

## Check normalized rates (Z-transformed)
#RateTable_Means_noCl_norm <- as.data.frame(lapply(RateTable_Medians[,2:9], z_transf))

#Normalize tree rate means
tree_norm<-tree
tree_norm@data <- tree@data %>% mutate(across(where(is.numeric) & contains("rate"), normalize_ln_minmax))

RateTable_MeansNorm <- get_clockrate_table_MrBayes(tree_norm, summary = "mean")

write.csv(RateTable_MeansNorm, file="RateTable_MeansNorm.csv", row.names = FALSE)

pl_rate_norm<-  ggplot(RateTable_MeansNorm, aes(x=rates5)) +
  geom_histogram(bins = 50, colour="black", fill="white") +
  geom_density(alpha=.2, fill="cyan") +
  geom_vline(aes(xintercept=mean(rates5)), color="red", linetype="dashed", linewidth=1) +
  labs(x = "Rates", y = "Density") +
  #xlim(c(-quantile(RateTable_MeansNorm$rates5, 0.95),quantile(RateTable_Means_noCl_norm$rates5, 0.95))) +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold", hjust = 0.5))
pl_rate_norm

ggplot2::ggsave	("rates5_Norm.pdf", width=10, height=5)




#### 3. Get clade Membership for each node in the phylogenetic tree 
#(necessary for downstream analyses)

## Plot tree node labels
tree_nodes<-ggtree::ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
                           branch.length="none", size = 0.05)+
  geom_tiplab(size=2, linesize = 0.01, color="black",  offset = 0.5)+
  geom_label(aes(label=node), size=3, color="purple", position = "dodge")
tree_nodes

## Save your plot to your working directory as a PDF
ggplot2::ggsave	("MCT_164t_nodes.pdf", width=14, height=20)


#Specify the ancestral nodes and their string names
ancestral_nodes <- list(
  Other_Lepidosauria = c(175,172),
  Gekkota = 192,
  Scincoidea = 209,
  Lacertoidea = "226 - 244",
  Amphisbaenia = 244,
  Anguiformes = 272,
  Iguania = 251,
  Caenophidia = 317,
  Early_Serpentes = "297 - 317"
)


#Generate the data frame
Nodes_Clade_Table <- clade_membership(
  tree = "./MCT_164t.tre", 
  ancestral_nodes = ancestral_nodes, 
  other_nodes_label = "Others")
#####

## Merge Rates Table with clade Membership table (use Nodes_Clade_Table from above)

RateTable_Means_Clades<- merge(RateTable_Means, Nodes_Clade_Table)
RateTable_MeansNorm_Clades<- merge(RateTable_MeansNorm, Nodes_Clade_Table)

write.csv(RateTable_Means_Clades, file="RateTable_Means_Clades.csv", row.names = FALSE)
write.csv(RateTable_MeansNorm_Clades, file="RateTable_MeansNorm_Clades.csv", row.names = FALSE)

#OR: Import rate table with clade membership (after new "clade" column added manually)
#RateTable_Means<- read.csv("RateTable_Means_Clades1.csv", header = TRUE)
#RateTable_MeansNorm<- read.csv("RateTable_MeansNorm_Clades1.csv", header = TRUE)



### 4. Get summary statistics table and plots

## Get summary statistics table for each clade by clock
clockrate_summary(RateTable_Means_Clades, "Sum_RateTable_Means_Clades.csv", digits=2)
clockrate_summary(RateTable_MeansNorm_Clades, "Sum_RateTable_MeansNorm_Clades.csv", digits=2)

### 5. Plot rates by clock partition and clade

#untransf rates
## Stacked plots with viridis color scale - Log10 Scale
Clock4P_164t<-clockrate_dens_plot(RateTable_Means_Clades, clock = 5:8, stack = TRUE, nrow = 4, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")+
  ggplot2::scale_x_log10()+
  labs(x = " Relative rate mean")
Clock4P_164t
ggplot2::ggsave	("Clock4P_164t_RatesByClade_Mean_P5-P8.pdf", width=6, height=12)


#Norm rates
## Stacked plots with viridis color scale
## Stacked plots with viridis color scale
Clock4P_164t_N<- clockrate_dens_plot(RateTable_MeansNorm_Clades, clock = 5:8, stack = TRUE, nrow = 4, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")+
  labs(x = " Relative rate mean (normalized)")
Clock4P_164t_N
ggplot2::ggsave	("Clock4P_164t_RatesbyClade_MeanNORM_P5-P8.pdf", width=4, height=10)










############################  EVOL RATE PLOTS ##############################


## Import summary tree produced by Mr. Bayes
tree<-treeio::read.mrbayes("Comb_164t_SFBD_noSA_ILN_4Mol4Mor_FixTree_MCT.t.con.tre") #edited as treeio misreads trees with multiple partitions within a single clock

#### Adjust MrB parameters names and class for data analysis
tree@data$`prob+-sd` <- as.factor(tree@data$`prob+-sd`)
tree@data <- tree@data %>% dplyr::mutate_if(is.character,as.numeric)
names(tree@data)  <- gsub("\\{|\\}", "", names(tree@data)) 
tidytree::get.fields(tree)


#Normalize tree rate means
#ln and min-max normalize rates
tree_norm<-tree
tree_norm@data <- tree@data %>% mutate(across(where(is.numeric) & contains("rate"), normalize_ln_minmax))


#Rect Node labels
ggtree(tree_norm, layout = "rectangular", ladderize=TRUE, right=TRUE, 
       branch.length="none", size = 0.05)+
  geom_tiplab(size=2.5, linesize = 0.01, color="black",  offset = 0.2)+
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIgrBrlens_Morph_median, 0.5)), vjust=-.5, size=1)+
  geom_label(aes(label=node), size=3, color="purple", position = "dodge")

ggsave("Comb_164t_SFBD_noSA_ILN_4Mol4Mor_FixTree_MCT.pdf", width=20, height=22)



########Age+Probs 

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
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4)


Pl1a <- revts(Pl1a)
Pl1a

ggsave("SFBD_noSA_ILN_4Mol4Mor_Age&Probs_Clades.pdf", width=20, height=25)


 
#####

#Ladder Tree - Rate log

plot <- ggtree(tree_norm, layout = "rectangular", ladderize=TRUE, right=TRUE,
             aes(color=log(rateIlnBrlens8_mean)), size=2)+
  geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic", offset=1)+
  geom_tiplab(aes(x=branch, label=round(log(rateIlnBrlens8_mean, 0.5))), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50", 
                         name = "Relative rate mean \n (normalized))")+
 #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  geom_nodelab(aes(x=branch, label=round(rateIlnBrlens8_mean, 0.5)), vjust=-.5, size=3, color="black")+
  coord_geo(xlim=c(-330,60), ylim=c(0,Ntip(tree)+1), expand=FALSE,
            dat = list("epochs", "periods"), abbrv = list(TRUE, FALSE), 
            skip = c("Pliocene", "Pleistocene", "Holocene","Quaternary"),
            pos = list("bottom", "bottom"),alpha = 1, height = unit(1.5, "line"),
            rot = 0, size = list(4,5), neg = TRUE) +
  #scale_x_continuous(position = "bottom", breaks=seq(-330,60,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  #ggtitle("Evol Rates") +
  #theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust=-2))+
  theme(legend.position = "inside",
        legend.position.inside=c(.15, y=.1),
        legend.title = element_text(size =14, face = "bold", hjust=0.5),
        legend.text=element_text(size=14), 
        legend.key.height=unit(1.3, "cm"),
        legend.key.width=unit(1.3, "cm"))+
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
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4)

plot <- revts(plot)
plot

ggsave("4Mol4Mor_164t_RateMean_Ln_(Rect)_p5.pdf", width=20, height=25) 





#Full TREE (Circular)- Rate Log - BL (species names)

Pl2a<-ggtree(tree_norm,layout="fan", open.angle=0, ladderize=TRUE, right=TRUE,
             aes(color=rateIlnBrlens5_mean),  size=2)+
  geom_tiplab(size=4,  color="black", fontface = "italic", offset=2,
              align=TRUE, linetype = "dotted", linesize = 0.001, alpha = 0.5)+
  #geom_tiplab(aes(x=branch, label=round(rateIlnBrlens5_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50", 
                         name = "Relative (Ln) rate mean /n Skull Surface")+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIlnBrlens5_mean, 0.5)), vjust=-.5, size=3, color="black")+
  scale_x_continuous(position = "bottom", breaks=seq(-360,30,20), labels=abs(seq(-330,60,20)))+
  theme_tree()+
  theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=8), 
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))+
    geom_cladelabel(node=167, label="Archos.", align=TRUE, 
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
               barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
               geom = "text", color = "black", fill = NA)+
      geom_strip(taxa1= 132, taxa2= 136, label="E. Snakes", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset = 180, offset.text=10,
               geom = "text", color = "black", fill = NA)+
      geom_strip(taxa1= 6, taxa2= 7, label="E. Lep.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 75, horizontal = FALSE, hjust='center', offset = 180, offset.text=15,
               geom = "text", color = "black", fill = NA) 



Pl2a
ggsave("4Mol4Mor_MeanRateLn_Tip+Clades_p5(Circ_BL).pdf", width=16, height=16)


#Full TREE (Circular) - Clade Names (BL - Time-Tree)
# Function to plot and save tree for a given partition

plot_clockrate_by_part_circ <- function(partition_number) {
  partition_label <- paste0("rateIlnBrlens", partition_number, "_mean")
  
  plot <- ggtree(tree_norm, layout = "circular", ladderize=TRUE, right=TRUE,
                 aes_string(color=partition_label), size=3) +
    scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50", 
                           name = "Relative rate mean \n (normalized)") +
    theme_tree() +
    theme(legend.position=c(.07, y=.2),
          legend.title = element_text(size = 14, face = "bold", hjust=0.5),
          legend.text=element_text(size=12), 
          legend.key.height=unit(1, "cm"),
          legend.key.width=unit(1, "cm")) +
    geom_cladelabel(node=167, label="Archos.", align=TRUE, 
                    barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=244, label="Amphisb.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 22, taxa2= 25, label="E. Squam.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset.text = 15,
               geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 142, taxa2= 149, label="'Haenophidia'", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 39, horizontal = FALSE, hjust='center', offset.text = 15,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 132, taxa2= 136, label="E. Serp.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset.text=15,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 6, taxa2= 7, label="E. Lep.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 75, horizontal = FALSE, hjust='center', offset.text = 15,
               geom = "text", color = "black", fill = NA) 

  
  file_name <- paste0("4Mol4Mor_RateNorm_fan_BL_Clades_p", partition_number, ".pdf")
  ggsave(file_name, plot, width=14, height=14)
  return(plot)
}

# Produce individual and combined plots
# Function to combine all plots into a master plot
combine_plots_circ <- function(partition_numbers) {
  plots <- list()
  for (partition in partition_numbers) {
    plot <- plot_clockrate_by_part_circ(partition)
    plot <- plot + ggtitle(paste0("p_", partition, " rates")) +
      theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust=-10))
    plots[[length(plots) + 1]] <- plot
  }
  master_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2)
  ggsave("Morpho_AllParts_RateMeanNorm_(Circ_BL).pdf", master_plot, width=28, height=14 * length(partition_numbers)/2)
}

#Run
# Example: use 1:8, 5:8, or c(1, 3, 5) for specific partitions
combine_plots_circ(5:8)



#Full TREE (Circular) - Clade Names (no BL)  
# Function to plot and save tree for a given partition

plot_clockrate_by_part_circ <- function(partition_number) {
  partition_label <- paste0("rateIlnBrlens", partition_number, "_mean")
  
  plot <- ggtree(tree_norm, layout = "circular", ladderize=TRUE, right=TRUE,
                 aes_string(color=partition_label), branch.length="none", size=3) +
    scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                          na.value = "grey50", 
                           name = "Relative rate mean \n (normalized)") +
    theme_tree() +
    theme(legend.position=c(.07, y=.2),
          legend.title = element_text(size = 14, face = "bold", hjust=0.5),
          legend.text=element_text(size=12), 
          legend.key.height=unit(1, "cm"),
          legend.key.width=unit(1, "cm")) +
    geom_cladelabel(node=167, label="Archos.", align=TRUE, 
                    barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=244, label="Amphisb.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 1,
                    geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 22, taxa2= 25, label="E. Squam.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset.text = 1.5,
               geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 142, taxa2= 149, label="'Haenophidia'", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 39, horizontal = FALSE, hjust='center', offset.text = 1,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 132, taxa2= 136, label="E. Serp.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset.text=1.5,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 6, taxa2= 7, label="E. Lep.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 75, horizontal = FALSE, hjust='center', offset.text = 1.5,
               geom = "text", color = "black", fill = NA) 
  
  file_name <- paste0("4Mol4Mor_RateNorm_fan_noBL_Clades_p", partition_number, ".pdf")
  ggsave(file_name, plot, width=14, height=14)
  return(plot)
}

# Produce individual and combined plots
# Function to combine all plots into a master plot
combine_plots_circ <- function(partition_numbers) {
  plots <- list()
  for (partition in partition_numbers) {
    plot <- plot_clockrate_by_part_circ(partition)
    plot <- plot + ggtitle(paste0("p_", partition, " rates")) +
      theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust=-10))
    plots[[length(plots) + 1]] <- plot
  }
  master_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2)
  ggsave("Morpho_AllParts_RateMeanNorm_(Circ_noBL).pdf", master_plot, width=28, height=14 * length(partition_numbers)/2)
}

#Run
# Example: use 1:8, 5:8, or c(1, 3, 5) for specific partitions
combine_plots_circ(5:8)






##################  Selection mode inference ##########

### 1. Import and transform table

## Import rate table with clade membership
RateTable_Means<- read.csv("RateTable_Means_Clades1.csv", header = TRUE)

## Transform table from wide to long format
RatesByClade <- clock_reshape(RateTable_Means)


### 2. Import combined log file from all runs.

## Import all log (.p) files from all runs and combine them, with burn-in = 25% and downsampling to 2.5k trees in each log file
Comb_posterior <- combine_log("./P_files",
                              burnin = 0.25, downsample = 2500)
#OR
Comb_posterior<-read.csv("Comb_posterior.csv")

plot_back_rates2(type = "MrBayes", Comb_posterior, 
                 trans = "none", size = 10, quantile = 0.95)
plot_back_rates2(type = "MrBayes", Comb_posterior, 
                 trans = "log", size = 10, quantile = 0.95)

#Normalize? If yes,
#Comb_posterior <- Comb_posterior %>% mutate(across(where(is.numeric) & contains("rate"), z_transf))

### 3. Pairwise t-tests of Rate values

## Get table of pairwise t-tests for difference between the posterior mean and the rate for each tree node
RateSign_tests<- get_pwt_rates_MrBayes(RateTable_Means, Comb_posterior)

## Export the tables
write.csv(Comb_posterior, file="Comb_posterior.csv", row.names = FALSE)
write.csv(RateSign_tests, file="RateSign_tests.csv", row.names = FALSE)


### 5. Plot selection gradient on the summary tree (width=15, height=20)

#Normalize tree rate means? If yes,
tree_norm<-tree
tree_norm@data <- tree@data %>% mutate(across(where(is.numeric) & contains("rate"), z_transf))

## Plot tree using various thresholds for each clock partition
# Circular (BL)_Clade names - One Partition
S1<- plot_treerates_sgn2(type = "MrBayes", 
                         tree, Comb_posterior, trans = "none",
                         layout = "circular",
                         clock = 1,               #Show rates for clock partition 1
                         summary = "mean",        #sets summary stats to get from summary tree nodes
                         branch_size = 1, show_tip_labels = FALSE,                    #sets size for tree elements
                         threshold = c("1 SD", "2 SD", "3 SD"))+                       #sets threshold for selection mode
  labs(title = "skull surface")+
  geom_cladelabel(node=167, label="Archos.", align=TRUE, 
                  barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=244, label="Amphisb.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                  geom = "text", color = "black", fill = NA) +
  geom_strip(taxa1= 22, taxa2= 25, label="E. Squam.", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset.text = 15,
             geom = "text", color = "black", fill = NA) +
  geom_strip(taxa1= 142, taxa2= 149, label="'Haenophidia'", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 39, horizontal = FALSE, hjust='center', offset.text = 15,
             geom = "text", color = "black", fill = NA)+
  geom_strip(taxa1= 132, taxa2= 136, label="E. Serp.", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset.text=15,
             geom = "text", color = "black", fill = NA)+
  geom_strip(taxa1= 6, taxa2= 7, label="E. Lep.", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 75, horizontal = FALSE, hjust='center', offset.text = 15,
             geom = "text", color = "black", fill = NA) 
S1

## Save your plot to your working directory as a PDF
ggsave("4Mol4Mor_RateSign_TEST.pdf", width=14, height=14)

# Circular (BL)_Clade names - Function - Apply to all partitions

plot_rate_sgn_bypart <- function(partition_number) {
  partition_label <- paste0("rateIlnBrlens", partition_number, "_mean")
  
  plot <- plot_treerates_sgn2(type = "MrBayes", 
                              tree, Comb_posterior, trans = "none",
                              layout = "circular",
                              clock = partition_number,               #Show rates for clock partition 1
                              summary = "mean",        #sets summary stats to get from summary tree nodes
                              branch_size = 2.5, show_tip_labels = FALSE,                    #sets size for tree elements
                              threshold = c("2 SD", "3 SD"))+                       #sets threshold for selection mode
    labs(title = "partition_label")+
    geom_cladelabel(node=167, label="Archos.", align=TRUE, 
                    barsize = 1.5, fontsize = 5.5, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=175, label="Sphenodontia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto", horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=192, label="Gekkota", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=209, label="Scincoidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=227, label="Teiioidea", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=241, label="Lacert.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=244, label="Amphisb.", align=TRUE,
                    barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=272, label="Anguiformes", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=252, label="Acrodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=260, label="Pleurodonta", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=291, label="Mosasauria", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=317, label="Caenophidia", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_cladelabel(node=304, label="Scolec.", align=TRUE,
                    barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 10,
                    geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 22, taxa2= 25, label="E. Squam.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 40, horizontal = FALSE, hjust='center', offset.text = 15,
               geom = "text", color = "black", fill = NA) +
    geom_strip(taxa1= 142, taxa2= 149, label="'Haenophidia'", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 39, horizontal = FALSE, hjust='center', offset.text = 15,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 132, taxa2= 136, label="E. Serp.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 64, horizontal = FALSE, hjust='center', offset.text=15,
               geom = "text", color = "black", fill = NA)+
    geom_strip(taxa1= 6, taxa2= 7, label="E. Lep.", align=TRUE,
               barsize = 1.5, fontsize = 6, angle = 75, horizontal = FALSE, hjust='center', offset.text = 15,
               geom = "text", color = "black", fill = NA) 
  
  file_name <- paste0("4Mol4Mor_RateSign_Clades_p", partition_number, ".pdf")
  ggsave(file_name, plot, width=14, height=14)
  return(plot)
}

# Produce individual and combined plots
# Function to combine all plots into a master plot
combine_plots_ratesgn <- function(partition_numbers) {
  plots <- list()
  for (partition in partition_numbers) {
    plot <- plot_rate_sgn_bypart(partition)
    plot <- plot + ggtitle(paste0("p_", partition, " rates")) +
      theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust=-10))
    plots[[length(plots) + 1]] <- plot
  }
  master_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2)
  ggsave("Morpho_AllParts_RateSign.pdf", master_plot, width=28, height=14 * length(partition_numbers)/2)
}

#Run
# Example: use 1:8, 5:8, or c(1, 3, 5) for specific partitions
combine_plots_ratesgn(5:8)



# Circular (BL)_Species - Function - Apply to all partitions

plot_rate_sgn_bypart_spp <- function(partition_number) {
  partition_label <- paste0("rateIlnBrlens", partition_number, "_mean")
  
  plot <- plot_treerates_sgn2(type = "MrBayes", 
                              tree, Comb_posterior, trans = "none",
                              layout = "circular",
                              clock = partition_number,               #Show rates for clock partition 1
                              summary = "mean",        #sets summary stats to get from summary tree nodes
                              branch_size = 2.5, show_tip_labels = TRUE, tip_size=4,                    #sets size for tree elements
                              threshold = c("2 SD", "3 SD"))+                       #sets threshold for selection mode
    labs(title = "partition_label")
  
  file_name <- paste0("4Mol4Mor_RateSign_Spp_p", partition_number, ".pdf")
  ggsave(file_name, plot, width=16, height=16)
  return(plot)
}

# Produce individual and combined plots
# Function to combine all plots into a master plot
combine_plots_ratesgn_spp <- function(partition_numbers) {
  plots <- list()
  for (partition in partition_numbers) {
    plot <- plot_rate_sgn_bypart_spp(partition)
    plot <- plot + ggtitle(paste0("p_", partition, " rates")) +
      theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust=+3))
    plots[[length(plots) + 1]] <- plot
  }
  master_plot <- cowplot::plot_grid(plotlist = plots, ncol = 2)
  ggsave("Morpho_AllParts_RateSign_Spp.pdf", master_plot, width=36, height=19 * length(partition_numbers)/2)
}

#Run
# Example: use 1:8, 5:8, or c(1, 3, 5) for specific partitions
combine_plots_ratesgn_spp(5:8)




######################### EXTRAS #######################################
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


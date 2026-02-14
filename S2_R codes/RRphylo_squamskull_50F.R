library(ape)
library(RRphylo)
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


setwd("D:/Programas/Rcodes/RRphylo/Squam_skull/PCs_95var/Squam_50F")

tree<- read.nexus("pruned_tree50.nex")
PCs_data <- read.csv("PCs-95var_50F.csv", header = TRUE)
row.names(PCs_data) <- PCs_data[,1]
PCs_data<- PCs_data[-c(1)]

# Perform RRphylo and search.trend
RRsquam<- RRphylo(tree,PCs_data)

RR_BL_tip<-RRsquam$tip.path #root-tip BL  (matrix)
RR_BL_intnode<-RRsquam$node.path #root-int node BL  (matrix)
RR_RateMult<- RRsquam$rates # rates at every node
colnames(RR_RateMult)<- "all_PCs"
RR_RatesbyPC <- RRsquam$multiple.rates #multiv rates at every node
RR_anc<- RRsquam$aces #ancestral character values at ancestral nodes
RR_predpheno<- RRsquam$predicted.phenotype #the vector of estimated tip values. It is a matrix in the case of multivariate


#Create rates data table
RR_all_rates<-data.frame(node_name = row.names(RR_RateMult),
                         RR_RateMult,
                         RR_RatesbyPC, row.names = NULL)

RR_all_rates<-RR_all_rates[-c(1),]#remove root stats (not in tree_gmm)

#Rect Node labels
ggtree(tree, layout = "rectangular", ladderize=TRUE, right=TRUE,
       branch.length="none", size = 0.05)+
  geom_tiplab(size=2.5, linesize = 0.01, color="black",  offset = 0.2)+
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  #geom_nodelab(aes(x=branch, label=round(rateIgrBrlens_Morph_median, 0.5)), vjust=-.5, size=1)+
  geom_label(aes(label=node), size=3, color="purple", position = "dodge")

ggsave("Tree_NodeNo.pdf", width=20, height=22)


######### Get NODE NUMBERS for int and tip nodes in the phylogenetic tree
#get tip names and equivalent tip node # (use node = root node number)

tip_names<-tips(tree, node=107,labels=TRUE)
tip_values<-tips(tree, node=107,labels=FALSE)

#combine them
translate_df <- data.frame(node_number = as.character(tip_values),
                           taxa_name  = as.character(tip_names),
                           stringsAsFactors = FALSE)

# Merge data table with the translation table
gmm_rates <- RR_all_rates %>%
  left_join(translate_df, by = c("node_name" = "taxa_name")) %>%
  mutate(node_name = ifelse(!is.na(node_number), as.character(node_number), node_name)) %>%
  select(-node_number)  # Drop the 'node_number' column after replacement
gmm_rates$node<- as.integer(gmm_rates$node_name)
gmm_rates<-gmm_rates[-c(1)] #drop redundant node_name

write.csv(gmm_rates, "./All_rates_PCs5_50F.csv", row.names = FALSE)


### Obtain relative and normalized rates
gmm_rates$AllPCs_RelRate <- gmm_rates$all_PCs/mean(gmm_rates$all_PCs)
gmm_rates$AllPCs_RelRateLn <- log(gmm_rates$AllPCs_RelRate)
gmm_rates$AllPCs_RelRateNorm <- normalize_ln_minmax(gmm_rates$AllPCs_RelRate)


#Write post summary stats
write.csv(gmm_rates, file = "./All_rates_PCs_50F_Final.csv", row.names = FALSE)






################## Evolutionary Rates Statistics and Plots ##############

#gmm_rates<-read.csv("All_rates_phyloPCs_50F_Final.csv")
#tree<- read.nexus("pruned_tree50.nex")

######### Get clade Membership for each node in the phylogenetic tree
#(necessary for downstream analyses)

# 1.Specify the ancestral nodes and their string names
ancestral_nodes <- list(
  Other_Lepidosauria = c(92, 89),
  Gekkota = 197,
  Scincoidea = 186,
  Lacertoidea = "171 - 182",
  Amphisbaenia = 182,
  Anguiformes = 132,
  Iguania = 116,
  Caenophidia = 150,
  Early_Serpentes = "144 - 150"
)


# Generate the data frame
Nodes_Clade_Table <- clade_membership(
  tree = "./pruned_tree50.nex",
  ancestral_nodes = ancestral_nodes,
  other_nodes_label = "Others")


##
gmm_rates_clades<- merge(gmm_rates, Nodes_Clade_Table)

#Write post summary stats
write.csv(gmm_rates_clades, file = "All_rates_Squam50F_stats_RR_Final_Clades.csv", row.names = FALSE)


############## Evolutionary Rates Statistics

#setwd("D:/Programas/Rcodes/RRphylo/Squam_skull/PCs_95var_Stand/Squam_50F")
#gmm_rates_clades<-read.csv("All_rates_Squam50F_stats_RR_Final_Clades.csv")

#Check rate dist
pl_rate<-  ggplot(gmm_rates_clades, aes(x=AllPCs_RelRate)) +
  geom_histogram(bins = 50, colour="black", fill="white") +
  geom_density(alpha=.2, fill="cyan") +
  geom_vline(aes(xintercept=mean(AllPCs_RelRate)), color="red", linetype="dashed", linewidth=1) +
  labs(x = "Rates", y = "Density") +
  #xlim(c(0,quantile(gmm_rates_clades$AllPCs_AllPCs_RelRate, 0.95))) +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold", hjust = 0.5))
pl_rate

ggplot2::ggsave	("rates_50F_RelRate_dist.pdf", width=10, height=5)


#Check normalized rate dist
pl_rate_norm<-  ggplot(gmm_rates_clades, aes(x=AllPCs_RelRateNorm)) +
  geom_histogram(bins = 50, colour="black", fill="white") +
  geom_density(alpha=.2, fill="cyan") +
  geom_vline(aes(xintercept=mean(AllPCs_RelRateNorm)), color="red", linetype="dashed", linewidth=1) +
  labs(x = "Rates", y = "Density") +
  #xlim(c(-quantile(gmm_rates_clades$AllPCs_RelRateNorm, 0.95),quantile(RateTable_Means_noCl_norm$AllPCs_RelRateNorm, 0.95))) +
  theme_bw() +
  theme(plot.title = element_text(size = 2, face = "bold", hjust = 0.5))
pl_rate_norm

ggplot2::ggsave	("rates_50F_AllPCs_RelRateNorm_dist.pdf", width=10, height=5)


### Get summary statistics table and plots

#convert gmm_data to EvoPhylo stly tables
RateTable_Means_Clades<-data.frame(nodes = gmm_rates_clades$node,
                                   rates = gmm_rates_clades$AllPCs_RelRate,
                                   clade = gmm_rates_clades$clade)
RateTable_MeansNorm_Clades<-data.frame(nodes = gmm_rates_clades$node,
                                       rates = gmm_rates_clades$AllPCs_RelRateNorm,
                                       clade = gmm_rates_clades$clade)

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
RR_106t_Rel<- clockrate_dens_plot(RateTable_Means_Clades, stack = TRUE, nrow = 2, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")+
  ggplot2::scale_x_log10()+
  labs(x = "Log10 (relative rate mean)")
RR_106t_Rel
#ggplot2::ggsave	("RR_106t_RatesByClade_SkullSurface.pdf", width=8, height=6)


#Norm rates
## Stacked plots with viridis color scale
RR_106t_N<- clockrate_dens_plot(RateTable_MeansNorm_Clades, stack = TRUE, nrow = 4, scales = "fixed")+
  ggplot2::scale_colour_viridis_d(option = "turbo")+
  ggplot2::scale_fill_viridis_d(option = "turbo")+
  labs(x = " Relative rate mean (normalized)")
RR_106t_N
ggplot2::ggsave	("RR_106t_Rates_MeanNorm_ByClade_SkullSurface.pdf", width=8, height=6)







############# Tree rate plots

#merge time tree and gmm rate & data stats

tree_gmm_data <- full_join(tree, gmm_rates, by = 'node')

#rel rate
p1 <- ggtree(tree_gmm_data, aes(color=AllPCs_RelRate), layout = 'circular',
             ladderize = FALSE, size=2) +
  scale_color_gradientn(colours=c('blue','cyan3','orange','red'),
                        transform = "log", na.value = "grey50",
                        name = "Relative rates)") +
  geom_tiplab(size=3, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
  theme_tree() +
  theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=12),
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))
p1
#ggsave("Tree106t_GMM_BL_Tips_RelRate(Circ).pdf", width=16, height=16)


#(circular) rel rate Ln  - BL & tip labels
p2 <- ggtree(tree_gmm_data, aes(color=AllPCs_RelRateNorm), layout = 'circular',
             ladderize = FALSE, size=3) +
  geom_tiplab(size=3.5, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50",
                        name = "Ln (relative rates)")+
  theme_tree() +
  theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=12),
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))
p2
ggsave("Tree106t_GMM_BL_Tips_RelRatelLn(Circ).pdf", width=16, height=16)

#(circular) rel rate normalized [-1,+1] - BL & tip labels
p3 <- ggtree(tree_gmm_data, aes(color=AllPCs_RelRateNorm), layout = 'circular',
             ladderize = FALSE, size=3) +
  geom_tiplab(size=3.5, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50",
                        name = "Relative rates \n (normalized)")+
  theme_tree() +
  theme(legend.position = "inside",
        legend.position.inside=c(.1, y=.1),
        legend.title = element_text(size = 12, face = "bold", hjust=0.5),
        legend.text=element_text(size=12),
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm"))
p3
ggsave("Tree106t_GMM_BL_Tips_RelRateNorm(Circ).pdf", width=16, height=16)



#(Ladder Tree) rel rate means normalized [-1,+1] - BL & tip labels
Pl1a<-ggtree(tree_gmm_data, layout = "rectangular", ladderize=TRUE, right=TRUE,
             aes(color=AllPCs_RelRateNorm), size=2)+
  geom_tiplab(size=4, linesize = 0.01, color="black", fontface = "italic", hjust = -.1)+
  geom_tiplab(aes(x=branch, label=round(AllPCs_RelRateNorm, 0.5)), vjust=-.5, size=3, color="black")+
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50",
                        name = "Relative rates \n (normalized)")+
  #geom_range(range='rateTK02Brlens_0.95HPD', color='red', alpha=.6, size=2) +
  geom_nodelab(aes(x=branch, label=round(AllPCs_RelRateNorm, 0.5)), vjust=-.5, size=3, color="black")+
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
  geom_cladelabel(node=197, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=186, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=172, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=182, label="Amph.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=132, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=129, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=117, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=169, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=161, label="'Haenophidia'", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4) +
  geom_cladelabel(node=150, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = 90, hjust='center',
                  geom = "text", color = "black", fill = NA, offset = 45, offset.text=4)


Pl1a <- revts(Pl1a)
Pl1a
ggsave("Tree106t_GMM_BL_RelRateNorm_Tip+Clade(Rect).pdf", width=25, height=25)



# Tree circular plot with mean rates - BL
pl1c <- ggtree(tree_gmm_data, layout = "circular", ladderize=TRUE, right=TRUE,
               aes(color=AllPCs_RelRateNorm), size=3) +
  scale_color_gradientn(colours=c('blue','cyan3','orange',"red"),
                        na.value = "grey50",
                        name = "Relative rates \n (normalized)")+
  theme_tree() +
  theme(legend.position=c(.07, y=.2),
        legend.title = element_text(size = 14, face = "bold", hjust=0.5),
        legend.text=element_text(size=12),
        legend.key.height=unit(1, "cm"),
        legend.key.width=unit(1, "cm")) +
  ggtitle(paste0("Skull Surface GMM\n (100% fossils) ")) +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5, vjust=-10))+
  geom_cladelabel(node=197, label="Gekkota", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=186, label="Scincoidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=172, label="Teiioidea", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=182, label="Amphisb.", align=TRUE,
                  barsize = 1.5, fontsize = 5.5, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=132, label="Anguiformes", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=129, label="Acrodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=117, label="Pleurodonta", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=169, label="Scolec.", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_cladelabel(node=150, label="Caenophidia", align=TRUE,
                  barsize = 1.5, fontsize = 6, angle = "auto",horizontal = FALSE, hjust='center', offset.text = 6,
                  geom = "text", color = "black", fill = NA) +
  geom_strip(taxa1= 99, taxa2= 20, label="'Haenophidia'", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 72, horizontal = FALSE, hjust='center', offset.text = 10,
             geom = "text", color = "black", fill = NA) +
  geom_strip(taxa1= 75, taxa2= 75, label="Mosasauria", align=TRUE,
             barsize = 1.5, fontsize = 6, angle = 10, horizontal = FALSE, hjust='center', offset.text = -50,
             geom = "text", color = "black", fill = NA)
pl1c

ggsave("Tree106t_GMM_RelRateNorm_Clade(Circ_BL).pdf", width=14, height=14)





#####################################
# Trends in pheno and rates (MUST HAVE FOSSILS IN THE TREE)

setwd("D:/Programas/Rcodes/RRphylo/Squam_skull/PCs_95var_Stand/Squam_50F")

#Import data
tree<- read.nexus("pruned_tree50.nex")
PCs_data <- read.csv("PCs-95var_50F.csv", header = TRUE)
row.names(PCs_data) <- PCs_data[,1]
PCs_data<- PCs_data[-c(1)]

# Perform RRphylo and search.trend
RRsquam<- RRphylo(tree,PCs_data)

##### Trends

Trend <- search.trend(RRsquam,PCs_data)
Trend$rate.regression
Trend$phenotypic.regression

# Visualize search.trend results
pTrend <-plotTrend(Trend)
pTrend

par(mfrow=c(1,2),mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
PC1<-pTrend$plotSTphen("Comp1") #PC1
PC1<-pTrend$plotSTrates(1) #PC1
PC2_phen<-pTrend$plotSTphen("Comp2") #PC1
PC2_r<-pTrend$plotSTrates(1) #PC1
PC3_phen<-pTrend$plotSTphen("Comp3") #PC1
PC3_r<-pTrend$plotSTrates(1) #PC1
PC4_phen<-pTrend$plotSTphen("Comp4") #PC1
PC4_r<-pTrend$plotSTrates(1) #PC1
PC5_phen<-pTrend$plotSTphen("Comp5") #PC1
PC5_r<-pTrend$plotSTrates(1) #PC1


# Perform search.trend setting Smilodontini and Pantherini as individual clades
STclades <-search.trend(RRsquam,PCs_data,node=c(129,154))

# Visualize search.trend results
pTrend2 <-plotTrend(STclades)

#trait
par(mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
pTrend2$plotSTphenNode("PC1",node=1:2)
#rate
pTrend2$plotSTratesNode("PC1",node=c(154,129)) # This is the same as indicating node= 2:1

## Customizing parameters
par(mar=c(4,3,3,1),mgp=c(1.5,0.5,0))
#trait
pTrend2$plotSTphenNode("PC2",node=1:2,
                       plot.args = list(pch.node=c(23,24),pch=1,col="gray70",cex=1.2),
                       lineTree.args = list(col="black",lwd=3),lineNode.args = list(lwd=5),
                       node.palette = c("orangered","chartreuse"))
#rate
pTrend2$plotSTratesNode("rate",node=c(154,129),
                        plot.args = list(pch.node=c(5,6),pch=16,col="gray70",cex=1.2,lwd=2),
                        lineTree.args = list(col="gold",lwd=3,lty=4),
                        lineNode.args = list(lwd=5),
                        node.palette = c("deeppink","cyan2"))





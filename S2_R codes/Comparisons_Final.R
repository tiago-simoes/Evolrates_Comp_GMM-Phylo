# Load the packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(ape)
library(treeio)


############### Method comparison: Compare within and among clades rates for clock vs PCm analyses  #############
setwd("D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons")

# Load the table_meth_norm from a CSV file
table_meth_norm <- read.csv("Table_All_SkullSurface_Norm_PC.csv")


summary(table_meth_norm)

#### Summary stats by tree size

Sum_mean_by_tree <-FSA::Summarize(mean ~ tree_size,
                                  data=table_meth_norm,
                                  digits=2)
metric <- c("mean")
Sum_mean_by_tree <- cbind(metric, Sum_mean_by_tree)
#Sum_mean_by_tree <-Sum_mean_by_tree[-c(11)]

Sum_SD_by_tree <-FSA::Summarize(sd ~ tree_size,
                                data=table_meth_norm,
                                digits=2)
metric <- c("SD")
Sum_SD_by_tree <- cbind(metric, Sum_SD_by_tree)
Sum_SD_by_tree<- Sum_SD_by_tree[-c(4)]

Sum_Metrics_byTree <- rbind (Sum_mean_by_tree, Sum_SD_by_tree)

write.csv(Sum_Metrics_byTree, file="PC_Sum_Metrics_byTree.csv", row.names = F)




############# Plots #################

table_meth_norm$tree_size<-factor(table_meth_norm$tree_size, levels=c("122t_(24%_fossils)", "106t_(14%_fossils)", "92t_(0%_fossils)"), ordered = T)
table_meth_norm$analysis<-factor(table_meth_norm$analysis, levels=c("Clocks", "BayesTraits", "RRphylo"), ordered = T)
table_meth_norm$tree_size
table_meth_norm$clade<-factor(table_meth_norm$clade, levels=c("Other_Lepidosauria", "Gekkota", "Scincoidea",
                                                              "Lacertoidea","Amphisbaenia","Iguania",
                                                              "Anguiformes","Early_Serpentes","Caenophidia","Others"), ordered = T)
table_meth_norm$clade



##### Tree size by all methods
# Hist of rate SD by Tree
P1a<-ggplot(table_meth_norm, aes(x = tree_size, y = mean, fill = analysis)) +
  geom_violin(alpha=0.5,  trim=FALSE, draw_quantiles = 0.5, size=0.5)+
  #geom_boxplot(alpha=0.8, width = 0.8, size=0.3, notch = F)+
  #stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  labs(title = "All methods",
       x = NULL,
       y = "Normalized rate means",
       fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust=0),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  scale_fill_discrete()
P1a
#ggsave("ALl_NormRate mean by Fossil Sampling  (all analyses).pdf", width=6, height=4)


P1b<-ggplot(table_meth_norm, aes(x = tree_size, y = sd, fill = analysis)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = NULL,
       x = NULL,
       y = "Normalized rate means SD",
       fill = NULL) +
  theme_minimal() +
  theme(#plot.title = element_text(size = 10, face = "bold", hjust = 0.5, vjust=0),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")+
  scale_fill_discrete()
P1b
#ggsave("ALl_NormRate SD by Fossil Sampling (all analyses).pdf", width=6, height=4)





# Hist of rate SD by Tree
P2a<-ggplot(table_meth_norm, aes(x = analysis, y = mean, fill = tree_size)) +
  geom_violin(alpha=0.5,  trim=FALSE, draw_quantiles = 0.5, size=0.5)+
  #geom_boxplot(alpha=0.8, width = 0.8, size=0.3, notch = F)+
  #stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  labs(title = NULL,
       x = NULL,
       y = "Normalized rate means",
       fill = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust=0),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  scale_fill_viridis_d()
P2a
#ggsave("ALl_NormRate mean by Analysis (all Fossil Samplings).pdf", width=6, height=4)


##### Tree size by all methods
P2b<-ggplot(table_meth_norm, aes(x = analysis, y = sd, fill = tree_size)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = NULL,
       x = NULL,
       y = "Normalized rate means SD",
       fill = NULL) +
  theme_minimal() +
  theme(#plot.title = element_text(size = 10, face = "plain", hjust = 0.5, vjust=0),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")+
  scale_fill_viridis_d()
P2b
#ggsave("ALl_NormRate SD by Analysis (all Fossil Samplings).pdf", width=6, height=4)



# Boxplot of rate mean by Clade and Analysis (ALL TREES)

P3a<-ggplot(table_meth_norm, aes(x = clade, y = mean, fill = analysis)) +
  geom_boxplot(alpha=0.8, width = 0.8, size=0.3, notch = F)+
  #stat_summary(fun=mean, geom="point", shape=16, size=2, color = "black")+
  labs(title = NULL,
       x = NULL,
       y = "Normalized rate means",
       fill = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  scale_fill_discrete()
P3a
#ggsave("PhyloPC_NormRate Mean by Clade.pdf", width=6, height=4)



P3b<-ggplot(table_meth_norm, aes(x = clade, y = sd, fill = analysis)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = NULL,
       x = NULL,
       y = "Normalized rate means SD",
       fill = NULL) +
  theme_minimal() +
  theme(#plot.title = element_text(size = 10, face = "plain", hjust = 0.5, vjust=0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom")+
  scale_fill_discrete()
P3b
#ggsave("PhyloPC_NormRate mean by Analysis (all Fossil Samplings).pdf", width=6, height=4)




########

#Combined
G1a<-plot_grid(P2a, P2b,
          labels = c('A', 'B'),
          align="hv", ncol = 1,
          rel_heights = c(1.5,1))
G1a
#ggsave("PhyloPC_Comp_byMethod.pdf", width=6, height=8)


G1b<-plot_grid(P3a, P3b,
          labels = c('A', 'B'),
          align="hv", ncol = 1,
          rel_heights = c(1.5,1))
G1b
#ggsave("PhyloPC_Comp_byMethod&Clade.pdf", width=6, height=8)

plot_grid(P2a, P3a, P2b, P3b,
          labels = c('A', 'B', 'C', 'D'),
          align="hv", ncol = 2,
          rel_heights = c(1.5,1),
          rel_widths = c(1,1.5))
ggsave("PC_Comp_byMethod_Clade_Sampling.pdf", width=9, height=7)





###############################  RATES & TIME BY METHOD #######################
#Ln and Min-Max Transf; range = [-1,+1]
normalize_ln_minmax <- function(x) {
  x_log<- log(x+ 1e-9)  # Adding a small constant to avoid log(0)
  x_normalized <-2 * ((x_log - min(x_log)) / (max(x_log) - min(x_log))) - 1
  return(x_normalized)
}


setwd("D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons")

#Tree import 122t
tree_gmm<- read.nexus("pruned_tree100.nex")

# BayesTraits

#Get rate table with node numbers and ages
#gmm_data_BT <- read.csv("D:/Programas/BayesTraits/Datasets/Squam(122t)/PCs_Stand_95var/delta/Squam_100F/results_plots/RateTable_MeansNorm_Clades.csv", header = TRUE)
#gmm_data_BT$node<-gmm_data_BT$nodes

gmm_data_BT <- read.csv("D:/Programas/BayesTraits/Datasets/Squam(122t)/PCs_Stand_95var/delta/Squam_100F/results_plots/Squam100F_stats_BT_PC_Final_Clades.csv", header = TRUE)

RateTable_BT <- data.frame(anc_node_BT = gmm_data_BT$ancNode,
                           desc_node_BT = gmm_data_BT$descNode,
                           node = gmm_data_BT$node,
                           branch_dur = gmm_data_BT$orgBL,
                           age_median = gmm_data_BT$end,
                           rate_BT = gmm_data_BT$RelmeanRateNorm,
                           clade_BT = gmm_data_BT$clade)

#RateTable_BT$age_median <- round(RateTable_BT$age_median, digits=2)
#RateTable_BT$branch_dur <- round(RateTable_BT$branch_dur, digits=2)
#write.csv(RateTable_BT, "D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons/RateTable_BT.csv", row.names = F)


# RRphylo
#Get rate table with node numbers and ages
gmm_data_RR <- read.csv("D:/Programas/Rcodes/RRphylo/Squam_skull/PCs_95var_Stand/Squam_100F/All_rates_Squam100F_stats_RR_Final_Clades.csv", header = TRUE)

RateTable_RR<- data.frame(anc_node_RR = gmm_data_RR$ancestral_node,
                          node = gmm_data_RR$node,
                          rate_RR = gmm_data_RR$AllPCs_RelRateNorm,
                          clade_RR = gmm_data_RR$clade)


# Clocks

# route 1
#Get rate table with node numbers and ages
clock_data <- read.csv("D:/Programas/Rcodes/PhyloParameters/Squamata_2024/4Mol4Mor_122t_100F/RateTable_MeansNorm_Clades.csv", header = TRUE)

#clock_data<- clock_data[!(row.names(clock_data) %in% "123"),] #remove root node (not available for other methods)
clock_data$node<-clock_data$nodes

tree_clock<-read.nexus("MCT_122t.tre")

tree_clock_data <- full_join(tree_clock, clock_data, by = 'node')

# Drop the row where node2 equals 123
tree_clock_data@extraInfo <- tree_clock_data@extraInfo %>% filter(node != 123)

RateTable_Clocks<- data.frame(anc_node_clock = tree_clock_data@phylo$edge[,1],
                              desc_node_clock = tree_clock_data@phylo$edge[,2],
                              branch_dur = tree_clock_data@phylo$edge.length,
                              node_clock = tree_clock_data@extraInfo$node,
                              rate_clock = tree_clock_data@extraInfo$rates5,
                              age_median = tree_clock_data@extraInfo$age_median,
                              clade_clock = tree_clock_data@extraInfo$clade,
                              row.names = NULL)

write.csv(RateTable_Clocks, "D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons/RateTable_Clocks.csv", row.names = F)


##############



#Combine Tables

RateTable_All_temp <- merge(RateTable_BT, RateTable_RR, by = "node")

write.csv(RateTable_All_temp, "D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons/RateTable_All_temp_PC.csv", row.names = F)

RateTable_All <- merge(RateTable_All_temp, RateTable_Clocks2, by = "age_median")

write.csv(RateTable_All, "D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons/RateTable_All.csv", row.names = F)


#######Plot Rates Corr

setwd("D:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons")

RateTable_All<- read.csv("RateTable_All_Final_PC.csv")
names(RateTable_All)


## Clocks vs BT
Lm1a<- lm(data=RateTable_All, rate_Clock ~ rate_BT)
summary(Lm1a)

Lm1a_r2 <- summary(Lm1a)$adj.r.squared
Lm1a_p_value <- summary(Lm1a)$coefficients[2, 4]  # p-value for rate_RR coefficient

R1a <- ggplot(RateTable_All, aes(y = rate_Clock, x = rate_BT)) +
  geom_point() +
  geom_smooth(method ="lm", se=TRUE)+
  scale_x_continuous()+
  scale_y_continuous()+
  labs(y = "Clock rates",
       x = "BT rates") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Adj R² = ", round(Lm1a_r2, 2),
                                                        "\nP-value = ", round(Lm1a_p_value, 4)),
                hjust = 1.1, vjust = 1.1, size = 5, color = "black")

R1a


## Clocks vs RR
Lm1b<- lm(data=RateTable_All, rate_Clock ~ rate_RR)
summary(Lm1b)

Lm1b_r2 <- summary(Lm1b)$adj.r.squared
Lm1b_p_value <- summary(Lm1b)$coefficients[2, 4]  # p-value for rate_RR coefficient

# Create the plot
R1b <- ggplot(RateTable_All, aes(y = rate_Clock, x = rate_RR)) +
  geom_point() +
  geom_smooth(method ="lm", se=TRUE)+
  scale_x_continuous()+
  scale_y_continuous()+
  labs(y = "Clock rates",
       x = "RRphylo rates") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Adj R² = ", round(Lm1b_r2, 2),
                                                   "\nP-value = ", round(Lm1b_p_value, 4)),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black")

R1b

## BT vs RR
Lm1c<- lm(data=RateTable_All, rate_BT ~ rate_RR)
summary(Lm1c)

Lm1c_r2 <- summary(Lm1c)$adj.r.squared
Lm1c_p_value <- summary(Lm1c)$coefficients[2, 4]  # p-value for rate_RR coefficient

R1c <- ggplot(RateTable_All, aes(y = rate_BT, x = rate_RR)) +
  geom_point() +
  geom_smooth(method ="lm", se=TRUE)+
  scale_x_continuous()+
  scale_y_continuous()+
  labs( y = "BT rates",
       x = "RRphylo rates") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Adj R² = ", round(Lm1c_r2, 2),
                                                   "\nP-value = ", round(Lm1c_p_value, 4)),
           hjust = 1.1, vjust = 1.1, size = 5, color = "black")

R1c


#Combine
plot_grid(R1a, R1b, R1c,
          labels = c('A', 'B', 'C'),
          align="hv", ncol = 3)
ggsave("Fig_Rate_LinReg_PC.pdf", width=12, height=4)


############### Body parts comparison: Compare within and among clades for clock rates    #############

setwd("E:/Programas/Rcodes/PhyloParameters/Squamata_2024/Comparisons")


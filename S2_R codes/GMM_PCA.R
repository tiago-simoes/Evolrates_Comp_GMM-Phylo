#########################
###GMM Skull squamates###
#########################

library(rlist) #prunning tree
library(geomorph) #GMM
library(ggplot2) #graphics
library(ggtree) #graphics with phylo
library(ggrepel) #graphics
library(svglite) #graphics
library(tidyverse) #graphics
library(phytools)
library(ape)
library(phylotate)
library(rlist)
library(devtools) #graphics
library(treeio)
library(geiger)
library(readxl) #import
library(xlsx) #export
library(ggtreeSpace)

#Setting directory (ctrl+shift+H) 

#GMM ####
#Importing TPS file:
squama = readland.tps(file = 'Squama_skull_10Phylo.TPS',
                      specID = "ID")

#Defining semi-landmarks:
sliders = rbind(define.sliders(c(3, 36:49, 4)),
                define.sliders(c(34, 50:71, 16)),
                define.sliders(c(35, 72:93, 17)),
                define.sliders(c(21, 94:115,23)),
                define.sliders(c(22, 116:137,24)))

#Performing GPA:
GPA_squama = gpagen(squama, 
                    curves = sliders, 
                    ProcD = T,
                    approxBE = T, Proj = T)
plot(GPA_squama)


#Checking for outlier:
outliers_squama=plotOutliers(GPA_squama$coords,
                                inspect.outliers = T)

#Exploring PCA without phylogeny####
coords_aligned_squama=two.d.array(GPA_squama$coords)

##Grids and exploring graph
PCA_squama=gm.prcomp(GPA_squama$coords)
PCA_plot_squama=plot(PCA_squama)
picknplot.shape(PCA_plot_squama) #take specimens in the PCA plot to see deformation

broken_stick_fun = function(x){
  eigenvalues = x$sdev^2
  total_length = sum(eigenvalues)
  expected_lengths = numeric(length(eigenvalues))
  for (i in seq_along(eigenvalues)) {
    expected_lengths[i]=(total_length/(length(eigenvalues)-i+1))*i
  }
  analysis_result <- data.frame(
    Principal_Component = seq_along(eigenvalues),
    Eigenvalue = eigenvalues,
    Expected_Length = expected_lengths,
    Significance = eigenvalues > expected_lengths
  )
  
  return(analysis_result)
}

BS_PCA = broken_stick_fun(PCA_squama)
print(BS_PCA)

BS_PCA = BS_PCA[1:30,]

#Ploting broken stick results:
Fig_BS_PCA = ggplot(BS_PCA, aes(x=Principal_Component, y=Eigenvalue)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  geom_line(aes(x=Principal_Component, y = Expected_Length), stat = 'identity', color = 'red', group = 1) + 
  geom_point(aes(x=Principal_Component, y = Expected_Length)) + 
  labs(title = "Broken stick of 30 first PCs", x = "Principal Component", y = "Eigenvalue") + 
  scale_y_continuous("Eigenvalue", sec.axis = sec_axis(~., name = "Expected Length"))
Fig_BS_PCA

#Exporting PCA data do .xlsx tables:
DF.coords_aligned_squama = as.data.frame(coords_aligned_squama)
write.xlsx(DF.coords_aligned_squama, file = "Coords_aligned_squama.xlsx",
           sheetName = 'coords_aligned_squama', col.names = T, row.names = T)

DF.PCA_squama = as.data.frame(PCA_squama$x)
write.xlsx(DF.PCA_squama, file = "PCA_squama.xlsx",
           sheetName = 'PCA_squama', col.names = T, row.names = T)

write.xlsx(PCA_squama$rotation, file = "loadings.xlsx", 
           col.names = T, row.names = T)

CS = as.data.frame(GPA_squama$Csize)
write.xlsx(CS, file = "CS_squama.xlsx",
           sheetName = 'CS', col.names = T, row.names = T)

#PCA graphics#####
Comp_PCA_squama <- read_excel("Compiled data_squama.xlsx", 
                                   sheet = "PCA")

hull_Comps12 <- Comp_PCA_squama %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))
color_phygroups <- read_excel("Classification.xlsx", 
                                                 sheet = "phygroups")

material = Comp_PCA_squama$Material
material = replace(material, material == 'fossil', 22)
material = replace(material, material == 'extant', 21)
material = as.numeric(material)

#Fig_PCA12 = 
  ggplot(data = Comp_PCA_squama, aes(x=Comp1, y= -Comp2, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps12, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                     linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material, size = 3, color = 'black',fill = Comp_PCA_squama$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  labs(x = 'PC 1 (24.92%)', y = 'PC 2 (15.37%)') + #choosing axes labs
  ggtitle("PCA (extant + 100% fossil)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
               force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
Fig_PCA12
#EXPORT 800 X 500 px

hull_Comps13 <- Comp_PCA_squama %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp3))

#Fig_PCA13 = 
ggplot(data = Comp_PCA_squama, aes(x=Comp1, y= Comp3, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps13, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                     linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material, size = 3, color = 'black',fill = Comp_PCA_squama$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  labs(x = 'PC 1 (24.92%)', y = 'PC 3 (8.04%)') + #choosing axes labs
  ggtitle("PCA (extant + 100% fossil)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
#EXPORT 800 X 500 px
Fig_PCA13

hull_Comps23 <- Comp_PCA_squama %>% group_by(Group.2) %>% 
  slice(chull(Comp2, Comp3))

#Fig_PCA23 = 
ggplot(data = Comp_PCA_squama, aes(x=-Comp2, y= Comp3, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps23, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                     linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material, size = 3, color = 'black',fill = Comp_PCA_squama$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  labs(x = 'PC 2 (15.37%)', y = 'PC 3 (8.03%)') + #choosing axes labs
  ggtitle("PCA (extant + 100% fossil)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                 force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
#EXPORT 800 X 500 px
  
###Ploting deformation in relation to the consensus:####
Links <- read.table("Links.txt", quote="\"", comment.char="")
Links = as.matrix(Links)

#Ploting min and max along PCs:
scoresPC1 = PCA_squama$x[,"Comp1"]
shapesPC1 = PCA_squama$shapes$shapes.comp1

plotRefToTarget(GPA_squama$consensus,
                shapesPC1$min,
                links = Links, method = 'vector')
plotRefToTarget(GPA_squama$consensus,
                shapesPC1$max,
                links = Links, method = 'vector')

scoresPC2 = PCA_squama$x[,"Comp2"]
shapesPC2 = PCA_squama$shapes$shapes.comp2
plotRefToTarget(GPA_squama$consensus,
                shapesPC2$min,
                links = Links, method = 'vector')
plotRefToTarget(GPA_squama$consensus,
                shapesPC2$max,
                links = Links, method = 'vector')

scoresPC3 = PCA_squama$x[,"Comp3"]
shapesPC3 = PCA_squama$shapes$shapes.comp3
plotRefToTarget(GPA_squama$consensus,
                shapesPC3$min,
                links = Links, method = 'vector')
plotRefToTarget(GPA_squama$consensus,
                shapesPC3$max,
                links = Links, method = 'vector')


####including phylogenetic data####
#Get original tree with a full set of taxa
#@S3 class tree
tree_full <- read.nexus("Comb_164t_SFBD_noSA_ILN_4Mol4Mor_FixTree_MCT.t.con.tre")
tree_full = ladderize(tree_full)

taxa_table <- read_excel("Classification.xlsx", 
                              sheet = "taxa_gmm")

keep_taxa <- taxa_table$name_final

#Drop the tips not included in the reduced tree or keep task a list
#@S3 class tree
pruned_tree <- ape::keep.tip(tree_full, keep_taxa)

plot.phylo(pruned_tree)

names_phylo = pruned_tree$tip.label
names_phylo = as.data.frame(names_phylo)

name.check(pruned_tree, coords_aligned_squama)

## Exporting the tree
write.tree(pruned_tree, file = "pruned_tree.tre") # saves as .tre 
writeNexus(pruned_tree, file = "pruned_tree.nex") # save as.nex
#dput(ttrees_equal, file = "pruned_tree, equal method.txt") # for keeping the root age (good for other analyses in R because it keeps the root age)


#Including phylogenetic data to PCA:####
### PaCA Alignment of data to physlogenetic signal ###
#PaCA (OLS)
#PACA_OLS_squama = gm.prcomp(coords_aligned_squama, phy = pruned_tree,
#                            align.to.phy = T, GLS = F) ### Attention! This is a PACA, not PCA (alignment of data to phylo signal rather than variance)
#summary(PACA_OLS_squama)
#plot(PACA_OLS_squama, phylo = TRUE, main = "PACA (OLS)")

#PaCA (PGLS) 
#PACA_GLS_squama = gm.prcomp(coords_aligned_squama, phy = pruned_tree,
#                            align.to.phy = T, GLS = T) ### Attention! This is a PACA, not PCA (alignment of data to phylo signal rather than variance)
#summary(PACA_GLS_squama)
#plot(PACA_GLS_squama, phylo = TRUE, main = "PACA (PGLS)")

### Phylo PCA ###
### Phylo PCA (PGLS) with phylo corrected (transformed) projection
phyPCA_squama = gm.prcomp(coords_aligned_squama, phy = pruned_tree,
                          GLS = TRUE, transform = TRUE, align.to.phy = F)
summary(phyPCA_squama)
plot(phyPCA_squama, phylo = TRUE, main = "phylo PCA")

BS_pPCA = broken_stick_fun(phyPCA_squama)
print(BS_pPCA)
BS_pPCA = BS_pPCA[1:30,]

#Ploting broken stick results:
Fig_BS_pPCA = ggplot(BS_pPCA, aes(x=Principal_Component, y=Eigenvalue)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  geom_line(aes(x=Principal_Component, y = Expected_Length), stat = 'identity', color = 'red', group = 1) + 
  geom_point(aes(x=Principal_Component, y = Expected_Length)) + 
  labs(title = "Broken stick of 30 first PCs", x = "Principal Component", y = "Eigenvalue") + 
  scale_y_continuous("Eigenvalue", sec.axis = sec_axis(~., name = "Expected Length"))
Fig_BS_pPCA

#Exporting data:
DF.phyPCA_squama = as.data.frame(phyPCA_squama$x)
write.xlsx(DF.phyPCA_squama, file = 'phyPCA_squama.xlsx',
           sheetName = 'phyPCA_squama', col.names = T, row.names = T)

ancestors.phyPCA_squama = as.data.frame(phyPCA_squama$ancestors)
write.xlsx(ancestors.phyPCA_squama, file = 'ancestors_phyPCA_squama.xlsx',
           sheetName = 'ancestors_phyPCA_squama', col.names = T, row.names = T)

write.xlsx(phyPCA_squama$rotation, file = "phyloadings.xlsx", 
           col.names = T, row.names = T)

Comp_phyloPCA_squama <- read_excel("Compiled data_squama.xlsx", 
                              sheet = "phyloPCA")
row.names(Comp_phyloPCA_squama) = Comp_phyloPCA_squama$File_Name

#Checking for phylogenetic signal PCs####
#Test for phylogenetic signal in shape (Kmult)
PS_K.shape <- physignal(A = coords_aligned_squama, phy = pruned_tree)
summary(PS_K.shape)
plot(PS_K.shape) #sign phylo sign as expected under BM (K=1), or even higher than by BM (K>1)
plot(PS_K.shape$PACA, phylo = TRUE)

# Phylogenetic signal profile (K stat by PC)
# important to see where the phylo signal is concentrated!
PS_K.shape$K.by.p # First PCs much higher signal than last PCs
plot(PS_K.shape$K.by.p)

#Test for phylogenetic signal in shape (Z-stat)
PS_Z.shape <- physignal.z(A = coords_aligned_squama, phy = pruned_tree,
                          lambda = "front")
summary(PS_Z.shape)

# Problem with ill-conditioned residual covariance matrix;  shaving last several PCS keeping most relevant ones (e.g., top 10)
PS_Z.shape <- physignal.z(A = coords_aligned_squama, phy = pruned_tree,
                          lambda = "front", PAC.no = 43) #93.59%
summary(PS_Z.shape)
plot(PS_Z.shape)
plot(PS_Z.shape$PACA, phylo = TRUE)

#Checking for phylogenetic signal CS####
#Test for phylogenetic signal in shape (Kmult)
coords.CS_squama <- read_excel("Compiled data_squama.xlsx", 
                               sheet = "coords_aligned and CS")

log.CS = as.matrix(coords.CS_squama$Log.CS)
row.names(log.CS) = coords.CS_squama$File_Name 

PS_K.CS <- physignal(A = log.CS, phy = pruned_tree)
summary(PS_K.CS)
plot(PS_K.CS) #sign phylo sign as expected under BM (K=1), or even higher than by BM (K>1)

PS_Z.CS <- physignal.z(A = log.CS, phy = pruned_tree,
                          lambda = "front")
summary(PS_Z.CS)
plot(PS_Z.CS)

#Performing contMap with PCs and CS:####
#see: http://blog.phytools.org/2021/11/showing-different-orders-or-other.html
plotTree(pruned_tree,ftype="i",fsize=0.8)
nodelabels(cex=0.6,frame="circle",bg="white") #checking the node numbers

color_phygroups_nodes <- read_excel("Classification.xlsx", 
                                    sheet = "node")

nodes = color_phygroups_nodes$Node #Nodes of target groups
clades = color_phygroups_nodes$Group
clade_colors = color_phygroups_nodes$Color

squama_groups=pruned_tree
for(i in 1:length(nodes)){
  squama_groups<-paintSubTree(squama_groups,node=nodes[i],
                              state=clades[i],anc.state=if(i==1) "NA" else NULL,
                              stem=T)
}

#Checking if the nodes and colors match
plot(squama_groups, colors <- setNames(c('grey', clade_colors),
                                       c('NA', clades)), ftype = 'i', fsize = 0.8)
nodelabels(cex=0.6, frame = 'circle', bg = 'white')

#Performing Contmap for pPC1:####
pComp1 = -phyPCA_squama$x[, 1]

fit_pPC1 <- phytools::fastAnc(pruned_tree, pComp1, vars=TRUE, CI=TRUE)

td_pPC1 <- data.frame(node = nodeid(pruned_tree, names(pComp1)),
                 trait = pComp1)
nd_pPC1 <- data.frame(node = names(fit_pPC1$ace), trait = fit_pPC1$ace)

d_pPC1 <- rbind(td_pPC1, nd_pPC1)
d_pPC1$node <- as.numeric(d_pPC1$node)
tree_pPC1 <- full_join(pruned_tree, d_pPC1, by = 'node')

#Ploting: 
Fig_contpPC1 =  ggtree(tree_pPC1, layout = 'circular', ladderize = T,
         right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
             continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 1 (14.28%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
    geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
    geom_cladelab(node = 133, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
    geom_cladelab(node = 145, label = 'Acrod.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#56AA21', barcolor = '#56AA21') + 
    geom_cladelab(node = 151, label = 'Neoang.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#0B6345', barcolor = '#0B6345') + 
    geom_cladelab(node = 159, label = 'Paleo.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#45AA88', barcolor = '#45AA88') + 
    geom_cladelab(node = 170, label = 'Caenophidia', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#2B386B', barcolor = '#2B386B') + 
    geom_cladelab(node = 189, label = 'Scolec.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
    geom_cladelab(node = 192, label = 'Teiioidea', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#F1CA3A', barcolor = '#F1CA3A') + 
    geom_cladelab(node = 205, label = 'Amphisb.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#8C373E', barcolor = '#8C373E') + 
    geom_cladelab(node = 209, label = 'Lacert.', align = T,
                  barsize = 1, fontsize = 4, angle = 60, horizontal = F,
                  hjust = 'center', offset.text = 40,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
    geom_cladelab(node = 210, label = 'Scincoidea', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#251EDE', barcolor = '#251EDE') + 
    geom_cladelab(node = 221, label = 'Gekkota', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#7DFF56', barcolor = '#7DFF56') + 
    geom_cladelab(node = 237, label = 'Sphen.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#715128', barcolor = '#715128') + 
    geom_cladelab(node = 242, label = 'Arch.', align = T,
                  barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                  hjust = 'center', offset.text = 15,
                  geom = 'text', color = 'black', fill = NA,
                  textcolor = '#EA4F0D', barcolor = '#EA4F0D') + 
    geom_strip(taxa1 = 100, taxa2 = 107, label = "'Haenop.'", align=TRUE,
               barsize = 1, fontsize = 4, angle = 60, horizontal = F,
               hjust = 'center', offset.text = 15,
               geom = 'text', color = '#5566AB', fill = NA) + 
    geom_strip(taxa1 = 93, taxa2 = 93, label = "Mosasauria", align=TRUE,
               barsize = 1, fontsize = 4, angle = 0, horizontal = F,
               hjust = 'center', offset.text = 10,
               geom = 'text', color = '#FE922A', fill = NA)

Fig_contpPC1
#save as PDF, in size 10x10 inches

#Performing Contmap for pPC2:####
pComp2 = -phyPCA_squama$x[, 2]

fit_pPC2 <- phytools::fastAnc(pruned_tree, pComp2, vars=TRUE, CI=TRUE)

td_pPC2 <- data.frame(node = nodeid(pruned_tree, names(pComp2)),
                      trait = pComp2)
nd_pPC2 <- data.frame(node = names(fit_pPC2$ace), trait = fit_pPC2$ace)

d_pPC2 <- rbind(td_pPC2, nd_pPC2)
d_pPC2$node <- as.numeric(d_pPC2$node)
tree_pPC2 <- full_join(pruned_tree, d_pPC2, by = 'node')

#Ploting: 
Fig_contpPC2 = ggtree(tree_pPC2, aes(color=trait), layout = 'circular', size = 1.5, 
                      ladderize = F, continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 2 (10.99%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 133, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 145, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 151, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 159, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 170, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 189, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 192, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') + 
  geom_cladelab(node = 205, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 209, label = 'Lacert.', align = T,
                barsize = 1, fontsize = 4, angle = 60, horizontal = F,
                hjust = 'center', offset.text = 40,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 210, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 221, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7DFF56', barcolor = '#7DFF56') + 
  geom_cladelab(node = 237, label = 'Sphen.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#715128', barcolor = '#715128') + 
  geom_cladelab(node = 242, label = 'Arch.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#EA4F0D', barcolor = '#EA4F0D') + 
  geom_strip(taxa1 = 100, taxa2 = 107, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 60, horizontal = F,
             hjust = 'center', offset.text = 15,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 93, taxa2 = 93, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)
Fig_contpPC2

#Performing Contmap for pPC3:####
pComp3 = -phyPCA_squama$x[, 3]

fit_pPC3 <- phytools::fastAnc(pruned_tree, pComp3, vars=TRUE, CI=TRUE)

td_pPC3 <- data.frame(node = nodeid(pruned_tree, names(pComp3)),
                      trait = pComp3)
nd_pPC3 <- data.frame(node = names(fit_pPC3$ace), trait = fit_pPC3$ace)

d_pPC3 <- rbind(td_pPC3, nd_pPC3)
d_pPC3$node <- as.numeric(d_pPC3$node)
tree_pPC3 <- full_join(pruned_tree, d_pPC3, by = 'node')

#Ploting: 
Fig_contpPC3 = ggtree(tree_pPC3, aes(color=trait), layout = 'circular', size = 1.5, 
                      ladderize = F, continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 3 (8.11%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 133, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 145, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 151, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 159, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 170, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 189, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 192, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') + 
  geom_cladelab(node = 205, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 209, label = 'Lacert.', align = T,
                barsize = 1, fontsize = 4, angle = 60, horizontal = F,
                hjust = 'center', offset.text = 40,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 210, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 221, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7DFF56', barcolor = '#7DFF56') + 
  geom_cladelab(node = 237, label = 'Sphen.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#715128', barcolor = '#715128') + 
  geom_cladelab(node = 242, label = 'Arch.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#EA4F0D', barcolor = '#EA4F0D') + 
  geom_strip(taxa1 = 100, taxa2 = 107, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 60, horizontal = F,
             hjust = 'center', offset.text = 15,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 93, taxa2 = 93, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)
Fig_contpPC3 

#Performing Contmap for log.CS:####

fit_CS <- phytools::fastAnc(pruned_tree, log.CS, vars=TRUE, CI=TRUE)

td_CS <- data.frame(node = nodeid(pruned_tree, names(log.CS)),
                      trait = log.CS)
nd_CS <- data.frame(node = names(fit_CS$ace), trait = fit_CS$ace)

d_CS <- rbind(td_CS, nd_CS)
d_CS$node <- as.numeric(d_CS$node)
tree_CS <- full_join(pruned_tree, d_CS, by = 'node')

#Ploting: 
Fig_contCS = ggtree(tree_CS, aes(color=trait), layout = 'circular', size = 1.5, 
                      ladderize = F, continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'Log10(CS)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 133, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 145, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 151, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 159, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 170, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 189, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 192, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') + 
  geom_cladelab(node = 205, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 209, label = 'Lacert.', align = T,
                barsize = 1, fontsize = 4, angle = 60, horizontal = F,
                hjust = 'center', offset.text = 40,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 210, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 221, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7DFF56', barcolor = '#7DFF56') + 
  geom_cladelab(node = 237, label = 'Sphen.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#715128', barcolor = '#715128') + 
  geom_cladelab(node = 242, label = 'Arch.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#EA4F0D', barcolor = '#EA4F0D') + 
  geom_strip(taxa1 = 100, taxa2 = 107, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 60, horizontal = F,
             hjust = 'center', offset.text = 15,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 93, taxa2 = 93, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)
Fig_contCS


plot_list(Fig_contpPC1, Fig_contpPC2, Fig_contpPC3, Fig_contCS, 
          ncol=2, nrow = 2, tag_levels='A')

#Ploting phylomorphospace####

#Plotting in phylomorphospace
taxon_names = data.frame(number_tip = c(1:122), names = pruned_tree[["tip.label"]])
taxon_names = taxon_names[order(taxon_names$names),]
taxon_names = data.frame(taxon_names, number = Comp_phyloPCA_squama$Number, color = Comp_phyloPCA_squama$Color)
taxon_names = taxon_names[order(taxon_names$number_tip),]

MphyPCA_squama12 = as.matrix(-phyPCA_squama$x[,1:2])
MphyPCA_squama13 = as.matrix(-phyPCA_squama$x[,c(1,3)])
MphyPCA_squama23 = as.matrix(-phyPCA_squama$x[,2:3])

# Flip the desired principal components (e.g., the first principal component)

#Ploting#### 

phylomorphospace(squama_groups, MphyPCA_squama12, colors = colors,
                ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                xlab = 'pPC 1 (14.28%)', ylab = 'pPC 2 (10.99%)', 
                lwd = 2)
#legend('bottomright', names(colors), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names$number, frame = 'none', col = taxon_names$color) #To add specimen number
#500 x 500 px

phylomorphospace(squama_groups, MphyPCA_squama13, colors = colors,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 1 (14.28%)', ylab = 'pPC 3 (8.11%)', 
                 lwd = 2)
#legend('bottomright', names(colors), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names$number, frame = 'none', col = taxon_names$color) #To add specimen number
#500 x 500 px

phylomorphospace(squama_groups, MphyPCA_squama23, colors = colors,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 2 (10.99%)', ylab = 'pPC 3 (8.11%)', 
                 lwd = 2)
#legend('bottomright', names(colors), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names$number, frame = 'none', col = taxon_names$color) #To add specimen number
#500 x 500 px

#Ploting phylomorphospace with convex hulls####
hull_phyComps12 <- Comp_phyloPCA_squama %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))

ggplot(data = Comp_phyloPCA_squama, aes(x=-Comp1, y=-Comp2, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps12, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                     linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material, size = 3, color = 'black',fill = Comp_phyloPCA_squama$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 1 (14.28%)', y = 'pPC 2 (10.99%)') + #choosing axes labs
  ggtitle("phyloPCA") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right') 

#EXPORT 800 X 500 px

hull_phyComps13 <- Comp_phyloPCA_squama %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp3))

ggplot(data = Comp_phyloPCA_squama, aes(x=-Comp1, y=-Comp3, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps13, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                        linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material, size = 3, color = 'black',fill = Comp_phyloPCA_squama$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 1 (14.28%)', y = 'pPC 3 (8.11%)') + #choosing axes labs
  ggtitle("phyloPCA") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right') 

#EXPORT 800 X 500 px

hull_phyComps23 <- Comp_phyloPCA_squama %>% group_by(Group.2) %>% 
  slice(chull(Comp2, Comp3))

ggplot(data = Comp_phyloPCA_squama, aes(x=-Comp2, y=-Comp3, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps23, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                        linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material, size = 3, color = 'black',fill = Comp_phyloPCA_squama$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 2 (10.99%)', y = 'pPC 3 (8.11%)') + #choosing axes labs
  ggtitle("phyloPCA") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right') 

#EXPORT 800 X 500 px

################################################################################
###Fossil subsampling 50%#### 
#Importing TPS file:
squama50 = readland.tps(file = 'Squama_skull_10Phylo_sub50.TPS',
                        specID = "ID")

#Performing GPA:
GPA_squama50 = gpagen(squama50, 
                      curves = sliders, 
                      ProcD = T,
                      approxBE = T, Proj = T)
plot(GPA_squama50)


#Checking for outlier:
outliers_squama50=plotOutliers(GPA_squama50$coords,
                               inspect.outliers = T)

#Exploring PCA without phylogeny####
coords_aligned_squama50=two.d.array(GPA_squama50$coords)

##Grids and exploring graph
PCA_squama50=gm.prcomp(GPA_squama50$coords)
PCA_plot_squama50=plot(PCA_squama50)

BS_squama50 = broken_stick_fun(PCA_squama50)
print(BS_squama50)
BS = BS_squama50[1:30,]

#Ploting broken stick results:
ggplot(BS, aes(x=Principal_Component, y=Eigenvalue)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  geom_line(aes(x=Principal_Component, y = Expected_Length), stat = 'identity', color = 'red', group = 1) + 
  geom_point(aes(x=Principal_Component, y = Expected_Length)) + 
  labs(title = "Broken stick of 30 first PCs", x = "Principal Component", y = "Eigenvalue") + 
  scale_y_continuous("Eigenvalue", sec.axis = sec_axis(~., name = "Expected Length"))

#Exporting PCA data do .xlsx tables:
DF.coords_aligned_squama50 = as.data.frame(coords_aligned_squama50)
write.xlsx(DF.coords_aligned_squama50, file = "Coords_aligned_squama50.xlsx",
           sheetName = 'coords_aligned_squama50', col.names = T, row.names = T)

DF.PCA_squama50 = as.data.frame(PCA_squama50$x)
write.xlsx(DF.PCA_squama50, file = "PCA_squama50.xlsx",
           sheetName = 'PCA_squama50', col.names = T, row.names = T)

write.xlsx(PCA_squama50$rotation, file = "loadings50.xlsx", 
           col.names = T, row.names = T)

CS50 = as.data.frame(GPA_squama50$Csize)
write.xlsx(CS50, file = "CS_squama50.xlsx",
           sheetName = 'CS50', col.names = T, row.names = T)

#PCA graphics#####
#Classifiers <- read.csv("C:/Users/Arthur/Arthur/GMM/TPS/Classifiers.csv", sep=";")
Comp_PCA_squama50 <- read_excel("Compiled data_squama.xlsx", 
                                sheet = "PCA_sub50")

hull_Comps12.50 <- Comp_PCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))

Comp_PCA_squama50 <- read_excel("Compiled data_squama.xlsx", 
                                sheet = "PCA_sub50")

hull_Comps12.50 <- Comp_PCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))

material50 = Comp_PCA_squama50$Material
material50 = replace(material50, material50 == 'fossil', 22)
material50 = replace(material50, material50 == 'extant', 21)
material50 = as.numeric(material50)

#Fig_PCA12.50 = 
ggplot(data = Comp_PCA_squama50, aes(x=Comp1, y= Comp2, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps12.50, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                        linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material50, size = 3, color = 'black',fill = Comp_PCA_squama50$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  labs(x = 'PC 1 (26.81%)', y = 'PC 2 (14.62%)') + #choosing axes labs
  ggtitle("PCA (extant + 50% fossil)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                 force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
Fig_PCA12.50
#EXPORT 800 X 500 px

hull_Comps13.50 <- Comp_PCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp3))

#Fig_PCA13.50 = 
  ggplot(data = Comp_PCA_squama50, aes(x=Comp1, y= -Comp3, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps13.50, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                        linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material50, size = 3, color = 'black',fill = Comp_PCA_squama50$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  labs(x = 'PC 1 (26.81%)', y = 'PC 3 (8.82%)') + #choosing axes labs
  ggtitle("PCA (extant + 50% fossil)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                 force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
#EXPORT 800 X 500 px
Fig_PCA13.50

hull_Comps23.50 <- Comp_PCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp2, Comp3))

#Fig_PCA23.50 = 
  ggplot(data = Comp_PCA_squama50, aes(x=Comp2, y= -Comp3, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps23.50, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                        linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material50, size = 3, color = 'black',fill = Comp_PCA_squama50$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  labs(x = 'PC 2 (14.62%)', y = 'PC 3 (8.82%)') + #choosing axes labs
  ggtitle("PCA (extant + 50% fossil)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                 force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
#EXPORT 800 X 500 px
Fig_PCA23.50

####including phylogenetic data####

taxa_table50 <- read_excel("Classification.xlsx", 
                           sheet = "taxa_gmm_sub50")
keep_taxa50 <- taxa_table50$name_final

#Drop the tips not included in the reduced tree or keep tak a list
#@S3 class tree
pruned_tree50 <- ape::keep.tip(tree_full, keep_taxa50)

plot.phylo(pruned_tree50)

names_phylo50 = pruned_tree50$tip.label
names_phylo50 = as.data.frame(names_phylo50)

name.check(pruned_tree50, coords_aligned_squama50)

## Exporting the tree
write.tree(pruned_tree50, file = "pruned_tree50.tre") # saves as .tre 
writeNexus(pruned_tree50, file = "pruned_tree50.nex") # save as.nex
#dput(ttrees_equal, file = "pruned_tree, equal method.txt") # for keeping the root age (good for other analyses in R because it keeps the root age)


#Including phylogenetic data to PCA:####
### PaCA Alignment of data to physlogenetic signal ###
#PaCA (OLS)
#PACA_OLS_squama50 = gm.prcomp(coords_aligned_squama50, phy = pruned_tree50,
#                            align.to.phy = T, GLS = F) ### Attention! This is a PACA, not PCA (alignment of data to phylo signal rather than variance)
#summary(PACA_OLS_squama50)
#plot(PACA_OLS_squama50, phylo = TRUE, main = "PACA (OLS)")

#PaCA (PGLS) 
#PACA_GLS_squama50 = gm.prcomp(coords_aligned_squama50, phy = pruned_tree50,
#                            align.to.phy = T, GLS = T) ### Attention! This is a PACA, not PCA (alignment of data to phylo signal rather than variance)
#summary(PACA_GLS_squama50)
#plot(PACA_GLS_squama50, phylo = TRUE, main = "PACA (PGLS)")

### Phylo PCA ###
### Phylo PCA (PGLS) with phylo corrected (transformed) projection
phyPCA_squama50 = gm.prcomp(coords_aligned_squama50, phy = pruned_tree50,
                          GLS = TRUE, transform = TRUE)
summary(phyPCA_squama50)
plot(phyPCA_squama50, phylo = TRUE, main = "phylo PCA (50% fossil)")

BS50 = broken_stick_fun(phyPCA_squama50)
print(BS50)
BS = BS50[1:30,]

#Ploting broken stick results:
ggplot(BS, aes(x=Principal_Component, y=Eigenvalue)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  geom_line(aes(x=Principal_Component, y = Expected_Length), stat = 'identity', color = 'red', group = 1) + 
  geom_point(aes(x=Principal_Component, y = Expected_Length)) + 
  labs(title = "Broken stick of 30 first PCs", x = "Principal Component", y = "Eigenvalue") + 
  scale_y_continuous("Eigenvalue", sec.axis = sec_axis(~., name = "Expected Length"))

#Exporting data:
DF.phyPCA_squama50 = as.data.frame(phyPCA_squama50$x)
write.xlsx(DF.phyPCA_squama50, file = 'phyPCA_squama50.xlsx',
           sheetName = 'phyPCA_squama50', col.names = T, row.names = T)

ancestors.phyPCA_squama50 = as.data.frame(phyPCA_squama50$ancestors)
write.xlsx(ancestors.phyPCA_squama50, file = 'ancestors_phyPCA_squama50.xlsx',
           sheetName = 'ancestors_phyPCA_squama50', col.names = T, row.names = T)

Comp_phyloPCA_squama50 <- read_excel("Compiled data_squama.xlsx", 
                                   sheet = "phyloPCA_sub50")
row.names(Comp_phyloPCA_squama50) = Comp_phyloPCA_squama50$File_name

#Checking for phylogenetic signal####
#Test for phylogenetic signal in shape (Kmult)
PS_K.shape50 <- physignal(A = coords_aligned_squama50, phy = pruned_tree50)
summary(PS_K.shape50)
plot(PS_K.shape50) #sign phylo sign as expected under BM (K=1), or even higher than by BM (K>1)
plot(PS_K.shape50$PACA, phylo = TRUE)

# Phylogenetic signal profile (K stat by PC)
# important to see where the phylo signal is concentrated!
PS_K.shape50$K.by.p # First PCs much higher signal than last PCs

#Test for phylogenetic signal in shape (Z-stat)
PS_Z.shape50 <- physignal.z(A = coords_aligned_squama50, phy = pruned_tree50,
                          lambda = "front")
summary(PS_Z.shape50)

#Test for phylogenetic signal in shape (Kmult)
# Problem with ill-conditioned residual covariance matrix;  shaving last several PCS keeping most relevant ones (e.g., top 10)
PS_Z.shape50 <- physignal.z(A = coords_aligned_squama50, phy = pruned_tree50,
                          lambda = "front", PAC.no = 35) #99.97%
summary(PS_Z.shape50)
plot(PS_Z.shape50)
plot(PS_Z.shape50$PACA, phylo = TRUE)

#Checking for phylogenetic signal in CS####
#Test for phylogenetic signal in shape (Kmult)
coords.CS_squama50 <- read_excel("Compiled data_squama.xlsx", 
                                 sheet = "coords_aligned and CS_sub50")

row.names(coords.CS_squama50) = coords.CS_squama50$File_Name
log.CSv50 = coords.CS_squama50$Log.CS
names(log.CSv50) = coords.CS_squama50$File_Name 

PS_K.CS50 <- physignal(A = log.CSv50, phy = pruned_tree50)
summary(PS_K.CS50)
plot(PS_K.CS50) #sign phylo sign as expected under BM (K=1), or even higher than by BM (K>1)
plot(PS_K.CS50$PACA, phylo = TRUE)

# Phylogenetic signal profile (K stat by PC)
# important to see where the phylo signal is concentrated!
PS_K.CS50$K.by.p # First PCs much higher signal than last PCs

#Test for phylogenetic signal in shape (Z-stat)
PS_Z.CS50 <- physignal.z(A = log.CSv50, phy = pruned_tree50,
                            lambda = "front")
summary(PS_Z.CS50)

summary(PS_Z.CS50)
plot(PS_Z.CS50)

#Performing Contmap for pPC1:####
pComp1_50 = phyPCA_squama50$x[, 1]

fit_pPC1_50 <- phytools::fastAnc(pruned_tree50, pComp1_50, vars=TRUE, CI=TRUE)

td_pPC1_50 <- data.frame(node = nodeid(pruned_tree50, names(pComp1_50)),
                         trait = pComp1_50)
nd_pPC1_50 <- data.frame(node = names(fit_pPC1_50$ace), trait = fit_pPC1_50$ace)

d_pPC1_50 <- rbind(td_pPC1_50, nd_pPC1_50)
d_pPC1_50$node <- as.numeric(d_pPC1_50$node)
tree_pPC1_50 <- full_join(pruned_tree50, d_pPC1_50, by = 'node')

#Ploting: 

Fig_contpPC1_50 =  ggtree(tree_pPC1_50, layout = 'circular', ladderize = T,
                          right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                          continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 1 (15.6%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 117, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 129, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 133, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 140, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 150, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 169, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 182, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 172, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 186, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 197, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 42, taxa2 = 42, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 3, taxa2 = 3, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = -10, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Arch.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 30,
             geom = 'text', color = '#EA4F0D', fill = NA) + 
  geom_strip(taxa1 = 84, taxa2 = 91, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 77, taxa2 = 77, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 15, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)

Fig_contpPC1_50
#save as PDF, in size 10x10 inches

#Performing Contmap for pPC2:####
pComp2_50 = phyPCA_squama50$x[, 2]

fit_pPC2_50 <- phytools::fastAnc(pruned_tree50, pComp2_50, vars=TRUE, CI=TRUE)

td_pPC2_50 <- data.frame(node = nodeid(pruned_tree50, names(pComp2_50)),
                         trait = pComp2_50)
nd_pPC2_50 <- data.frame(node = names(fit_pPC2_50$ace), trait = fit_pPC2_50$ace)

d_pPC2_50 <- rbind(td_pPC2_50, nd_pPC2_50)
d_pPC2_50$node <- as.numeric(d_pPC2_50$node)
tree_pPC2_50 <- full_join(pruned_tree50, d_pPC2_50, by = 'node')

#Ploting: 

Fig_contpPC2_50 =  ggtree(tree_pPC2_50, layout = 'circular', ladderize = T,
                          right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                          continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 2 (11.82%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 117, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 129, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 133, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 140, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 150, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 169, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 182, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 172, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 186, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 197, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 42, taxa2 = 42, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 3, taxa2 = 3, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = -10, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Arch.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 30,
             geom = 'text', color = '#EA4F0D', fill = NA) + 
  geom_strip(taxa1 = 84, taxa2 = 91, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 77, taxa2 = 77, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 15, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)

Fig_contpPC2_50
#save as PDF, in size 10x10 inches

#Performing Contmap for pPC3:####
pComp3_50 = phyPCA_squama50$x[, 3]

fit_pPC3_50 <- phytools::fastAnc(pruned_tree50, pComp3_50, vars=TRUE, CI=TRUE)

td_pPC3_50 <- data.frame(node = nodeid(pruned_tree50, names(pComp3_50)),
                         trait = pComp3_50)
nd_pPC3_50 <- data.frame(node = names(fit_pPC3_50$ace), trait = fit_pPC3_50$ace)

d_pPC3_50 <- rbind(td_pPC3_50, nd_pPC3_50)
d_pPC3_50$node <- as.numeric(d_pPC3_50$node)
tree_pPC3_50 <- full_join(pruned_tree50, d_pPC3_50, by = 'node')

#Ploting: 

Fig_contpPC3_50 =  ggtree(tree_pPC3_50, layout = 'circular', ladderize = T,
                          right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                          continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 3 (7.89%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 117, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 129, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 133, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 140, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 150, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 169, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 182, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 172, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 186, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 197, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 42, taxa2 = 42, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 3, taxa2 = 3, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = -10, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Arch.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 30,
             geom = 'text', color = '#EA4F0D', fill = NA) + 
  geom_strip(taxa1 = 84, taxa2 = 91, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 77, taxa2 = 77, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 15, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)

Fig_contpPC3_50
#save as PDF, in size 10x10 inches

#Performing Contmap for pPC1:####
fit_CS50 <- phytools::fastAnc(pruned_tree50, log.CSv50, vars=TRUE, CI=TRUE)

td_CS50 <- data.frame(node = nodeid(pruned_tree50, names(log.CSv50)),
                      trait = log.CSv50)
nd_CS50 <- data.frame(node = names(fit_CS50$ace), trait = fit_CS50$ace)

d_CS50 <- rbind(td_CS50, nd_CS50)
d_CS50$node <- as.numeric(d_CS50$node)
tree_CS50 <- full_join(pruned_tree50, d_CS50, by = 'node')

#Ploting: 

Fig_contCS_50 =  ggtree(tree_CS50, layout = 'circular', ladderize = T,
                        right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                        continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'Log10(CS)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 117, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 129, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 133, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 140, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 150, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 169, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 182, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 172, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 186, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 197, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 42, taxa2 = 42, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 3, taxa2 = 3, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = -10, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Arch.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 30,
             geom = 'text', color = '#EA4F0D', fill = NA) + 
  geom_strip(taxa1 = 84, taxa2 = 91, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA) + 
  geom_strip(taxa1 = 77, taxa2 = 77, label = "Mosasauria", align=TRUE,
             barsize = 1, fontsize = 4, angle = 15, horizontal = F,
             hjust = 'center', offset.text = 10,
             geom = 'text', color = '#FE922A', fill = NA)

Fig_contCS_50
#save as PDF, in size 10x10 inches

#Ploting phylomorphospace####
plotTree(pruned_tree50,ftype="i",fsize=0.8)
nodelabels(cex=0.6,frame="circle",bg="white") #checking the node numbers

color_phygroups_nodes50 <- read_excel("Classification.xlsx", 
                                    sheet = "node_50")

nodes = color_phygroups_nodes50$Node #Nodes of target groups
clades = color_phygroups_nodes50$Group
clade_colors = color_phygroups_nodes50$Color

squama_groups50=pruned_tree50
for(i in 1:length(nodes)){
  squama_groups50<-paintSubTree(squama_groups50,node=nodes[i],
                              state=clades[i],anc.state=if(i==1) "NA" else NULL,
                              stem=T)
}

#Checking if the nodes and colors match
plot(squama_groups50, colors50 <- setNames(c('grey', clade_colors),
                                         c('NA', clades)), ftype = 'i', fsize = 0.8)
nodelabels(cex=0.6, frame = 'circle', bg = 'white')

#Plotting in phylomorphospace
taxon_names50 = data.frame(number_tip = c(1:106), names = pruned_tree50[["tip.label"]])
taxon_names50 = taxon_names50[order(taxon_names50$names),]
taxon_names50 = data.frame(taxon_names50, number = Comp_phyloPCA_squama50$Number, color = Comp_phyloPCA_squama50$Color)
taxon_names50 = taxon_names50[order(taxon_names50$number_tip),]

MphyPCA_squama12_50 = as.matrix(phyPCA_squama50$x[,1:2])
MphyPCA_squama13_50 = as.matrix(phyPCA_squama50$x[,c(1,3)])
MphyPCA_squama23_50 = as.matrix(phyPCA_squama50$x[,2:3])

phylomorphospace(squama_groups50, MphyPCA_squama12_50, colors = colors50,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 1 (15.6%)', ylab = 'pPC 2 (11.82%)', 
                 lwd = 2)
#legend('bottomright', names(colors50), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree50$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names50$number, frame = 'none', col = taxon_names50$color) #To add specimen number
#500 x 500 px

phylomorphospace(squama_groups50, MphyPCA_squama13_50, colors = colors50,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 1 (15.6%)', ylab = 'pPC 3 (7.89%)', 
                 lwd = 2)
#legend('bottomright', names(colors50), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree50$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names50$number, frame = 'none', col = taxon_names50$color) #To add specimen number
#500 x 500 px

phylomorphospace(squama_groups50, MphyPCA_squama23_50, colors = colors50,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 2 (11.82%)', ylab = 'pPC 3 (7.89%)', 
                 lwd = 2)
#legend('bottomright', names(colors50), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree50$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names50$number, frame = 'none', col = taxon_names50$color) #To add specimen number
#500 x 500 px

#Ploting phylomorphospace with convex hulls####
hull_phyComps12.50 <- Comp_phyloPCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))

ggplot(data = Comp_phyloPCA_squama50, aes(x=Comp1, y=Comp2, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps12.50, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                        linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material50, size = 3, color = 'black',fill = Comp_phyloPCA_squama50$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 1 (15.6%)', y = 'pPC 2 (11.82%)') + #choosing axes labs
  ggtitle("phyloPCA 50") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right') 

#EXPORT 800 X 500 px

hull_phyComps13.50 <- Comp_phyloPCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp3))

ggplot(data = Comp_phyloPCA_squama50, aes(x=Comp1, y=Comp3, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps13.50, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                           linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material50, size = 3, color = 'black',fill = Comp_phyloPCA_squama50$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 1 (15.6%)', y = 'pPC 3 (7.89%)') + #choosing axes labs
  ggtitle("phyloPCA 50") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')

#EXPORT 800 X 500 px

hull_phyComps23.50 <- Comp_phyloPCA_squama50 %>% group_by(Group.2) %>% 
  slice(chull(Comp2, Comp3))

ggplot(data = Comp_phyloPCA_squama50, aes(x=Comp2, y=Comp3, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps23.50, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                           linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material50, size = 3, color = 'black',fill = Comp_phyloPCA_squama50$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 2 (11.82%)', y = 'pPC 3 (7.89%)') + #choosing axes labs
  ggtitle("phyloPCA 50") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')

#EXPORT 800 X 500 px

################################################################################
###Fossil subsampling no fossil#### 
#Importing TPS file:
squama0 = readland.tps(file = 'Squama_skull_10Phylo_sub0.TPS', #no sinaeo
                        specID = "ID")

#Performing GPA:
GPA_squama0 = gpagen(squama0, 
                      curves = sliders, 
                      ProcD = T,
                      approxBE = T, Proj = T)
plot(GPA_squama0)


#Checking for outlier:
outliers_squama0=plotOutliers(GPA_squama0$coords,
                               inspect.outliers = T)

#Exploring PCA without phylogeny####
coords_aligned_squama0=two.d.array(GPA_squama0$coords)

##Grids and exploring graph
PCA_squama0=gm.prcomp(GPA_squama0$coords)
PCA_plot_squama0=plot(PCA_squama0)

BS_squama0 = broken_stick_fun(PCA_squama0)
print(BS_squama0)
BS = BS_squama0[1:30,]

#Ploting broken stick results:
ggplot(BS, aes(x=Principal_Component, y=Eigenvalue)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  geom_line(aes(x=Principal_Component, y = Expected_Length), stat = 'identity', color = 'red', group = 1) + 
  geom_point(aes(x=Principal_Component, y = Expected_Length)) + 
  labs(title = "Broken stick of 30 first PCs", x = "Principal Component", y = "Eigenvalue") + 
  scale_y_continuous("Eigenvalue", sec.axis = sec_axis(~., name = "Expected Length"))

#Exporting PCA data do .xlsx tables:
DF.coords_aligned_squama0 = as.data.frame(coords_aligned_squama0)
write.xlsx(DF.coords_aligned_squama0, file = "Coords_aligned_squama_sub0.xlsx",
           sheetName = 'coords_aligned_squama0', col.names = T, row.names = T)

DF.PCA_squama0 = as.data.frame(PCA_squama0$x)
write.xlsx(DF.PCA_squama0, file = "PCA_squama0.xlsx",
           sheetName = 'PCA_squama0', col.names = T, row.names = T)

write.xlsx(PCA_squama0$rotation, file = "loadings0.xlsx", 
           col.names = T, row.names = T)

CS0 = as.data.frame(GPA_squama0$Csize)
write.xlsx(CS0, file = "CS_squama0.xlsx",
           sheetName = 'CS0', col.names = T, row.names = T)

#PCA graphics#####
#Classifiers <- read.csv("C:/Users/Arthur/Arthur/GMM/TPS/Classifiers.csv", sep=";")
Comp_PCA_squama0 <- read_excel("Compiled data_squama.xlsx", 
                               sheet = "PCA_sub0")

hull_Comps12.0 <- Comp_PCA_squama0 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))

color_phygroups0 <- read_excel("Classification.xlsx", 
                               sheet = "phygroups0")

material0 = Comp_PCA_squama0$Material
material0 = replace(material0, material0 == 'fossil', 22)
material0 = replace(material0, material0 == 'extant', 21)
material0 = as.numeric(material0)

#Fig_PCA12.0 = 
  ggplot(data = Comp_PCA_squama0, aes(x=Comp1, y= Comp2, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps12.0, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                       linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material0, size = 3, color = 'black',fill = Comp_PCA_squama0$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups0$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups0$Color) + #choosing point colors manually
  labs(x = 'PC 1 (27.35%)', y = 'PC 2 (14.24%)') + #choosing axes labs
  ggtitle("PCA (extant)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                 force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
Fig_PCA12.0
#EXPORT 800 X 500 px

hull_Comps13.0 <- Comp_PCA_squama0 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp3))

#Fig_PCA13.0 = 
  ggplot(data = Comp_PCA_squama0, aes(x=Comp1, y= Comp3, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps13.0, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                       linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material0, size = 3, color = 'black',fill = Comp_PCA_squama0$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups0$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups0$Color) + #choosing point colors manually
  labs(x = 'PC 1 (27.35%)', y = 'PC 3 (9.2%)') + #choosing axes labs
  ggtitle("PCA (extant)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                 force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
Fig_PCA13.0
#EXPORT 800 X 500 px

hull_Comps23.0 <- Comp_PCA_squama0 %>% group_by(Group.2) %>% 
  slice(chull(Comp2, Comp3))

#Fig_PCA23.0 = 
  ggplot(data = Comp_PCA_squama0, aes(x=Comp2, y= Comp3, colour = Group.2)) + #basic arguments
  geom_polygon(data = hull_Comps23.0, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                       linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = material0, size = 3, color = 'black',fill = Comp_PCA_squama0$Color, aes(color = Group.2)) + #specifying general parameters of points
  scale_fill_manual(values = color_phygroups0$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups0$Color) + #choosing point colors manually
  labs(x = 'PC 2 (14.24%)', y = 'PC 3 (9.2%)') + #choosing axes labs
  ggtitle("PCA (extant)") + #adding title
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')
Fig_PCA23.0
#EXPORT 800 X 500 px

####including phylogenetic data####

taxa_table0 <- read_excel("Classification.xlsx", 
                           sheet = "taxa_gmm_sub0")
keep_taxa0 <- taxa_table0$name_final

#Drop the tips not included in the reduced tree or keep tak a list
#@S3 class tree
pruned_tree0 <- ape::keep.tip(tree_full, keep_taxa0)

plot.phylo(pruned_tree0)

names_phylo0 = pruned_tree0$tip.label
names_phylo0 = as.data.frame(names_phylo0)

name.check(pruned_tree0, coords_aligned_squama0)

## Exporting the tree
write.tree(pruned_tree0, file = "pruned_tree0.tre") # saves as .tre 
writeNexus(pruned_tree0, file = "pruned_tree0.nex") # save as.nex
#dput(ttrees_equal, file = "pruned_tree, equal method.txt") # for keeping the root age (good for other analyses in R because it keeps the root age)

#Including phylogenetic data to PCA:####
### PaCA Alignment of data to physlogenetic signal ###
#PaCA (OLS)
#PACA_OLS_squama0 = gm.prcomp(coords_aligned_squama0, phy = pruned_tree0,
#                            align.to.phy = T, GLS = F) ### Attention! This is a PACA, not PCA (alignment of data to phylo signal rather than variance)
#summary(PACA_OLS_squama0)
#plot(PACA_OLS_squama0, phylo = TRUE, main = "PACA (OLS)")

#PaCA (PGLS) 
#PACA_GLS_squama0 = gm.prcomp(coords_aligned_squama0, phy = pruned_tree0,
#                            align.to.phy = T, GLS = T) ### Attention! This is a PACA, not PCA (alignment of data to phylo signal rather than variance)
#summary(PACA_GLS_squama0)
#plot(PACA_GLS_squama0, phylo = TRUE, main = "PACA (PGLS)")

### Phylo PCA ###
### Phylo PCA (PGLS) with phylo corrected (transformed) projection
phyPCA_squama0 = gm.prcomp(coords_aligned_squama0, phy = pruned_tree0,
                            GLS = TRUE, transform = TRUE)
summary(phyPCA_squama0)
plot(phyPCA_squama0, phylo = TRUE, main = "phylo PCA (no fossil)")

BS0 = broken_stick_fun(phyPCA_squama0)
print(BS0)
BS = BS0[1:30,]

#Ploting broken stick results:
ggplot(BS, aes(x=Principal_Component, y=Eigenvalue)) + 
  geom_bar(stat = 'identity', width = 0.5) + 
  geom_line(aes(x=Principal_Component, y = Expected_Length), stat = 'identity', color = 'red', group = 1) + 
  geom_point(aes(x=Principal_Component, y = Expected_Length)) + 
  labs(title = "Broken stick of 30 first PCs", x = "Principal Component", y = "Eigenvalue") + 
  scale_y_continuous("Eigenvalue", sec.axis = sec_axis(~., name = "Expected Length"))

#Exporting data:
DF.phyPCA_squama0 = as.data.frame(phyPCA_squama0$x)
write.xlsx(DF.phyPCA_squama0, file = 'phyPCA_squama0.xlsx',
           sheetName = 'phyPCA_squama0', col.names = T, row.names = T)

write.xlsx(phyPCA_squama0$rotation, file = 'phyloadings0.xlsx',
           col.names = T, row.names = T)

ancestors.phyPCA_squama0 = as.data.frame(phyPCA_squama0$ancestors)
write.xlsx(ancestors.phyPCA_squama0, file = 'ancestors_phyPCA_squama0.xlsx',
           sheetName = 'ancestors_phyPCA_squama0', col.names = T, row.names = T)

Comp_phyloPCA_squama0 <- read_excel("Compiled data_squama.xlsx", 
                                     sheet = "phyloPCA_sub0")
row.names(Comp_phyloPCA_squama0) = Comp_phyloPCA_squama0$File_name

#Checking for phylogenetic signal####
#Test for phylogenetic signal in shape (Kmult)
PS_K.shape0 <- physignal(A = coords_aligned_squama0, phy = pruned_tree0)
summary(PS_K.shape0)
plot(PS_K.shape0) #sign phylo sign as expected under BM (K=1), or even higher than by BM (K>1)
plot(PS_K.shape0$PACA, phylo = TRUE)

# Phylogenetic signal profile (K stat by PC)
# important to see where the phylo signal is concentrated!
PS_K.shape0$K.by.p # First PCs much higher signal than last PCs

#Test for phylogenetic signal in shape (Z-stat)
PS_Z.shape0 <- physignal.z(A = coords_aligned_squama0, phy = pruned_tree0,
                            lambda = "front")
summary(PS_Z.shape0)

# Problem with ill-conditioned residual covariance matrix;  shaving last several PCS keeping most relevant ones (e.g., top 10)
PS_Z.shape0 <- physignal.z(A = coords_aligned_squama0, phy = pruned_tree0,
                         lambda = "front", PAC.no = 90) #99.83%
summary(PS_Z.shape0)
plot(PS_Z.shape0)
plot(PS_Z.shape0$PACA, phylo = TRUE)

#Test for phylogenetic signal in shape (Kmult)
coords.CS_squama0 <- read_excel("Compiled data_squama.xlsx", 
                                 sheet = "coords_aligned and CS_sub0")

row.names(coords.CS_squama0) = coords.CS_squama0$File_Name
log.CSv0 = coords.CS_squama0$Log.CS
names(log.CSv0) = coords.CS_squama0$File_Name 

PS_K.CS0 <- physignal(A = log.CSv0, phy = pruned_tree0)
summary(PS_K.CS0)
plot(PS_K.CS0) #sign phylo sign as expected under BM (K=1), or even higher than by BM (K>1)
plot(PS_K.CS0$PACA, phylo = TRUE)

# Phylogenetic signal profile (K stat by PC)
# important to see where the phylo signal is concentrated!
PS_K.CS0$K.by.p # First PCs much higher signal than last PCs

#Test for phylogenetic signal in shape (Z-stat)
PS_Z.CS0 <- physignal.z(A = log.CSv0, phy = pruned_tree0,
                         lambda = "front")
summary(PS_Z.CS0)

# Problem with ill-conditioned residual covariance matrix;  shaving last several PCS keeping most relevant ones (e.g., top 10)
PS_Z.CS0 <- physignal.z(A = log.CSv0, phy = pruned_tree0,
                         lambda = "front", PAC.no = 1) #99.83%
summary(PS_Z.CS0)
plot(PS_Z.CS50)
plot(PS_Z.CS50$PACA, phylo = TRUE)

#Performing contMap with PCs and CS:####
#see: http://blog.phytools.org/2021/11/showing-different-orders-or-other.html
plotTree(pruned_tree0,ftype="i",fsize=0.8)
nodelabels(cex=0.6,frame="circle",bg="white") #checking the node numbers

color_phygroups_nodes0 <- read_excel("Classification.xlsx", 
                                    sheet = "node_0")

nodes0 = color_phygroups_nodes0$Node #Nodes of target groups
clades0 = color_phygroups_nodes0$Group
clade_colors0 = color_phygroups_nodes0$Color

squama_groups0=pruned_tree0
for(i in 1:length(nodes0)){
  squama_groups0<-paintSubTree(squama_groups0,node=nodes0[i],
                              state=clades0[i],anc.state=if(i==1) "NA" else NULL,
                              stem=T)
}

#Checking if the nodes and colors match
plot(squama_groups0, colors0 <- setNames(c('grey', clade_colors0),
                                       c('NA', clades0)), ftype = 'i', fsize = 0.8)
nodelabels0(cex=0.6, frame = 'circle', bg = 'white')

#Performing Contmap for pPC1:####
pComp1_0 = phyPCA_squama0$x[, 1]

fit_pPC1_0 <- phytools::fastAnc(pruned_tree0, pComp1_0, vars=TRUE, CI=TRUE)

td_pPC1_0 <- data.frame(node = nodeid(pruned_tree0, names(pComp1_0)),
                        trait = pComp1_0)
nd_pPC1_0 <- data.frame(node = names(fit_pPC1_0$ace), trait = fit_pPC1_0$ace)

d_pPC1_0 <- rbind(td_pPC1_0, nd_pPC1_0)
d_pPC1_0$node <- as.numeric(d_pPC1_0$node)
tree_pPC1_0 <- full_join(pruned_tree0, d_pPC1_0, by = 'node')

#Ploting: 

Fig_contpPC1_0 =  ggtree(tree_pPC1_0, layout = 'circular', ladderize = F,
                         right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                         continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 1 (15.2%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 101, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 112, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 115, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 121, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 127, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 146, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 156, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 149, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 160, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 171, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 35, taxa2 = 35, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 70, taxa2 = 77, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 80, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA)


Fig_contpPC1_0
#save as PDF, in size 10x10 inches

#Performing Contmap for pPC2:####
pComp2_0 = phyPCA_squama0$x[, 2]

fit_pPC2_0 <- phytools::fastAnc(pruned_tree0, pComp2_0, vars=TRUE, CI=TRUE)

td_pPC2_0 <- data.frame(node = nodeid(pruned_tree0, names(pComp2_0)),
                        trait = pComp2_0)
nd_pPC2_0 <- data.frame(node = names(fit_pPC2_0$ace), trait = fit_pPC2_0$ace)

d_pPC2_0 <- rbind(td_pPC2_0, nd_pPC2_0)
d_pPC2_0$node <- as.numeric(d_pPC2_0$node)
tree_pPC2_0 <- full_join(pruned_tree0, d_pPC2_0, by = 'node')

#Ploting: 

Fig_contpPC2_0 =  ggtree(tree_pPC2_0, layout = 'circular', ladderize = F,
                         right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                         continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 2 (12.75%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 101, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 112, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 115, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 121, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 127, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 146, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 156, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 149, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 160, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 171, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 35, taxa2 = 35, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 70, taxa2 = 77, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 80, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA)


Fig_contpPC2_0
#save as PDF, in size 10x10 inches

#Performing Contmap for pPC3:####
pComp3_0 = phyPCA_squama0$x[, 3]

fit_pPC3_0 <- phytools::fastAnc(pruned_tree0, pComp3_0, vars=TRUE, CI=TRUE)

td_pPC3_0 <- data.frame(node = nodeid(pruned_tree0, names(pComp3_0)),
                        trait = pComp3_0)
nd_pPC3_0 <- data.frame(node = names(fit_pPC3_0$ace), trait = fit_pPC3_0$ace)

d_pPC3_0 <- rbind(td_pPC3_0, nd_pPC3_0)
d_pPC3_0$node <- as.numeric(d_pPC3_0$node)
tree_pPC3_0 <- full_join(pruned_tree0, d_pPC3_0, by = 'node')

#Ploting: 

Fig_contpPC3_0 =  ggtree(tree_pPC3_0, layout = 'circular', ladderize = F,
                         right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                         continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'pPC 3 (7.82%)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 101, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 112, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 115, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 121, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 127, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 146, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 156, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 149, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 160, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 171, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 35, taxa2 = 35, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 70, taxa2 = 77, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 80, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA)


Fig_contpPC3_0
#save as PDF, in size 10x10 inches

#Ploting CS:

fit_CS0 <- phytools::fastAnc(pruned_tree0, log.CSv0, vars=TRUE, CI=TRUE)

td_CS0 <- data.frame(node = nodeid(pruned_tree0, names(log.CSv0)),
                     trait = log.CSv0)
nd_CS0 <- data.frame(node = names(fit_CS0$ace), trait = fit_CS0$ace)

d_CS0 <- rbind(td_CS0, nd_CS0)
d_CS0$node <- as.numeric(d_CS0$node)
tree_CS0 <- full_join(pruned_tree0, d_CS0, by = 'node')
#Performing Contmap for CS#####

#Ploting: 

Fig_contCS_0 =  ggtree(tree_CS0, layout = 'circular', ladderize = F,
                       right = T, aes(color=trait), branch.length = 'branch.length', size = 1.5,
                       continuous = 'colour') +
  scale_color_gradientn(colours=c("#440154", "#31688e", "#35b779", "#fde725"), name = 'Log10(CS)') + 
  theme_tree() +
  theme(legend.position = c(x = .07, y = .2),
        legend.title = element_text(size = 10, face = 'bold', hjust = 0.5),
        legend.text=element_text(size=12), 
        legend.key.height=unit(.4, "cm"),
        legend.key.width=unit(.4, "cm")) +
  geom_treescale(fontsize=5, linesize=1.5, offset=2, width = 30, offset.label = 2) +
  geom_cladelab(node = 101, label = 'Pleurodonta', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#558835', barcolor = '#558835') + 
  geom_cladelab(node = 112, label = 'Acrod.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#56AA21', barcolor = '#56AA21') + 
  geom_cladelab(node = 115, label = 'Neoang.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#0B6345', barcolor = '#0B6345') + 
  geom_cladelab(node = 121, label = 'Paleo.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#45AA88', barcolor = '#45AA88') + 
  geom_cladelab(node = 127, label = 'Caenophidia', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#2B386B', barcolor = '#2B386B') + 
  geom_cladelab(node = 146, label = 'Scolec.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#7F8DC5', barcolor = '#7F8DC5') + 
  geom_cladelab(node = 156, label = 'Amphisb.', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#8C373E', barcolor = '#8C373E') + 
  geom_cladelab(node = 149, label = 'Teiioidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#F1CA3A', barcolor = '#F1CA3A') +
  geom_cladelab(node = 160, label = 'Scincoidea', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#251EDE', barcolor = '#251EDE') + 
  geom_cladelab(node = 171, label = 'Gekkota', align = T,
                barsize = 1, fontsize = 4, angle = 'auto', horizontal = F,
                hjust = 'center', offset.text = 15,
                geom = 'text', color = 'black', fill = NA,
                textcolor = '#9E8925', barcolor = '#9E8925') + 
  geom_strip(taxa1 = 35, taxa2 = 35, label = "Lacert.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 70, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#F1CA3A', fill = NA) + 
  geom_strip(taxa1 = 1, taxa2 = 1, label = "Sphen.", align=TRUE,
             barsize = 1, fontsize = 4, angle = 0, horizontal = F,
             hjust = 'center', offset.text = 40,
             geom = 'text', color = '#715128', fill = NA) + 
  geom_strip(taxa1 = 70, taxa2 = 77, label = "'Haenop.'", align=TRUE,
             barsize = 1, fontsize = 4, angle = 80, horizontal = F,
             hjust = 'center', offset.text = 20,
             geom = 'text', color = '#5566AB', fill = NA)


Fig_contCS_0
#save as PDF, in size 10x10 inches

#Ploting phylomorphospace####

#Plotting in phylomorphospace
taxon_names0 = data.frame(number_tip = c(1:92), names = pruned_tree0[["tip.label"]])
taxon_names0 = taxon_names0[order(taxon_names0$names),]
taxon_names0 = data.frame(taxon_names0, number = Comp_phyloPCA_squama0$Number, color = Comp_phyloPCA_squama0$Color)
taxon_names0 = taxon_names0[order(taxon_names0$number_tip),]

MphyPCA_squama12_0 = as.matrix(phyPCA_squama0$x[,1:2])
MphyPCA_squama13_0 = as.matrix(phyPCA_squama0$x[,c(1,3)])
MphyPCA_squama13_0[,2] = -MphyPCA_squama13_0[,2]
MphyPCA_squama23_0 = as.matrix(phyPCA_squama0$x[,2:3])
MphyPCA_squama23_0[,2] = -MphyPCA_squama23_0[,2]

phylomorphospace(squama_groups0, MphyPCA_squama12_0, colors = colors0,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 1 (15.2%)', ylab = 'pPC 2 (12.75%)', 
                 lwd = 2)
#legend('bottomright', names(colors0), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree0$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names0$number, frame = 'none', col = taxon_names0$color) #To add specimen number
#500 x 500 px

phylomorphospace(squama_groups0, MphyPCA_squama13_0, colors = colors0,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 1 (15.2%)', ylab = 'pPC 3 (7.82%)', 
                 lwd = 2)
#legend('bottomright', names(colors0), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree0$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names0$number, frame = 'none', col = taxon_names0$color) #To add specimen number
#500 x 500 px

phylomorphospace(squama_groups0, MphyPCA_squama23_0, colors = colors0,
                 ftype = 'off', bty = 'l', node.size=c(1,1.5), node.by.map=T, pch = NULL, 
                 xlab = 'pPC 2 (12.75%)', ylab = 'pPC 3 (7.82%)', 
                 lwd = 2)
#legend('bottomright', names(colors0), pch = 21,
#       pt.bg = colors, pt.cex = 1.5, bg = NULL, box.lwd = NULL)
#title(main = "Phylomorphospace of phyloPCA")
nodelabels(text = pruned_tree0$node.label, frame = 'none', col = 'black')
tiplabels(text = taxon_names0$number, frame = 'none', col = taxon_names0$color) #To add specimen number
#500 x 500 px

#Ploting phylomorphospace with convex hulls####
hull_phyComps12.0 <- Comp_phyloPCA_squama0 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp2))

ggplot(data = Comp_phyloPCA_squama0, aes(x=Comp1, y=Comp2, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups0$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups0$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps12.0, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                           linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = 21, size = 3, color = 'black',fill = Comp_phyloPCA_squama0$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 1 (15.2%)', y = 'pPC 2 (12.75%)') + #choosing axes labs
  ggtitle("phyloPCA 0") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right') 

#EXPORT 800 X 500 px

hull_phyComps13.0 <- Comp_phyloPCA_squama0 %>% group_by(Group.2) %>% 
  slice(chull(Comp1, Comp3))

ggplot(data = Comp_phyloPCA_squama0, aes(x=Comp1, y=-Comp3, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups0$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups0$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps13.0, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                          linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = 21, size = 3, color = 'black',fill = Comp_phyloPCA_squama0$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 1 (15.2%)', y = 'pPC 3 (7.82%)') + #choosing axes labs
  ggtitle("phyloPCA 0") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right') 

#EXPORT 800 X 500 px

hull_phyComps23.0 <- Comp_phyloPCA_squama0 %>% group_by(Group.2) %>% 
  slice(chull(Comp2, Comp3))

ggplot(data = Comp_phyloPCA_squama0, aes(x=Comp2, y=-Comp3, colour = Group.2)) + #basic arguments
  scale_fill_manual(values = color_phygroups0$Color) + #filling inside convex hulls manually
  scale_color_manual(values = color_phygroups0$Color) + #choosing point colors manually
  geom_polygon(data = hull_phyComps23.0, alpha = 0.4, aes(color = Group.2, fill = Group.2, 
                                                          linewidth = NULL, linetype = NULL)) + #adding convex hulls
  geom_point(shape = 21, size = 3, color = 'black',fill = Comp_phyloPCA_squama0$Color, aes(color = Group.2)) + #specifying general parameters of points
  geom_text_repel(min.segment.length = 0.5, seed = 50, max.overlaps = 20, #to avoid overlap of the numbers
                  force = 15, force_pull = 5, aes(label = Number)) + 
  labs(x = 'pPC 2 (12.75%)', y = 'pPC 3 (7.82%)') + #choosing axes labs
  ggtitle("phyloPCA 0") + #adding title
  theme_classic() + #defining theme
  theme(axis.text = element_text(size = 10, color = 'black'), #defining txt elements
        axis.line = element_line(size = 0.5), #axis line width
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.position = 'right')

#EXPORT 800 X 500 px
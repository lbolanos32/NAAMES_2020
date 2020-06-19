#NAAMESNTWK.R contains the commands neccesary to reproduce the network analysis performed on the data of the four NAAMES cruises. 

#Authors: Luis M. Bolaños and Stephen Giovannoni

library("phyloseq")
library(ggnet)
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("dplyr")
library("phangorn")
library("data.table")
library("cowplot")
library(grid)
library(gridExtra) 
library("ggdendro")
library(igraph)
library("microbiome") # data analysis and visualisation
library("microbiomeutilities") # some utility tools
library("ggpubr") # publication quality figures, based on ggplot2
library("DT") # interactive tables in html and markdown
library(SpiecEasi)
library(intergraph)
library(network)
library("qgraph")
library("sna")

setwd("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/GithubDocs/Tables/")

#Photosynthetic fraction
count_tab <- read.table("Photlink_final_collapsedVR2.otu", header=T, row.names=1, check.names=F)
sample_info_tab <- read.table("Merge_envFileDNA.txt", header=T, row.names=1, check.names=F, sep ="\t")
tax_tab <- as.matrix(read.table("Photlink_final_collapsedVR.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

#List to remove
contam_asvs <-c("N287","N370","N504","N932","N1195","N1575","N2862","N2934","N3307","N3544","N8532","N24639","N25035","N25138","N25145","N25317","N25556","N26238")

####Prune them from the main Phyloseq object

asv_tab_no_contam <- count_tab[!row.names(count_tab) %in% contam_asvs, ]# making new OTU table

  
asv_tax_no_contam <- tax_tab[!row.names(tax_tab) %in% contam_asvs, ] # making new taxonomy table

sample_info_tab$depth <- factor(sample_info_tab$depth,levels = c("5", "25", "50", "75", "100","150","200","300"))

OTU = otu_table(asv_tab_no_contam, taxa_are_rows = TRUE)
TAX = tax_table(asv_tax_no_contam)
SAM = sample_data(sample_info_tab)


physeqphot = phyloseq(OTU,TAX,SAM)

euph = get_variable(physeqphot, "depth") %in% c("5", "25", "50", "75", "100")
sample_data(physeqphot)$euph <- factor(euph)
phyeuph<-subset_samples(physeqphot, euph %in% TRUE)

st = get_variable(phyeuph, "Type") %in% "Cast"
sample_data(phyeuph)$st <- factor(st)

physt<-subset_samples(phyeuph, st %in% TRUE)
phystV1filt= filter_taxa(physt, function(x) sum(x > 2) > (0.00075*length(x)), TRUE)

phyeuphminN1N2 = prune_samples(sample_sums(phystV1filt) > 1600,  phystV1filt)

#Heterotrophic fraction

count_tabH <- read.table("Hetlink_final_collapsedV2.otu", header=T, row.names=1, check.names=F)
sample_info_tabH <- read.table("Merge_envFileDNA.txt", header=T, row.names=1, check.names=F, sep ="\t")
tax_tabH <- as.matrix(read.table("HetlinkTax_202_11.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

asv_tab_no_contam_het <- count_tabH [!row.names(count_tabH ) %in% contam_asvs, ]# making new OTU table

  
asv_tax_no_contam_het <- tax_tabH[!row.names(tax_tabH) %in% contam_asvs, ] # making new taxonomy table


sample_info_tabH$depth <- factor(sample_info_tabH$depth,levels = c("5", "25", "50", "75", "100","150","200","300"))

OTU_H = otu_table(asv_tab_no_contam_het , taxa_are_rows = TRUE)
TAX_H = tax_table(asv_tax_no_contam_het)
SAM_H = sample_data(sample_info_tabH)

physeqhet = phyloseq(OTU_H,TAX_H,SAM_H)

het = get_variable(physeqhet, "depth") %in% c("5", "25", "50", "75", "100")
sample_data(physeqhet)$het <- factor(het)
phyeuphet<-subset_samples(physeqhet, het %in% TRUE)

sth = get_variable(phyeuphet, "Type") %in% "Cast"

sample_data(phyeuphet)$sth <- factor(sth)

physth<-subset_samples(phyeuphet, sth %in% TRUE)
phystV1filth1<-filter_taxa(physth, function(x) sum(x > 1) > (0.000075*length(x)), TRUE)
phystV1filth <- prune_samples(sample_sums(phystV1filth1)>15000, phystV1filth1)

#####Remove non-abundant taxa to 

phystV1filth2<-filter_taxa(phystV1filth, function(x) sum(x > 6) > (0.4*length(x)), TRUE) #Remove taxa not seen more than 6 times in at least 40% of the samples
phystV1filt2<-filter_taxa(phystV1filt, function(x) sum(x > 6) > (0.1*length(x)), TRUE) #Remove taxa not seen more than 6 times in at least 10% of the samples


merged_All<-merge_phyloseq(phystV1filt2,phystV1filth2)

merged_All10k <- prune_samples(sample_sums(merged_All)>10000, merged_All)

mycolorsphl <- scale_color_manual(values = c("Chlorophyta"="forestgreen","Cyanobacteria"="brown","Stramenopiles"="blue","Cryptophyta"="yellow3","Haptophyta"="darkgoldenrod3","Rhodophyta"="chocolate3","Proteobacteria"="lightcoral","Chloroflexi"="lightskyblue","SAR406_clade"="firebrick1","Bacteroidetes"="#7FC97F","Actinobacteria"="purple","Rappemonad"="gold4","new_euk_C"="maroon1","Verrucomicrobia"="deeppink1","Eusiphoniidae"="black","k__Unclassified"="gray34","k__Bacteria"="cornsilk4", "Gemmatimonadetes"="darkslateblue", "Planctomycetes"="slategray4"))

######NAAMES_ALL: Merge Photo and Het to start the analysis

merged_All10k.f1<- microbiomeutilities::format_to_besthit(merged_All10k) #change to  besthit
otu.merged_All10k <- t(otu_table(merged_All10k.f1)@.Data) #transform 
tax.merged_All10k <- as.data.frame(tax_table(merged_All10k.f1)@.Data) #extract the taxa table from phyloseq object

#set.seed(1777)
#net.N1 <- spiec.easi(otu.merged_All10k, method='mb', icov.select.params=list(rep.num=50)) #run spiec-easi covariance method. 

#SPIEC.EASI covariance analysis is a computing exhaustive process. For the first time, user can run. We save the .rds of this analysis and load it. (Provided in GitHub)


net.N1 <- readRDS("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/analysis/netNAAMESall.rds")

merged_All10k.c <- symBeta(getOptBeta(net.N1)) #output the covariance matrix

colnames(merged_All10k.c) <- rownames(merged_All10k.c) <- colnames(otu.merged_All10k) # add to the matrix sample names and taxa names

vsizeall <- log2(apply(otu.merged_All10k, 2, mean)) #estimate the log2 of the mean of each ASV so it can be the size of the node 

region.ig.merged_All10k <- graph.adjacency(merged_All10k.c, mode='undirected', add.rownames = TRUE, weighted = TRUE) #Create an igraph object 

coords.N1fdr = layout_with_fr(region.ig.merged_All10k) #Create the coordinates of the igraph 
E(region.ig.merged_All10k)[weight > 0]$color<-"steelblue" #set positivevalues to orange color for the edge
E(region.ig.merged_All10k)[weight < 0]$color<-"orange" #set negative values to orange color for the edge

merged_All10k.net <- asNetwork(region.ig.merged_All10k) #create a network object from the igraph 

network::set.edge.attribute(merged_All10k.net, "color", ifelse(merged_All10k.net %e% "weight" > 0, "steelblue", "orange")) #set the colors again for the network
phylaall <- map_levels(colnames(otu.merged_All10k), from = "best_hit", to = "Phylum", tax_table(merged_All10k.f1)) #extract phylum to be added  
merged_All10k.net %v% "Phylum" <- phylaall #add the phyla information 
merged_All10k.net %v% "nodesize" <- vsizeall #add the node information 

#P will contain the network plot with positive and negative associations. For modularity we will consider a positive and negative assocation as equal contributors to the module not only positives. 

p <- ggnet2(merged_All10k.net, legend.size = 12, node.color = "Phylum", 
            label = FALSE, node.size = "nodesize", 
            label.size = 2, edge.color = "color",alpha = 0.75) + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolorsphl + ggtitle("All samples from 5 to 100m")+theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

#MODULES#

#Negative associations turn them into absolute values 

merged_All10k.abs<-abs(merged_All10k.c)
all10k.ig.abs <- graph.adjacency(merged_All10k.abs, mode='undirected', add.rownames = TRUE, weighted = TRUE)

modules =cluster_fast_greedy(all10k.ig.abs)

V(all10k.ig.abs)$module=modules$membership

#IGRAPH clustering fast greedy, groups: 14, mod: 0.63

merged_All10k.abs.net <- asNetwork(all10k.ig.abs)
merged_All10k.abs.net %v% "Phylum" <- phylaall #add the phyla information
merged_All10k.abs.net %v% "nodesize" <- vsizeall  #add the phyla information 

x = gplot.layout.fruchtermanreingold(merged_All10k.abs.net, NULL)
merged_All10k.abs.net %v% "x" = x[, 1]
merged_All10k.abs.net %v% "y" = x[, 2]

#IGRAPH clustering fast greedy, groups: 13, mod: 0.57

a<-ggnet2(merged_All10k.abs.net, legend.size = 12, node.color = "Phylum",mode = c("x", "y"), 
            label = FALSE, node.size = "nodesize", 
            label.size = 2,alpha = 0.85) + guides(color=guide_legend(title="Phylum"), size = FALSE) + mycolorsphl +theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) + theme(legend.position = "none")

b<-ggnet2(merged_All10k.abs.net, legend.size = 12, node.color = "module",mode = c("x", "y"), 
            label = FALSE, node.size = "nodesize", 
            label.size = 2,alpha = 0.85) + guides(color=guide_legend(title="Module"), size = FALSE) + theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))+scale_color_manual(values=c(c("purple", "coral4", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red", "goldenrod", "black", "orange", "firebrick", "cyan"))) + theme(legend.position = "none")

#b will contain the same plot as p, but nodes colored by module

######PLOT B and A together SUPPLEMENTARY FIGURE #######

#No titles or legends
#svg("ntwg_woleg.svg", width=12)
plot_grid(b,a, nrow=1,labels = c('a', 'b'))
dev.off()

#with legends 
#svg("ntwg_leg.svg", width=12)
plot_grid(b,a, nrow=1,labels = c('a', 'b'))
dev.off()


######NETWORK ANALYSIS#######

#USE igraph all10k.ig.abs
deg <- igraph::degree(all10k.ig.abs, mode = "all")
# Calculate  degree distribution
deg.dist <- degree_distribution(all10k.ig.abs, mode = "all", cumulative = F)

# plot degree distribution

#the degree of a node in a network is the number of connections it has to other nodes and the degree distribution is the probability distribution of these degrees over the whole network

#La centralidad de los nodos en la red se mide de dos formas: closeness y betweenness

clos <- igraph::closeness(all10k.ig.abs, mode = "all")
betw <- igraph::betweenness(all10k.ig.abs, v = V(all10k.ig.abs))

centralityPlot(all10k.ig.abs, include = c("Betweenness", "Closeness", "Degree")) + 
  theme(axis.text.y = element_blank())

#Keystone spp -> Node Degree + Node Centrality

# Sort nodes from greater to lesser degree
deg_sort <- sort(deg, decreasing = TRUE)
#head(deg_sort)

#Plot an histogram to visualize the degree values of the nodes  
# Transform 'deg_sort' to a 'data frame'
deg_df <- as.data.frame(deg_sort)
# Graficar histograma
# También usamos la función `geom_vline()` del paquete ggplot2 para dibujar una línea discontinua para marcar la media o promedio de los valores de degree

# Estimate betweenness
bn <- igraph::betweenness(all10k.ig.abs)
# Ordenar nodos según betweenness de mayor a menor
bn_sort <- sort(bn, decreasing = TRUE)

#Plot an histograma to visualize  betweenness values of the nodes 
# Transform 'bn_sort' to 'data frame'
bn_df <- as.data.frame(bn_sort)
# Plot histogram
# También usamos la función `geom_vline()` del paquete ggplot2 para dibujar una línea discontinua para marcar la media o promedio de los valores de betweenness

#Personal note: Un monton de mis nodos ~80 tienen cero

#Degree > 15

#Betweenness > 1250

deg_df$TaxID <- row.names(deg_df) 
head(deg_df)

#HACER UNA TABLA QUE INCLUYA MAS informacion para poderla graficar



####Submodule STATS#####

#We want number of nodes, clustering coefficient, modularity index and average degree en este debe de estar la pertenencia de los modules

#merged_All10k.abs.net#

#COLORS for all 

color_mod<-c("Chlorophyta"="forestgreen","Cyanobacteria"="brown","Stramenopiles"="blue","Cryptophyta"="yellow3","Haptophyta"="darkgoldenrod3","Rhodophyta"="chocolate3","Proteobacteria"="lightcoral","Chloroflexi"="lightskyblue","SAR406_clade"="firebrick1","Bacteroidetes"="#7FC97F","Actinobacteria"="purple","Rappemonad"="gold4","new_euk_C"="maroon1","Verrucomicrobia"="deeppink1","Eusiphoniidae"="black","k__Unclassified"="gray34","k__Bacteria"="cornsilk4", "Gemmatimonadetes"="darkslateblue", "Planctomycetes"="slategray4")



###EJEMPLO con uno####

m1<-V(all10k.ig.abs)[modules$membership == 1]
m1_subnet <- induced_subgraph(all10k.ig.abs, m1)
transitivity(m1_subnet, type =  "global")  #Global, in this case module, clustering coefficient
unique(ave(igraph::degree(m1_subnet, mode = "all"))) #Average degree
length(igraph::degree(m1_subnet, mode = "all")) #Number of nodes

#Loop to get all this data and transferred to the excel sup 1 

#Number of phytoplankton 
mod_phy_table #tiene todos por phylum de ahi contar Chlorophyta Cryptophyta Cyanobacteria Haptophyta Rhodophyta Stramenopiles

rowSums(mod_phy_table[,c(4,5,6,8,11,13)])

#Allmod1 Allmod10 Allmod11 Allmod12 Allmod13 Allmod14  Allmod2  Allmod3  Allmod4  Allmod5  Allmod6  Allmod7  Allmod8  Allmod9 
#      0        1        3        5        0        0       24        1       35        0        2       33        1       23

####Submodule PLOTS from the region.ig.merged_All10k##### Individual plots and easthatetics later 

#stablish the color palette for the nodes
mycolorsphl_ntw <- c("Chlorophyta"="forestgreen","Cyanobacteria"="brown","Stramenopiles"="blue","Cryptophyta"="yellow3","Haptophyta"="darkgoldenrod3","Rhodophyta"="chocolate3","Proteobacteria"="lightcoral","Chloroflexi"="lightskyblue","SAR406_clade"="firebrick1","Bacteroidetes"="#7FC97F","Actinobacteria"="purple","Rappemonad"="gold4","new_euk_C"="maroon1","Verrucomicrobia"="deeppink1","Eusiphoniidae"="black","k__Unclassified"="gray34","k__Bacteria"="cornsilk4", "Gemmatimonadetes"="darkslateblue", "Planctomycetes"="slategray4")

V(region.ig.merged_All10k)$module=modules$membership

m1 <- V(region.ig.merged_All10k)[modules$membership==1]
m1_subnet <- induced_subgraph(region.ig.merged_All10k, m1)

m1_names<-V(m1_subnet)$name
m1_taxa<-tax.merged_All10k[tax.merged_All10k$best_hit %in% m1_names,]

# Creamos el objeto de clase network (input necesario para ggnet2)

m1_class <- as_adjacency_matrix(m1_subnet, type = "both")# Extraer matriz de adyacencia
m1_class <- network(as.matrix(m1_class),  matrix.type = "adjacency", directed = F) # Generar objeto de clase 'network'
m1_ntwk<-ggnet2(m1_class, mode = "circle", color = m1_taxa$Phylum, palette =mycolorsphl_ntw, label = TRUE, label.size = 3,edge.color = "black", size = "degree")

#svg("mod1.svg", width=9)
m1_ntwk
#dev.off()

##BARPLOTS FOR THE MODULES GRAPH######

otu_modules<-as.data.frame(otu_table(MeanModules))
tax_modules<-as.data.frame(tax_table(MeanModules)[,1])
rownames(otu_modules)<-tax_modules$module

otu_modules_t<-(t(otu_modules))

md_to_add<-as.data.frame(sample_data(MeanModules))[,c(1,2,3,7,12,16,21,34,40,43,45,47,79,51,52,84,99,102,127)] #Cruise Station depth Temperature name

md_to_add$name<-row.names(md_to_add)

modules_array<-cbind(otu_modules_t, md_to_add)

mod_melt<-melt(modules_array,id.vars=c("Cruise","position","Station","Latitude","depth","Subregion","name"), measure.vars = c("Allmod6","Allmod2","Allmod12","Allmod9","Allmod7","Allmod1","Allmod11","Allmod3","Allmod5","Allmod10","Allmod8","Allmod14","Allmod4","Allmod13"))

mod_melt$variable<-factor(mod_melt$variable, levels = c("Allmod1","Allmod2","Allmod3","Allmod4","Allmod5","Allmod6","Allmod7","Allmod8","Allmod9","Allmod10","Allmod11","Allmod12","Allmod13","Allmod14"))

#svg("barplot_modules_allrel.svg", width=12)
ggplot(mod_melt, aes(x = position, y = value, fill = variable)) + geom_bar(stat = "identity",width=.85)+ scale_fill_manual(values=c("mediumpurple1", "coral4", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "red", "goldenrod", "black", "orange", "firebrick", "cyan")) + theme_bw()+ ylab("Relative contribution [%]") +theme(strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=12)) +facet_grid(~Cruise,scales = "free_x",space = "free_x")

#All_modules

#svg("hm_modules.svg", width=12)
ggplot(mod_melt, aes(x = position, y = variable)) + geom_tile(aes(fill = value),colour = "white")+scale_fill_gradient(low = "aliceblue", high = "firebrick4")+theme(strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=12))+ ylab("Module") +facet_grid(~Cruise,scales = "free_x",space = "free_x")


#svg("hm_modules_depth.svg", width=8)
ggplot(mod_melt, aes(x = depth, y = variable)) + geom_tile(aes(fill = value),colour = "white")+scale_fill_gradient(low = "aliceblue", high = "firebrick4")+theme(strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=12))+ ylab("Module") +facet_grid(~depth,scales = "free_x",space = "free_x")

#svg("hm_modules_region.svg", width=8)
ggplot(mod_melt, aes(x = Subregion, y = variable)) + geom_tile(aes(fill = value),colour = "white")+scale_fill_gradient(low = "aliceblue", high = "firebrick4")+theme(strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=12))+ ylab("Module") +facet_grid(~Subregion,scales = "free_x",space = "free_x")
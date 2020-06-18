
#OrdinationNAAMES.R script to generate the Figure 2 of the Seasonality of the microbial community composition in the western North Atlantic manuscript 

#Authors: Luis M. Bolaños and Stephen Giovannoni

#Load libraries
library("phyloseq")
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

setwd("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/GithubDocs/Tables/")

####CREATE A PHOTOSYNTETIC phyloseq object####

count_tab <- read.table("Photlink_final_collapsedVR2.otu", header=T, row.names=1, check.names=F)
sample_info_tab <- read.table("/Users/luisbolanos/Documents/OSU_postdoc/NAAMES/MergeAll/FINAL/Merge_envFileDNA.txt", header=T, row.names=1, check.names=F, sep ="\t")
tax_tab <- as.matrix(read.table("Taxa_custom.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

sample_info_tab$depth <- factor(sample_info_tab$depth,levels = c("5", "25", "50", "75", "100","150","200","300"))

OTU = otu_table(count_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
SAM = sample_data(sample_info_tab)

physeqphot = phyloseq(OTU,TAX,SAM)
sample_data(physeqphot)$id <- paste(sample_data(physeqphot)$Cruise, sample_data(physeqphot)$Station,sample_data(physeqphot)$depth)

euph = get_variable(physeqphot, "depth") %in% c("5", "25", "50", "75", "100")
sample_data(physeqphot)$euph <- factor(euph)
phyeuph<-subset_samples(physeqphot, euph %in% TRUE)


st = get_variable(phyeuph, "Type") %in% "Cast"
sample_data(phyeuph)$st <- factor(st)

physt<-subset_samples(phyeuph, st %in% TRUE)

##Get total reads for photo from this physt
##a<-colSums(otu_table(physt))
##write.table(a,file= "totalPhytoreads", quote=FALSE, sep = "\t")

phystV1filt= filter_taxa(physt, function(x) sum(x > 2) > (0.00075*length(x)), TRUE) ###### Reducing table from 2815 to 1369 ASVs

phyeuphminN1N2 = prune_samples(sample_sums(phystV1filt) > 1600,  phystV1filt)

bray_not_na <- phyloseq::distance(physeq = phyeuphminN1N2, method = "bray")
cap_ord <- ordinate(physeq = phyeuphminN1N2, method = "MDS", distance = bray_not_na)

###PLOT the photosynthetic ordination and save it on "p"###

p <-  plot_ordination(phyeuphminN1N2, cap_ord , color="Subregion", shape="depth")+aes(shape = depth) + scale_shape_manual(values=c(15,16,17,18,1,5,6,7)) + geom_point(aes(colour = Subregion), alpha = 0.9, size = 4) +geom_text(aes(label=Station),hjust=-.5, vjust=-.5)+
 geom_point(colour = "grey90", size = 1.5)+
 scale_color_manual(values =c("red","blue","chocolate","cyan","#a65628","#ffae19", "red",  "#1919ff", "darkorchid3", "magenta","black","cyan")) +theme(axis.text.x=element_text(angle=0,size=12), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=12))+theme_bw()+facet_wrap(~Cruise, scales="free")+theme_minimal() + theme(axis.line=element_line(),axis.text = element_text( size = 12,face = "bold" ),axis.title = element_text( size = 14, face = "bold" ),strip.text = element_text(size=16),legend.text=element_text(size=12)) + scale_x_continuous(limits=c(-.6,.6)) + scale_y_continuous(limits=c(-.6,.6))


####CREATE A Heterotrophic phyloseq object####

count_tabH <- read.table("Hetlink_final_collapsedV2.otu", header=T, row.names=1, check.names=F)
sample_info_tabH <- read.table("Merge_envFileDNA.txt", header=T, row.names=1, check.names=F, sep ="\t")
tax_tabH <- as.matrix(read.table("HetlinkTax_202_11.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

sample_info_tabH$depth <- factor(sample_info_tabH$depth,levels = c("5", "25", "50", "75", "100","150","200","300"))

OTU_H = otu_table(count_tabH, taxa_are_rows = TRUE)
TAX_H = tax_table(tax_tabH)
SAM_H = sample_data(sample_info_tabH)

physeqhet = phyloseq(OTU_H,TAX_H,SAM_H)

het = get_variable(physeqhet, "depth") %in% c("5", "25", "50", "75", "100","150","200","300")
sample_data(physeqhet)$het <- factor(het)
phyeuphet<-subset_samples(physeqhet, het %in% TRUE)

sth = get_variable(phyeuphet, "Type") %in% "Cast"

sample_data(phyeuphet)$sth <- factor(sth)

physth<-subset_samples(phyeuphet, sth %in% TRUE)

phystV1filth1<-filter_taxa(physth, function(x) sum(x > 1) > (0.000075*length(x)), TRUE)

#Remove the three samples below 15000 reads after the above QC

phystV1filth <- prune_samples(sample_sums(phystV1filth1)>15000, phystV1filth1)

bray_not_nah <- phyloseq::distance(physeq = phystV1filth, method = "bray")
cap_ordh <- ordinate(physeq = phystV1filth, method = "MDS", distance = bray_not_nah)

d <-  plot_ordination(phystV1filth, cap_ordh , color="Subregion", shape="depth")+aes(shape = depth) + scale_shape_manual(values=c(15,16,17,18,1,5,6,7)) + geom_point(aes(colour = Subregion), alpha = 0.9, size = 4) +
geom_point(colour = "grey90", size = 1.5)+ scale_color_manual(values =c("red","blue","chocolate","cyan","#a65628","#ffae19", "red",  "#1919ff", "darkorchid3", "magenta","black","cyan")) +theme(axis.text.x=element_text(angle=0,size=8), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=6))+theme_bw()+facet_wrap(~Cruise, scales="free")+theme_minimal() + theme(axis.line=element_line(),axis.text = element_text( size = 12,face = "bold" ),axis.title = element_text( size = 14, face = "bold" ),strip.text = element_text(size=16),legend.text=element_text(size=12)) + scale_x_continuous(limits=c(-.5,.45)) + scale_y_continuous(limits=c(-.4,.45))

#Plot both ordinations in the same file##
plot_grid(p,d, nrow=1)
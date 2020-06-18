#HetCCNAAMES.R script to generate the Figure 4 Heterotrophic Bacteria in the Community Composition of the Seasonality of the microbial community composition in the western North Atlantic manuscript 

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

#We want to estimate the heterotrophic 

setwd("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/GithubDocs/Tables/")

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

##Get total reads for het from this physt
##b<-colSums(otu_table(physth))
##write.table(b,file= "totalHetreads", quote=FALSE, sep = "\t") #No hice ningún tipo de filtrado para poner el total de reads >S

phystV1filth1<-filter_taxa(physth, function(x) sum(x > 1) > (0.000075*length(x)), TRUE)

#Remove the three samples below 15000 reads after the above QC

phystV1filth <- prune_samples(sample_sums(phystV1filth1)>15000, phystV1filth1)

glomV1filth<-tax_glom(phystV1filth, taxrank="Order")

glomV1filthperc<- transform_sample_counts(glomV1filth, function(x) x/sum(x))

A1<-otu_table(glomV1filthperc) # sacar la tabla para manipularla 
B1<-tax_table(glomV1filthperc)
C1<- sample_data(glomV1filthperc)

dim(A1)
[1] 196 245

major_taxa_proportions_tab_for_plot <- A1

major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(major_taxa_proportions_tab_for_plot)

A1<-data.frame(A1)
temp_filt_major_taxa_proportions_tab <- data.frame(A1[apply(A1, 1, max) > .01, ])

Nlessthan1 <- 1-colSums(temp_filt_major_taxa_proportions_tab)
df<-rbind(temp_filt_major_taxa_proportions_tab, "N1k"=Nlessthan1)

Tax43<-row.names(temp_filt_major_taxa_proportions_tab) 
prune_taxa(Tax43, glomV1filthperc)
phy43<-prune_taxa(Tax43, glomV1filthperc)
Tax43_table<-as.table(tax_table(phy43))

taxa<- colnames(Tax43_table)

Tax44<-(t(data.frame(N1k= c("Order < 1% abund", "Order < 1% abund", "Order < 1% abund", "Order < 1% abund","Order < 1% abund" ,"NA"))))

Tax44_table<-rbind(Tax43_table,Tax44)

OTU_H2 = otu_table(df, taxa_are_rows = TRUE)
TAX_H2 = tax_table(Tax44_table)
rownames(C1) <- gsub("-", ".", rownames(C1))

phy_sum = phyloseq(OTU_H2,TAX_H2,C1)

#Hetortrophic fraction

y4 <- psmelt(phy_sum)
y4$Order <- as.character(y4$Order)#convert to character
y4$Order[y4$Abundance < 0.01] <- "Order < 1% abund" #rename Order with < 1% abundance
y4$Order[y4$Order == "Order_Incertae_Sedis"]<-"Other"
names_ordenados<-c("SAR202_clade1","SAR202_clade2","SAR202_clade3","SAR202_clade4","SAR202_clade5","SAR202_clade","Flavobacteriales","Cytophagales","Sphingobacteriales","Gemmatimonadetes:BD2-11","Nitrospirales","Planctomycetales","OM190","SAR11_clade","Rhodobacterales","Rhodospirillales","Rickettsiales","OCS116_clade","Rhizobiales","Methylophilales","Nitrosomonadales","SAR324_clade(Marine_group_B)","Desulfobacterales","Bdellovibrionales","Oceanospirillales","Cellvibrionales","Thiotrichales","Salinisphaerales","Alteromonadales","Xanthomonadales","E01-9C-26_marine_group","KI89A_clade","SAR406_clade","Puniceicoccales","Verrucomicrobiales","Verrucomicrobia:Arctic97B-4","Verrucomicrobia:MB11C04","Acidobacteria_s6","Acidimicrobiales","NA","Order < 1% abund","Other")

rev_names_ordenados<-rev(names_ordenados)
y4$Order <- factor(y4$Order,levels = rev_names_ordenados)
#order incertea_seds change to other)
#colourCount = length(unique(y4$Order))
ordercolors<-c("SAR202_clade1"="#386CB0","SAR202_clade2"="#A6CEE3","SAR202_clade3"="aquamarine","SAR202_clade4"="darkslateblue","SAR202_clade5"="darkcyan","SAR202_clade"="darkblue","Flavobacteriales"="#7FC97F","Cytophagales"="#1B9E77","Sphingobacteriales"="#66A61E","Gemmatimonadetes:BD2-11"="lightcoral","Nitrospirales"="hotpink3","Planctomycetales"="thistle2","OM190"="thistle","SAR11_clade"="lightgoldenrod","Rhodobacterales"="gold3","Rhodospirillales"="yellow2","Rickettsiales"="khaki3","OCS116_clade"="khaki1","Rhizobiales"="orange","Methylophilales"="yellowgreen","Nitrosomonadales"="lightgoldenrod3","SAR324_clade(Marine_group_B)"="#B15928","Desulfobacterales"="lightsalmon3","Bdellovibrionales"="tan3","Oceanospirillales"="darkolivegreen","Cellvibrionales"="tomato2","Thiotrichales"="wheat","Salinisphaerales"="snow3","Alteromonadales"="orange2","Xanthomonadales"="darkorange1","E01-9C-26_marine_group"="tomato3","KI89A_clade"="maroon4","SAR406_clade"="palevioletred3","Puniceicoccales"="indianred2","Verrucomicrobiales"="firebrick2","Verrucomicrobia:Arctic97B-4"="red3","Verrucomicrobia:MB11C04"="orangered4","Acidobacteria_s6"="#B3DE69","Acidimicrobiales"="#80B1D3","NA"="black","Order < 1% abund"="lightcyan4","Other"="gray")

#Format orden de los Orders 

n <- 42
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))

ordercolors<-col_vector[1:42]

a<-ggplot(y4, aes(x = position, y = Abundance, fill = Order)) + geom_bar(stat = "identity",width=.85)+ scale_fill_manual(values = ordercolors) + theme_bw()+ ylab("Relative contribution [%]") +theme(legend.key.size = unit(1,"line"), strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=10)) +facet_grid(~Cruise,scales = "free_x",space = "free_x")+guides(fill=guide_legend(ncol=1))+xlab(NULL)

Bactplot<- a +theme(legend.position="none") #without legend

B<-ggplot(y4, aes(x = position, y = BactProd)) + geom_point() + geom_line()+ylab("Bact Prod [pmol/L/h]") +xlab(NULL)+facet_grid(~Cruise,scales = "free_x",space = "free_x")+theme_bw() +theme(strip.background = element_blank(),axis.text.y=element_text(size=8),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

C<-ggplot(y4, aes(x = position, y = BactAbund)) + geom_point() + geom_line()+ylab("Bac Ab [10^8 cells/L]") +xlab(NULL)+facet_grid(~Cruise,scales = "free_x",space = "free_x")+theme_bw() +theme(strip.background = element_blank(),axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

D<- ggplot(y4, aes(x = position, y = O2_Winkler)) + geom_point()+ geom_line() +ylab("Chl a [mg/m^3]") +xlab(NULL) +facet_grid(~Cruise,scales = "free_x",space = "free_x")+ylab("Oxygen [mg/L]")+theme_bw()+theme(strip.background = element_blank(),axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

plot_grid(B,C,D,Bactplot, align = "v", nrow = 4, rel_heights = c(1/6, 1/6,1/6, 1/2))
#SAR11NAAMES.R script to generate the Figure 5 SAR11 Bacteria in the Community Composition of the Seasonality of the microbial community composition in the western North Atlantic manuscript 

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


############Create phystV1filth (as in HetCCNAAMES.R) ############### 

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


############SAR11############### 
# Extract SAR11 from phystV1filth

SAR11_physeq = subset_taxa(phystV1filth, Order=="SAR11_clade")
Otu11<-otu_table(SAR11_physeq)
Env11<-sample_data(SAR11_physeq)

Taxsar <-as.matrix(read.table("SAR11V2.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

Tax11<-tax_table(Taxsar)

#Objeto Phyloseq
PHYSAR11<-phyloseq(Otu11,Env11,Tax11)

glomSar11<-tax_glom(PHYSAR11, taxrank="phylotype")

MeanStV2rel<-transform_sample_counts(glomSar11, function(x){x / sum(x)})
SAR11Subp<-as.data.frame(tax_table(MeanStV2rel)[,6])
ASV_frame11<-as.data.frame(otu_table(MeanStV2rel))
ASV_frame11[ "Taxa" ] <- SAR11Subp[,1]

ASV_frameSAR11_2 <- ASV_frame11[,-243]
rownames(ASV_frameSAR11_2) <- ASV_frame11[,243]

ASV_frw<-t(ASV_frameSAR11_2)
md_to_addSAR11<-as.data.frame(sample_data(MeanStV2rel))[,c(1,2,12,21)] #Cruise Station depth Temperature name
md_to_addSAR11$name<-row.names(md_to_addSAR11)

final_sar11<-cbind(ASV_frw,md_to_addSAR11)


fn2a_melt<-melt(final_sar11,id.vars=c("Cruise","position","depth","Temperature","name"), measure.vars = c("I","Ia.1","Ia.3","Ia.4","Ib","Ib.1","Ib.2","Ic.1","Ic.2","Id","uncharact1","II","IIa","IIa.A","IIa.B","IIb","clade N2","III","IIIa","IV","sar11"))

coloresbarplot = c("I"="lightgreen","Ia.1"="darkseagreen","Ia.3"="darkolivegreen4","Ia.4"="darkkhaki","Ib"="darkgoldenrod","Ib.1"="darkorange","Ib.2"="coral","Ic.1"="mediumpurple1","Ic.2"="palevioletred","Id"="rosybrown1","uncharact1"="wheat2","II"="slategray2","IIa"="skyblue","IIa.A"="midnightblue","IIa.B"="royalblue","IIb"="lightseagreen","clade N2"="gray34","III"="lightpink4","IIIa"="orangered4","IV"="darkslategray4","sar11"="black")


plota<-ggplot(fn2a_melt, aes(x = position, y = value, fill = factor(variable, levels=c("I","Ia.1","Ia.3","Ia.4","Ib","Ib.1","Ib.2","Ic.1","Ic.2","Id","uncharact1","II","IIa","IIa.A","IIa.B","IIb","clade N2","III","IIIa","IV","sar11")))) + geom_bar(stat = "identity",width=.85)+ scale_fill_manual(values = coloresbarplot) + theme_bw()+ ylab("Relative contribution [%]") +theme(strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=12)) +facet_grid(~Cruise,scales = "free_x",space = "free_x")

plot1<-plota+theme(legend.position="none")

#svg("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/NAAMES_SAR11/SAR11_withLeg2.svg", width=12, height=6)
plota
#dev.off()

#svg("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/NAAMES_SAR11/SAR11_woLeg.svg", width=12, height=6)
plot1
#dev.off()

#Chao1 comparison

set.seed(717)
phyeuphminrarH = rarefy_even_depth(phystV1filth, sample.size = 6740)

dfindexesH<- estimate_richness(phyeuphminrarH, measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))

#Df with indexes as column, we need to add two more columns linking region and cruise 

df_sampledatararH<-data.frame(sample_data(phyeuphminrarH),stringsAsFactors=FALSE)


#row.names don’t match, change "-" from dfindexes to "."

rownames(df_sampledatararH) <- gsub("-", ".", rownames(df_sampledatararH))

#Now merge columns region and cruise

dfindexesH$cruise = df_sampledatararH[match(row.names(dfindexesH), row.names(df_sampledatararH)),"Cruise"]
dfindexesH$date = df_sampledatararH[match(row.names(dfindexesH), row.names(df_sampledatararH)),"Date"]
dfindexesH$depth = df_sampledatararH[match(row.names(dfindexesH), row.names(df_sampledatararH)),"depth"]
dfindexesH$Station = df_sampledatararH[match(row.names(dfindexesH), row.names(df_sampledatararH)),"Station"]
dfindexesH$position = df_sampledatararH[match(row.names(dfindexesH), row.names(df_sampledatararH)),"position"]

##############Agregar a dfindexesH SAR11 chao1 

#First rarefy to the min -> 6740
set.seed(717)
PHYSAR11rar = rarefy_even_depth(PHYSAR11, sample.size = 6740)

dfindexesSAR<- estimate_richness(PHYSAR11rar, measures=c("Observed", "Chao1"))

rownames(dfindexesSAR) <- gsub("-", ".", rownames(dfindexesSAR))

dfindexesH$SAR11Observed = dfindexesSAR[match(row.names(dfindexesH), row.names(dfindexesSAR)),"Observed"]
dfindexesH$SAR11Chao1 = dfindexesSAR[match(row.names(dfindexesH), row.names(dfindexesSAR)),"Chao1"]

Index<-ggplot(dfindexesH, aes(x = position, y = Chao1)) + geom_point(alpha=0.25)+ geom_line(alpha=0.25)+geom_line(aes(y = SAR11Chao1, colour = "red")) + geom_point(aes(y = SAR11Chao1, colour = "red")) +ylab("Chao1") +xlab(NULL) +facet_grid(~cruise,scales = "free_x",space = "free_x")+theme_bw()+theme(strip.background = element_blank(),axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

IndexWOleg<-Index+theme(legend.position="none")

#svg("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_annual/NAAMES_SAR11/SAR11_fig5corrected.svg", width=12,height=7)
plot_grid(IndexWOleg,plot1, align = "v", nrow = 2, rel_heights = c(1/4, 3/4))
#dev.off()

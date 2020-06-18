
#PhytoCCNAAMES.R script to generate the Figure 3 of the Seasonality of the microbial community composition in the western North Atlantic manuscript 

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

count_tab <- read.table("Photlink_final_collapsedVR2.otu", header=T, row.names=1, check.names=F)
sample_info_tab <- read.table(Merge_envFileDNA.txt", header=T, row.names=1, check.names=F, sep ="\t")
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

glomV1filt<-tax_glom(phyeuphminN1N2, taxrank="Taxa")


colSums(otu_table(phyeuphminN1N2))

MeanStV2rel<-transform_sample_counts(glomV1filt, function(x){x / sum(x)})

######### Build the bar plots ############


#taxaSubp<-as.data.frame(tax_table(MeanStV2rel)[,2])
#ASV_frame<-as.data.frame(otu_table(MeanStV2rel))
#ASV_frame[ "Taxa" ] <- taxaSubp[,1]

#dim(ASV_frame)
#[1] 22 134
#ASV_frame2 <- ASV_frame[,-134]
#rownames(ASV_frame2) <- ASV_frame[,134]


#ASV_frw<-t(ASV_frame2)
#md_to_add<-as.data.frame(sample_data(MeanStV2rel))[,c(1,2,11,20)] #Cruise Station depth Temperature name

#md_to_add$name<-row.names(md_to_add)

#final_2a<-cbind(ASV_frw,md_to_add)

#write.table(final_2a,file= "barplots_addposition_min1600.txt", quote=FALSE, sep = "\t")

#I added manually the stats of the run so we can get everything together from the provided table 

final_2a<-read.table("barplots_addposition_min1600.txt",header=T, row.names=1, check.names=F, sep = "\t")

fn2a_melt<-melt(final_2a,id.vars=c("Cruise","Station","depth","Temperature","name", "position","Chao1Phyto","Totalphyto","TotalHet","Totalreads","phytoperc","hetperc","TChl_a"), measure.vars = c("Diatoms","Bolidophyceae","Dictyochophyceae","Pelagophyceae","Chrysophyceae","Micromonas","Bathycoccus","OstreococcusI","OstreococcusII","Cryptophyceae","Prymnesiophyceae","Rappemonad", "ASV357","Other plastid","ProchlorococcusHLI","ProchlorococcusHLII","ProchlorococcusLLI","SynechococcusI","SynechococcusIV","SynechococcusII","Other Cyanobacteria","Not assigned"))

coloresbarplot = c("Diatoms"="blue","Bolidophyceae"="cadetblue","Dictyochophyceae"="lightskyblue","Pelagophyceae"="aquamarine","Chrysophyceae"="turquoise","Prymnesiophyceae"="darkgoldenrod3 ","Rappemonad"=	"gold2","Cryptophyceae"="coral3" ,"Micromonas"="forestgreen","Bathycoccus"="limegreen","OstreococcusII"="olivedrab","OstreococcusI"="palegreen4","PrasinophyceaeI"="greenyellow","Other plastid"="lightgreen","ASV357"="lemonchiffon3","ProchlorococcusHLI"="lightcoral","ProchlorococcusHLII"="hotpink","ProchlorococcusLLI"="maroon1","SynechococcusI"="blueviolet","SynechococcusII"="mediumpurple","SynechococcusIV"="plum3","Other Cyanobacteria"="mediumvioletred","others"="cornsilk4","Not assigned"="gray34")

plota<-ggplot(fn2a_melt, aes(x = position, y = value, fill = variable)) + geom_bar(stat = "identity",width=.85)+ scale_fill_manual(values = coloresbarplot) + theme_bw()+ ylab("Relative contribution [%]") +theme(strip.background = element_blank(),strip.text.x =element_blank(),axis.text.y=element_text(size=12), axis.text.x=element_text(size=12,angle = 60, hjust = 1, vjust=.5), text = element_text(size=12),strip.text = element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=12)) +facet_grid(~Cruise,scales = "free_x",space = "free_x")

plot1<-plota+theme(legend.position="none")

d<- ggplot(fn2a_melt, aes(x = position, y = .1)) + geom_tile(aes(fill = Temperature),colour = "white")+scale_fill_gradient(low = "blue", high = "firebrick4")+ theme_bw()+theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ ylab("Chao1") +facet_grid(~Cruise,scales = "free_x",space = "free_x")

d2<-d+theme(legend.position="none")

Bottom_p<-plot_grid(d2,plot1, align = "v", nrow = 2, rel_heights = c(1/6, 5/6))

#Add other panels

#Chao1
a<-ggplot(fn2a_melt, aes(x = position, y = Chao1Phyto)) + geom_point()+ theme_bw()+theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ ylab("Chao1") +facet_grid(~Cruise,scales = "free_x",space = "free_x")

#b chloroplast perce
b<-ggplot(fn2a_melt, aes(x = position, y = phytoperc)) + geom_point()+ theme_bw()+theme(strip.background = element_blank(),strip.text.x = element_blank(), axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Rel. contribution[%]") +facet_grid(~Cruise,scales = "free_x",space = "free_x")

#c chla
c<-ggplot(fn2a_melt, aes(x = position, y = TChl_a)) + geom_point()+ theme_bw() +theme(strip.background = element_blank(),axis.text.y=element_text(size=12),axis.title.y=element_text(size=10),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Chl a [mg/m^3]") +facet_grid(~Cruise,scales = "free_x",space = "free_x")

plot_grid(c,b,a,Bottom_p, nrow = 4, rel_heights = c(1/6, 1/6,1/6, 1/2))

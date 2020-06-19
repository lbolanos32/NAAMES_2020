setwd("/Users/luisbolanos/Documents/OSU_postdoc/NAAMES/Vits_Chris")
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
library("gplots")
library("VennDiagram")
library("UpSetR")
library("cowplot")
library("ggrepel")
library(data.table)
#library("vegan")
library(decontam)

setwd("/Users/luisbolanos/Documents/MisDrafts/InProgress/NAAMES_vits_CSetal/GithubChris/")

##########FINAL ANALYSIS ##############
count_tab <- read.table("vits.otu", header=T, row.names=1, check.names=F)
sample_info_tab <- read.table("vits.env", header=T, row.names=1, check.names=F, sep ="\t")
tax_tab <- as.matrix(read.table("vitsf.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

#######CleanContamin

OTU = otu_table(count_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
SAM = sample_data(sample_info_tab)

List to remove
contam_asvs <-c("N287","N370","N504","N932","N1195","N1575","N2862","N2934","N3307","N3544","N8532","N24639","N25035","N25138","N25145","N25317","N25556","N26238")

####Prune them from the VITS object

asv_tab_no_contam <- count_tab[!row.names(count_tab) %in% contam_asvs, ]# making new OTU table
  
asv_tax_no_contam <- tax_tab[!row.names(tax_tab) %in% contam_asvs, ] # making new taxonomy table

#factorizing environmental to remove the SLB contaminant sample 

sample_info_tab$depth <- factor(sample_info_tab$depth,levels = c("5", "25", "75", "100","150","300","1000","2600"))
sample_info_tab$station<- factor(sample_info_tab$station,levels = c("1", "2", "2.1", "3","4"))

#Remove 1k,2.6k,2.8k and SLBneg

asv_tab_no_contam1<-asv_tab_no_contam[,-c(12,13,14,25)] 


OTUc = otu_table(asv_tab_no_contam1, taxa_are_rows = TRUE)
TAXc = tax_table(asv_tax_no_contam)
SAMc= sample_data(sample_info_tab)

VITSc<-phyloseq(OTUc,TAXc,SAMc) #Clean phyloseq object 

VITScf= filter_taxa(VITSc, function(x) sum(x > 2) > (0.0015*length(x)), TRUE) #clean and filter

#VITScfrel<-transform_sample_counts(VITScf, function(x){x / sum(x)})


####################Constrained ordination###################

#CAP# Cada que se crea un phyloseq object format 

###Preliminary info:
### min(colSums(otu_table(VITScf)))
####[1] 40762

### min(colSums(otu_table(BACT)))
####[1] 30772

### min(colSums(otu_table(PHY)))
####[1] 619 N4S4-300_S72, 1746 N4S3-300_S62. #a(probar con 1740) y luego todo. 

set.seed(717)
phyeuphminrarVITScf = rarefy_even_depth(VITScf, sample.size = 40500)
#phyeuphminrarVITSc = rarefy_even_depth(BACT, sample.size = 30700)
#phyeuphminrarVITSc = rarefy_even_depth(PHY, sample.size = 1740)

phyna <- phyeuphminrarVITScf %>%
subset_samples(
!is.na(DOC)& #Chris picks
!is.na(Fluorescence)& #Chris picks
!is.na(PO4)& #Chris picks
!is.na(Density00)& #Chris picks
!is.na(BactAbund)& #Chris picks
!is.na(BactProd)& #Chris picks
!is.na(AmMP)& #Chris picks
!is.na(B1)& #Chris picks
!is.na(cHET)& #Chris picks
!is.na(HET)& #Chris picks
!is.na(HMP) #Chris picks
)
#!is.na(BactAbund)&
#!is.na(DOC)&
#!is.na(Fluorescence)&
#!is.na(N_N)&
#!is.na(NH4)&
#!is.na(NO2)&
#!is.na(Oxygen_uM)&
#!is.na(PO4)&
#!is.na(Sal11)&
#!is.na(SiO4)&
#!is.na(Temperature)&
#!is.na(AmMP)&
#!is.na(B1)&
#!is.na(cHET)&
#!is.na(HET)&
#!is.na(HMP)
#)



bray_not_na <- phyloseq::distance(physeq = phyna, method = "bray")
 
cap_ord <- ordinate(
physeq = phyna, 
 method = "CAP",
 distance = bray_not_na,
 formula = ~ BactAbund+BactProd+DOC+Fluorescence+PO4+Density00+AmMP+B1+cHET+HET+HMP) #New
# formula = ~ BactAbund+DOC+Fluorescence+N_N+NH4+NO2+Oxygen_uM+PO4+Sal11+SiO4+Temperature+AmMP+B1+cHET+HET+HMP)
 
#PLOT
cap_plot <- plot_ordination(
physeq = phyna, 
ordination = cap_ord, 
color = "station", 
 axes = c(1,2)
 ) + 
 aes(shape = depth) + scale_shape_manual(values=c(15,16,17,18,1,5,6,7)) + 
 geom_point(aes(colour = station), alpha = 0.9, size = 5) + 
 geom_point(colour = "grey90", size = 1.5)+ 
 scale_color_manual(values=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00", "#FFFF33")) + theme(axis.text.x=element_text(angle=0,size=14), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=14))+theme_bw()


#Arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
    yend = CAP2, 
    x = 0, 
    y = 0, 
    shape = NULL, 
    color = NULL, 
    label = labels)

label_map <- aes(x = 1.3 * CAP1, 
    y = 1.3 * CAP2, 
    shape = NULL, 
    color = NULL, 
    label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .6, 
    data = arrowdf, 
    color = "gray48", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )


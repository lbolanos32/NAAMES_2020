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



##########FINAL ANALYSIS ##############
count_tab <- read.table("vits.otu", header=T, row.names=1, check.names=F)
sample_info_tab <- read.table("vits.env", header=T, row.names=1, check.names=F, sep ="\t")
tax_tab <- as.matrix(read.table("vitsf.tax", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))

#######CleanContamin

OTU = otu_table(count_tab, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
SAM = sample_data(sample_info_tab)

#List to remove
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

Phylumcolors<-c("Proteobacteria"="cadetblue2","Chloroflexi"="seagreen3","Marinimicrobia_(SAR406_clade)"="brown3","Bacteroidetes"="mediumorchid4","Verrucomicrobia"="green","Actinobacteria"="darkgoldenrod3","Gemmatimonadetes"="mediumorchid4","Acidobacteria"="lightsalmon","Planctomycetes"="pink","Deinococcus-Thermus"="black","Nitrospirae"="gray32","Spirochaetae"="cyan","Gracilibacteria"="orange","SHA-109"="brown","Tenericutes"="forestgreen3","TM6"="dark olive","Lentisphaerae"="coral2","Parcubacteria"="darkcyan","Saccharibacteria"="chocolate4","Cyanobacteria"="bisque","Firmicutes"="azure3","Microgenomates"="gold4","Hydrogenedentes"="red","Omnitrophica"="firebrick4","Chlorophyta"="forestgreen","Stramenopiles"="darkgoldenrod3","Cryptophyta"="mediumpurple3","Prymnesiophyceae"="lemonchiffon3","NA"="gray32")

#STATION 1#

N4S1 = get_variable(VITScf, "names") %in% c("N4S1_5","N4S1_25","N4S1_75")
sample_data(VITScf)$N4S1 <- factor(N4S1)
phystN4S1<-subset_samples(VITScf, N4S1 %in% TRUE)

##pattern 1 station 1 type A##

diagdds_pol = phyloseq_to_deseq2(phystN4S1, ~H_AmMP) #este pattern aplica para H_HET, H_HMP
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resA = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01
fold_a <- resA[which(abs(resA$log2FoldChange)>2), ]sigtab_a = fold_a [which(fold_a$padj < alpha), ]
sigtab_a = cbind(as(sigtab_a, "data.frame"), as(tax_table(phystN4S1)[rownames(sigtab_a), ], "matrix"))
sigtab_names_a<-row.names(sigtab_a)
sigtab_a$names<-row.names(sigtab_a) 
sigtab_a$pattern<-"a" 
differentVITScA <- subset(otu_table(phystN4S1), rownames(otu_table(phystN4S1)) %in% sigtab_names_a)

#A with legend
plotA<-ggplot(sigtab_a, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_a))) + geom_point(data = sigtab_a, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=14))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")


##pattern 2 station 1 type B##
diagdds_pol = phyloseq_to_deseq2(phystN4S1, ~H_B1) #este pattern aplica para H_cHET tmb
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resB = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_b <- resB[which(abs(resB$log2FoldChange)>2), ]sigtab_b = fold_b[which(fold_b$padj < alpha), ]
sigtab_b = cbind(as(sigtab_b, "data.frame"), as(tax_table(phystN4S1)[rownames(sigtab_b), ], "matrix"))
sigtab_names_b<-row.names(sigtab_b)
sigtab_b$names<-row.names(sigtab_b)
sigtab_b$pattern<-"b" 
differentVITScB <- subset(otu_table(phystN4S1), rownames(otu_table(phystN4S1)) %in% sigtab_names_b)

#B with legend
plotB<-ggplot(sigtab_b, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_b))) + geom_point(data = sigtab_b, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=14))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")



#STATION 2#

N4S2 = get_variable(VITScf, "names") %in% c("N4S2_5","N4S2_25","N4S2_75")
sample_data(VITScf)$N4S2 <- factor(N4S2)
phystN4S2<-subset_samples(VITScf, N4S2 %in% TRUE)

##pattern 1 station 2 type C##
diagdds_pol = phyloseq_to_deseq2(phystN4S2, ~H_AmMP) #este pattern aplica para H_HET, H_HMP
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resC = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_c <- resC[which(abs(resC$log2FoldChange)>2), ] #LOG2sigtab_c = fold_c[which(fold_c$padj < alpha), ] #padj
sigtab_c$names<-row.names(sigtab_c) 
sigtab_c = cbind(as(sigtab_c, "data.frame"), as(tax_table(phystN4S2)[rownames(sigtab_c), ], "matrix"))
sigtab_naames_c<-row.names(sigtab_c)
sigtab_c$pattern<-"c" 
differentVITScC <- subset(otu_table(phystN4S2), rownames(otu_table(phystN4S2)) %in% sigtab_naames_c)
plotC<-ggplot(sigtab_c, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_c))) + geom_point(data = sigtab_c, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red") ##NO HUBO DIFERENCIA 

##pattern 2 station 2 type D##
diagdds_pol = phyloseq_to_deseq2(phystN4S2, ~H_B1) 
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resD = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_d <- resD[which(abs(resD$log2FoldChange)>2), ]sigtab_d = fold_d[which(fold_d$padj < alpha), ]
sigtab_d = cbind(as(sigtab_d, "data.frame"), as(tax_table(phystN4S2)[rownames(sigtab_d), ], "matrix"))
sigtab_naames_d<-row.names(sigtab_d)
sigtab_d$names<-row.names(sigtab_d)
sigtab_d$pattern<-"d" 
differentVITScD <- subset(otu_table(phystN4S2), rownames(otu_table(phystN4S2)) %in% sigtab_naames_d)
plotD<-ggplot(sigtab_d, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_d))) + geom_point(data = sigtab_d, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")

##pattern 3 station 2 type E###
diagdds_pol = phyloseq_to_deseq2(phystN4S2, ~H_cHET) 
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resE = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_e <- resE[which(abs(resE$log2FoldChange)>2), ]sigtab_e = fold_e[which(fold_e$padj < alpha), ]
sigtab_e = cbind(as(sigtab_e, "data.frame"), as(tax_table(phystN4S2)[rownames(sigtab_e), ], "matrix"))
sigtab_naames_e<-row.names(sigtab_e)
sigtab_e$names<-row.names(sigtab_e) 
sigtab_e$pattern<-"e" 
differentVITScE <- subset(otu_table(phystN4S2), rownames(otu_table(phystN4S2)) %in% sigtab_naames_e)
plotE<-ggplot(sigtab_e, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_e))) + geom_point(data = sigtab_e, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red") ####NO significative differences

F is station 2.1 not included here 


#STATION 3#

N4S3 = get_variable(VITScf, "names") %in% c("N4S3_5","N4S3_25","N4S3_75")
sample_data(VITScf)$N4S3 <- factor(N4S3)
phystN4S3<-subset_samples(VITScf, N4S3 %in% TRUE)

##pattern 1 station 3 type G##

diagdds_pol = phyloseq_to_deseq2(phystN4S3, ~H_B1) 
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resG = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_g <- resG[which(abs(resG$log2FoldChange)>2), ]sigtab_g = fold_g[which(fold_g$padj < alpha), ]
sigtab_g = cbind(as(sigtab_g, "data.frame"), as(tax_table(phystN4S3)[rownames(sigtab_g), ], "matrix"))
sigtab_naames_g<-row.names(sigtab_g)
sigtab_g$names<-row.names(sigtab_g) 
sigtab_g$pattern<-"g" 
differentVITScG <- subset(otu_table(phystN4S3), rownames(otu_table(phystN4S3)) %in% sigtab_naames_g)
plotG<-ggplot(sigtab_g, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_g))) + geom_point(data = sigtab_g, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")

##pattern 2 station 3 type H##

diagdds_pol = phyloseq_to_deseq2(phystN4S3, ~H_HET) 
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resH = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_h <- resH[which(abs(resH$log2FoldChange)>2), ]sigtab_h = fold_h[which(fold_h$padj < alpha), ]
sigtab_h = cbind(as(sigtab_h, "data.frame"), as(tax_table(phystN4S3)[rownames(sigtab_h), ], "matrix"))
sigtab_naames_h<-row.names(sigtab_h)
sigtab_h$names<-row.names(sigtab_h)
sigtab_h$pattern<-"h"
differentVITScH <- subset(otu_table(phystN4S3), rownames(otu_table(phystN4S3)) %in% sigtab_naames_h)
plotH<-ggplot(sigtab_h, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_h))) + geom_point(data = sigtab_h, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")

#STATION 4#	 

N4S4 = get_variable(VITScf, "names") %in% c("N4S4_5","N4S4_25","N4S4_75")
sample_data(VITScf)$N4S4 <- factor(N4S4)
phystN4S4<-subset_samples(VITScf, N4S4 %in% TRUE)

##pattern 1 station 4 type I##
diagdds_pol = phyloseq_to_deseq2(phystN4S4, ~H_HET) 
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resI = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_i <- resI[which(abs(resI$log2FoldChange)>2), ]sigtab_i = fold_i[which(fold_i$padj < alpha), ]
sigtab_i = cbind(as(sigtab_i, "data.frame"), as(tax_table(phystN4S4)[rownames(sigtab_i), ], "matrix"))
sigtab_naames_i<-row.names(sigtab_i)
sigtab_i$names<-row.names(sigtab_i)
sigtab_i$pattern<-"i"
differentVITScI <- subset(otu_table(phystN4S4), rownames(otu_table(phystN4S4)) %in% sigtab_naames_i)
plotI<-ggplot(sigtab_i, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_i))) + geom_point(data = sigtab_i, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")

##pattern 2 station 4 type J##
diagdds_pol = phyloseq_to_deseq2(phystN4S4, ~H_cHET) 
diagdds_pol = DESeq(diagdds_pol, test="Wald", fitType="parametric")
resJ = results(diagdds_pol, cooksCutoff = FALSE)alpha = 0.01fold_j <- resJ[which(abs(resJ$log2FoldChange)>2), ]sigtab_j = fold_j[which(fold_j$padj < alpha), ]
sigtab_j = cbind(as(sigtab_j, "data.frame"), as(tax_table(phystN4S4)[rownames(sigtab_j), ], "matrix"))
sigtab_naames_j<-row.names(sigtab_j)
sigtab_j$names<-row.names(sigtab_j) 
differentVITScJ <- subset(otu_table(phystN4S4), rownames(otu_table(phystN4S4)) %in% sigtab_naames_j)
plotJ<-ggplot(sigtab_j, aes(x=Class, y=log2FoldChange, color=Phylum, label=rownames(sigtab_j))) + geom_point(data = sigtab_j, aes(x = Class, y = log2FoldChange, size = baseMean), position = position_jitter(width = 0.3, height = 0.15))+theme(text = element_text(size=12), axis.text.x = element_text(angle = -45, hjust = 0, vjust=0.5),legend.text=element_text(size=12))+scale_color_manual(values =Phylumcolors)+geom_hline(yintercept=0)+geom_hline(yintercept=c(-5,5),linetype="dashed", color = "red")

#########################Vitamins B1#####################

##BIND ONLY differences in sigtabB (ST1), sigtabD (ST2)##

sigtab_fulldifB1<-data.frame(rbind(sigtab_b,sigtab_d)) 

sigtab_fulldifB1<-rbind(sigtab_b,sigtab_d) ##CREATE A MASTER LIST OF THE DIFFERENTIAL ASVs

unique_sigtab_B1<-sigtab_fulldifB1[!duplicated(sigtab_fulldifB1[,'names']),] ##CREATE A unique list by names ASVs

uniq_names_extract<-unique_sigtab_B1$names #### step two to get the list of uniques

diffB<- subset(resB, rownames(resB) %in% uniq_names_extract) #B1
new_diffB<-data.frame(diffB[,c(1,2,6)])
new_diffB['Station']='Station 1'
new_diffB$names <- rownames(new_diffB)

diffD<- subset(resD, rownames(resD) %in% uniq_names_extract) #B1
new_diffD<-data.frame(diffD[,c(1,2,6)])
new_diffD['Station']='Station 2'
new_diffD$names <- rownames(new_diffD)


mergeB1<-rbind(new_diffB,new_diffD)
rownames(mergeB1) <- c()

non_dup_namesB1<-mergeB1[!duplicated(mergeB1[,'names']),]$names #SacarTaxa de mergB1 no duplicada
taxB1<-tax_table(VITSc)#Extraer de physeq
mergeB1_taxa<-data.frame(taxB1[row.names(taxB1) %in% non_dup_namesB1, ])
mergeB1_taxa$taxname<-paste(mergeB1_taxa$Order) #mergeB1_taxa$Family,row.names(mergeB1_taxa), sep = ":")
mergeB1_taxa$names<-row.names(mergeB1_taxa) 

mergeB1_taxa$taxname[mergeB1_taxa$names == "N853"] <- "Diatom:NA"
mergeB1_taxa$taxname[mergeB1_taxa$names == "N24786"] <- "Bacteria:NA"
mergeB1_taxa$taxname[mergeB1_taxa$names == "N688"] <- "SAR406"
mergeB1_taxa$taxname[mergeB1_taxa$names == "N1599"] <- "Bacteria:NA"
taxa_desorder<-mergeB1_taxa[,c(7,8)]
taxa_order<-taxa_desorder[order(taxa_desorder$taxname),]
taxa_prueba<- rev(taxa_order$taxname)
for_format<-rev(taxa_order$names)
mergeB1$names <- factor(mergeB1$names, levels =for_format )

mergeB1$log2FoldChange<-mergeB1$log2FoldChang * -1 ## set the log2foldchange backwards to take as reference enrichment in station 2

library(scales)

p<- ggplot(mergeB1, aes(x = Station, y = names)) + geom_tile(data =subset(mergeB1, padj<0.01),aes(x = Station, y = names, fill = log2FoldChange))+scale_fill_gradient2(low = muted("blue"), mid = "white",high = muted("red"), midpoint = 0)+geom_text(data =subset(mergeB1, baseMean>100), aes(label=ifelse(padj<0.01, "*","")))+ylab("Order") +theme_bw()
p+scale_y_discrete(labels=format(taxa_prueba))+ theme(axis.text.y = element_text(size=16,color="black"), axis.text.x=element_text(size=18,color="black"),axis.title.y = element_text(color="black", size=18))

#svg("B1_diffabV2.svg", height=12, width=9)
p+scale_y_discrete(labels=format(taxa_prueba))+ theme(axis.text.y = element_text(size=16,color="black"), axis.text.x=element_text(size=18,color="black"),axis.title.y = element_text(color="black", size=18))
#dev.off()

#README SAR11 Phylogenetic tree used as DB/frame to place the NAAMES V1_V2 Sequences

As part of the NAAMES study, in particular of the manuscript "Seasonality of the microbial community composition in the western North Atlantic", we provide a DB to accurately clasiffied SAR11 16S rRNA amplicons.

In this directory we provide the elements to recreate the DB through the setupdb.pl command in Phyloassigner. 

The phylogenetic tree was build based on ~full-length 16S rRNA genes. Therefore it is suitable to run any variable region the user decided to analyze in their high-throughput amplicon experiment.  

Disclaimer: It has been consistently observed that SAR11 16S rRNA region V1-V2 resolves accurately the phylogenetic position within the clade. V4 and V9 regions tend to be more conserved, lacking the necessary phylogenetic signal to resolve into exterior terminal nodes. 

Also we provide the directory ready to run. Feel free to download it and try it.

Phyloassigner can be downloaded from this link: https://www.awi.de/en/science/special-groups/scientific-computing/bioinformatics/software.html

Instructions:

1) Download the directory: Sar11CoreV3_aln_inpclean1_10.phyloassignerdb

2) Install Phyloassigner

3) Run phyloassigner as follows:

perl phyloassigner.pl --hmmerdir /PATH/TO/HHMER/ --pplacerdir PATH/TO/PPLACER/ -o PATH/TO/DIRECTORY/OUTPUT_TO_BE_STORED /PATH/TO/Sar11CoreV3_aln_inpclean1_10.phyloassignerdb YourSAR11Sequences.fasta

4) Viola!




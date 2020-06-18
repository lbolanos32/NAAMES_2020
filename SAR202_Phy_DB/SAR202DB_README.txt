#README SAR202 Phylogenetic tree used as DB/frame to place the NAAMES V1_V2 Sequences

This Phyloassigner DB was built based on the 16S rRNA tree phylogenetic shown on: SAR202 Genomes from the Dark Ocean Predict Pathways for the Oxidation of Recalcitrant Dissolved Organic Matter 
Zachary Landry, Brandon K. Swan, Gerhard J. Herndl, Ramunas Stepanauskas, Stephen J. Giovannoni

Citation:
Please cite Landry et al., 2017 and Bolaños & Giovannoni (in prep, to be updated)

We provide the directory ready to run. Feel free to download it and try it.

Disclaimer:
This DB was made using the ~full-length 16S rRNA of the corresponding organisms. Therefore it can be used to classify any amplicon region or any 16S rRNA fragment retrieved from SAGS, MAGS, etc. It has only be used to classify V1-V2 16S rRNA region, which outputs accurate taxonomic classification.

Phyloassigner can be downloaded from this link: https://www.awi.de/en/science/special-groups/scientific-computing/bioinformatics/software.html

Instructions:

1) Download the directory: 

2) Install Phyloassigner

3) Run phyloassigner as follows:

perl phyloassigner.pl --hmmerdir /PATH/TO/HHMER/ --pplacerdir PATH/TO/PPLACER/ -o PATH/TO/DIRECTORY/OUTPUT_TO_BE_STORED /PATH/TO/ss_aln.phyloassignerdb YourSAR202Sequences.fasta

4) Viola!

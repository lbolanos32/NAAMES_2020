library(dada2)

#filter
path <- "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2”

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(220,190), 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

saveRDS(errF, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/errF.rds")
saveRDS(errR, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/errR.rds")

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

saveRDS(dadaFs, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/dadaFs_N.rds")
saveRDS(dadaRs, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/dadaRs_N.rds")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

saveRDS(mergers, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/mergers.rds")
seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab)

dim(seqtab.nochim)

saveRDS(seqtab.nochim, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/seqtab.nochim.rds")

taxa <- assignTaxonomy(seqtab.nochim, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/BIOS/ProcessSeqs/SEQ1pr/silva_nr_v123_train_set.fa")
unname(head(taxa))

saveRDS(taxa, "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/taxa.rds")

write.table(cbind(t(seqtab.nochim) , taxa), "/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/seqtab-nochimtaxa.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(taxa,"/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/N1N2/taxa.txt“, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE )
uniquesToFasta(seqtab.nochim, fout="/nfs0/Giovannoni_Lab/workspaces/bolanosl/NAAMES/DADA2/N1N2/seqtab-nochim.reprs.fna", ids=colnames(seqtab.nochim))

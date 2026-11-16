# Project: Tri5MB R2 + BarleyEndos + HackFlex test
# Title: DADA2 pipeline Fungi 
# Author: Briana Whitaker
# Collaborators: Guixia Hao, Nick Rhoades, Nate Tyler

# Fungi data from NextSeq P2 run (~92GB)   # took 4 hrs to run main code on 40 threads, ~1.5hrs to run for-loop seq count at beginning
################################################################################

set.seed(233)
require(dada2)
setwd("/project/mpm_fusarium_genomics/whitaker/niko19/fungi")
FWD <- "GTGARTCATCGAATCTTTG" #fITS7
REV <- "TCCTSCGCTTATTGATATGC" #ITS4ngs
################################################################################

path <- "/project/mpm_fusarium_genomics/whitaker/niko19/fungi"  
#list.files(path) #654 files
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))

#test
#fnFs <- fnFs[1:10]
#fnRs <- fnRs[1:10]

length(fnFs); length(fnRs)  #327 files each
get.Sample.ID <- function(fname) strsplit(basename(fname), "_")[[1]][1]  
ids.dat <- as.data.frame(cbind("num" = unname(sapply(fnFs, get.Sample.ID)), fnFs, fnRs))
Sample.ID <- ids.dat$num
#sort(Sample.ID)
################################################################################

#print("loop count seqs")
#i = 1 #useful for checking
#fwdSeqs <- list()
#revSeqs <- list()
#for (i in 1:length(fnFs)) {
#fwdSeqs[[i]] <- length(sapply(fnFs[i], getSequences))
#revSeqs[[i]] <- length(sapply(fnRs[i], getSequences))
#}
#identical(c(unlist(fwdSeqs)),c(unlist(revSeqs))) #TRUE
#SeqsOrig.df <- data.frame(SampleID = c(basename(fnFs)) ,
#         OrigSeqsFwd = c(unlist(fwdSeqs)),  OrigSeqsRev = c(unlist(revSeqs)))
#rownames(SeqsOrig.df) <- SeqsOrig.df$SampleID
#SeqsOrig.df <- SeqsOrig.df[,-1]
#write.csv(SeqsOrig.df, "/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_TrackSequences_PriorFiltering.csv")
SeqsOrig.df <- read.csv("/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_TrackSequences_PriorFiltering.csv", row.names = 1)
################################################################################

print("first filter")
fnFs.filtN <- file.path(path, "filtN-", basename(fnFs)) 
fnRs.filtN <- file.path(path, "filtN-", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread=TRUE, verbose=TRUE)
################################################################################

print("cutadapt step")
#module load cutadapt/4.0
cutadapt <- "/software/el9/apps/cutadapt/4.0/bin/cutadapt" #change to location on your machine
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt-")
if (!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for (i in seq_along(fnFs)) {
 system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-m", 50, "-e", 0.1,
             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                   fnFs.filtN[i],     fnRs.filtN[i] ))   }
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))
################################################################################

print("quality figures")
pdf("/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung-SequenceQuality.pdf", width=16/2.54, height=16/2.54)
plotQualityProfile(cutFs[1:4])
plotQualityProfile(cutRs[1:4])
dev.off()
################################################################################

print("second filter")
filtFs <- file.path(path, "filtered-final", basename(cutFs))
filtRs <- file.path(path, "filtered-final", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                truncQ = 2, minLen = 50, rm.phix = TRUE, multithread=TRUE, verbose=TRUE)
save(out, file="/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_FilterOut.RData")
#load("/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_FilterOut.RData")
write.csv(out, "/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung-SecondFilteringStep.csv")
names(filtFs) <- Sample.ID
names(filtRs) <- Sample.ID
################################################################################

print("learn errors")
not.lost <- out[,"reads.out"] > 0  
print("No. of samples not.lost and out")
length(not.lost); dim(out)[1]  
filtFs <- filtFs[not.lost]
filtRs <- filtRs[not.lost]
errF <- learnErrors(filtFs, nbases = 1e8, multithread=TRUE, randomize=TRUE) 
errR <- learnErrors(filtRs, nbases = 1e8, multithread=TRUE, randomize=TRUE)
save(errF, file="/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_errF.RData")
save(errR, file="/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_errR.RData")
pdf("/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung-ErrorLearning.pdf", width=16/2.54, height=16/2.54)
plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE)
dev.off()
################################################################################

print("derep-denoise-mergers loop")

#make lists to hold the loop output
mergers <- vector("list", length(Sample.ID))
names(mergers) <- Sample.ID
ddF <- vector("list", length(Sample.ID))
names(ddF) <- Sample.ID
ddR <- vector("list", length(Sample.ID))
names(ddR) <- Sample.ID

for(i in Sample.ID) {
  cat("Processing:", i, "\n")
    derepF <- derepFastq(filtFs[[i]])
    dadaF <- dada(derepF, err=errF, multithread=TRUE)
    ddF[[i]] <- dadaF
    derepR <- derepFastq(filtRs[[i]])
    dadaR <- dada(derepR, err=errR, multithread=TRUE)
    ddR[[i]] <- dadaR
    merger <- mergePairs(ddF[[i]], derepF, ddR[[i]], derepR)
    mergers[[i]] <- merger
}
rm(derepF); rm(derepR)
print("make seqtab")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
save(seqtab, file="/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_seqtab.RData") 
print("merged seqtab sequence length and median")
table(nchar(getSequences(seqtab)))
med.seqtab <- median(nchar(getSequences(seqtab))); med.seqtab
print("remove chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
print("frequency of chimeric sequences")
sum(seqtab.nochim)/sum(seqtab)
print("dim of seqtab.nochim")
dim(seqtab.nochim)
save(seqtab.nochim, file="/project/mpm_fusarium_genomics/whitaker/niko19/R2Fung_seqtab.nochim.RData")
print("seqtab.nochim sequence length and median")
table(nchar(getSequences(seqtab.nochim)))
med.seqtab.nochim <- median(nchar(getSequences(seqtab.nochim))); med.seqtab.nochim
################################################################################


sessionInfo()
#end
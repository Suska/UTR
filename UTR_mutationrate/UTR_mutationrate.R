#Dr. Susanne Weller
#13/09/2014

#####################################
#UTR - abnormal mutation rate finder#
#####################################

#The aim is to find UTRs that are significantly more mutated compared to a baseline.
#Task 1: Find UTRs that are not related to cancer (via filter of xx) 
#Task 2: Use these to create a baseline mutation rate (with range)??(in mutations/bp)
# -> correct number of mutations by length of UTR
# -> for within patient mutation rate as well as over all patients
# -> correct this for gene expression (similar to MutSigCV)
#Task 3: Compare these to actual mutation rate of UTRs (again corrected for their expression)
#-------------------------------------------------------------------------------------------------------------------------------------
cll3 <- read.table("/home/suska/work/07_UTR/UTR/trunk/CLL3primeutr.txt", sep="\t", header=TRUE, strip.white = TRUE, na.strings=c("", "Sequence_Unavailable" ))
cll5 <- read.table("/home/suska/work/07_UTR/UTR/trunk/CLL5primeutr.txt", sep="\t", header=TRUE, strip.white = TRUE, na.strings=c("", "Sequence_Unavailable" ))

# count length of sequence and check
cll3$UCSClength <-nchar(as.character(cll3$X3_Prime_UTR_Sequence_UCSC))
cll3$UCSC_Seq_Length == cll3$UCSC_Seq_Length

cll5$UCSClength <-nchar(as.character(cll5$X5_Prime_UTR_Sequence_UCSC))
cll5$UCSC_Seq_Length == cll5$UCSC_Seq_Length

# reduce dataframe
cll3$X3_Prime_UTR_Sequence_BioMart <- NULL
cll3$X3_Prime_UTR_Sequence_UCSC <- NULL

cll5$X5_Prime_UTR_Sequence_BioMart <- NULL
cll5$X5_Prime_UTR_Sequence_UCSC <- NULL

# calculate number of mutations and number of samples affected per bp
cll3$mutations_per_bp <- cll3$Mutations_Number / cll3$UCSClength
cll3$samplesaffected_per_bp <- cll3$No_Samples_Affected / cll3$UCSClength

cll5$mutations_per_bp <- cll5$Mutations_Number / cll5$UCSClength
cll5$samplesaffected_per_bp <- cll5$No_Samples_Affected / cll5$UCSClength
#---------------------------------------------------------------------------------------------------------------------------------------
#Task 1:Find UTRs that are not related to cancer
cancerdirect <- read.table("/home/suska/work/07_UTR/UTR/trunk/UTR_noncancergenes/UTR_Genes_Cancer_Direct_Effect_Excluded.txt",
                           sep="\t", header=TRUE, quote="")
UTRfinder <- function(x){
  if(grepl("3'UTR", x)==TRUE){
  return("3UTR")
  }
  if(grepl("5'UTR", x)==TRUE){
  return("5UTR")  
  }
}
cancerdirect$UTR <- unlist(lapply(cancerdirect$Gene.Region, UTRfinder))

#separate UTRs:
cancerdirect3 <- subset(cancerdirect, cancerdirect$UTR=="3UTR")
#write.csv( cancerdirect3,"/home/suska/work/07_UTR/UTR/trunk/cancerdirect3prime.txt")
#write.csv( mergecancerdirect3,"/home/suska/work/07_UTR/UTR/trunk/cancerdirect3prime_merge.txt")
cancerdirect5 <- subset(cancerdirect, cancerdirect$UTR=="5UTR")

#get UTR length for cancer-unrelated UTRs:
mergecancerdirect3 <- merge(cancerdirect3, cll3, by.x="Gene.Symbol", by.y="Gene_ID")
mergecancerdirect5 <- merge(cancerdirect5, cll5, by.x="Gene.Symbol", by.y="RefSeq_Transcript_ID")

#get gene expression value from MutSigCV results:
mutsigcv <-read.table("/home/suska/work/02_CLLpilot/CLLpilot/trunk/MutSigCV_results/CLL_A1_VC1_VF02_S42_mutsigcv.sig_genes.txt",
           sep="\t", header=TRUE)
mergecancerdirect3 <- merge(mergecancerdirect3, mutsigcv, by.x="Gene.Symbol", by.y="gene")
mergecancerdirect5 <- merge(mergecancerdirect5, mutsigcv, by.x="Gene.Symbol", by.y="gene")

#clean up for clarity:
cancerdirect3 <- subset(mergecancerdirect3, ,c(Gene.Symbol, Chromosome, Position, Variation.Type,
                                              Case.Samples.With.Variant,
                                              UTR,
                                              Mutations_Number,
                                              No_Samples_Affected,
                                              UCSClength,
                                              mutations_per_bp,
                                              samplesaffected_per_bp,
                                              expr))
cancerdirect5 <- subset(mergecancerdirect5, ,c(Gene.Symbol, Chromosome, Position, Variation.Type,
                                               Case.Samples.With.Variant,
                                               UTR,
                                               Mutations_Number,
                                               No_Samples_Affected,
                                               UCSClength,
                                               mutations_per_bp,
                                               samplesaffected_per_bp,
                                               expr))

#-------------------------------------------------------------------------------------------------------------
##Task 2: Use these to create a baseline mutation rate (with range)??(in mutations/bp)
# -> correct number of mutations by length of UTR
# -> for within patient mutation rate as well as over all patients
# -> correct this for gene expression (similar to MutSigCV)

#correct by length of utr:
cancerdirect3$mutations_per_bp_utr <- cancerdirect3$Case.Samples.With.Variant / cancerdirect3$UCSClength
hist(cancerdirect3$mutations_per_bp_utr, breaks=20)
boxplot(cancerdirect3$mutations_per_bp_utr)
par(mfrow=c(1,1))
boxplot(cll3$mutations_per_bp, ylim=c(0, 0.015), main="all utr mutations")
boxplot(cancerdirect3$mutations_per_bp_utr, ylim=c(0, 0.015),main="cancer direct excluded")
summary(cancerdirect3$mutations_per_bp_utr)
summary(cll3$mutations_per_bp)

#correct for gene expression:
summary(cancerdirect3$expr)
hist(cancerdirect3$expr, breaks=40)

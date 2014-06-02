#Dr. Susanne Weller
#28/05/2014

# First look at UTRs in CLL pilot project
# - Get length of UTR
# - correct no of mutations for length of UTR
# - Combine 5` and 3` results

setwd("~/work/07_UTR/UTR/trunk")

# load both data tables from Adam
cll3 <- read.table("CLL3primeutr.txt", sep="\t", header=TRUE, strip.white = TRUE, na.strings=c("", "Sequence_Unavailable" ))
cll5 <- read.table("CLL5primeutr.txt", sep="\t", header=TRUE, strip.white = TRUE, na.strings=c("", "Sequence_Unavailable" ))

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

# plot these against each other

plot(cll3$mutations_per_bp, cll3$samplesaffected_per_bp)
points(cll5$mutations_per_bp, cll5$samplesaffected_per_bp, col="red", pch=4)

# take out everything below 0.002

subcll3 <- subset(cll3, cll3$mutations_per_bp > 0.002)
subcll5 <- subset(cll5, cll5$mutations_per_bp > 0.002)

# plot with ggplot
library(ggplot2)
ggplot(cll3, aes(x = cll3$mutations_per_bp, y = cll3$samplesaffected_per_bp)) + 
  geom_point() + 
  geom_text(aes(label = cll3$Gene_ID), vjust = 0.9)+
  scale_x_log10()+
  scale_y_log10()
  geom_line(data = river, aes(y = lat))ggplot


ggplot(mapping=aes(x=mutations_per_bp, y=samplesaffected_per_bp)) + 
  geom_point(data=subcll3, aes(x=subcll3$mutations_per_bp, y=subcll3$samplesaffected_per_bp), col="darkblue") + 
  geom_text(data=subcll3, aes(label = subcll3$Gene_ID), vjust = 0.9, hjust = 1.4, size=3, col="darkblue")+             
  geom_point(data=subcll5, aes(x=subcll5$mutations_per_bp, y=subcll5$samplesaffected_per_bp), pch=4, col="darkred")+
  geom_text(data=subcll5, aes(label = subcll5$RefSeq_Transcript_ID), vjust = 0.9, hjust = -0.1, size=3, col="darkred")+
  scale_x_log10()+
  scale_y_log10()

#Using multiplot function
source("multiplotfunction.R")

plot_3prime <- ggplot(mapping=aes(x=mutations_per_bp, y=samplesaffected_per_bp)) + 
  geom_point(data=subcll3, aes(x=subcll3$mutations_per_bp, y=subcll3$samplesaffected_per_bp), col="darkblue") + 
  geom_text(data=subcll3, aes(label = subcll3$Gene_ID), vjust = 0.9, hjust = 1.4, size=2.5, col="darkblue")+             
  scale_x_log10(limits = c(0.002,0.015))+
  scale_y_log10()
  

plot_5prime <- ggplot(mapping=aes(x=mutations_per_bp, y=samplesaffected_per_bp)) + 
  geom_point(data=subcll5, aes(x=subcll5$mutations_per_bp, y=subcll5$samplesaffected_per_bp), col="darkred")+
  geom_text(data=subcll5, aes(label = subcll5$RefSeq_Transcript_ID), vjust = 0.9, 
            hjust = -0.1, size=2.5, col="darkred", position=position_jitter(h=0.01,w=0.01))+
  scale_x_log10(limits = c(0.002,0.015))+
  scale_y_log10()


multiplot(plot_3prime, plot_5prime, cols=1)


BioMartdetector <- function(old, new){
if (old != new){
    return(c(old, new, (old-new)))
}
}

BioMarttest1 <- BioMartdetector(5, 8)


BioMartquality <- sapply(cll3,MARGIN=1, FUN=BioMartdetector, old=cll3$BioMart_Seq_Length, new=cll3$Biomartlength)
BioMartquality <- as.data.frame(apply(cll3, MARGIN=1, function (x, y) BioMartdetector(cll3$BioMart_Seq_Length, cll3$Biomartlength)))



BioMartquality <- do.call(BioMartdetector, old="cll3$BioMart_Seq_Length", new="cll3$Biomartlength")
BioMartquality <- do.call(BioMartdetector, old=1, new=4, cll3 )

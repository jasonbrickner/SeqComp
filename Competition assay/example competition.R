#This script finds the peak height of T and G from input ab files, calculates a ratio of the strains and normalizes based on the initial ratio 
#and corrects for initial densities not equal to 1.
library(plyr)
library(tidyverse)  
library(sangerseqR)
library(stringr)

#draw peak heights and calculate normalized ratios for each strain, condition and time point.
#Step 1: collect data for WT cells act vs react
fnames <- as.character(list.files("~/Downloads/Competition assay/example competition raw data/"))
WT_LUT <- read.csv("~/Downloads/Competition assay/LUT for example competition.csv")
WT_LUT$file <- paste(WT_LUT$file_number, ".ab1", sep = "")
SNPdf <- data.frame(file = character(),
                      A = numeric(),
                      C = numeric(),
                      CtoA_ratio = numeric())

for(i in 1:length(fnames)){
    SNPdf[i,1] <- fnames[i]
    sub_LUT <- subset(WT_LUT, file == fnames[i])
    df <- read.abif(paste("~/Downloads/Competition assay/example competition raw data/", fnames[i], sep = ""))@data
    seq <- sangerseq(read.abif(paste("~/Downloads/Competition assay/example competition raw data/", fnames[i], sep = "")))
    position <- as.numeric(str_locate_all(pattern ="GGGTTTTCCC", primarySeq(seq, string = "TRUE"))[[1]][1,2])+1
    window <- round((df$PLOC.2[position]-df$PLOC.2[position-1])/4, 0)
    tracestart <- df$PLOC.2[(position)]-window
    traceend <- df$PLOC.2[(position)]+window
    SNPdf[i, 2] <- max(seq@traceMatrix[c(tracestart:traceend), 1])
    SNPdf[i, 3] <- max(seq@traceMatrix[c(tracestart:traceend), 2])
    SNPdf[i, 4] <- round(SNPdf[i, 3]/SNPdf[i, 2], 2)
}
  SNP_final <- merge(WT_LUT[,-c(1:2)], SNPdf, by = "file")[,-1]
  SNP_sum1 <- ddply(SNP_final, .(A_strain, C_strain, time, bio_rep), summarize,
                   tech_ratio = mean(CtoA_ratio))
  zeros <- subset(SNP_sum1, time == 0)[,-3] %>% dplyr::rename(zero = tech_ratio)
  SNP_norm <- merge(SNP_sum1, zeros, by = c("A_strain", "C_strain", "bio_rep")) %>% dplyr::mutate(normalized_ratio = tech_ratio/zero)

  SNP_sum2 <- ddply(SNP_norm, .(A_strain, C_strain, time), summarize,
                    mean_ratio = mean(normalized_ratio),
                    sd_ratio = sd(normalized_ratio),
                    se_ratio = sd_ratio/sqrt(length(bio_rep)))
  
SNP_sum2$ID <- paste(SNP_sum2$C_strain, SNP_sum2$A_strain, SNP_sum2$time, sep = "_")

  


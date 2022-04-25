#This script finds the peak height of T and G from input ab files, calculates a ratio of the strains and normalizes based on the initial ratio 
#and corrects for initial densities not equal to 1. The script also plots these ratios to generate Figure 1D.
#Required packages: 
#stringr (CRAN) and sangerseqR (https://bioconductor.org/packages/release/bioc/html/sangerseqR.html)
#BiocManager::install("sangerseqR")

library(plyr)
library(tidyverse)  
library(sangerseqR)
library(stringr)

#draw peak heights and calculate normalized ratios for each strain, condition and time point.
#Step 1: collect data for WT cells act vs react
fnames <- as.character(list.files("~/Box/ Jason/papers in preparation/Bethany Memory/Figure 1/Fitness pRS306 WT Standard Curve/Raw Data/"))
std_LUT <- read.csv("~/Box/ Jason/papers in preparation/Bethany Memory/Figure 1/Fitness pRS306 WT Standard Curve/Scripts/StdCurve LUT for paper.csv")
std_LUT$file <- paste(std_LUT$filename, ".ab1", sep = "")
SNPdf <- data.frame(file = character(),
                      A = numeric(),
                      C = numeric(),
                      measured_percent_A = numeric())

for(i in 1:length(fnames)){
    SNPdf[i,1] <- fnames[i]
    sub_LUT <- subset(std_LUT, file == fnames[i])
    df <- read.abif(paste("~/Box/ Jason/papers in preparation/Bethany Memory/Figure 1/Fitness pRS306 WT Standard Curve/Raw Data/", fnames[i], sep = ""))@data
    seq <- sangerseq(read.abif(paste("~/Box/ Jason/papers in preparation/Bethany Memory/Figure 1/Fitness pRS306 WT Standard Curve/Raw Data/", fnames[i], sep = "")))
    position <- as.numeric(str_locate_all(pattern ="GGGTTTTCCC", primarySeq(seq, string = "TRUE"))[[1]][1,2])+1
    window <- round((df$PLOC.2[position]-df$PLOC.2[position-1])/4, 0)
    tracestart <- df$PLOC.2[(position)]-window
    traceend <- df$PLOC.2[(position)]+window
    SNPdf[i, 2] <- max(seq@traceMatrix[c(tracestart:traceend), 1])
    SNPdf[i, 3] <- max(seq@traceMatrix[c(tracestart:traceend), 2])
    SNPdf[i, 4] <- 100*round(SNPdf[i, 2]/(SNPdf[i, 3]+SNPdf[i, 2]), 2)
}
  SNP_final <- merge(std_LUT[,-1], SNPdf, by = "file")[,-1]
  SNP_sum <- ddply(SNP_final, .(percent_A), summarize,
                   mean_measured = mean(measured_percent_A),
                   sd_measured = sd(measured_percent_A))
  
  sc <- lm(SNP_sum$mean_measured~SNP_sum$percent_A)
  SNP_sum$fit <- sc$coefficients[1]+SNP_sum$percent_A*sc$coefficients[2]
  Rsq <- (cor(SNP_sum$mean_measured, SNP_sum$percent_A))^2
  
standard_cuvrve <-   ggplot()+
  geom_line(data = SNP_sum, aes(x = percent_A, y = fit), col = "blue", linetype = 2, size = 1.5, alpha = 0.4)+
  geom_errorbar(data = SNP_sum, aes(x = percent_A,ymax = mean_measured+sd_measured, ymin = mean_measured-sd_measured), col = "#333333", size = .55, width = 3, alpha = 0.7)+
  geom_point(data = SNP_sum, aes(x = percent_A, y= mean_measured), col = "black", fill = "white", shape = 21, size = 3, alpha = 0.5)+
  theme_bw()+
  ylab("measured percent A")+
  xlab("predicted percent A")+
  theme(text = element_text(family = "Arial", size = 12, color = "#333333"),
          axis.text = element_text(family = "Arial", size = 12, color = "#333333"),
          legend.position = "none")

ggsave("~/Box/ Jason/papers in preparation/Bethany Memory/Figure 1/panels/standard curve.png", standard_cuvrve, width = 3, height = 3)

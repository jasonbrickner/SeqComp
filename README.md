# SeqComp
Sanger sequencing-based quantification of relative abundance of two strains with a single SNP.

These two R scripts will analyze the accompanying ab files to either generate a standard curve (from mixing strains in defined ratios before preparing DNA, PCR amplifying the relevant region of the genome and then subjecting to Sanger sequencing) or generate relative abundances from a competition experiment.  Downloading the folder to your /Downloads/ folder, they should be directed to be able to run.  If you move the folder to another location, update the setwd command to reflect the new location.

The script searches for a polymorphism that was introduced into the pRS303 plasmid to change an A within the multiple cloning site to a C, which is downstream of the sequence GGGTTTTCCC.  To search for a different polymorphism, update this search sequence on line 25.  Here is how you can quantify the abundance of each of the four possible bases on line 29 and 30:
A:  SNPdf[i, 2] <- max(seq@traceMatrix[c(tracestart:traceend), 1])
C:  SNPdf[i, 2] <- max(seq@traceMatrix[c(tracestart:traceend), 2])
G:  SNPdf[i, 2] <- max(seq@traceMatrix[c(tracestart:traceend), 3])
T:  SNPdf[i, 2] <- max(seq@traceMatrix[c(tracestart:traceend), 4])

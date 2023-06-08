#=====rxCOV based metric for establishing analyte fidelity=====
#Shikhar Uttam
#Dept of Computational and Systems Biology
#UPMC Hillman Cancer Center
#University of Pittsburgh
rm(list=ls());
graphics.off()
library(ggplot2);
library(ggpubr)
library(reshape2)
library(umap)
library(openxlsx)
library(Matrix)
library(glmnet)
library(splitstackshape)
library(zoo)

load('./rxCOVdata.RData')
#=======================
#Noise Estimation
#======================
absdiff.MSD <- abs(df.one - df.two);

# ##Calculate phase diff averages and SD##
phase2diff.av.MSD <- apply(absdiff.MSD, 2, mean);
phase2diff.sd.MSD <- apply(absdiff.MSD, 2, sd);

#============================
#Differential signal Estimation
#============================
##Mean 
#AC
AC.av <- apply(AC.visit, 2, mean);
#UC
UC.av <- apply(UC.visit, 2, mean);
#AL
AL.av <- apply(AL.visit, 2, mean);
#UL
UL.av <- apply(UL.visit, 2, mean);
#L
L.av <- apply(rbind(AL.av,UL.av), 2, mean);

##Standard deviation
#AC
AC.sd <- apply(AC.visit, 2, sd);
#UC
UC.sd <- apply(UC.visit, 2, sd);
#AL
AL.sd <- apply(AL.visit, 2, sd);
#UL
UL.sd <- apply(UL.visit, 2, sd);
#L
L.sd <- apply(rbind(AL.av,UL.av), 2, sd);

##Take the average differences between groups##
ACUCdiff.av <- abs(AC.av - UC.av)#/(AC.av+UC.av);
ACALdiff.av <- abs(AC.av - AL.av)#/(AC.av+AL.av);
ACULdiff.av <- abs(AC.av - UL.av)#/(AC.av+UL.av);
UCALdiff.av <- abs(UC.av - AL.av)#/(UC.av+AL.av)
ULALdiff.av <- abs(UL.av - AL.av)#/(UL.av+AL.av)
UCULdiff.av <- abs(UC.av - UL.av)#/(UC.av+UL.av)
ACLdiff.av <- abs(AC.av - L.av)#/(AC.av+L.av)
UCLdiff.av <- abs(UC.av - L.av)#/(UC.av+L.av)



##Take the standard deviations between groups##
ACUCdiff.sd <- sqrt(AC.sd^2 + UC.sd^2)
ACALdiff.sd <- sqrt(AC.sd^2 + AL.sd^2)
ACULdiff.sd <- sqrt(AC.sd^2 + UL.sd^2)
UCALdiff.sd <- sqrt(UC.sd^2 + AL.sd^2)
ULALdiff.sd <- sqrt(UL.sd^2 + AL.sd^2)
UCULdiff.sd <- sqrt(UC.sd^2 + UL.sd^2)
ACLdiff.sd <- sqrt(AC.sd^2 + L.sd^2)
UCLdiff.sd <- sqrt(UC.sd^2 + L.sd^2)


#=====coefficient of variation====#
diff.cv.df <- as.data.frame(rbind(ACUCdiff.av*ACUCdiff.sd,
                                   ACALdiff.av*ACALdiff.sd,
                                   ACULdiff.av* ACULdiff.sd,
                                   UCALdiff.av*UCALdiff.sd,
                                   ULALdiff.av*ULALdiff.sd,
                                   UCULdiff.av*UCULdiff.sd,
                                   ACLdiff.av*ACLdiff.sd,
                                   UCLdiff.av*UCLdiff.sd));
rownames(diff.cv.df) <- c('AC vs UC', 'AC vs AL','AC vs UL','UC vs AL','UL vs AL','UC vs UL','AC vs L','UC vs L');


variability.cv <- c(phase2diff.av.MSD*phase2diff.sd.MSD); #++++++++ variability +++++++++
variability.cv.df <- as.data.frame(rbind(variability.cv,
                                         variability.cv,
                                         variability.cv,
                                         variability.cv,
                                         variability.cv,
                                         variability.cv,
                                         variability.cv,
                                         variability.cv))

rownames(variability.cv.df) <- c('AC vs UC', 'AC vs AL','AC vs UL','UC vs AL','UL vs AL','UC vs UL','AC vs L','UC vs L');

ratio.cv.df <- data.frame(log2((diff.cv.df+.Machine$double.eps)/(variability.cv.df+.Machine$double.eps)))
ratio.cv.df <- ratio.cv.df[c(6,5,2),]

#
rxCOV <- as.numeric(ratio.cv.df[1,])
par(mar = c(6.1, 4.1, 4.1, 2.1))
rxCOV.indicator <- as.numeric(rxCOV > 0) + 1;
cols <- c("red","blue")
barplot(rxCOV, names.arg = colnames(ratio.cv.df),las=3,col = cols[rxCOV.indicator], cex.names = 0.8,font.axis=2,
        ylab='rxCOV',font.lab=2)
abline(h=0,col='black')

rxCOV <- as.numeric(ratio.cv.df[2,])
par(mar = c(6.1, 4.1, 4.1, 2.1))
rxCOV.indicator <- as.numeric(rxCOV > 0) + 1;
cols <- c("red","blue")
barplot(rxCOV, names.arg = colnames(ratio.cv.df),las=3,col = cols[rxCOV.indicator], cex.names = 0.8,font.axis=2,
        ylab='rxCOV',font.lab=2)
abline(h=0,col='black')

rxCOV <- as.numeric(ratio.cv.df[3,])
par(mar = c(6.1, 4.1, 4.1, 2.1))
rxCOV.indicator <- as.numeric(rxCOV > 0) + 1;
cols <- c("red","blue")
barplot(rxCOV, names.arg = colnames(ratio.cv.df),las=3,col = cols[rxCOV.indicator], cex.names = 0.8,font.axis=2,
        ylab='rxCOVe',font.lab=2)
abline(h=0,col='black')

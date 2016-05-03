#This software is a rerelease of CyberTOR_WY.R module of Microarray Expression
#Potential Suite with alterations in the input and ouput format in order to
#improve transition between modules of eMESS.

#CyberT algorithm from:
#CyberT for use with spotted microaray data based on the derivation in
#"Improved Statistical Inference from DNA Microarray Data Using Analysis of
#Variance and Bayesian Statistical Framework" A. D. Long et al, J. Bio. Chem.,
#276(22), pp 19937-19944, 2001.
#Resampling-based FDR control of Benjamini-Hochberg (third alternative, p371)
#for use with microarray data based on the derivation in "Identifying
#differentially expressed genes using false discovery rate controlling
#procedures", A. Reiner et al, Bioinformatics, 19(3), pp 368-375, 2003.




#No part of this text or software may be used or reproduced in any matter
#whatsoever without written permission, except in the case of brief quotations
#embodied in critical articles of or for reviews.  For information address
#Brian W. Kirk (bwkbem@yahoo.com).
#Software Copyright (C) October 2002  Brian W. Kirk
#Form/Text Copyright (C) October 2002  Brian W. Kirk
#Released for use with R v1.5.1 (C)

#Copyright (C) October 2002  Brian W. Kirk
#Released for R v1.5.1 (C)

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

###############################################################################
#                                FUNCTIONS                                    #
###############################################################################

library(MASS);

#FUNCTION THAT CALCULATES FOR A VECTOR OF SIGNALS A ROBUST SUMMARY USING THE
#HUBER M-ESTIMATOR AND THE MAD ESTIMATOR FOR SCALE.  RESULTS ARE COMBINED INTO
#A DATA.FRAME, WITH INCOMPATIBLE SETS (IE WITH ONLY 1 ENTRY) AMENDED TO THE END
#OF THE DATA.FRAME.

  fhuber.summarize <- function(SignalVec=vec, Filter=filter,
                               NameIdentifier=name, HuberConf=k, HuberTol=tol)
  {
   #Prime Header of Output Table
   huber.report.template <- c("ProbeID", "mu", "s", "probes");
   
   #Determine the Number of Probes Present for each Probeset 
   probes <- NameIdentifier[Filter];
   probesetf <- factor(probes);
   probesets <- levels(probesetf);

   probe.number <- as.vector(tapply(probes, probesetf, length));
   
   #probe.number.table <- data.frame(cbind(probesets, probe.number))

   #Select Porbesets With More Than 1 Probe For Huber Analysis
   huber.ok <- probesets[which(as.numeric(
                                     as.vector(probe.number)) > 1)];

   #Apply Estimators To Qualified Probesets
   if (length(huber.ok) > 0) {
     for (i in (1:length(huber.ok))) {

       #Test That More Than Half of the Probes Are Not Identical
       test.probe <- SignalVec[as.logical(
                        match(NameIdentifier, huber.ok[i], nomatch=0)*Filter)];
     
       test.probef <- factor(test.probe);
       test.probes <- levels(test.probef);
       level.lengths <- tapply(test.probe, test.probef, length);

       continue <- TRUE;
       for (j in 1:length(test.probes)) {
         if ((length(test.probe) - level.lengths[j]) < 0.5*length(test.probe)){
           continue <- FALSE
         }
       }

       #Process Accordingly Based on the Previous Test
       if (continue) {
         xtemp.matrix <-  as.matrix(huber(test.probe, k=HuberConf,
                                          tol=HuberTol));
         huber.report.template <- cbind(huber.report.template,
                        c(as.character(huber.ok[i]),
                          as.numeric(xtemp.matrix[1,1]),
                          as.numeric(xtemp.matrix[2,1]),
                          as.numeric(as.vector(
                              probe.number[as.logical(match(
                              probesets, huber.ok[i], nomatch=0))]))));
       } else {
         mu <- median(test.probe);
         scale <- sd(test.probe);

         huber.report.template <- cbind(huber.report.template,
                        c(as.character(huber.ok[i]),
                          mu,
                          scale,
                          as.numeric(as.vector(
                              probe.number[as.logical(match(
                              probesets, huber.ok[i], nomatch=0))]))));
   
       } 
     }
   }
   
   #Organize Porbesets That Did Not Qualify For Estimator Analysis
   huber.lost <- probesets[which(as.numeric(
                                    as.vector(probe.number)) <= 1)];
   if (length(huber.lost) > 0) {
     for (i in (1:length(huber.lost))) {
       huber.report.template <- cbind(huber.report.template,
                            c(as.character(huber.lost[i]),
                              as.numeric(SignalVec
                                          [as.logical(match(NameIdentifier,
                                           huber.lost[i], nomatch=0)
                                           * Filter)]), NA, 1))
     }
   }

   #Combine Both Probeset Types Into A Single Report
   huber.report.template2 <- t(huber.report.template);
   huber.report.template3 <-
                        huber.report.template2[2:(length(probesets)+1),];
   huber.sort.filter <- order(as.character(huber.report.template3[,1]));
   
   return(data.frame(cbind(
                  ProbeID=
                    as.character(huber.report.template3[huber.sort.filter, 1]),
                  mu=
                    as.numeric(huber.report.template3[huber.sort.filter, 2]),
                  s=
                    as.numeric(huber.report.template3[huber.sort.filter, 3]),
                  probes=
                    as.numeric(
                    as.vector(huber.report.template3[huber.sort.filter, 4])))))
  }


###############################################################################
#                           CONSTANT DECLARATION                              #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x, y, or z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANENCY

                 ############################################
                 #######  DEFINE PERMANENT CONSTANTS  #######
                 ############################################

#Define Output Dirctory with Prefix

if (MatchedArrays) {
    OutputPrefix <- paste(OutputDirectory, "m", OutputFile,"A",
                          as.character(MinArrays), "P", as.character(MinProbes),
                          sep="");
} else {
    OutputPrefix <- paste(OutputDirectory, "u", OutputFile,"A",
                          as.character(MinArrays), "P", as.character(MinProbes),
                          sep="");
}




#DEFINE GRAPHICS WINDOW SIZE
WinWidth <- 8.5;           
WinHeight <- 11;


             #####################################################
             ##########  DEFINE TEMPORARY CONSTANTS  #############
             #####################################################


#Column position in LNProbesetDatawBGSubtract.RData
xIDPos <- 1;
MuFoldPos <- 2;
xMuAmpPos <- 3;
xProbeNumberPos <- 4;


###############################################################################
#                             VARIABLE DECLARATION                            #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x", "y", or "z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANCY



#Load LNProbesetData for each array into a list of data.frames
TotArray <- length(ProbesetFiles);
ArrayList <- vector("list", TotArray);
m <- 0;
for (Array in ProbesetFiles) {
  m <- m + 1;
  zNormalizedData <- paste(InputDirectory, Array, sep="");
  ArrayList[[m]] <- data.frame(read.delim(zNormalizedData, skip=HeaderLength));
}

#Concatenate respective vectors from each element of the ArrayList
IDVec <- as.character(as.vector(
                             ArrayList[[1]][,xIDPos]));
LogFoldVec <- as.numeric(as.vector(
                          ArrayList[[1]][,MuFoldPos]));
zRMSIntensityVec <- as.numeric(as.vector(
                               ArrayList[[1]][,xMuAmpPos]));
zProbeNumberVec  <- as.numeric(as.vector(
                               ArrayList[[1]][,xProbeNumberPos]));

ArrayID <- rep(1, length(IDVec));

for (n in 2:TotArray) {
  
  tempIDVec <- as.character(as.vector(
                                 ArrayList[[n]][,xIDPos]));
  
  IDVec <- c(IDVec, tempIDVec);
  
  LogFoldVec <- c(LogFoldVec,
                   as.numeric(as.vector(
                              ArrayList[[n]][,MuFoldPos])));
  
  zRMSIntensityVec <- c(zRMSIntensityVec,
                        as.numeric(as.vector(
                                   ArrayList[[n]][,xMuAmpPos])));
  
  zProbeNumberVec  <- c(zProbeNumberVec,
                       as.numeric(as.vector(
                                  ArrayList[[n]][,xProbeNumberPos])));
  
  ArrayID <- c(ArrayID, rep(n, length(tempIDVec)));
}

#Filter for data that passes the minimum number of probes required
ProbeFilter <- as.logical(zProbeNumberVec >= MinProbes);

IDVec <- IDVec[ProbeFilter];
LogFoldVec <- LogFoldVec[ProbeFilter];
zRMSIntensityVec <- zRMSIntensityVec[ProbeFilter];
zProbeNumberVec <- zProbeNumberVec[ProbeFilter];
ArrayID <- ArrayID[ProbeFilter];
  
#Filter for matched pairs ArrayIDs 1&2 3&4 5&6 7&8

if (MatchedArrays) {
    MatchArrayFilter <- 0;
    for (n in 2*(1:(TotArray/2))) {
        MatchArrayFilter <- c(MatchArrayFilter,
                        as.logical(match(IDVec[which(ArrayID==n-1)],
                                      IDVec[which(ArrayID==n)], nomatch=0)));

        MatchArrayFilter <- c(MatchArrayFilter,
                        as.logical(match(IDVec[which(ArrayID==n)],
                                      IDVec[which(ArrayID==n-1)], nomatch=0)));

    }

    MatchArrayFilter <- as.logical(MatchArrayFilter[2:length(MatchArrayFilter)]);

    IDVec <- IDVec[MatchArrayFilter];
    LogFoldVec <- LogFoldVec[MatchArrayFilter];
    zRMSIntensityVec <- zRMSIntensityVec[MatchArrayFilter];
    zProbeNumberVec <- zProbeNumberVec[MatchArrayFilter];
    ArrayID <- ArrayID[MatchArrayFilter];

}



#Input Ordered Genes from Annotated Genome
GenomeDataIn <- paste(GenomeDirectory, GenomeFile, sep="");
xGenomeData <- data.frame(read.delim(GenomeDataIn,
                                       skip=GenomeHeaderLength));
GenomeGenes <- as.character(as.vector(xGenomeData[,GenomeIDPos]));

rm(list=ls(pat="^x"));
###############################################################################
#                               MAIN BODY                                     #
###############################################################################

                 #############################################
                 ####     Apply huberM/MAD Estimators     #### 
                 #############################################

DummyFilter <- rep(TRUE, length(IDVec));

zLogFoldHuberReport <- fhuber.summarize(LogFoldVec, DummyFilter,
                                IDVec, HuberConf, HuberTol);

zRMSHuberReport <- fhuber.summarize(zRMSIntensityVec, DummyFilter,
                                IDVec, HuberConf, HuberTol);


#Match Huber Reports for Intensity data and LogFold Change Data
#(ie if gene is present in only one, then remove it)
zRMSHuberSetFilter <- !as.logical(
                             match(zRMSHuberReport[,3], c(NA, 0), nomatch=0) |
                             as.numeric(as.vector(zRMSHuberReport[,4])) <
                                        MinArrays);

zLogFoldHuberSetFilter <- !as.logical(
                         match(zLogFoldHuberReport[,3], c(NA, 0), nomatch=0) |
                         as.numeric(as.vector(zLogFoldHuberReport[,4])) <
                                    MinArrays);

zTotHuberSetFilter <- as.logical(zRMSHuberSetFilter * zLogFoldHuberSetFilter);

HuberGenes <- as.character(zLogFoldHuberReport[zTotHuberSetFilter, 1]);
            
Intensity <- as.numeric(as.vector(zRMSHuberReport[zTotHuberSetFilter, 2]));

LogFoldChange <- as.numeric(as.vector(
                                  zLogFoldHuberReport[zTotHuberSetFilter, 2]));

LogFoldStndDev <- as.numeric(as.vector(
                                  zLogFoldHuberReport[zTotHuberSetFilter, 3]));

TotArrays <- as.numeric(as.vector(zRMSHuberReport[zTotHuberSetFilter, 4]));


if (MatchedArrays) {

    TotBiosamples <- TotArrays/2;

} else {

    TotBiosamples <- 0;
    for (i in 1:length(HuberGenes)) {
        
        BiosampleFilter <- as.logical(match(IDVec, HuberGenes[i], nomatch=0));
        PickedArrayID <- ArrayID[BiosampleFilter];

        count <- 0;
        Flag1 <- FALSE;
        Flag2 <- FALSE;
        Flag3 <- FALSE;
        Flag4 <- FALSE;
        for (j in 1:length(PickedArrayID)) {

            if (!Flag1 & (PickedArrayID[j] == 1|PickedArrayID[j] == 2)) { 
                count <- count + 1;
                Flag1 <- TRUE;
            }    
            if (!Flag2 & (PickedArrayID[j] == 3|PickedArrayID[j] == 4)) {
                count <- count + 1;
                Flag2 <- TRUE;
            }
            if (!Flag3 & (PickedArrayID[j] == 5|PickedArrayID[j] == 6)) {
                count <- count + 1;
                Flag3 <- TRUE;
            } 
            if (!Flag4 & (PickedArrayID[j] == 7|PickedArrayID[j] == 8)) {
                count <- count + 1;
                Flag4 <- TRUE;
            }

        }

        TotBiosamples <- c(TotBiosamples, count);
    }
    TotBiosamples <- TotBiosamples[2:length(TotBiosamples)];
}


StndErr <- LogFoldStndDev/sqrt(TotBiosamples);

### Remove non Gene Probes and unrepresented genes in annotation
GenomeFilter <- as.logical(as.vector(
                             match(HuberGenes, GenomeGenes, nomatch=0)));

HuberGenes <- HuberGenes[GenomeFilter];
            
Intensity <- Intensity[GenomeFilter];

LogFoldChange <- LogFoldChange[GenomeFilter];

LogFoldStndDev <- LogFoldStndDev[GenomeFilter];

TotArrays <- TotArrays[GenomeFilter];

TotBiosamples <- TotBiosamples[GenomeFilter];

StndErr <- StndErr[GenomeFilter];


rm(list=ls(pat="^z"));
rm(list=ls(pat="^f"));
rm(list=ls(pat="^g"));
      
                 ##############################################
                 ####    CALCULATE ADJ STANDRD DEVIATION   #### 
                 ##############################################

# LOWESS fit data and organize for extraction of LOWESS values
xBGSigmaLowess <- lowess(Intensity, LogFoldStndDev, f=LowessfParam);

BGStDev <- 0;
xLowessFilter <- order(Intensity);
for (j in 1:length(xLowessFilter)) {
  tempBGStDev <- xBGSigmaLowess$y[as.logical(
                                      match(xLowessFilter, j, nomatch=0))]
  BGStDev <- c(BGStDev, tempBGStDev);
}
                                
BGStDev <- BGStDev[2:(length(BGStDev))]


#Calculate Adjusted Stndard Err
#nu is defined to result in a constant lambda for all genes
nu <- 0;
xAdjustedVar <- 0;
for (k in 1:length(HuberGenes)) {
  if (TotBiosamples[k] >= Lambda) {
    tempnu <- 0;
    } else {
      tempnu <- Lambda - TotBiosamples[k];
    }
  xSigmaSquared <- (tempnu * BGStDev[k]^2 +
       (TotBiosamples[k]-1) * LogFoldStndDev[k]^2)/
           (tempnu + TotBiosamples[k]-2);
    
  xAdjustedVar <- c(xAdjustedVar, xSigmaSquared);
  nu <- c(nu, tempnu);
}
      
xAdjustedVar <- xAdjustedVar[2:(length(xAdjustedVar))];
nu <- nu[2:length(nu)];
      
AdjustedStndDev <- sqrt(xAdjustedVar);
AdjustedStndErr <- sqrt(xAdjustedVar/Lambda);

rm(list=ls(pat="^temp"));
      
              #### Plot Fits Before and After Adjustment  ####
if (PlotResults) {
  
  #Plot to file
  postscript(file=paste(OutputPrefix, "CyberT", "L",
             as.character(Lambda), ".eps", sep=""),  horizontal=FALSE,
             onefile=TRUE);

  chh <-par()$cxy[2];
  chw <-par()$cxy[2];
  par(mar=c(4,5,3,2), omi=c(0.5,0.5,0.5,0.5));

  par(mfrow=c(2,1));

  y0 <- as.numeric(as.vector(LogFoldStndDev));
  ymax <- max(y0);
  x0 <- as.numeric(as.vector(Intensity));
 
  plot(y0 ~ x0, type="p",  main="StndDev", 
     xlab="A", ylab = "StndDev \nLogFoldProbeSets",
     cex= 1.0, xlim=c(0, 15),
     ylim=c(0, ymax),
     pch=16);
  lines(xBGSigmaLowess$x, xBGSigmaLowess$y, col=2);

  y0 <- AdjustedStndDev;
  x0 <- as.numeric(as.vector(Intensity));

  #Best fit after adjustment
  xCorrectedLowess <- lowess(Intensity, AdjustedStndDev, f=LowessfParam);
  
  plot(y0 ~ x0, type="p",  main="Adj Stnd Dev", 
     xlab="A", ylab = "AdjStndDev \nLogFoldProbeSets",
     cex= 1.0, xlim=c(0, 15),
     ylim=c(0, ymax),
     pch=16);
  lines(xBGSigmaLowess$x, xBGSigmaLowess$y, col=2);
  lines(xCorrectedLowess$x, xCorrectedLowess$y, col=3);
 
  dev.off()
}
rm(list=ls(pat="^x"));
rm(list=ls(pat="^y"));


           ###################################################
           ####    PERFORM t-test AND CALCULATE pVALUES   #### 
           ###################################################


#Calculate t-Test statistic using a pairwise t-Test with log fold data
tStat <- abs(LogFoldChange/AdjustedStndErr);

#Calculate pValue for tStat with AdjStndError
tStatFilter <- order(tStat, decreasing=TRUE);
RankedtStat <- tStat[tStatFilter];
RankedTotBiosamples <- TotBiosamples[tStatFilter];

DegreesOfFreedom <- Lambda - 2;
pValue <- 2 * (1 - pt(RankedtStat, DegreesOfFreedom));

#Calculate pValue for tStat with StndError
tStatNoCyberT <- abs(LogFoldChange/StndErr);  #For no CyberT
RankedtStatNoCyberT <- tStatNoCyberT[tStatFilter]; #For no CyberT

DOFnoCyberT <- RankedTotBiosamples - 1;   #For no CyberT
pValueNoCyberT <- 2 * (1 - pt(RankedtStatNoCyberT, DOFnoCyberT));


#################Calculate Multiplicity p-values for CyberT###############

xNumOfHyp <- length(RankedtStat);
Rank <- 1:xNumOfHyp;

FreqpValue <- Rank/xNumOfHyp;

BHpValue <- pValue * xNumOfHyp / Rank;
BHpValue[which(BHpValue > 1)] <- 1;
for (i in Rank) {
  BHpValue[i] <- min(BHpValue[i:xNumOfHyp]);   
}

PKu <- 1-(2 * (1 - pt(sum(RankedtStat), DegreesOfFreedom)));
Pg <- 1 - pValue;
PhenompValue <- PKu - Pg;

BonpValue <- xNumOfHyp * pValue;
BonpValue[which(BonpValue > 1)] <- 1;

SidakpValue <- 1 - Pg^xNumOfHyp;

HolmpValue <- (xNumOfHyp - Rank + 1) * pValue;
HolmpValue[which(HolmpValue > 1)] <- 1;
for (i in xNumOfHyp:1) {
  HolmpValue[i] <- max(HolmpValue[1:i]);   
}

HochpValue <- (xNumOfHyp - Rank + 1) * pValue;
HochpValue[which(HochpValue > 1)] <- 1;
for (i in Rank) {
  HochpValue[i] <- min(HochpValue[i:xNumOfHyp]);   
}

######################### Permutation Methods ############################
              #####    Perform Permutation Resampling  ####
#All data will be referenced based on relative rank determined by t-test
#statistic

if (MultCorrPerm) {

    if (CyberTOFF) {
        SampletStat <- tStatNoCyberT;
    } else {
        SampletStat <- tStat;
    }
    
    BStat <- rep(0, length(SampletStat));
    WYBStat <- rep(0, length(SampletStat));
    B <- TotalPermutations;

    IDVecFilter <- as.logical(match(IDVec, HuberGenes, nomatch=0));
    PerIDVec <- IDVec[IDVecFilter];
    PerIDVecf <- factor(PerIDVec);

    PerLogFoldVector <- LogFoldVec[IDVecFilter];
    PerArrayID <- ArrayID[IDVecFilter];

    xNumOfHyp <- length(RankedtStat);
    for (b in 1:B) {

        if (BootstrapSample) {
            PertStat <- sample(SampletStat, replace=TRUE);
        } else {
            PertStat <- sample(SampletStat);
        }
        
        for (k in 1:xNumOfHyp) {
            bStatk <- rep(0, xNumOfHyp);
            tStatk <- rep(RankedtStat[k], xNumOfHyp);
            bStatk[which(PertStat[1:xNumOfHyp] >= tStatk)] <- 1;
            BStatk <- sum(bStatk);
            BStat[k] <- BStat[k] + BStatk;
            
        }

        UStat <- rep(0, xNumOfHyp);
        UStat[xNumOfHyp] <- PertStat[xNumOfHyp];
        WYBStat[xNumOfHyp] <- WYBStat[xNumOfHyp] + 1;
        for (k in (xNumOfHyp-1):1) {
            UStat[k] <- max(UStat[k+1], PertStat[k]);
            if (UStat[k] >= RankedtStat[k]) {
                WYBStat[k] <- WYBStat[k] + 1
            }
        }   
    }

    WYpValue <- WYBStat / B;

    Rank <- 1:xNumOfHyp;
    EstpValue <-  BStat/(xNumOfHyp * B);
    
    BHpValuePerm <- EstpValue * xNumOfHyp / Rank;
    BHpValuePerm[which(BHpValuePerm > 1)] <- 1;
    for (i in Rank) {
        BHpValuePerm[i] <- min(BHpValuePerm[i:xNumOfHyp]);   
    }   
} else {
    xNumOfHyp <- length(RankedtStat);
    BHpValuePerm <- rep("NA", xNumOfHyp);
    WYpValue <- rep("NA", xNumOfHyp);
    TotalPermutations <- "NA";
    BootstrapSample <- "NA";
}


rm(ArrayList)
rm(list=ls(pat="^x"));
rm(list=ls(pat="^temp"));
rm(list=ls(pat="^Per"));
rm(list=ls(pat="^Next"));

                       ############################
                       #####  OUTPUT RESULTS  #####
                       ############################
      


#Combine vectors into a table for output
CyberTBHTable <- cbind(Rank, ID=HuberGenes[tStatFilter],
                       LogFoldChange=LogFoldChange[tStatFilter],
                       RMSIntensity=Intensity[tStatFilter], 
                       RankedTotBiosamples,
                       TotArrays=TotArrays[tStatFilter],
                       RankedtStat, RankedtStatNoCyberT,  
                       pValue, PhenompValue, pValueNoCyberT, 
                       FreqpValue, BHpValue, BonpValue, SidakpValue,
                       HolmpValue, HochpValue, BHpValuePerm, WYpValue,
                       AdjStndErr=AdjustedStndErr[tStatFilter],
                       StndErr=StndErr[tStatFilter],
                       AdjStndDev=AdjustedStndDev[tStatFilter],
                       BGStDev=BGStDev[tStatFilter],
                       StDev=LogFoldStndDev[tStatFilter]);



#NEED TO ADD MATCHEDARRAYS FLAG!!!!!!!!!!!!!!!!!!!!!!!!!!
FilterHeader1 <- paste("Input Directory = ", InputDirectory, sep="");
FilterHeader2 <- c("NormFiles = ", paste(as.character(ProbesetFiles),
                                             sep=""));
FilterHeader3 <- date();
FilterHeader4 <- paste("Lambda = ", as.character(Lambda), sep="");
FilterHeader5 <- c(paste("MinArrays = ", as.character(MinArrays), sep=""),
                   paste("MinProbes = ", as.character(MinProbes), sep=""));
FilterHeader6 <- paste("MatchedArrays = ", as.character(MatchedArrays), sep="");
FilterHeader7 <- c(paste("MultCorrPerm = ", as.character(MultCorrPerm), sep=""),
             paste("CyberTOFF = ", as.character(CyberTOFF), sep=""),
             paste("BootstrapSample = ", as.character(BootstrapSample), sep=""),
             paste("TotalPermutations = ", as.character(TotalPermutations),
                         sep=""));
FilterHeader8 <- c(paste("LowessfParam = ", as.character(LowessfParam),
                         sep=""),
                   paste("HuberConf = ", as.character(HuberConf), sep=""),
                   paste("HuberTol = ", as.character(HuberTol), sep=""));

FilterHeaderTable <-  c("Rank", "ID", "LogFoldChange", "RMSIntensity",
                        "TotBiosamples", "TotArrays", "tStatCyberT",
                        "tStatNoCyberT","pValue", "PhenompValue",
                        "pValueNoCyberT","FreqpValue", "BHpValue", "BonpValue",
                        "SidakpValue", "HolmpValue", "HochpValue",
                        "BHpValuePerm", "WYpValuePerm",
                        "AdjStndErr", "StndErr",
                        "AdjStndDev", "BGStndDev","StndDev");

cat(FilterHeader1, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\n");
cat(FilterHeader2, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader3, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader4, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader5, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader6, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader7, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader8, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeaderTable, file=paste(OutputPrefix, "L",
                         as.character(Lambda), ".CyTDS", sep=""), sep="\t",
                         append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
             ".CyTDS", sep=""), sep="\n", append=TRUE);

write.table(CyberTBHTable, file=paste(OutputPrefix, "L",
            as.character(Lambda), ".CyTDS", sep=""), quote=FALSE, sep="\t",
            col.names=FALSE, row.names=FALSE, append=TRUE);


#### NEED TO WORK ON THIS OUTPUT FILE!!!!!!!!!!!!!!!!!

##################
##  Genes ordered to genome for GenomeCrawling
##################
RankedGeneNames <- HuberGenes[tStatFilter];
#RankedtStat <- tStat[tStatFilter];
#RankedtStatNoCyberT <- tStatNoCyberT[tStatFilter]; 
RankedLogFoldChange <- LogFoldChange[tStatFilter];
RankedIntensity <- Intensity[tStatFilter];
#RankedTotBiosamples <-  TotBiosamples[tStatFilter];
RankedTotArrays <- TotArrays[tStatFilter];
#pValue is already ranked  
#pValueNoCyberT is already ranked  

GenetStat <- rep(0, length(GenomeGenes));
GenetStatNoCyberT <- rep(0, length(GenomeGenes));
GeneLogFoldChange <- rep(0, length(GenomeGenes));
GeneIntensity <- rep(0, length(GenomeGenes));
GeneTotBiosamples <- rep(0, length(GenomeGenes));
GeneTotArrays <- rep(0, length(GenomeGenes));
GenepValue <- rep(1, length(GenomeGenes))
GenepValueNoCyberT <- rep(1, length(GenomeGenes));

GenePos <- 1:length(GenomeGenes);
GeneNames <- GenomeGenes;


#### Output

for (i in 1:length(RankedGeneNames)) {

    GeneFilter <- as.logical(match(GeneNames, RankedGeneNames[i], nomatch=0));

    GenetStat[GeneFilter] <- RankedtStat[i];
    GenetStatNoCyberT[GeneFilter] <- RankedtStatNoCyberT[i] 
    GeneLogFoldChange[GeneFilter] <- RankedLogFoldChange[i];
    GeneIntensity[GeneFilter] <- RankedIntensity[i];
    GeneTotBiosamples[GeneFilter] <- RankedTotBiosamples[i];
    GeneTotArrays[GeneFilter] <- RankedTotArrays[i];
    GenepValue[GeneFilter] <- pValue[i];
    GenepValueNoCyberT[GeneFilter] <- pValueNoCyberT[i];

}

FilterHeaderTable <-  c("ID", "GenomePos", "tStat", "tStatNoCyberT", "pValue",
                        "pValueNoCyberT", "LogFoldChng", "Intenisty",
                        "BioSamples", "TotArrays");

FilterHeader1 <- paste("Input Directory = ", InputDirectory, sep="");
FilterHeader2 <- c("NormFiles = ", paste(as.character(ProbesetFiles),
                                             sep=""));
FilterHeader3 <- date();
FilterHeader4 <- paste("Lambda = ", as.character(Lambda), sep="");
FilterHeader5 <- c(paste("MinArrays = ", as.character(MinArrays), sep=""),
                   paste("MinProbes = ", as.character(MinProbes), sep=""));
FilterHeader6 <- paste("MatchedArrays = ", as.character(MatchedArrays), sep="");
FilterHeader7 <- paste("LowessfParam = ", as.character(LowessfParam),
                         sep="");          

FilterHeader8 <- c(paste("HuberConf = ", as.character(HuberConf), sep=""),
                   paste("HuberTol = ", as.character(HuberTol), sep=""));

SingleGeneTable <- cbind(GeneNames, GenePos, GenetStat, GenetStatNoCyberT,
                         GenepValue, GenepValueNoCyberT, GeneLogFoldChange,
                         GeneIntensity, GeneTotBiosamples, GeneTotArrays);

cat(FilterHeader1, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\n");
cat(FilterHeader2, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader3, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader4, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader5, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader6, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader7, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);
cat(FilterHeader8, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTDS", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTDS", sep=""), sep="\n", append=TRUE);

cat(FilterHeaderTable, file=paste(OutputPrefix, "L",
                         as.character(Lambda), ".GeCyTDS", sep=""),
                         sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
             ".GeCyTDS", sep=""), sep="\n", append=TRUE);

write.table(SingleGeneTable, file=paste(OutputPrefix, "L",
            as.character(Lambda), ".GeCyTDS", sep=""), quote=FALSE,
            sep="\t", col.names=FALSE, row.names=FALSE, append=TRUE);

rm(SingleGeneTable);

#rm(list=ls());


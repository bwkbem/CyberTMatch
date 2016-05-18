
MinProbes <- 2;
MinArrays <- 4;  

MatchedArrays <- FALSE; # See Input section of this file more detail

################ Permutation Parameters ##############

MultCorrPerm <- FALSE;  #Long Run Times!!!!  Performs WY and BH with resampling
BootstrapSample <- TRUE;
TotalPermutations <- 2000;

######### Statistical Parameters ####################
#PARAMETERS FOR LOWESS
LowessfParam <- 0.4;

#PARAMETERS FOR HUBER/MAD ESTIMATOR
HuberConf <- 1.345; #(ie 95% Confidence Limit)
HuberTol <- 1e-6;


            ################# INPUT / OUTPUT  #################

InputDirectory <- "/home/anon/Projects/Data/SPyM1/Lowess/";

# Probsets are matched to biosamples in order of c(1, 1, 2, 2, 3, 3, 4, 4).
# Matching algorithm is based on this if not this pattern set MatchedArrays <- FALSE
ProbesetFiles <- c("S730.LowPset", "F730.LowPset",
                   "S806.LowPset", "F806.LowPset",
                   "S812.LowPset", "F812.LowPset",
                   "S924.LowPset", "F924.LowPset");
HeaderLength <- 9;  #DO NOT EDIT

GenomeDirectory <- "/home/anon/Projects/Data/SPyM1/Genome/"

GenomeFile <- "SPyM1Table.txt";
GenomeIDPos <- 6;   #DO NOT EDIT
GenomeHeaderLength <- 0;   #DO NOT EDIT

OutputDirectory <- "/home/anon/Projects/Data/SPyM1/Development/CyberT/Unmatched/";



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

    TotBioSamples <- TotArrays/2;

} else {

    TotBioSamples <- 0;
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

        TotBioSamples <- c(TotBioSamples, count);
    }
    TotBioSamples <- TotBioSamples[2:length(TotBioSamples)];
}


StndErr <- LogFoldStndDev/sqrt(TotBioSamples);

### Remove non Gene Probes and unrepresented genes in annotation
GenomeFilter <- as.logical(as.vector(
                             match(HuberGenes, GenomeGenes, nomatch=0)));

HuberGenes <- HuberGenes[GenomeFilter];
            
Intensity <- Intensity[GenomeFilter];

LogFoldChange <- LogFoldChange[GenomeFilter];

LogFoldStndDev <- LogFoldStndDev[GenomeFilter];

TotArrays <- TotArrays[GenomeFilter];

TotBioSamples <- TotBioSamples[GenomeFilter];

StndErr <- StndErr[GenomeFilter];


rm(list=ls(pat="^z"));
rm(list=ls(pat="^f"));
rm(list=ls(pat="^g"));



           ###################################################
           ####    PERFORM t-test AND CALCULATE pVALUES   #### 
           ###################################################

#Calculate pValue for tStat with StndError
tStatNoCyberT <- abs(LogFoldChange/StndErr);  #For no CyberT
tStatFilterNoCyberT <- order(tStatNoCyberT, decreasing=TRUE);
RankedtStatNoCyberT <- tStatNoCyberT[tStatFilterNoCyberT]; 
RankedTotBioSamples <- TotBioSamples[tStatFilterNoCyberT];
RankedGenes <- HuberGenes[tStatFilterNoCyberT];

DOFnoCyberT <- RankedTotBioSamples - 1;   #For no CyberT
pValueNoCyberT <- 2 * (1 - pt(RankedtStatNoCyberT, DOFnoCyberT));

############## Calculate Multiplicity p-values for CyberT ###############

######################### Permutation Methods ############################
              #####    Perform Permutation Resampling  ####
#All data will be referenced based on relative rank determined by t-test
#statistic

    if (MultCorrPerm) {

    PermFilter1 <- as.logical(match(RankedTotBioSamples, 2, nomatch=0));
    PermFilter2 <- as.logical(match(RankedTotBioSamples, 3, nomatch=0));
    PermFilter3 <- as.logical(match(RankedTotBioSamples, 4, nomatch=0));

    SampletStat1 <- RankedtStatNoCyberT[PermFilter1];
    SampletStat2 <- RankedtStatNoCyberT[PermFilter2];
    SampletStat3 <- RankedtStatNoCyberT[PermFilter3];


        
#       SampletStat <- tStatNoCyberT;
  
        BStat1 <- rep(0, length(SampletStat1));
        BStat2 <- rep(0, length(SampletStat2));
        BStat3 <- rep(0, length(SampletStat3));

        WYBStat1 <- rep(0, length(SampletStat1));
        WYBStat2 <- rep(0, length(SampletStat2));
        WYBStat3 <- rep(0, length(SampletStat3));

        xNumOfHyp1 <- length(SampletStat1);
        xNumOfHyp2 <- length(SampletStat2);
        xNumOfHyp3 <- length(SampletStat3);

        B <- TotalPermutations;

        #####  Sample 1
    
        for (b in 1:B) {

            if (BootstrapSample) {
                PertStat1 <- sample(SampletStat1, replace=TRUE);
            } else {
                PertStat1 <- sample(SampletStat1);
            }
        
            for (k in 1:xNumOfHyp1) {
                bStatk <- rep(0, xNumOfHyp1);
                tStatk <- rep(SampletStat1[k], xNumOfHyp1);
                bStatk[which(PertStat1[1:xNumOfHyp1] >= tStatk)] <- 1;
                BStatk <- sum(bStatk);
                BStat1[k] <- BStat1[k] + BStatk;
            
            }

            UStat <- rep(0, xNumOfHyp1);
            UStat[xNumOfHyp1] <- PertStat1[xNumOfHyp1];
            WYBStat1[xNumOfHyp1] <- WYBStat1[xNumOfHyp1] + 1;
            for (k in (xNumOfHyp1-1):1) {
                UStat[k] <- max(UStat[k+1], PertStat1[k]);
                if (UStat[k] >= SampletStat1[k]) {
                    WYBStat1[k] <- WYBStat1[k] + 1
                }
            }   
        }

        WYpValue1 <- WYBStat1 / B;

        Rank1 <- 1:xNumOfHyp1;
        EstpValue1 <-  BStat1/(xNumOfHyp1 * B);
    
        BHpValuePerm1 <- EstpValue1 * xNumOfHyp1 / Rank1;
        BHpValuePerm1[which(BHpValuePerm1 > 1)] <- 1;
        for (i in Rank1) {
            BHpValuePerm1[i] <- min(BHpValuePerm1[i:xNumOfHyp1]);   
        }

        #### Sample 2

        for (b in 1:B) {

            if (BootstrapSample) {
                PertStat2 <- sample(SampletStat2, replace=TRUE);
            } else {
                PertStat2 <- sample(SampletStat2);
            }
        
            for (k in 1:xNumOfHyp2) {
                bStatk <- rep(0, xNumOfHyp2);
                tStatk <- rep(SampletStat2[k], xNumOfHyp2);
                bStatk[which(PertStat2[1:xNumOfHyp2] >= tStatk)] <- 1;
                BStatk <- sum(bStatk);
                BStat2[k] <- BStat2[k] + BStatk;
            
            }

            UStat <- rep(0, xNumOfHyp2);
            UStat[xNumOfHyp2] <- PertStat2[xNumOfHyp2];
            WYBStat2[xNumOfHyp2] <- WYBStat2[xNumOfHyp2] + 1;
            for (k in (xNumOfHyp2-1):1) {
                UStat[k] <- max(UStat[k+1], PertStat2[k]);
                if (UStat[k] >= SampletStat2[k]) {
                    WYBStat2[k] <- WYBStat2[k] + 1
                }
            }   
        }

        WYpValue2 <- WYBStat2 / B;

        Rank2 <- 1:xNumOfHyp2;
        EstpValue2 <-  BStat2/(xNumOfHyp2 * B);
    
        BHpValuePerm2 <- EstpValue2 * xNumOfHyp2 / Rank2;
        BHpValuePerm2[which(BHpValuePerm2 > 1)] <- 1;
        for (i in Rank2) {
            BHpValuePerm2[i] <- min(BHpValuePerm2[i:xNumOfHyp2]);   
        }

       ###### Sample 3

       for (b in 1:B) {

            if (BootstrapSample) {
                PertStat3 <- sample(SampletStat3, replace=TRUE);
            } else {
                PertStat3 <- sample(SampletStat3);
            }
        
            for (k in 1:xNumOfHyp3) {
                bStatk <- rep(0, xNumOfHyp3);
                tStatk <- rep(SampletStat3[k], xNumOfHyp3);
                bStatk[which(PertStat3[1:xNumOfHyp3] >= tStatk)] <- 1;
                BStatk <- sum(bStatk);
                BStat3[k] <- BStat3[k] + BStatk;
            
            }

            UStat <- rep(0, xNumOfHyp1);
            UStat[xNumOfHyp3] <- PertStat3[xNumOfHyp3];
            WYBStat3[xNumOfHyp3] <- WYBStat3[xNumOfHyp3] + 1;
            for (k in (xNumOfHyp3-1):1) {
                UStat[k] <- max(UStat[k+1], PertStat3[k]);
                if (UStat[k] >= SampletStat3[k]) {
                    WYBStat3[k] <- WYBStat3[k] + 1
                }
            }   
        }

        WYpValue3 <- WYBStat3 / B;

        Rank3 <- 1:xNumOfHyp3;
        EstpValue3 <-  BStat3/(xNumOfHyp3 * B);
    
        BHpValuePerm3 <- EstpValue3 * xNumOfHyp3 / Rank3;
        BHpValuePerm3[which(BHpValuePerm3 > 1)] <- 1;
        for (i in Rank3) {
            BHpValuePerm3[i] <- min(BHpValuePerm3[i:xNumOfHyp3]);   
        }


    
    # Combine all three BHpvalue1-3 sets ordered to ranked tStat

    BHpValuePerm <- rep(NA, length(RankedGenes));
    BHGenes1 <- RankedGenes[PermFilter1];
    BHGenes2 <- RankedGenes[PermFilter2];
    BHGenes3 <- RankedGenes[PermFilter3];

    for (i in 1:length(BHGenes1)) {
        BHpValuePerm[as.logical(match(RankedGenes, BHGenes1[i],
                                               nomatch=0))] <- BHpValuePerm1[i];
    }

    for (i in 1:length(BHGenes2)) {
        BHpValuePerm[as.logical(match(RankedGenes, BHGenes2[i],
                                               nomatch=0))] <- BHpValuePerm2[i];
    }

    for (i in 1:length(BHGenes3)) {
        BHpValuePerm[as.logical(match(RankedGenes, BHGenes3[i],
                                               nomatch=0))] <- BHpValuePerm3[i];
    }

    # Combine all three WYpvalue1-3 sets ordered to ranked tStat

    WYpValuePerm <- rep(NA, length(RankedGenes));
    WYGenes1 <- RankedGenes[PermFilter1];
    WYGenes2 <- RankedGenes[PermFilter2];
    WYGenes3 <- RankedGenes[PermFilter3];

    for (i in 1:length(WYGenes1)) {
        WYpValuePerm[as.logical(match(RankedGenes, WYGenes1[i],
                                               nomatch=0))] <- WYpValue1[i];
    }

    for (i in 1:length(WYGenes2)) {
        WYpValuePerm[as.logical(match(RankedGenes, WYGenes2[i],
                                               nomatch=0))] <- WYpValue2[i];
    }

    for (i in 1:length(WYGenes3)) {
        WYpValuePerm[as.logical(match(RankedGenes, WYGenes3[i],
                                               nomatch=0))] <- WYpValue3[i];
    }

    TBHpValuePerm <- BHpValuePerm;
    TWYpValuePerm <- WYpValuePerm;
}

########################### End of Permutation Methods ########################

xNumOfHyp <- length(RankedtStatNoCyberT);
Rank <- 1:xNumOfHyp;

FreqpValue <- Rank/xNumOfHyp;

BHpValue <- pValueNoCyberT * xNumOfHyp / Rank;
BHpValue[which(BHpValue > 1)] <- 1;
for (i in Rank) {
    BHpValue[i] <- min(BHpValue[i:xNumOfHyp]);   
}

#### BH sorted pValues #####
#sortFilter <- order(pValueNoCyberT);
#sortpValueNoCyberT <- pValueNoCyberT[sortFilter];
#sortGenes <- RankedGenes[sortFilter];

#sortBHpValue <- sortpValueNoCyberT * xNumOfHyp / Rank;
#sortBHpValue[which(sortBHpValue > 1)] <- 1;
#for (i in Rank) {
#    sortBHpValue[i] <- min(sortBHpValue[i:xNumOfHyp]);   
#}

#sBHpValue <- rep(NA, length(RankedGenes));
#for (i in 1:length(sortGenes)) {
#    sBHpValue[as.logical(match(RankedGenes, sortGenes[i],
#                                               nomatch=0))] <- sortBHpValue[i];
#}

#### BH performed for each degree of freedom separately

BHFilter1 <- as.logical(match(RankedTotBioSamples, 2, nomatch=0));
BHFilter2 <- as.logical(match(RankedTotBioSamples, 3, nomatch=0));
BHFilter3 <- as.logical(match(RankedTotBioSamples, 4, nomatch=0));

BHpValueIn1 <- pValueNoCyberT[BHFilter1];
BHpValueIn2 <- pValueNoCyberT[BHFilter2];
BHpValueIn3 <- pValueNoCyberT[BHFilter3];

xBHNumOfHyp1 <- length(BHpValueIn1);
BHRank1 <- 1:xBHNumOfHyp1;
BHpValue1 <- BHpValueIn1 * xBHNumOfHyp1 / BHRank1;
BHpValue1[which(BHpValue1 > 1)] <- 1;
for (i in BHRank1) {
    BHpValue1[i] <- min(BHpValue1[i:xBHNumOfHyp1]);
}

xBHNumOfHyp2 <- length(BHpValueIn2);
BHRank2 <- 1:xBHNumOfHyp2;
BHpValue2 <- BHpValueIn2 * xBHNumOfHyp2 / BHRank2;
BHpValue2[which(BHpValue1 > 1)] <- 1;
for (i in BHRank2) {
    BHpValue2[i] <- min(BHpValue2[i:xBHNumOfHyp2]);
}

xBHNumOfHyp3 <- length(BHpValueIn3);
BHRank3 <- 1:xBHNumOfHyp3;
BHpValue3 <- BHpValueIn3 * xBHNumOfHyp3 / BHRank3;
BHpValue3[which(BHpValue1 > 1)] <- 1;
for (i in BHRank3) {
    BHpValue3[i] <- min(BHpValue3[i:xBHNumOfHyp3]);
}

# Combine all three BHpvalue1-3 sets ordered to ranked tStat

BHpValueMixed <- rep(NA, length(RankedGenes));
BHGenes1 <- RankedGenes[BHFilter1];
BHGenes2 <- RankedGenes[BHFilter2];
BHGenes3 <- RankedGenes[BHFilter3];

for (i in 1:length(BHGenes1)) {
    BHpValueMixed[as.logical(match(RankedGenes, BHGenes1[i],
                                               nomatch=0))] <- BHpValue1[i];
}

for (i in 1:length(BHGenes2)) {
    BHpValueMixed[as.logical(match(RankedGenes, BHGenes2[i],
                                               nomatch=0))] <- BHpValue2[i];
}

for (i in 1:length(BHGenes3)) {
    BHpValueMixed[as.logical(match(RankedGenes, BHGenes3[i],
                                               nomatch=0))] <- BHpValue3[i];
}

#####


#PKu <- 1-(2 * (1 - pt(sum(RankedtStatNoCyberT), DegreesOfFreedom)));
Pg <- 1 - pValueNoCyberT;
DeleuzepValue <- pValueNoCyberT;

BonpValue <- xNumOfHyp * pValueNoCyberT;
BonpValue[which(BonpValue > 1)] <- 1;

SidakpValue <- 1 - Pg^xNumOfHyp;

HolmpValue <- (xNumOfHyp - Rank + 1) * pValueNoCyberT;
HolmpValue[which(HolmpValue > 1)] <- 1;
for (i in xNumOfHyp:1) {
    HolmpValue[i] <- max(HolmpValue[1:i]);   
}

HochpValue <- (xNumOfHyp - Rank + 1) * pValueNoCyberT;
HochpValue[which(HochpValue > 1)] <- 1;
for (i in Rank) {
    HochpValue[i] <- min(HochpValue[i:xNumOfHyp]);   
}

TLambda <- "null"
TRankedtStat <-  RankedtStatNoCyberT;
TBHpValue <- BHpValue;
TBHpValueMixed <- BHpValueMixed;
TFreqpValue <- FreqpValue;
TDeleuzepValue <- DeleuzepValue;
TSidakpValue <- SidakpValue;
TBonpValue <- BonpValue;
THolmpValue <- HolmpValue;

TGenes <- RankedGenes;
TtStat <- RankedtStatNoCyberT;
TTotBioSamples <- RankedTotBioSamples;

for (Lambda in c(4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 25, 30)) {
               #(2, 3, 4, 5, 6, 7,  8,  9, 10, 11, 12, 13, 14
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
        if (TotBioSamples[k] >= Lambda) {
            tempnu <- 0;
        } else {
            tempnu <- Lambda - TotBioSamples[k];
        }
        xSigmaSquared <- (tempnu * BGStDev[k]^2 +
                              (TotBioSamples[k]-1) * LogFoldStndDev[k]^2)/
                              (tempnu + TotBioSamples[k]-2);
    
        xAdjustedVar <- c(xAdjustedVar, xSigmaSquared);
        nu <- c(nu, tempnu);
    }
      
    xAdjustedVar <- xAdjustedVar[2:(length(xAdjustedVar))];
    nu <- nu[2:length(nu)];
      
    AdjustedStndDev <- sqrt(xAdjustedVar);
    AdjustedStndErr <- sqrt(xAdjustedVar/Lambda);

    rm(list=ls(pat="^temp"));


    #Calculate t-Test statistic using a pairwise t-Test with log fold data
    tStat <- abs(LogFoldChange/AdjustedStndErr);

    #Calculate pValue for tStat with AdjStndError
    tStatFilter <- order(tStat, decreasing=TRUE);
    RankedtStat <- tStat[tStatFilter];
    RankedTotBioSamples <- TotBioSamples[tStatFilter];
    RankedGenes <- HuberGenes[tStatFilter];

    DegreesOfFreedom <- Lambda - 2;
    pValue <- 2 * (1 - pt(RankedtStat, DegreesOfFreedom));


    #############Calculate Multiplicity p-values for CyberT###############


######################### Permutation Methods ############################
              #####    Perform Permutation Resampling  ####
#All data will be referenced based on relative rank determined by t-test
#statistic

    if (MultCorrPerm) {

    PermFilter1 <- as.logical(match(RankedTotBioSamples, 2, nomatch=0));
    PermFilter2 <- as.logical(match(RankedTotBioSamples, 3, nomatch=0));
    PermFilter3 <- as.logical(match(RankedTotBioSamples, 4, nomatch=0));

    SampletStat1 <- RankedtStat[PermFilter1];
    SampletStat2 <- RankedtStat[PermFilter2];
    SampletStat3 <- RankedtStat[PermFilter3];


        
#       SampletStat <- tStatNoCyberT;
  
        BStat1 <- rep(0, length(SampletStat1));
        BStat2 <- rep(0, length(SampletStat2));
        BStat3 <- rep(0, length(SampletStat3));

        WYBStat1 <- rep(0, length(SampletStat1));
        WYBStat2 <- rep(0, length(SampletStat2));
        WYBStat3 <- rep(0, length(SampletStat3));

        xNumOfHyp1 <- length(SampletStat1);
        xNumOfHyp2 <- length(SampletStat2);
        xNumOfHyp3 <- length(SampletStat3);

        B <- TotalPermutations;

        #####  Sample 1
    
        for (b in 1:B) {

            if (BootstrapSample) {
                PertStat1 <- sample(SampletStat1, replace=TRUE);
            } else {
                PertStat1 <- sample(SampletStat1);
            }
        
            for (k in 1:xNumOfHyp1) {
                bStatk <- rep(0, xNumOfHyp1);
                tStatk <- rep(SampletStat1[k], xNumOfHyp1);
                bStatk[which(PertStat1[1:xNumOfHyp1] >= tStatk)] <- 1;
                BStatk <- sum(bStatk);
                BStat1[k] <- BStat1[k] + BStatk;
            
            }

            UStat <- rep(0, xNumOfHyp1);
            UStat[xNumOfHyp1] <- PertStat1[xNumOfHyp1];
            WYBStat1[xNumOfHyp1] <- WYBStat1[xNumOfHyp1] + 1;
            for (k in (xNumOfHyp1-1):1) {
                UStat[k] <- max(UStat[k+1], PertStat1[k]);
                if (UStat[k] >= SampletStat1[k]) {
                    WYBStat1[k] <- WYBStat1[k] + 1
                }
            }   
        }

        WYpValue1 <- WYBStat1 / B;

        Rank1 <- 1:xNumOfHyp1;
        EstpValue1 <-  BStat1/(xNumOfHyp1 * B);
    
        BHpValuePerm1 <- EstpValue1 * xNumOfHyp1 / Rank1;
        BHpValuePerm1[which(BHpValuePerm1 > 1)] <- 1;
        for (i in Rank1) {
            BHpValuePerm1[i] <- min(BHpValuePerm1[i:xNumOfHyp1]);   
        }

        #### Sample 2

        for (b in 1:B) {

            if (BootstrapSample) {
                PertStat2 <- sample(SampletStat2, replace=TRUE);
            } else {
                PertStat2 <- sample(SampletStat2);
            }
        
            for (k in 1:xNumOfHyp2) {
                bStatk <- rep(0, xNumOfHyp2);
                tStatk <- rep(SampletStat2[k], xNumOfHyp2);
                bStatk[which(PertStat2[1:xNumOfHyp2] >= tStatk)] <- 1;
                BStatk <- sum(bStatk);
                BStat2[k] <- BStat2[k] + BStatk;
            
            }

            UStat <- rep(0, xNumOfHyp2);
            UStat[xNumOfHyp2] <- PertStat2[xNumOfHyp2];
            WYBStat2[xNumOfHyp2] <- WYBStat2[xNumOfHyp2] + 1;
            for (k in (xNumOfHyp2-1):1) {
                UStat[k] <- max(UStat[k+1], PertStat2[k]);
                if (UStat[k] >= SampletStat2[k]) {
                    WYBStat2[k] <- WYBStat2[k] + 1
                }
            }   
        }

        WYpValue2 <- WYBStat2 / B;

        Rank2 <- 1:xNumOfHyp2;
        EstpValue2 <-  BStat2/(xNumOfHyp2 * B);
    
        BHpValuePerm2 <- EstpValue2 * xNumOfHyp2 / Rank2;
        BHpValuePerm2[which(BHpValuePerm2 > 1)] <- 1;
        for (i in Rank2) {
            BHpValuePerm2[i] <- min(BHpValuePerm2[i:xNumOfHyp2]);   
        }

       ###### Sample 3

       for (b in 1:B) {

            if (BootstrapSample) {
                PertStat3 <- sample(SampletStat3, replace=TRUE);
            } else {
                PertStat3 <- sample(SampletStat3);
            }
        
            for (k in 1:xNumOfHyp3) {
                bStatk <- rep(0, xNumOfHyp3);
                tStatk <- rep(SampletStat3[k], xNumOfHyp3);
                bStatk[which(PertStat3[1:xNumOfHyp3] >= tStatk)] <- 1;
                BStatk <- sum(bStatk);
                BStat3[k] <- BStat3[k] + BStatk;
            
            }

            UStat <- rep(0, xNumOfHyp1);
            UStat[xNumOfHyp3] <- PertStat3[xNumOfHyp3];
            WYBStat3[xNumOfHyp3] <- WYBStat3[xNumOfHyp3] + 1;
            for (k in (xNumOfHyp3-1):1) {
                UStat[k] <- max(UStat[k+1], PertStat3[k]);
                if (UStat[k] >= SampletStat3[k]) {
                    WYBStat3[k] <- WYBStat3[k] + 1
                }
            }   
        }

        WYpValue3 <- WYBStat3 / B;

        Rank3 <- 1:xNumOfHyp3;
        EstpValue3 <-  BStat3/(xNumOfHyp3 * B);
    
        BHpValuePerm3 <- EstpValue3 * xNumOfHyp3 / Rank3;
        BHpValuePerm3[which(BHpValuePerm3 > 1)] <- 1;
        for (i in Rank3) {
            BHpValuePerm3[i] <- min(BHpValuePerm3[i:xNumOfHyp3]);   
        }


    
    # Combine all three BHpvalue1-3 sets ordered to ranked tStat

    BHpValuePerm <- rep(NA, length(RankedGenes));
    BHGenes1 <- RankedGenes[PermFilter1];
    BHGenes2 <- RankedGenes[PermFilter2];
    BHGenes3 <- RankedGenes[PermFilter3];

    for (i in 1:length(BHGenes1)) {
        BHpValuePerm[as.logical(match(RankedGenes, BHGenes1[i],
                                               nomatch=0))] <- BHpValuePerm1[i];
    }

    for (i in 1:length(BHGenes2)) {
        BHpValuePerm[as.logical(match(RankedGenes, BHGenes2[i],
                                               nomatch=0))] <- BHpValuePerm2[i];
    }

    for (i in 1:length(BHGenes3)) {
        BHpValuePerm[as.logical(match(RankedGenes, BHGenes3[i],
                                               nomatch=0))] <- BHpValuePerm3[i];
    }

    # Combine all three WYpvalue1-3 sets ordered to ranked tStat

    WYpValuePerm <- rep(NA, length(RankedGenes));
    WYGenes1 <- RankedGenes[PermFilter1];
    WYGenes2 <- RankedGenes[PermFilter2];
    WYGenes3 <- RankedGenes[PermFilter3];

    for (i in 1:length(WYGenes1)) {
        WYpValuePerm[as.logical(match(RankedGenes, WYGenes1[i],
                                               nomatch=0))] <- WYpValue1[i];
    }

    for (i in 1:length(WYGenes2)) {
        WYpValuePerm[as.logical(match(RankedGenes, WYGenes2[i],
                                               nomatch=0))] <- WYpValue2[i];
    }

    for (i in 1:length(WYGenes3)) {
        WYpValuePerm[as.logical(match(RankedGenes, WYGenes3[i],
                                               nomatch=0))] <- WYpValue3[i];
    }


    TBHpValuePerm <- cbind(TBHpValuePerm, BHpValuePerm);
    TWYpValuePerm <- cbind(TWYpValuePerm, WYpValuePerm);
}
    
########################### End of Permutation Methods ########################

    
    xNumOfHyp <- length(RankedtStat);
    Rank <- 1:xNumOfHyp;

    FreqpValue <- Rank/xNumOfHyp;

    BHpValue <- pValue * xNumOfHyp / Rank;
    BHpValue[which(BHpValue > 1)] <- 1;
    for (i in Rank) {
        BHpValue[i] <- min(BHpValue[i:xNumOfHyp]);   
    }

    #### BH performed for each degree of freedom separately

    BHFilter1 <- as.logical(match(RankedTotBioSamples, 2, nomatch=0));
    BHFilter2 <- as.logical(match(RankedTotBioSamples, 3, nomatch=0));
    BHFilter3 <- as.logical(match(RankedTotBioSamples, 4, nomatch=0));

    BHpValueIn1 <- pValue[BHFilter1];
    BHpValueIn2 <- pValue[BHFilter2];
    BHpValueIn3 <- pValue[BHFilter3];

    xBHNumOfHyp1 <- length(BHpValueIn1);
    BHRank1 <- 1:xBHNumOfHyp1;
    BHpValue1 <- BHpValueIn1 * xBHNumOfHyp1 / BHRank1;
    BHpValue1[which(BHpValue1 > 1)] <- 1;
    for (i in BHRank1) {
        BHpValue1[i] <- min(BHpValue1[i:xBHNumOfHyp1]);
    }

    xBHNumOfHyp2 <- length(BHpValueIn2);
    BHRank2 <- 1:xBHNumOfHyp2;
    BHpValue2 <- BHpValueIn2 * xBHNumOfHyp2 / BHRank2;
    BHpValue2[which(BHpValue1 > 1)] <- 1;
    for (i in BHRank2) {
        BHpValue2[i] <- min(BHpValue2[i:xBHNumOfHyp2]);
    }

    xBHNumOfHyp3 <- length(BHpValueIn3);
    BHRank3 <- 1:xBHNumOfHyp3;
    BHpValue3 <- BHpValueIn3 * xBHNumOfHyp3 / BHRank3;
    BHpValue3[which(BHpValue1 > 1)] <- 1;
    for (i in BHRank3) {
        BHpValue3[i] <- min(BHpValue3[i:xBHNumOfHyp3]);
    }

    # Combine all three BHpvalue1-3 sets ordered to ranked tStat

    BHpValueMixed <- rep(NA, length(RankedGenes));
    BHGenes1 <- RankedGenes[BHFilter1];
    BHGenes2 <- RankedGenes[BHFilter2];
    BHGenes3 <- RankedGenes[BHFilter3];

    for (i in 1:length(BHGenes1)) {
        BHpValueMixed[as.logical(match(RankedGenes, BHGenes1[i],
                                               nomatch=0))] <- BHpValue1[i];
    }

    for (i in 1:length(BHGenes2)) {
        BHpValueMixed[as.logical(match(RankedGenes, BHGenes2[i],
                                               nomatch=0))] <- BHpValue2[i];
    }

    for (i in 1:length(BHGenes3)) {
        BHpValueMixed[as.logical(match(RankedGenes, BHGenes3[i],
                                               nomatch=0))] <- BHpValue3[i];
    }


#####

    

    
    PKu <- 1-(2 * (1 - pt(sum(RankedtStat), DegreesOfFreedom)));
    Pg <- 1 - pValue;
    DeleuzepValue <- PKu - Pg;

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


    TLambda <- c(TLambda, Lambda); 
    TRankedtStat <-  cbind(TRankedtStat, RankedtStat);
    TBHpValue <- cbind(TBHpValue, BHpValue);
    TBHpValueMixed <- cbind(TBHpValueMixed, BHpValueMixed);
    TFreqpValue <- cbind(TFreqpValue, FreqpValue);
    TDeleuzepValue <- cbind(TDeleuzepValue, DeleuzepValue);
    TSidakpValue <- cbind(TSidakpValue, SidakpValue);
    TBonpValue <- cbind(TBonpValue, BonpValue);
    THolmpValue <- cbind(THolmpValue, HolmpValue);

    TGenes <- cbind(TGenes, RankedGenes);
    TtStat <- cbind(TtStat, RankedtStat);
    TTotBioSamples <- cbind(TTotBioSamples, RankedTotBioSamples);

} 


rm(ArrayList)
rm(list=ls(pat="^x"));
rm(list=ls(pat="^temp"));
rm(list=ls(pat="^Per"));
rm(list=ls(pat="^Next"));


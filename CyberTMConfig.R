

Lambda <- 6;

MinProbes <- 2;
MinArrays <- 4;  

MatchedArrays <- FALSE; # See Input section of this file more detail

#########  Permutation resampling Methods ONLY ##############
MultCorrPerm <- FALSE;  #Long Run Times!!!!  Performs WY and BH with resampling
BootstrapSample <- TRUE;
CyberTOFF <- FALSE;  #For permutation methods only
TotalPermutations <- 2000;


######### Statistical Parameters ####################
#PARAMETERS FOR LOWESS
LowessfParam <- 0.4;

#PARAMETERS FOR HUBER/MAD ESTIMATOR
HuberConf <- 1.345; #(ie 95% Confidence Limit)
HuberTol <- 1e-6;


            ################# INPUT / OUTPUT  #################

PlotResults <- TRUE;

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



OutputDirectory <- "/home/anon/Projects/Data/SPyM1/Development/";
OutputFile <- "G8b";



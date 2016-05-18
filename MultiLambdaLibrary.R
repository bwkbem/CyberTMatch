                       ############################
                       #####  OUTPUT RESULTS  #####
                       ############################
      

#DEFINE GRAPHICS WINDOW SIZE
WinWidth <- 8.5;           
WinHeight <- 11;



#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1logpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TDeleuzepValue[,2]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TDeleuzepValue[,3]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="blue");

y <- log10(TDeleuzepValue[,4]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="green");


legend(0, -2, c("Deleuzian Null", "Lambda=4", "Lambda=5", "Lambda=6"), cex=2, pch=1, col=c("black", "red", "blue", "green"));

dev.off();

####

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TDeleuzepValue[,3]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="blue");

legend(0, -2, c("Deleuzian Null", "Lambda=5"), cex=2, pch=1, col=c("black", "blue"));

dev.off();

######

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TDeleuzepValue[,4]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="green");

legend(0, -2, c("Deleuzian Null", "Lambda=6"), cex=2, pch=1, col=c("black", "green"));

dev.off();


#####

x11(width=11, height=8);

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,2]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TBHpValue[,3]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="blue");

y <- log10(TBHpValue[,4]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="green");


legend(0, -2, c("BH Null", "Lambda=4", "Lambda=5", "Lambda=6"), cex=2, pch=1, col=c("black", "red", "blue", "green"));


dev.off();


###########

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TBHpValue[,2]);
x <- log10(length(TRankedtStat[,2]):1);

points(x, y, cex=1, pch=1, col="orange");

y <- log10(TBHpValue[,3]);
x <- log10(length(TRankedtStat[,2]):1);

points(x, y, cex=1, pch=1, col="yellow");

y <- log10(TBHpValue[,4]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="green");

y <- log10(TBHpValue[,5]);
x <- log10(length(TRankedtStat[,2]):1);

points(x, y, cex=1, pch=1, col="blue");

y <- log10(TBHpValue[,6]);
x <- log10(length(TRankedtStat[,2]):1);

points(x, y, cex=1, pch=1, col="purple");


legend(0, -2, c("Deleuze Null", "BH Null", "BH Lambda=4", "BH Lambda=5", "BH Lambda=6", "BH Lambda=7", "BH Lambda=8"), cex=2, pch=1, col=c("black", "red", "orange", "yellow", "green", "blue", "purple"));


dev.off();

################################################################################



postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1alogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

#x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TBHpValue[,11]);
x <- log10(length(TRankedtStat[,11]):1);

points(x, y, cex=1, pch=1, col="orange");

y <- log10(TBHpValue[,13]);
x <- log10(length(TRankedtStat[,13]):1);

points(x, y, cex=1, pch=1, col="green");

y <- log10(TBHpValue[,15]);
x <- log10(length(TRankedtStat[,15]):1);
points(x, y, cex=1, pch=1, col="blue");

legend(0, -2, c("Deleuzian Null", "BH Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1blogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

#x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(2.6, max(x)),
                       ylim=c(-4.2, -0.2),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TBHpValue[,11]);
x <- log10(length(TRankedtStat[,11]):1);

points(x, y, cex=1, pch=1, col="orange");

y <- log10(TBHpValue[,13]);
x <- log10(length(TRankedtStat[,13]):1);

points(x, y, cex=1, pch=1, col="green");

y <- log10(TBHpValue[,15]);
x <- log10(length(TRankedtStat[,15]):1);
points(x, y, cex=1, pch=1, col="blue");

legend(2.6, -2, c("Deleuzian Null", "BH Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig2aBHNull.eps", sep=""),  horizontal=TRUE, onefile=TRUE);

#x11(width=11, height=8);

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(2.5, max(x)),
                       ylim=c(-0.9, -0.4),
     cex=2, pch=1, col="red");


dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig2bBHNull.eps", sep=""),  horizontal=TRUE, onefile=TRUE);

#x11(width=11, height=8);

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(2.5, max(x)),
                       ylim=c(-0.9, -0.4),
     cex=2, pch=1, col="red");

abline(v=2.6031444, lty=2, lwd=2);
abline(v=2.6919651, lty=2, lwd=2);
abline(v=2.7993405, lty=2, lwd=2);

dev.off();

#######



#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1alogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(-4.2, 0),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TBHpValue[,11]);
x <- log10(length(TRankedtStat[,11]):1);

points(x, y, cex=1, pch=1, col="orange");

y <- log10(TFreqpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="yellow");

y <- log10(TBHpValue[,13]);
x <- log10(length(TRankedtStat[,13]):1);

points(x, y, cex=1, pch=1, col="green");

y <- log10(TBHpValue[,15]);
x <- log10(length(TRankedtStat[,15]):1);
points(x, y, cex=1, pch=1, col="blue");

legend(0, -2, c("PS Null", "BH Null", "Freq Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "yellow", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
#dev.off();


#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1blogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(2.6, max(x)),
                       ylim=c(-4.2, -0.2),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="red");

y <- log10(TBHpValue[,11]);
x <- log10(length(TRankedtStat[,11]):1);

points(x, y, cex=1, pch=1, col="orange");

y <- log10(TFreqpValue[,1]);
x <- log10(length(TRankedtStat[,1]):1);

points(x, y, cex=1, pch=1, col="yellow");

y <- log10(TBHpValue[,13]);
x <- log10(length(TRankedtStat[,13]):1);

points(x, y, cex=1, pch=1, col="green");

y <- log10(TBHpValue[,15]);
x <- log10(length(TRankedtStat[,15]):1);
points(x, y, cex=1, pch=1, col="blue");

legend(2.6, -2, c("PS Null", "BH Null", "yellow", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
#dev.off();


#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1alogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- TDeleuzepValue[,1];
x <- length(TRankedtStat[,1]):1;

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(0, max(y)),
     cex=1, pch=1, col="black");

y <- TBHpValue[,1];
x <- length(TRankedtStat[,1]):1;

points(x, y, cex=1, pch=1, col="red");

y <- TBHpValue[,11];
x <- length(TRankedtStat[,11]):1;

points(x, y, cex=1, pch=1, col="orange");

y <- TFreqpValue[,1];
x <- length(TRankedtStat[,1]):1;

points(x, y, cex=1, pch=1, col="yellow");

y <- TBHpValue[,13];
x <- length(TRankedtStat[,13]):1;

points(x, y, cex=1, pch=1, col="green");

y <- TBHpValue[,15];
x <- length(TRankedtStat[,15]):1;
points(x, y, cex=1, pch=1, col="blue");

legend(500, 1, c("PS Null", "BH Null", "Freq Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "yellow", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
#dev.off();


#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1blogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- TDeleuzepValue[,1];
x <- length(TRankedtStat[,1]):1;

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(300, max(x)),
                       ylim=c(0, 0.6),
     cex=1, pch=1, col="black");

y <- TBHpValue[,1];
x <- length(TRankedtStat[,1]):1;

points(x, y, cex=1, pch=1, col="red");

y <- TBHpValue[,11];
x <- length(TRankedtStat[,11]):1;

points(x, y, cex=1, pch=1, col="orange");

y <- TFreqpValue[,1];
x <- length(TRankedtStat[,1]):1;

points(x, y, cex=1, pch=1, col="yellow");

y <- TBHpValue[,13];
x <- length(TRankedtStat[,13]):1;

points(x, y, cex=1, pch=1, col="green");

y <- TBHpValue[,15];
x <- length(TRankedtStat[,15]):1;
points(x, y, cex=1, pch=1, col="blue");

legend(550, 0.6, c("PS Null", "BH Null", "Freq Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "yellow", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
#dev.off();


#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1alogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);


x11(width=11, height=8);

y <- TDeleuzepValue[,1];
x <- TRankedtStat[,1];

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(0, max(y)),
     cex=1, pch=1, col="black");

y <- TBHpValue[,1];
x <- TRankedtStat[,1];

points(x, y, cex=1, pch=1, col="red");

y <- TBHpValue[,11];
x <- TRankedtStat[,11];

points(x, y, cex=1, pch=1, col="orange");

y <- TFreqpValue[,1];
x <- TRankedtStat[,1];

points(x, y, cex=1, pch=1, col="yellow");

y <- TBHpValue[,13];
x <- TRankedtStat[,13];

points(x, y, cex=1, pch=1, col="green");

y <- TBHpValue[,15];
x <- TRankedtStat[,15];
points(x, y, cex=1, pch=1, col="blue");

legend(10, 1, c("PS Null", "BH Null", "Freq Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "yellow", "orange", "green", "blue"));

#c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 24, 29, 34, 39)
#dev.off();


#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1blogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- TDeleuzepValue[,1];
x <- TRankedtStat[,1];

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(2, max(x)),
                       ylim=c(0, 0.6),
     cex=1, pch=1, col="black");

y <- TBHpValue[,1];
x <- TRankedtStat[,1];

points(x, y, cex=1, pch=1, col="red");

y <- TBHpValue[,11];
x <- TRankedtStat[,11];

points(x, y, cex=1, pch=1, col="orange");

y <- TFreqpValue[,1];
x <- TRankedtStat[,1];

points(x, y, cex=1, pch=1, col="yellow");

y <- TBHpValue[,13];
x <- TRankedtStat[,13];

points(x, y, cex=1, pch=1, col="green");

y <- TBHpValue[,15];
x <- TRankedtStat[,15];
points(x, y, cex=1, pch=1, col="blue");

legend(4, 0.6, c("PS Null", "BH Null", "Freq Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "yellow", "orange", "green", "blue"));

#dev.off();

#postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/Fig1blogpvslogTBHNu.eps", sep=""),  horizontal=TRUE, onefile=TRUE);
 

x11(width=11, height=8);

y <- log10(TDeleuzepValue[,1]);
x <- log10(TRankedtStat[,1]);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(0, 1),
                       ylim=c(-1, 0),
     cex=1, pch=1, col="black");

y <- log10(TBHpValue[,1]);
x <- log10(TRankedtStat[,1]);

points(x, y, cex=1, pch=1, col="red");

y <- TBHpValue[,11];
x <- TRankedtStat[,11];

points(x, y, cex=1, pch=1, col="orange");

y <- TFreqpValue[,1];
x <- TRankedtStat[,1];

points(x, y, cex=1, pch=1, col="yellow");

y <- TBHpValue[,13];
x <- TRankedtStat[,13];

points(x, y, cex=1, pch=1, col="green");

y <- TBHpValue[,15];
x <- TRankedtStat[,15];
points(x, y, cex=1, pch=1, col="blue");

legend(4, 0.6, c("PS Null", "BH Null", "Freq Null", "BH nu=9", "BH nu=11", "BH nu=13"), cex=2, pch=1, col=c("black", "red", "yellow", "orange", "green", "blue"));

#dev.off();

            
postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/pvsTMicrodof3L14.eps", sep=""),  horizontal=TRUE, onefile=TRUE);


#x11(width=11, height=8);

Column <- 12;

y <- TSidakpValue[,Column];
x <- length(TRankedtStat[,Column]):1;

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(min(TDeleuzepValue[,Column]),
                              max(TDeleuzepValue[,Column])),
     cex=2, pch=1, col="red");


y <- TFreqpValue[,Column];
x <- length(TRankedtStat[,Column]):1;

points(x, y, cex=2, pch=1, col="blue");


y <- TBHpValue[,Column];
x <- length(TRankedtStat[,Column]):1;

points(x, y, cex=2, pch=1, col="green");

y <- TDeleuzepValue[,Column];
x <- length(TRankedtStat[,1]):1;

points(x, y, cex=2, pch=1, col="black");


legend(7, 1, c("Deleuze", "BH", "Freq", "Sidak"), cex=2, pch=1, col=c("black", "green", "blue", "red"));

dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/pvsTMicrodof3", as.character(Lambda), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);


#x11(width=11, height=8);


y <- SidakpValue;
x <- RankedtStat;

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(min(DeleuzepValue), max(y)),
     cex=2, pch=1, col="red");


y <- FreqpValue;
x <- RankedtStat;

points(x, y, cex=2, pch=1, col="blue");


y <- BHpValue;
x <- RankedtStat;

points(x, y, cex=2, pch=1, col="green");

y <- DeleuzepValue;
x <- RankedtStat;

points(x, y, cex=2, pch=1, col="black");


legend(12.6, 1, c("Deleuze", "BH", "Freq", "Sidak"), cex=2, pch=1, col=c("black", "green", "blue", "red"));

dev.off()



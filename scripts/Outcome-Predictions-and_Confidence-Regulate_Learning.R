
library(R.matlab)
library(languageR)
library(lme4)
library(MASS)
library(ggplot2) #cookbook for R/ Graphs
library(hexbin)
library(memisc)
library(reshape)
library(reshape2) #melt and cast -> restructure and aggregate data
library(data.table)
library(coin) #for permutation tests
library(psych)
library(doBy)
library(heplots)
library(plyr) #necessary for ddply
library(matrixStats) 
library(foreign) 
library(Hmisc)
library(lmerTest)
library (stringr)
library("piecewiseSEM")
library(gdata)
library(effects)
library(gridExtra)
library(ggExtra)
library(simr)
library(sjPlot)

### load logfiles  
setwd("/Users/Romy/Dropbox (Brown)/Confidence/logfiles")

## get all relevant files
fl <- list.files()

for (i in 1:length (fl))
{ vp   = (substring2(fl[i], 1,2))
  
  vpn  <- rep(vp, 250)
  
  tmp2 <- (read.table(fl[i] , header = TRUE, sep = "", #dec = ",", bad... kills everything!
                   fill = TRUE))
  tmp <- cbind(vpn,tmp2)
  if (i==1)
  {a1 <- tmp}
  else {a1 <- rbind(a1,tmp)}
}


## data transformations
a1 <- subset(a1, vpn !="04" & vpn !="08") # exclude participants who did the task incorrectly
str(a1)
a1$vpn <- factor(a1$vpn)
a1$Trial  <- as.numeric(a1$Trial)
a1$RT  <- as.numeric(a1$RT)
a1$Difference  <- as.numeric(a1$Difference) # objective error
a1$sDifference  <- as.numeric(a1$sDifference) # received feedback (known error)
a1$Direction  <- as.numeric(a1$Direction) # scaled prediction
a1$osDifference <-  a1$sDifference 
a1$oDirection  <- a1$Direction 
a1$Confidence  <- as.numeric(a1$Confidence)
a1$dssDifference <- round(a1$sDifference/4) # Feedback scaled to same resolution as prediction
a1$acE <- abs(a1$Direction - a1$dssDifference) *32 # SPE computed based on same scale actual and predicted outcome 
a1$sPE <- (abs(a1$Direction) - abs(a1$dssDifference)) *32 # RPE computed based on same scale actual and predicted outcome 
a1$absE <- abs(a1$Difference) # objective error magnitude
a1$sDifference  <- a1$dssDifference * 32 # get same scale for prediction, feedback and actual performance error
a1$Direction  <- a1$Direction*32 # get same scale for prediction, feedback and actual performance error
a1$asDifference  <- abs(a1$sDifference) # Error Magnitude irrespective of error direction
a1$aDirection  <- abs(a1$Direction) # same for predicted error


### add blocks

a1$block  <- a1$Trial
a1$block[a1$Trial<51]  <- 1
a1$block[a1$Trial>50 &a1$Trial<101 ]  <- 2
a1$block[a1$Trial>100 &a1$Trial<151 ]  <- 3
a1$block[a1$Trial>150 &a1$Trial<201 ]  <- 4
a1$block[a1$Trial>200]  <- 5

a1$fblock  <- as.factor(a1$block)
contrasts(a1$fblock) <- contr.sdif(5)

# rescale variables to appoximately similar ranges (e.g. time measures to seconds instead of ms)
a1$sTrial <- (a1$Trial-1)/100
a1$oConfidence <- a1$Confidence/31 ## Confidence 
a1$assDifference <- a1$asDifference/1000
a1$ssDifference <- a1$sDifference/1000
a1$sDirection <- a1$Direction/1000
a1$asDirection <- a1$aDirection/1000
a1$sacE <- a1$acE/1000 # this is SPE
a1$ssacE <- scale(a1$sacE, center=TRUE, scale=FALSE) # center OPE
a1$ssPE <- a1$sPE/1000 # this is RPE - I'm a terrible person
a1$block.c <-scale(a1$block, center=TRUE, scale=FALSE)/2 # center and scale block
a1$sassDifference <- scale(a1$assDifference, center=TRUE, scale=FALSE) # center error magnitude


a1$f3Conf  <- cut(a1$oConfidence,
                  quantile(a1$oConfidence, seq(0, 1, 1/3), na.rm=TRUE), include.lowest=T,
                  labels=c("low","medium","high"))

## for Fig 4 
a1$f5Conf  <- cut(a1$oConfidence,
                  quantile(a1$oConfidence, seq(0, 1, 1/5), na.rm=TRUE), include.lowest=T,
                  labels=c("verylow","low","medium","high","veryhigh"))

a1[a1$RT>6000,4:34] <-  NA_real_ # set values of variables for outliers to NAs
#writeMat('/Volumes/daten/romy/Confidence/Export/behav_n.mat',a1=a1) # export behavior matrix

###### Import ERPs

# Peak to Peak FRN at FCz
input_file = '/Users/Romy/Dropbox (Brown)/Confidence/Export/P2PFCZ.mat'
P2PFCZin = readMat(input_file)
a1$P2PFCZ  <- P2PFCZin$P2PFCZ

# P3a
input_file = '/Users/Romy/Dropbox (Brown)/Confidence/Export/P3.mat'
P3 = readMat(input_file)

P3 = P3$P3 #--> da muss man die Variable angeben, als die man in Matlab gespeichert hat
P3 = as.data.frame(P3)
colnames(P3)[c( 1, 2, 3, 4, 5, 
                6, 7, 8, 9, 10, 
                11, 12, 13, 14, 15, 
                16, 17, 18, 19, 20, 
                21, 22, 23, 24, 25, 
                26, 27, 28, 29, 30,
                31, 32, 33, 34, 35,
                36, 37, 38, 39, 40,
                41, 42, 43, 44, 45,
                46, 47, 48, 49, 50,
                51, 52, 53, 54, 55,
                56, 57, 58, 59, 60,
                61, 62, 63, 64, 65)]=c('FP1P3', 'FPzP3', 'FP2P3', 'AF3P3', 
                                       'AFzP3', 'AF4P3', 'F7P3', 'F5P3', 'F3P3', 
                                       'F1P3','FzP3', 'F2P3', 'F4P3', 'F6P3', 
                                       'F8P3', 'FT9P3', 'FT7P3', 'FC5P3', 'FC3P3', 
                                       'FC1P3', 'FCzP3', 'FC2P3', 'FC4P3', 'FC6P3', 
                                       'FT8P3', 'FT10P3', 'T7P3', 'C5P3', 'C3P3', 
                                       'C1P3','C2P3', 'C4P3', 'C6P3', 'T8P3', 
                                       'TP9P3', 'TP7P3', 'CP5P3', 'CP3P3', 'CP1P3', 
                                       'CPzP3','CP2P3', 'CP4P3', 'CP6P3', 'TP8P3', 
                                       'TP10P3', 'P7P3', 'P5P3', 'P3P3', 'P1P3', 
                                       'PzP3', 'P2P3', 'P4P3', 'P6P3', 'P8P3', 
                                       'PO3P3', 'POzP3', 'PO4P3', 'O1P3', 'OzP3', 
                                       'O2P3', 'LO1P3', 'IO1P3', 'IO2P3', 'LO2P3', 
                                       'CzP3')
a1$P3 <- apply(cbind (P3$F1P3, P3$FzP3, P3$F2P3,P3$FC1P3,P3$FCzP3, P3$FC2P3, P3$C1P3, P3$CzP3, P3$C2P3), 1, mean)# 

# P3b
input_file = '/Users/Romy/Dropbox (Brown)/Confidence/Export/P3b.mat'
P3b = readMat(input_file)

P3b = P3b$P3b #--> da muss man die Variable angeben, als die man in Matlab gespeichert hat
P3b = as.data.frame(P3b)
colnames(P3b)[c( 1, 2, 3, 4, 5, 
                 6, 7, 8, 9, 10, 
                 11, 12, 13, 14, 15, 
                 16, 17, 18, 19, 20, 
                 21, 22, 23, 24, 25, 
                 26, 27, 28, 29, 30,
                 31, 32, 33, 34, 35,
                 36, 37, 38, 39, 40,
                 41, 42, 43, 44, 45,
                 46, 47, 48, 49, 50,
                 51, 52, 53, 54, 55,
                 56, 57, 58, 59, 60,
                 61, 62, 63, 64, 65)]=c('FP1P3b', 'FPzP3b', 'FP2P3b', 'AF3P3b', 
                                        'AFzP3b', 'AF4P3b', 'F7P3b', 'F5P3b', 'F3P3b', 
                                        'F1P3b','FzP3b', 'F2P3b', 'F4P3b', 'F6P3b', 
                                        'F8P3b', 'FT9P3b', 'FT7P3b', 'FC5P3b', 'FC3P3b', 
                                        'FC1P3b', 'FCzP3b', 'FC2P3b', 'FC4P3b', 'FC6P3b', 
                                        'FT8P3b', 'FT10P3b', 'T7P3b', 'C5P3b', 'C3P3b', 
                                        'C1P3b','C2P3b', 'C4P3b', 'C6P3b', 'T8P3b', 
                                        'TP9P3b', 'TP7P3b', 'CP5P3b', 'CP3P3b', 'CP1P3b', 
                                        'CPzP3b','CP2P3b', 'CP4P3b', 'CP6P3b', 'TP8P3b', 
                                        'TP10P3b', 'P7P3b', 'P5P3b', 'P3P3b', 'P1P3b', 
                                        'PzP3b', 'P2P3b', 'P4P3b', 'P6P3b', 'P8P3b', 
                                        'PO3P3b', 'POzP3b', 'PO4P3b', 'O1P3b', 'OzP3b', 
                                        'O2P3b', 'LO1P3b', 'IO1P3b', 'IO2P3b', 'LO2P3b', 
                                        'CzP3b')
a1$P3b <- apply(cbind (P3b$CP1P3b, P3b$CPzP3b, P3b$CP2P3b,P3b$P1P3b,P3b$PzP3b, P3b$P2P3b,  P3b$PO3P3b,  P3b$POzP3b,  P3b$PO4P3b), 1, mean)# 

## quick-correct here
as <- a1 # safe for latter matching with EEG data frames

## continue here after quick-correct
a1 <- a1[!is.na(a1$Difference),] # omit NAs

### plot Difference by Block
library(gridExtra)


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p2  <- ggplot(a1, aes(x=Direction, y=sDifference, group=f3Conf, color=f3Conf, shape=f3Conf,fill=f3Conf)) +
  geom_point(size=1, alpha=0.2) +scale_color_manual(name="Confidence", values=cbPalette) +scale_fill_manual(name="Confidence", values=cbPalette) +  # Use hollow circles , method=lm
  geom_smooth(method=lm, size = 1,  inherit.aes = TRUE) + ggtitle("B")+ scale_shape_manual(name="Confidence",values = c(15, 19, 17) ) + #+stat_smooth(aes(group = 1), formula = y ~ x, se = TRUE, colour = "#666666") 
  theme_bw()+xlab("predicted Error [ms]") + ylab("actual Error [ms]")+theme(plot.title = element_text(hjust = 0))+ geom_vline(xintercept = 0, linetype="dotted")+ geom_hline(yintercept=0, linetype="dotted")+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position="bottom")

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/Distribution_PE_AE_Conf.pdf", width = 5, height = 5)#, units = 'cm', res = 200, compression = 'lzw'
ggMarginal(p2,groupColour = TRUE, groupFill = TRUE)
dev.off()

###### CHANGE ACROSS BLOCKS!!!!! FIGURE S 2
pbeff  <- ggplot(a1, aes(x=Direction, y=sDifference, group=f3Conf, color=f3Conf, shape=f3Conf,fill=f3Conf)) + facet_wrap(~ fblock, nrow=1)+
  geom_point(size=1, alpha=0.2) +scale_color_manual(name="Confidence", values=cbPalette) +scale_fill_manual(name="Confidence", values=cbPalette) +  # Use hollow circles , method=lm
  geom_smooth(method=lm, size = 1,  inherit.aes = TRUE) + ggtitle("C")+ scale_shape_manual(name="Confidence",values = c(15, 19, 17) ) + #+stat_smooth(aes(group = 1), formula = y ~ x, se = TRUE, colour = "#666666") 
  theme_bw()+xlab("predicted Error [ms]") + ylab("actual Error [ms]")+theme(plot.title = element_text(hjust = 0))+ geom_vline(xintercept = 0, linetype="dotted")+ geom_hline(yintercept=0, linetype="dotted")+ theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ theme(legend.position="bottom")+theme(strip.background = element_rect(colour="white",fill="white"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/Outcome_Prediction_block.pdf", width = 10, height = 5)#, units = 'cm', res = 200, compression = 'lzw'
pbeff
dev.off()

### This is for the individual differences 
### compute correlation coefficients (Conf OPE) and average performance by subject (and also block)

for (i in levels(a1$vpn) ) # if I say print(i) it gives me the subject ids     {print(j)} returns trial numbers
{ print(i)
  a1$avE[a1$vpn==i]  <- mean(a1$absE[a1$vpn==i], na.rm=TRUE)/1000
  a1$avOPE[a1$vpn==i]  <- mean(a1$acE[a1$vpn==i], na.rm=TRUE)/1000
  a1$avConf[a1$vpn==i]  <- mean(a1$oConfidence[a1$vpn==i], na.rm=TRUE)
  for (j in levels(a1$fblock))
    {print(i) 
    print(j)
    a1$avEb[a1$vpn==i& a1$fblock==j]  <- mean(a1$absE[a1$vpn==i& a1$fblock==j], na.rm=TRUE)/1000
  }
  outnpb <-partial.r(a1[a1$vpn==i,c("oConfidence", "acE","avEb")], c(1,2), c(3)) #,c(19, 9,36)
  a1$corcoefslmpb[a1$vpn==i] <- outnpb[2,1] 
  
}

a1$savE  <- scale(a1$avE, center=TRUE, scale=FALSE)
a1$savEb  <- scale(a1$avEb, center=TRUE, scale=FALSE)
a1$savEOPE  <- scale(a1$avOPE, center=TRUE, scale=FALSE)
a1$savConf  <- scale(a1$avConf, center=TRUE, scale=FALSE)/10
## for inset plot 
subcorrs <- ddply(a1, .(vpn), summarise,
              indavE      = mean(avE),
              indcorrEb   = mean(corcoefslmpb))


## average error magnitude as a function of Confidence calibration part 1 of figure 2 (inset)
quantile(subcorrs$indcorrEb*-1, seq(0, 1, 1/3))
pcdb <- ggplot(data=subcorrs, aes(x=indcorrEb*-1, y=indavE))+stat_smooth(aes(group = 1), method=lm, formula = y ~ x, se = TRUE, colour = "black")+ geom_vline(xintercept=0, linetype="dotted")+ geom_point(aes(colour = cut(indcorrEb*-1, c(-Inf, 0.104, 0.201, Inf)))) +theme_bw(12)+ #facet_wrap(~ssConfidence, nrow=1)+#geom_ribbon(data=IA, aes(x=block.c, max = upper, min = lower,alpha=0.1, fill=sDirection),show.legend =FALSE)
  xlab("confidence accuracy") + ylab("average error magnitude [s]") + theme( panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none")+
  scale_color_manual(name = "Confidence accuracy",
                     values = c("(-Inf,0.104]" = "#999999",
                                "(0.104,0.201]" =  "#E69F00",
                                "(0.201, Inf]" = "#56B4E9"),
                     labels = c("low", "medium", "high"))+ theme(plot.title = element_text(hjust = 0))

####################### LMMs######################
######################################################

##################### Performance Prediction !!!

### Prediction

print (summary(Emod0 <- lmer(sDifference~ (sDirection*block.c)*oConfidence+(oConfidence*sDirection+block.c|vpn), a1, 
                             REML=FALSE))) 

sv1_max <- svd(getME(Emod0, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)  ### final model


## Prediction error; decreases significantly across blocks and from block 1 to 2
# linear trend
print (summary(OPEmod0s <- lmer(acE~ block.c+(block.c|vpn), a1, 
                                REML=FALSE))) 

sv1_max <- svd(getME(OPEmod0s, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

## by block comparison
detach("package:lmerTest", unload = TRUE)

print (summary(OPEmod0s <- lmer(acE~ fblock+(block.c|vpn), a1,
                                REML=FALSE))) 

sv1_max <- svd(getME(OPEmod0s, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

tab_model(OPEmod0s)

#################### Confidence Calibration model ############
a1$ConfCali <- scale(a1$corcoefslmpb*-1, scale = FALSE, center= TRUE)

print (summary(Esmod0lMetaf <- lmer(log(absE+1)~ ((ConfCali)*sTrial + I(sTrial*sTrial))+(sTrial |vpn), a1, 
                                    REML=FALSE))) 

sv1_max <- svd(getME(Esmod0lMetaf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1) 


### figure 2D
eff_df <- Effect(c("sTrial","ConfCali"), Esmod0lMetaf, xlevels=list(ConfCali =c(-0.13, 0, 0.13), sTrial =seq(min(a1$sTrial), max(a1$sTrial), 0.01)))
IA <- as.data.frame(eff_df)

IA$ConfCali <- as.factor(IA$ConfCali)
IA$ConfCali<- revalue(IA$ConfCali, c("-0.13"="low", "0"="medium", "0.13"= "high"))
IA$sTrial <- (IA$sTrial+min(a1$sTrial))*100+1
library(RColorBrewer)
darkcolsg <- brewer.pal(8, "Greens")
darkcolssubg <-  darkcolsg[c(4:8)]

p1 <- ggplot(data=IA, aes(x=sTrial, y=fit , color= ConfCali )) +scale_colour_manual(name="Confidence\nCalibration", values=c("#CC79A7", "#FFEC8B", "#79CC9E")) +theme_bw(12)+ geom_ribbon(data=IA, aes(x=sTrial, max = fit + se, min = fit- se, fill = ConfCali),alpha=0.1, inherit.aes = FALSE)+ geom_line(size= 1.5)+ #ylim(4.3, 5.2)+
  scale_fill_manual(name="Confidence\nCalibration", values=c("#CC79A7", "#FFEC8B", "#79CC9E"))+ xlab("Trial") + ylab("log Error Magnitude") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position=c(0.25, 0.2),legend.background=element_blank()) #+#c(0.8, 0.2)

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/ConfCaliNewTrial.pdf", width = 5, height = 5)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()


##### ERPs ####
# FRN 

print (summary(FRNmod0cf <- lmer(P2PFCZ~ block.c*(oConfidence*(ssPE+ssacE +sassDifference))+(sassDifference+block.c|vpn), a1, 
                               REML=FALSE))) 

sv1_max <- svd(getME(FRNmod0cf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1) 

print (summary(FRNmod1cf <- lmer(P2PFCZ~ oConfidence+(ssPE + ssacE+ sassDifference)+block.c+(sassDifference+block.c|vpn), a1, 
                               REML=FALSE))) 

sv1_max <- svd(getME(FRNmod1cf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1) 

anova(FRNmod1cf, FRNmod0cf)

## Fig 3D 
eff_df <- Effect(c("ssPE"), FRNmod1cf)
IA <- as.data.frame(eff_df)
p1 <-(ggplot(data=IA, aes(x=ssPE, y=fit))+ geom_line()+theme_bw(12)+ geom_ribbon(data=IA, aes(x=ssPE, max = upper, min = lower,alpha=0.1), fill="#CCCCCC",show.legend =FALSE, inherit.aes = TRUE) +
    xlab("RPE") + ylab("FRN amplitude [µV]") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/RPE_FRN.pdf", width =3, height = 3)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()



# P3a

print (summary(P3mod0mboNcf <- lmer(P3~  block.c*((oConfidence)*(ssacE+ssPE+sassDifference))+(ssacE|vpn), a1, #[a1$sPE>0,]
                                  REML=FALSE)))
sv1_max <- svd(getME(P3mod0mboNcf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

print (summary(P3mod1mboNcf <- lmer(P3~  (oConfidence)+block.c*(ssacE)+ssPE+sassDifference+(ssacE|vpn), a1, #[a1$sPE>0,]
                                  REML=FALSE)))
sv1_max <- svd(getME(P3mod1mboNcf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

anova(P3mod1mboNcf, P3mod0mboNcf)

print (summary(P3mod2mboNcf <- lmer(P3~  (oConfidence)+ block.c+(ssacE)+ssPE+sassDifference+(ssacE|vpn), a1, #[a1$sPE>0,]
                                    REML=FALSE)))
sv1_max <- svd(getME(P3mod2mboNcf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

anova(P3mod2mboNcf, P3mod1mboNcf)
anova(P3mod2mboNcf, P3mod0mboNcf)

### Fig 3J:

eff_df <- Effect(c("ssacE"), P3mod2mboNcf, xlevels=list(ssacE=seq(min(a1$ssacE), max(a1$ssacE), 0.1)))
IA <- as.data.frame(eff_df)#cbind(eff_df$x, eff_df$fit, eff_df$lower, eff_df$upper)
IA$ssacE <- IA$ssacE - min(a1$ssacE)
p1 <-(ggplot(data=IA, aes(x=ssacE, y=fit))+ geom_line()+theme_bw(12)+ geom_ribbon(data=IA, aes(x=ssacE, max = upper, min = lower,alpha=0.1), fill="#CCCCCC",show.legend =FALSE, inherit.aes = TRUE) +
        xlab("OPE") + ylab("P3a amplitude [µV]") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/OPE_P3a.pdf", width =3, height = 3)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()

## Fig 3 K:
eff_df <- Effect(c("oConfidence"), P3mod2mboNcf, xlevels=list(oConfidence=seq(min(a1$oConfidence), max(a1$oConfidence), 0.1)))
IA <- as.data.frame(eff_df)#cbind(eff_df$x, eff_df$fit, eff_df$lower, eff_df$upper)
IA$oConfidence <- IA$oConfidence - min(a1$oConfidence)
p1 <-(ggplot(data=IA, aes(x=oConfidence, y=fit))+ geom_line()+theme_bw(12)+ geom_ribbon(data=IA, aes(x=oConfidence, max = upper, min = lower,alpha=0.1), fill="#CCCCCC",show.legend =FALSE, inherit.aes = TRUE) +
        xlab("Confidence") + ylab("P3a amplitude [µV]") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/Confidence_P3a.pdf", width =3, height = 3)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()



#### P3b

print (summary(P3bmod0mboNcf <- lmer(P3b~ ( block.c*(oConfidence)*(ssacE+ssPE+sassDifference))+(ssacE+ssPE+oConfidence|vpn), a1, #[a1$sPE>0,][!a1$block==1,]
                                   REML=FALSE)))
sv1_max <- svd(getME(P3bmod0mboNcf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

print (summary(P3bmod1mboNcf <- lmer(P3b~  block.c*(oConfidence*ssPE+ ssacE)+sassDifference+(ssacE+ssPE+oConfidence|vpn), a1, #[a1$sPE>0,]
                                   REML=FALSE)))
sv1_max <- svd(getME(P3bmod1mboNcf, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

anova(P3bmod0mboNcf, P3bmod1mboNcf)


## Fig 4 C-F
eff_df <- Effect(c("ssacE"), P3bmod1mboNcf, xlevels=list(ssacE=seq(min(a1$ssacE), max(a1$ssacE), 0.1)))
IA <- as.data.frame(eff_df)#cbind(eff_df$x, eff_df$fit, eff_df$lower, eff_df$upper)
IA$ssacE <- IA$ssacE - min(a1$ssacE)
p1 <-(ggplot(data=IA, aes(x=ssacE, y=fit))+ geom_line()+theme_bw(12)+ geom_ribbon(data=IA, aes(x=ssacE, max = upper, min = lower,alpha=0.1), fill="#CCCCCC",show.legend =FALSE, inherit.aes = TRUE) +
        xlab("OPE") + ylab("P3b amplitude [µV]") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/OPE_P3b.pdf", width =3, height = 3)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()

eff_df <- Effect(c("ssPE"), P3bmod1mboNcf, xlevels=list(ssPE=seq(min(a1$ssPE), max(a1$ssPE), 0.1)))
IA <- as.data.frame(eff_df)#cbind(eff_df$x, eff_df$fit, eff_df$lower, eff_df$upper)
#IA$oConfidence <- IA$oConfidence - min(a1$oConfidence)
p1 <-(ggplot(data=IA, aes(x=ssPE, y=fit))+ geom_line()+theme_bw(12)+ geom_ribbon(data=IA, aes(x=ssPE, max = upper, min = lower,alpha=0.1), fill="#CCCCCC",show.legend =FALSE, inherit.aes = TRUE) +
        xlab("RPE") + ylab("P3b amplitude [µV]") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/RPE_P3b.pdf", width =3, height = 3)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()

eff_df <- Effect(c("sassDifference"), P3bmod1mboNcf, xlevels=list(sassDifference=seq(min(a1$sassDifference), max(a1$sassDifference), 0.1)))
IA <- as.data.frame(eff_df)#cbind(eff_df$x, eff_df$fit, eff_df$lower, eff_df$upper)
IA$sassDifference <- IA$sassDifference - min(a1$sassDifference)
p1 <-(ggplot(data=IA, aes(x=sassDifference, y=fit))+ geom_line()+theme_bw(12)+ geom_ribbon(data=IA, aes(x=sassDifference, max = upper, min = lower,alpha=0.1), fill="#CCCCCC",show.legend =FALSE, inherit.aes = TRUE) +
        xlab("Absolute Error") + ylab("P3b amplitude [µV]") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="none"))

pdf("/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Science/MS/Figures/EM_P3b.pdf", width =3, height = 3)#, units = 'cm', res = 200, compression = 'lzw'
p1
dev.off()



### post hoc tests for RPE by block by Confidence interaction


contrasts(a1$f3Conf) <- contr.sdif(3)
contrasts(a1$f5Conf) <- contr.sdif(5)

detach("package:lmerTest", unload = TRUE)
print (summary(P3bmod0mboN <- lmer(P3b~ fblock/( f3Conf/(ssPE))+ssacE+ssacE:block.c+sassDifference+(ssacE|vpn), a1, #[a1$sPE>0,]
                                   REML=FALSE)))
sv1_max <- svd(getME(P3bmod0mboN, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)

tab_model(P3bmod0mboN)

## for significance within 5-tiles for figure
print (summary(P3bmod0mboN <- lmer(P3b~ fblock/( f5Conf/(ssPE))+ssacE+ssacE:block.c+sassDifference+(ssacE+ssPE+oConfidence|vpn), a1, #[a1$sPE>0,]
                                   REML=FALSE)))
sv1_max <- svd(getME(P3bmod0mboN, "Tlist")[[1]]) 
sv1_max$d
round(sv1_max$d^2/sum(sv1_max$d^2)*100, 1)


###################### Generating effect topos for plotting ########################
library(lmerTest)
ap <- as

ap  <- cbind(as, P3)
#ap$ssConfidence  <- scale(ap$ssConfidence, scale=TRUE)
plotEmat <- as.data.frame(matrix(, nrow = 65, ncol=6))
plotTmat <- as.data.frame(matrix(, nrow = 65, ncol=6))

j=1
for (i in (ncol(as)+1):(ncol(ap))){
ap$DV  <- ap[,i]
Plotmod <- lmer(DV~ (oConfidence)+block.c+(ssacE)+ssPE+sassDifference+(ssacE|vpn), ap, #[a1$sPE>0,]
                REML=FALSE)
k <- summary(Plotmod)
if (i ==(ncol(as)+1)){
  f <- dimnames(k$coefficients)[[1]]
  colnames(plotEmat) <-  f
  colnames(plotTmat) <-  f}
  plotEmat[j,] <- k$coefficients[,1]
  plotTmat[j,] <- k$coefficients[,4]
j=j+1
}

# colnames(plotEmat) <- sedit(colnames(plotEmat), ":", "by")
# colnames(plotTmat) <- sedit(colnames(plotEmat), ":", "by")
colnames(plotEmat) <- sedit(colnames(plotEmat), "block.c", "block")
colnames(plotTmat) <- sedit(colnames(plotEmat),  "block.c", "block")

writeMat('/Users/Romy/Dropbox (Brown)/Confidence/Export/plotes.mat',Emat=plotEmat)
writeMat('/Users/Romy/Dropbox (Brown)/Confidence/Export/plotts.mat',Tmat=plotTmat)


ap  <- cbind(as, P3b)
#ap$ssConfidence  <- scale(ap$ssConfidence, scale=TRUE)
plotEmat2 <- as.data.frame(matrix(, nrow = 65, ncol=11))
plotTmat2 <- as.data.frame(matrix(, nrow = 65, ncol=11))

j=1
for (i in (ncol(as)+1):(ncol(ap))){
  ap$DV  <- ap[,i]
  Plotmod <- lmer(DV~ block.c*(oConfidence*ssPE+ ssacE)+sassDifference+(ssacE|vpn), ap, #[a1$sPE>0,]
                  REML=FALSE)
  k <- summary(Plotmod)
  if (i ==(ncol(as)+1)){
    f <- dimnames(k$coefficients)[[1]]
    colnames(plotEmat2) <-  f
    colnames(plotTmat2) <-  f}
  plotEmat2[j,] <- k$coefficients[,1]
  plotTmat2[j,] <- k$coefficients[,4]
  j=j+1
}

colnames(plotEmat2) <- sedit(colnames(plotEmat2), ":", "by")
colnames(plotTmat2) <- sedit(colnames(plotEmat2), ":", "by")
colnames(plotEmat2) <- sedit(colnames(plotEmat2), "block.c", "block")
colnames(plotTmat2) <- sedit(colnames(plotEmat2),  "block.c", "block")

writeMat('/Users/Romy/Dropbox (Brown)/Confidence/Export/plotes2.mat',Emat2=plotEmat2)
writeMat('/Users/Romy/Dropbox (Brown)/Confidence/Export/plotts2.mat',Tmat2=plotTmat2)


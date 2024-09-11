# clean current work space
rm(list=ls(all=T))
options(stringsAsFactors = F)   # no automatic data transformation
options("scipen" = 100, "digits" = 4) # suppress math annotation

library(descr)
library(ggplot2)
library(survival)
library(survminer)
library(lubridate)
library(ggpubr)
library(dplyr)

setwd('~/Siontis research/Amyloid_NSVT/')
ansvt <- read.csv('Copy of NSVT_WS.csv', na.strings = c("", "NA"))
ansvt <- ansvt[1:217,]
View(ansvt)

freq(ansvt$VT_VF_WS)
vtvf <- subset(x = ansvt,subset = ansvt$VT_VF_WS == 1)
vtvfsub <- vtvf[,c(2,3,16,18:21,25:28,31,51,61,63:66)]
View(vtvfsub)
View(vtvf)

# Total events per person-year
py <- sum(ansvt$Fig.2) / 365.25
py
11/py * 100
217 * 27 /12

poisson.test(x = 11, T = 602.2,conf.level = 0.95)


# Total followup time in months
summary(ansvt$Fig.3b)/30.437

# Time index to holter
summary(ansvt$Dxtoindex) / 30.437


freq(ansvt$Dxtoindex)

# Of the 82 who had holter first...
negindex <- subset(x = ansvt$Dxtoindex, subset = ansvt$Dxtoindex < 0)
str(negindex)
summary(negindex) / 30.437

freq(ansvt$Dxtoindex)

# in months
summary(negindex) / 30.437
# Number whose holter first was >1 month before dx
a <- negindex / 30.437 < -1
sum(a)

# Time index to dx for table 1
summary(nsvt$Dxtoindex) / 30.437
summary(nonsvt$Dxtoindex) / 30.437

a <- nsvt$Dxtoindex / 30.437
b <- nsvt$Dxtoindex / 30.437
t.test(a,b)

## for double checking dx to index ##
test <- subset(x = ansvt, subset = ansvt$Dxtoindex < -99)
View(test)


test[,c(2,28,29)]

vtvf <- subset(x = ansvt,subset = ansvt$VT_VF_WS == 1)
vtvf$LGE.extensive
vtvf$LGE2nd

# Key variaVE..# Key variables
  freq(ansvt$NSVT) # 116 had NSVT
  freq(ansvt$Device) # Prior device
  freq(ansvt$ICD.baseline) # 6 had a device at baseline

  # Key variables of devices
  freq(ansvt$ICD..all.) # 21 had an ICD
  freq(ansvt$ICD.FU) # 15 got ICD in f/u
  freq(ansvt$Loop.recorder) # 15 had a loop recorder
  freq(ansvt$pacemaker..all.) # 59 had a device of some sort
  freq(ansvt$PPM.baseline) # 17 had PPM
  freq(ansvt$PPM.during.FU) # 41 got a ppm in f/u

  # Figure 1 - variables not normalized
  freq(ansvt$Runs)
  freq(ansvt$Longest)

  # outcomes
  freq(ansvt$Mortality) # 46 died
  freq(ansvt$vt_vf) # 18 had a VT or VF
  freq(ansvt$vt_vf.1) # 15 had VT or VF in repeat chart review
  freq(ansvt$Composite) # 18 had one, 35 had both
  freq(ansvt$VT.VF.and.Mortality) # 53 met the composite endpoint

  # Key dates 
  freq(ansvt$diagnosis.date.Index) # months from diagnosis to index
  freq(ansvt$Last.FU..primary.)

  
#####
# Column 16 is NSVT
# 17 is VT/VF - old version with 18
# 59 is VT/VF - new version with 15
# Mortality is 51
# index date is 28
# days to primary outcome is 64
# Days to secondary outcome is 65
  
# Review table 1
  nsvt <- subset(x = ansvt, subset = ansvt$NSVT == 1)
  nonsvt <- subset(x = ansvt, subset = ansvt$NSVT == 0)
  dim(nsvt)
  dim(nonsvt)
  
        ### Replace variable names so don't have to copy/paste ###
        ## NSVT = VE >= 0.5%
        ## No NSVT = VE < 0.5%
        freq(ansvt$VE..)
        ve <- ansvt$VE..
        output <- rep(0,217)
        for (i in 1:217) {
          if (ve[i] >= 0.5) {
            output[i] <- 1
          } else {
            output[i] <- 0
          }
        }
        freq(output)
        nsvt <- ansvt[output == 1,]
        nonsvt <- ansvt[output == 0,]
        dim(nsvt)
        dim(nonsvt)
        
        
        
        ### Replace variable names for EF so don't have to copy/paste
        ## NSVT = EF </= 35 ##
        ## No NSVT = EF > 35 ##
        ef35 <- ansvt$EF
        output <- rep(0,217)
        freq(ansvt$EF) # Zero NAs
        
        for (i in 1:217) {
          if (ef35[i] <= 35) {
            output[i] <- 1
          } else {
            output[i] <- 0
          }
        }
        ef35
        output
        nsvt <- ansvt[output == 1,]
        nonsvt <- ansvt[output == 0,]
        dim(nsvt)
        dim(nonsvt)
        
        ### Replace variable names for EF so don't have to copy/paste
        ## NSVT = Mayo stage 2/3 ##
        ## No NSVT = Mayo stage 1 ##
        ms <- ansvt$Mayo.stage
        output <- rep(0,217)

        for (i in 1:217) {
          if (is.na(ms[i])) {
            output[i] <- 2
          } else if (ms[i] >= 2) {
            output[i] <- 1
          } else {
            output[i] <- 0
          }
        }
        nsvt <- ansvt[output == 1,]
        nonsvt <- ansvt[output == 0,]
        dim(nsvt)      
        dim(nonsvt)  

        
    # NSVT
        freq(ansvt$NSVT)
        freq(nsvt$NSVT)
        freq(nonsvt$NSVT)
        
        # By VE
        nsvtve <- matrix(c(83, 33, 28, 73), byrow = T, ncol = 2, nrow = 2)
        colnames(nsvtve) <- c("NSVT", "No NSVT")
        rownames(nsvtve) <- c("Male sex", "Female")
        nsvtve
        chisq.test(nsvtve)
        fisher.test(nsvtve)
        
        # By EF
        nsvtef <- matrix(c(15, 101, 9, 92), byrow = T, ncol = 2, nrow = 2)
        colnames(nsvtef) <- c("NSVT", "No NSVT")
        rownames(nsvtef) <- c("Male sex", "Female")
        nsvtef
        chisq.test(nsvtef)
        
        # By ms
        nsvtms <- matrix(c(62, 47, 43, 42), byrow = T, ncol = 2, nrow = 2)
        colnames(nsvtms) <- c("NSVT", "No NSVT")
        rownames(nsvtms) <- c("Male sex", "Female")
        nsvtms
        chisq.test(nsvtms)
            

    # Age
  summary(nsvt$Age)
  freq(nsvt$Age)
  summary(nonsvt$Age)  
  summary(ansvt$Age)
  wilcox.test(x = nsvt$Age, y = nonsvt$Age)
  ks.test(x = nsvt$Age, y = nonsvt$Age)
  t.test(x = nsvt$Age,y = nonsvt$Age)

  # Sex
  freq(nsvt$Sex)
  freq(nonsvt$Sex)
  freq(ansvt$Sex)
  
    # NSVT
    sex <- matrix(c(113, 92, 3, 9), byrow = T, ncol = 2, nrow = 2)
    colnames(sex) <- c("NSVT", "No NSVT")
    rownames(sex) <- c("Male sex", "Female")
    sex
    chisq.test(sex)
    fisher.test(sex)
    
    # VE
    sex <- matrix(c(108, 97, 3, 9), byrow = T, ncol = 2, nrow = 2)
    colnames(sex) <- c("NSVT", "No NSVT")
    rownames(sex) <- c("Male sex", "Female")
    sex
    chisq.test(sex)
    fisher.test(sex)
    
    # EF
    sex <- matrix(c(22, 183, 2, 10), byrow = T, ncol = 2, nrow = 2)
    colnames(sex) <- c("NSVT", "No NSVT")
    rownames(sex) <- c("Male sex", "Female")
    sex
    chisq.test(sex)
    fisher.test(sex)
    
    # MS
    sex <- matrix(c(99, 94, 6, 5), byrow = T, ncol = 2, nrow = 2)
    colnames(sex) <- c("NSVT", "No NSVT")
    rownames(sex) <- c("Male sex", "Female")
    sex
    chisq.test(sex)
    fisher.test(sex)
  
  # Time from dx - not done yet
  ansvt$Dxtoindex
  
  freq(nsvt$Dxtoindex)
  summary(nsvt$Dxtoindex) / 30.437
  summary(nonsvt$Dxtoindex) / 30.437
  summary(ansvt$Dxtoindex) / 30.437
  ks.test(x = nsvt$Dxtoindex, y = nonsvt$Dxtoindex)
  t.test(x = nsvt$Dxtoindex, y = nonsvt$Dxtoindex)
  
  # NYHA
  summary(nsvt$NYHA)
  summary(nonsvt$NYHA)
  summary(ansvt$NYHA)
  t.test(x = nsvt$NYHA, y = nonsvt$NYHA)
  
  # NYHA 3 or 4
  freq(nsvt$NYHA)
  freq(nonsvt$NYHA)
  freq(ansvt$NYHA)
  
    # NSVT
    nyha34 <- matrix(c(40, 33, 76, 68), byrow = T, ncol = 2, nrow = 2)
    colnames(nyha34) <- c("NSVT", "No NSVT")
    rownames(nyha34) <- c("3/4", "1/2")
    nyha34
    chisq.test(nyha34)
    fisher.test(nyha34)
    
    # VE
    nyha34 <- matrix(c(36, 37, 75, 69), byrow = T, ncol = 2, nrow = 2)
    colnames(nyha34) <- c("NSVT", "No NSVT")
    rownames(nyha34) <- c("3/4", "1/2")
    nyha34
    chisq.test(nyha34)
    fisher.test(nyha34)
    
    # ef
    nyha34 <- matrix(c(13, 60, 11, 133), byrow = T, ncol = 2, nrow = 2)
    colnames(nyha34) <- c("NSVT", "No NSVT")
    rownames(nyha34) <- c("3/4", "1/2")
    nyha34
    chisq.test(nyha34)
    fisher.test(nyha34)
    
    # ms
    nyha34 <- matrix(c(47, 58, 22, 77), byrow = T, ncol = 2, nrow = 2)
    colnames(nyha34) <- c("NSVT", "No NSVT")
    rownames(nyha34) <- c("3/4", "1/2")
    nyha34
    chisq.test(nyha34)
    fisher.test(nyha34)
  
  
  # Atrial fibrillation
  freq(nsvt$Atrial.fibrillation)
  freq(nonsvt$Atrial.fibrillation)
  freq(ansvt$Atrial.fibrillation)
  
    # NSVT
    afib <- matrix(c(84, 75, 32, 26), byrow = T, ncol = 2, nrow = 2)
    colnames(afib) <- c("NSVT", "No NSVT")
    rownames(afib) <- c("afib", "no afib")
    afib
    chisq.test(afib)
    fisher.test(afib)
    
    # VE
    afib <- matrix(c(82, 77, 29, 29), byrow = T, ncol = 2, nrow = 2)
    colnames(afib) <- c("NSVT", "No NSVT")
    rownames(afib) <- c("afib", "no afib")
    afib
    chisq.test(afib)
    fisher.test(afib)
    
    # EF
    afib <- matrix(c(21, 138, 3, 55), byrow = T, ncol = 2, nrow = 2)
    colnames(afib) <- c("NSVT", "No NSVT")
    rownames(afib) <- c("afib", "no afib")
    afib
    chisq.test(afib)
    fisher.test(afib)
    
    # ms
    afib <- matrix(c(73, 77, 32, 22), byrow = T, ncol = 2, nrow = 2)
    colnames(afib) <- c("NSVT", "No NSVT")
    rownames(afib) <- c("afib", "no afib")
    afib
    chisq.test(afib)
    fisher.test(afib)
  
  # ICD at baseline
  freq(nsvt$ICD.baseline)
  freq(nonsvt$ICD.baseline)
  freq(ansvt$ICD.baseline)
  
    # NSVT
    icd <- matrix(c(7, 0, 109, 101), byrow = T, ncol = 2, nrow = 2)
    colnames(icd) <- c("NSVT", "No NSVT")
    rownames(icd) <- c("ICD", "no ICD")
    icd
    chisq.test(icd)
    fisher.test(icd)
    
    # VE
    icd <- matrix(c(6, 1, 105, 105), byrow = T, ncol = 2, nrow = 2)
    colnames(icd) <- c("NSVT", "No NSVT")
    rownames(icd) <- c("ICD", "no ICD")
    icd
    chisq.test(icd)
    fisher.test(icd)
    
    # ef
    icd <- matrix(c(2, 5, 22, 188), byrow = T, ncol = 2, nrow = 2)
    colnames(icd) <- c("NSVT", "No NSVT")
    rownames(icd) <- c("ICD", "no ICD")
    icd
    chisq.test(icd)
    fisher.test(icd)
    
    # ms
    icd <- matrix(c(4, 3, 101, 96), byrow = T, ncol = 2, nrow = 2)
    colnames(icd) <- c("NSVT", "No NSVT")
    rownames(icd) <- c("ICD", "no ICD")
    icd
    chisq.test(icd)
    fisher.test(icd)
  
  # PPM baseline
  freq(nsvt$PPM.baseline)
  freq(nonsvt$PPM.baseline)
  freq(ansvt$PPM.baseline)
  
    # NSVT
    ppm <- matrix(c(12, 5, 104, 96), byrow = T, ncol = 2, nrow = 2)
    colnames(ppm) <- c("NSVT", "No NSVT")
    rownames(ppm) <- c("PPM", "no PPM")
    ppm
    chisq.test(ppm)
    fisher.test(ppm)
  
    # VE  
    ppm <- matrix(c(12, 5, 99, 101), byrow = T, ncol = 2, nrow = 2)
    colnames(ppm) <- c("NSVT", "No NSVT")
    rownames(ppm) <- c("PPM", "no PPM")
    ppm
    chisq.test(ppm)
    fisher.test(ppm)
    
    # EF
    ppm <- matrix(c(3, 14, 21, 179), byrow = T, ncol = 2, nrow = 2)
    colnames(ppm) <- c("NSVT", "No NSVT")
    rownames(ppm) <- c("PPM", "no PPM")
    ppm
    chisq.test(ppm)
    fisher.test(ppm)
    
    # ms
    ppm <- matrix(c(11, 5, 94, 94), byrow = T, ncol = 2, nrow = 2)
    colnames(ppm) <- c("NSVT", "No NSVT")
    rownames(ppm) <- c("PPM", "no PPM")
    ppm
    chisq.test(ppm)
    fisher.test(ppm)
  
  
  # Duration of monitoring
  summary(nsvt$Duration.of.monitoring..hours.)
  summary(nonsvt$Duration.of.monitoring..hours.)
  summary(ansvt$Duration.of.monitoring..hours.)
  t.test(nsvt$Duration.of.monitoring..hours., nonsvt$Duration.of.monitoring..hours.)
  
  freq(nsvt$Duration.of.monitoring..hours.)
  freq(nonsvt$Duration.of.monitoring..hours.)
  freq(ansvt$Duration.of.monitoring..hours.)
  
  # NSVT runs / 24h
  summary(nsvt$Runs.per.monitoring.period)
  summary(nonsvt$Runs.per.monitoring.period)
  summary(ansvt$Runs.per.monitoring.period)
  t.test(x = nsvt$Runs.per.monitoring.period, y = nonsvt$Runs.per.monitoring.period)
  
  # NSVT runs / 1h
  summary(nsvt$Runs.normalized)
  summary(nonsvt$Runs.normalized)
  
  # Longest run
  summary(nsvt$Longest)
  summary(nonsvt$Longest)
  summary(ansvt$Longest)
  freq(nsvt$Longest)
  freq(nonsvt$Longest)
  t.test(x = nsvt$Longest, y = nonsvt$Longest)
  
  # VE
  summary(nsvt$VE..)
  summary(nonsvt$VE..)
  ummary(ansvt$VE..)
  freq(nsvt$VE..)
  freq(nonsvt$VE..)
  ks.test(x = nsvt$VE..,y = nonsvt$VE..)
  t.test(nsvt$VE.., nonsvt$VE..)
  wilcox.test(x = nsvt$VE..,y = nonsvt$VE..)
  
  # EF
  summary(nsvt$EF)
  summary(nonsvt$EF)
  summary(ansvt$EF)
  ks.test(x = nsvt$EF,y = nonsvt$EF)
  t.test(x = nsvt$EF, y = nonsvt$EF)
  
    # NSVT
    sum(nsvt$EF <= 35)
    sum(nonsvt$EF <= 35)
    ef35 <- matrix(c(15, 9, 101, 92), byrow = T, ncol = 2, nrow = 2)
    colnames(ef35) <- c("NSVT", "No NSVT")
    rownames(ef35) <- c("eflow", "no eflow")
    ef35
    chisq.test(ef35)
    fisher.test(ef35)
    
    # VE
    ef35 <- matrix(c(19, 5, 92, 101), byrow = T, ncol = 2, nrow = 2)
    colnames(ef35) <- c("NSVT", "No NSVT")
    rownames(ef35) <- c("eflow", "no eflow")
    ef35
    chisq.test(ef35)
    fisher.test(ef35)
    
    # ms
    ef35 <- matrix(c(16, 7, 90, 92), byrow = T, ncol = 2, nrow = 2)
    colnames(ef35) <- c("NSVT", "No NSVT")
    rownames(ef35) <- c("eflow", "no eflow")
    ef35
    chisq.test(ef35)
    fisher.test(ef35)
  
  # Septal thickness
  summary(nsvt$IVST)
  summary(nonsvt$IVST)
  summary(ansvt$IVST)
  ks.test(x = nsvt$IVST, y = nonsvt$IVST)
  t.test(nsvt$IVST, nonsvt$IVST)

  # Posterior wall thickness
  summary(nsvt$PWT)
  summary(nonsvt$PWT)
  summary(ansvt$PWT)
  ks.test(x = nsvt$PWT,y = nonsvt$PWT)
  t.test(nsvt$PWT, nonsvt$PWT)
  
  # Left atrial volume index
  summary(nsvt$LAVI)
  summary(nonsvt$LAVI)
  summary(ansvt$LAVI)
  ks.test(x = nsvt$LAVI,y =  nonsvt$LAVI)
  t.test(nsvt$LAVI, nonsvt$LAVI)
  
  # Ventricular ectopy
  summary(nsvt$VE..)
  summary(nonsvt$VE..)
  summary(ansvt$VE..)
  t.test(nsvt$VE.., nonsvt$VE..)
  
  # Troponin collected
  freq(nsvt$Hs.cTnT)
  freq(nonsvt$Hs.cTnT)
  nsvt$troponin
  
    # NSVT
    tnt <- matrix(c(109, 96, 7, 5), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt) <- c("NSVT", "No NSVT")
    rownames(tnt) <- c("tnt", "no tnt")
    tnt
    chisq.test(tnt)
    
    # VE
    tnt <- matrix(c(104, 101, 7, 5), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt) <- c("NSVT", "No NSVT")
    rownames(tnt) <- c("tnt", "no tnt")
    tnt
    chisq.test(tnt)
    
    # EF
    tnt <- matrix(c(23, 183, 1, 10), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt) <- c("NSVT", "No NSVT")
    rownames(tnt) <- c("tnt", "no tnt")
    tnt
    chisq.test(tnt)
    
    # ms
    tnt <- matrix(c(105, 99, 0, 0), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt) <- c("NSVT", "No NSVT")
    rownames(tnt) <- c("tnt", "no tnt")
    tnt
    chisq.test(tnt)
    
    
  # Troponin value
  summary(nsvt$Hs.cTnT)
  summary(nonsvt$Hs.cTnT)
  summary(ansvt$Hs.cTnT)
  ks.test(x = nsvt$Hs.cTnT,y = nonsvt$Hs.cTnT)
  t.test(nsvt$Hs.cTnT, nonsvt$Hs.cTnT)

  
  # Troponin greater than 5
  a <- nsvt$Hs.cTnT[!is.na(nsvt$Hs.cTnT)]
  sum(a>0.05)
  b <- nonsvt$Hs.cTnT[!is.na(nonsvt$Hs.cTnT)]
  sum(b > 0.05)
  c <- ansvt$Hs.cTnT[!is.na(ansvt$Hs.cTnT)]
  sum(c > 0.05)
  
    # NSVT
    tnt15 <- matrix(c(49, 31, 60, 65), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt15) <- c("NSVT", "No NSVT")
    rownames(tnt15) <- c("tnt > 15", "tnt < 15")
    tnt15
    chisq.test(tnt15)
    fisher.test(tnt15)
    
    # VE
    tnt15 <- matrix(c(51, 29, 53, 72), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt15) <- c("NSVT", "No NSVT")
    rownames(tnt15) <- c("tnt > 15", "tnt < 15")
    tnt15
    chisq.test(tnt15)
    fisher.test(tnt15)
    
    # EF
    tnt15 <- matrix(c(11, 69, 12, 113), byrow = T, ncol = 2, nrow = 2)
    colnames(tnt15) <- c("NSVT", "No NSVT")
    rownames(tnt15) <- c("tnt > 15", "tnt < 15")
    tnt15
    chisq.test(tnt15)
    fisher.test(tnt15)
  
  
  # BNP collected
  freq(nsvt$NT.proBNP)
  freq(nonsvt$NT.proBNP)
  
    # NSVT
    bnp <- matrix(c(116, 100, 0, 1), byrow = T, ncol = 2, nrow = 2)
    colnames(bnp) <- c("NSVT", "No NSVT")
    rownames(bnp) <- c("bnp", "no bnp")
    bnp
    chisq.test(bnp)
    fisher.test(bnp)
    
    # VE
    bnp <- matrix(c(116, 105, 0, 1), byrow = T, ncol = 2, nrow = 2)
    colnames(bnp) <- c("NSVT", "No NSVT")
    rownames(bnp) <- c("bnp", "no bnp")
    bnp
    chisq.test(bnp)
    fisher.test(bnp)
    
    # EF
    bnp <- matrix(c(24, 192, 0, 1), byrow = T, ncol = 2, nrow = 2)
    colnames(bnp) <- c("NSVT", "No NSVT")
    rownames(bnp) <- c("bnp", "no bnp")
    bnp
    chisq.test(bnp)
    fisher.test(bnp)
    
  # BNP value
  freq(nsvt$NT.proBNP)
  freq(nonsvt$NT.proBNP)
  freq(ansvt$NT.proBNP)
  summary(nsvt$NT.proBNP)
  summary(nonsvt$NT.proBNP)
  summary(ansvt$NT.proBNP)
  ks.test(x = nsvt$NT.proBNP,y = nonsvt$NT.proBNP)
  t.test(nsvt$NT.proBNP, nonsvt$NT.proBNP)
    
  a <- nsvt$NT.proBNP[!is.na(nsvt$NT.proBNP)]
  sum(a>3000)
  b <- nonsvt$NT.proBNP[!is.na(nonsvt$NT.proBNP)]
  sum(b > 3000)
  c <- ansvt$NT.proBNP[!is.na(ansvt$NT.proBNP)]
  sum(c > 3000)
  
    # NSVT
    bnp3 <- matrix(c(36, 26, 80, 74), byrow = T, ncol = 2, nrow = 2)
    colnames(bnp3) <- c("NSVT", "No NSVT")
    rownames(bnp3) <- c("bnp > 3000", "bnp < 15")
    bnp3
    chisq.test(bnp3)
    fisher.test(bnp3)
    
    # VE
    bnp3 <- matrix(c(34, 28, 77, 77), byrow = T, ncol = 2, nrow = 2)
    colnames(bnp3) <- c("NSVT", "No NSVT")
    rownames(bnp3) <- c("bnp > 3000", "bnp < 15")
    bnp3
    chisq.test(bnp3)
    fisher.test(bnp3)
    
    # EF
    bnp3 <- matrix(c(12, 50, 12, 142), byrow = T, ncol = 2, nrow = 2)
    colnames(bnp3) <- c("NSVT", "No NSVT")
    rownames(bnp3) <- c("bnp > 3000", "bnp < 15")
    bnp3
    chisq.test(bnp3)
    fisher.test(bnp3)
  
  # Mayo stage
    freq(nsvt$Mayo.stage)
    freq(nonsvt$Mayo.stage)
    freq(ansvt$Mayo.stage)
  
  # Mayo staging available
    # NSVT
    ms <- matrix(c(109, 95, 7, 6), byrow = T, ncol = 2, nrow = 2)
    colnames(ms) <- c("NSVT", "No NSVT")
    rownames(ms) <- c("ms", "no ms")
    ms
    chisq.test(ms)
    fisher.test(ms)
    
    # VE
    ms <- matrix(c(104, 100, 7, 6), byrow = T, ncol = 2, nrow = 2)
    colnames(ms) <- c("NSVT", "No NSVT")
    rownames(ms) <- c("ms", "no ms")
    ms
    chisq.test(ms)
    fisher.test(ms)
    
    # EF
    ms <- matrix(c(23, 180, 1, 12), byrow = T, ncol = 2, nrow = 2)
    colnames(ms) <- c("NSVT", "No NSVT")
    rownames(ms) <- c("ms", "no ms")
    ms
    chisq.test(ms)
    fisher.test(ms)
    
    
  # Mayo stage discrete
    # NSVT
    ms <- matrix(c(47, 52, 37, 29, 25, 14), byrow = T, ncol = 2, nrow = 3)
    colnames(ms) <- c("NSVT", "No NSVT")
    rownames(ms) <- c("ms 1", "ms 2", "ms 3")
    ms
    chisq.test(ms)
    fisher.test(ms)
    ms[,1] / 109*100
    
    # VE
    ms <- matrix(c(43, 56, 38, 28, 23, 16), byrow = T, ncol = 2, nrow = 3)
    colnames(ms) <- c("NSVT", "No NSVT")
    rownames(ms) <- c("ms 1", "ms 2", "ms 3")
    ms
    chisq.test(ms)
    fisher.test(ms)
    ms[,1] / 109*100
    
    # EF
    ms <- matrix(c(7, 92, 8, 58, 8, 31), byrow = T, ncol = 2, nrow = 3)
    colnames(ms) <- c("NSVT", "No NSVT")
    rownames(ms) <- c("ms 1", "ms 2", "ms 3")
    ms
    chisq.test(ms)
    fisher.test(ms)
    ms[,1] / 109*100
    
    # Mayo stage continuous - don't use
    freq(nsvt$Mayo.stage)
    freq(nonsvt$Mayo.stage)
    summary(nsvt$Mayo.stage)
    summary(nonsvt$Mayo.stage)
    summary(ansvt$Mayo.stage)
    ks.test(x = nsvt$Mayo.stage,y = nonsvt$Mayo.stage)
    t.test(x = nsvt$Mayo.stage,y = nonsvt$Mayo.stage)
  
  # Beta blocker
  freq(nsvt$BB)
  freq(nonsvt$BB)
  freq(ansvt$BB)
  
    # NSVT
    bb <- matrix(c(75, 51, 41, 50), byrow = T, ncol = 2, nrow = 2)
    colnames(bb) <- c("NSVT", "No NSVT")
    rownames(bb) <- c("bb", "no bb")
    bb
    chisq.test(bb)
    fisher.test(bb)
    
    # VE
    bb <- matrix(c(67, 59, 44, 47), byrow = T, ncol = 2, nrow = 2)
    colnames(bb) <- c("NSVT", "No NSVT")
    rownames(bb) <- c("bb", "no bb")
    bb
    chisq.test(bb)
    fisher.test(bb)
    
    # EF
    bb <- matrix(c(16, 110, 8, 83), byrow = T, ncol = 2, nrow = 2)
    colnames(bb) <- c("NSVT", "No NSVT")
    rownames(bb) <- c("bb", "no bb")
    bb
    chisq.test(bb)
    fisher.test(bb)
    
    # MS
    bb <- matrix(c(59, 59, 46, 40), byrow = T, ncol = 2, nrow = 2)
    colnames(bb) <- c("NSVT", "No NSVT")
    rownames(bb) <- c("bb", "no bb")
    bb
    chisq.test(bb)
    fisher.test(bb)
  
  # Calcium channel blocker
  freq(nsvt$CCB)
  freq(nonsvt$CCB)
  freq(ansvt$CCB)
  
    # NSVT
    ccb <- matrix(c(24, 23, 92, 78), byrow = T, ncol = 2, nrow = 2)
    colnames(ccb) <- c("NSVT", "No NSVT")
    rownames(ccb) <- c("ccb", "no ccb")
    ccb
    chisq.test(ccb)
    
    # VE
    ccb <- matrix(c(18, 29, 93, 77), byrow = T, ncol = 2, nrow = 2)
    colnames(ccb) <- c("NSVT", "No NSVT")
    rownames(ccb) <- c("ccb", "no ccb")
    ccb
    chisq.test(ccb)
    
    # EF
    ccb <- matrix(c(3, 44, 21, 149), byrow = T, ncol = 2, nrow = 2)
    colnames(ccb) <- c("NSVT", "No NSVT")
    rownames(ccb) <- c("ccb", "no ccb")
    ccb
    chisq.test(ccb)
    
    # ms
    ccb <- matrix(c(20, 26, 85, 73), byrow = T, ncol = 2, nrow = 2)
    colnames(ccb) <- c("NSVT", "No NSVT")
    rownames(ccb) <- c("ccb", "no ccb")
    ccb
    chisq.test(ccb)
    
  # Antiarrhythmic
    freq(ansvt$Antiarrhythmic....name)
    freq(nsvt$Antiarrhythmic....name)
    freq(nonsvt$Antiarrhythmic....name)
    
    aa <- matrix(c(14, 26, 102, 75), byrow = T, ncol = 2, nrow = 2)
    colnames(aa) <- c("NSVT", "No NSVT")
    rownames(aa) <- c("aa", "no aa")
    aa
    chisq.test(aa)
    
    # VE
    aa <- matrix(c(17, 23, 94, 83 ), byrow = T, ncol = 2, nrow = 2)
    colnames(aa) <- c("VE", "No VE")
    rownames(aa) <- c("aa", "no aa")
    aa
    chisq.test(aa)
    
    # EF
    aa <- matrix(c(4, 36, 20, 257 ), byrow = T, ncol = 2, nrow = 2)
    colnames(aa) <- c("VE", "No VE")
    rownames(aa) <- c("aa", "no aa")
    aa
    chisq.test(aa)
  
    # ms
    aa <- matrix(c(14, 23, 91, 76 ), byrow = T, ncol = 2, nrow = 2)
    colnames(aa) <- c("VE", "No VE")
    rownames(aa) <- c("aa", "no aa")
    aa
    chisq.test(aa)  
      
  # Tafamadis
  freq(nsvt$Tafamidis.baseline)
  freq(nonsvt$Tafamidis.baseline)
  freq(ansvt$Tafamidis.baseline)

    # NSVT
    taf <- matrix(c(33, 23, 83, 78), byrow = T, ncol = 2, nrow = 2)
    colnames(taf) <- c("NSVT", "No NSVT")
    rownames(taf) <- c("taf", "no taf")
    taf
    chisq.test(taf)
    
    # VE
    taf <- matrix(c(28, 23, 77, 76), byrow = T, ncol = 2, nrow = 2)
    colnames(taf) <- c("NSVT", "No NSVT")
    rownames(taf) <- c("taf", "no taf")
    taf
    chisq.test(taf)
    
    # EF
    taf <- matrix(c(6, 50, 18, 143), byrow = T, ncol = 2, nrow = 2)
    colnames(taf) <- c("NSVT", "No NSVT")
    rownames(taf) <- c("taf", "no taf")
    taf
    chisq.test(taf)
    
  # CMR obtained
    freq(ansvt$Cardiac.MR.obtained.1...yes..0...no)
    freq(nsvt$Cardiac.MR.obtained.1...yes..0...no)
    freq(nonsvt$Cardiac.MR.obtained.1...yes..0...no)
    
    cmr <- matrix(c(10, 83, 14, 110), byrow = T, ncol = 2, nrow = 2)
    colnames(cmr) <- c("NSVT", "No NSVT")
    rownames(cmr) <- c("cmr", "no cmr")
    cmr
    chisq.test(cmr)
    
    # VE
    cmr <- matrix(c(45, 48, 66, 58), byrow = T, ncol = 2, nrow = 2)
    colnames(cmr) <- c("NSVT", "No NSVT")
    rownames(cmr) <- c("cmr", "no cmr")
    cmr
    chisq.test(cmr)
    
    # VE
    cmr <- matrix(c(41, 43, 64, 56), byrow = T, ncol = 2, nrow = 2)
    colnames(cmr) <- c("NSVT", "No NSVT")
    rownames(cmr) <- c("cmr", "no cmr")
    cmr
    chisq.test(cmr)
    
  
  # LGE present
    freq(ansvt$LGE.present..1...yes..0...no.)
    freq(nsvt$LGE.present..1...yes..0...no.)
    freq(nonsvt$LGE.present..1...yes..0...no.)
    
    
    
  # LGE qual
    showmatrix <- lge[,c(86, 88, 89, 90, 91, 131, 133)]
    View(showmatrix)
    
    lge <- subset(ansvt,subset = ansvt$Cardiac.MR.obtained.1...yes..0...no == 1)
    subendo <- ifelse(lge$Subendo.LGE == 1 & lge$Midwall.LGE == 0 & lge$Subepi.LGE == 0,1,0)
    sum(subendo)
    transmural <- ifelse(lge$Subendo.LGE == 1 & lge$Subepi.LGE == 1,1,0)
    sum(transmural)
    freq(lge$LGE2nd)
    
    # NSVT
    lgensvt <- subset(nsvt,subset = nsvt$Cardiac.MR.obtained.1...yes..0...no == 1)
    subendo <- ifelse(lgensvt$Subendo.LGE == 1 & lgensvt$Midwall.LGE == 0 & lgensvt$Subepi.LGE == 0,1,0)
    sum(subendo)
    subendo
    transmural <- ifelse(lgensvt$Subendo.LGE == 1 & lgensvt$Subepi.LGE == 1,1,0)
    sum(transmural)
    transmural
    freq(lgensvt$LGE2nd)
    
    # No NSVT
    lgenonsvt <- subset(nonsvt,subset = nonsvt$Cardiac.MR.obtained.1...yes..0...no == 1)
    subendo <- ifelse(lgenonsvt$Subendo.LGE == 1 & lgenonsvt$Midwall.LGE == 0 & lgenonsvt$Subepi.LGE == 0,1,0)
    sum(subendo)
    subendo
    transmural <- ifelse(lgenonsvt$Subendo.LGE == 1 & lgenonsvt$Subepi.LGE == 1,1,0)
    sum(transmural)
    transmural
    freq(lgenonsvt$LGE2nd)
    
    # Stats
      lgepat <- matrix(c(16, 13, 16, 12, 17, 19), byrow = T, ncol = 2, nrow = 3)
      colnames(lgepat) <- c("NSVT", "No NSVT")
      rownames(lgepat) <- c("subendo", "transmural", "other")
      lgepat
      chisq.test(lgepat)
    
      # VE
      lgepat <- matrix(c(17, 12, 15, 13, 13, 23), byrow = T, ncol = 2, nrow = 3)
      colnames(lgepat) <- c("NSVT", "No NSVT")
      rownames(lgepat) <- c("subendo", "transmural", "other")
      lgepat
      chisq.test(lgepat)
    
      # EF
      lgepat <- matrix(c(2, 27, 4, 24, 4, 32), byrow = T, ncol = 2, nrow = 3)
      colnames(lgepat) <- c("NSVT", "No NSVT")
      rownames(lgepat) <- c("subendo", "transmural", "other")
      lgepat
      chisq.test(lgepat)
    
      # Stage
      lgepat <- matrix(c(17, 9, 13, 12, 11, 22), byrow = T, ncol = 2, nrow = 3)
      colnames(lgepat) <- c("NSVT", "No NSVT")
      rownames(lgepat) <- c("subendo", "transmural", "other")
      lgepat
      chisq.test(lgepat)
      
    # LGE extensive
      freq(ansvt$LGE2nd)
      freq(nsvt$LGE2nd)
      freq(nonsvt$LGE2nd)
      
      # By NSVT
      lgep <- matrix(c(1, 2, 12, 7, 36, 35), byrow = T, ncol = 2, nrow = 3)
      colnames(lgep) <- c("NSVT", "No NSVT")
      rownames(lgep) <- c("None", "Non-extensive", "extensive")
      lgep
      chisq.test(lgep)
      
      # VE
      lgep <- matrix(c(0, 3, 11, 8, 34, 37), byrow = T, ncol = 2, nrow = 3)
      colnames(lgep) <- c("NSVT", "No NSVT")
      rownames(lgep) <- c("None", "Non-extensive", "extensive")
      lgep
      chisq.test(lgep)
      
      # EF
      lgep <- matrix(c(2, 1, 1, 18, 7, 64), byrow = T, ncol = 2, nrow = 3)
      colnames(lgep) <- c("NSVT", "No NSVT")
      rownames(lgep) <- c("None", "Non-extensive", "extensive")
      lgep
      chisq.test(lgep)
      
      # MS
      lgep <- matrix(c(1, 2, 6, 11, 34, 30), byrow = T, ncol = 2, nrow = 3)
      colnames(lgep) <- c("NSVT", "No NSVT")
      rownames(lgep) <- c("None", "Non-extensive", "extensive")
      lgep
      chisq.test(lgep)
      
      
  # Sustained VT
  freq(nsvt$VT.or.VF.WS)
  freq(nonsvt$VT.or.VF.WS)
  freq(ansvt$VT.or.VF.WS)
  
  # Sustained VT
    svt <- matrix(c(10, 1, 106, 100), byrow = T, ncol = 2, nrow = 2)
    colnames(svt) <- c("NSVT", "No NSVT")
    rownames(svt) <- c("svt", "no svt")
    svt
    chisq.test(svt)
    fisher.test(svt)

  # VF
    vf <- matrix(c(0, 0, 116, 101), byrow = T, ncol = 2, nrow = 2)
    colnames(vf) <- c("NSVT", "No NSVT")
    rownames(vf) <- c("vf", "no vf")
    vf
    chisq.test(vf)
    fisher.test(vf)
  
  # Mortality
  freq(nsvt$Mortality)
  freq(nonsvt$Mortality)
  freq(ansvt$Mortality)  
  died <- matrix(c(2, 0, 114, 101), byrow = T, ncol = 2, nrow = 2)
  colnames(died) <- c("NSVT", "No NSVT")
  rownames(died) <- c("died", "alive")
  died
  chisq.test(died)
  fisher.test(died)

  # Total followup
  summary(ansvt$Fig.3b) / 30.347
  
  negativedxtoindex <- subset(ansvt,subset = ansvt$Dxtoindex < 0)
  View(negativedxtoindex)
  dim(negativedxtoindex)
  summary(negativedxtoindex$Dxtoindex)
  
  nsvt <- subset(x = ansvt, subset = ansvt$NSVT == 1)
  nonsvt <- subset(x = ansvt, subset = ansvt$NSVT == 0)
  dim(nsvt)
  dim(nonsvt)
  
  summary(nsvt$Fig.3b) / 30.437
  summary(nonsvt$Fig.3b) / 30.437
  summary(ansvt$Fig.3b) / 30.437
  ks.test(nsvt$Fig.3b, nonsvt$Fig.3b)
  t.test(nsvt$Fig.3b, nonsvt$Fig.3b)
  
  
  
  
  nsvt <- subset(x = ansvt, subset = ansvt$NSVT == 1)
  nonsvt <- subset(x = ansvt, subset = ansvt$NSVT == 0)
  
  
  # ICDs placed
  ansvt$ICD.FU
  icdput <- subset(x = ansvt, subset = ansvt$ICD..all. == 1 & is.na(ansvt$ICD.baseline))
  dim(icdput)
    nsvticdput <- subset(x = nsvt, subset = nsvt$ICD..all. == 1 & is.na(nsvt$ICD.baseline))
    dim(nsvticdput)
    nonsvticdput <- subset(x = nonsvt, subset = nonsvt$ICD..all. == 1 & is.na(nonsvt$ICD.baseline))
    dim(nonsvticdput)
  icdp <- matrix(c(11,4, 206,97), byrow = T, ncol = 2, nrow = 2)
  colnames(icdp) <- c("NSVT", "No NSVT")
  rownames(icdp) <- c("icdp", "no icdp")
  icdp
  chisq.test(icdp)
  fisher.test(icdp)
  fisher.test(icdp)

  # PPM placed
  ppmput <- subset(x=ansvt, subset = ansvt$pacemaker..all. 
                   != 0 & ansvt$PPM.baseline == 0 & ansvt$ICD.FU == 0)
  dim(ppmput)
  freq(ppmput$pacemaker..all.)
  nsvtppmput <- subset(x=nsvt, subset = nsvt$pacemaker..all. 
                       != 0 & nsvt$PPM.baseline == 0 & nsvt$ICD.FU == 0)
  dim(nsvtppmput)
  nonsvtppmput <- subset(x=nonsvt, subset = nonsvt$pacemaker..all. 
                         != 0 & nonsvt$PPM.baseline == 0 & nonsvt$ICD.FU == 0)
  dim(nonsvtppmput)
  ppmp <- matrix(c(11, 4, 206,97), byrow = T, ncol = 2, nrow = 2)
  colnames(ppmp) <- c("NSVT", "No NSVT")
  rownames(ppmp) <- c("ppmp", "no ppmp")
  ppmp
  chisq.test(ppmp)
  fisher.test(ppmp)

  
  
  # loop placed
  loop <- subset(x = ansvt, subset = ansvt$Loop.recorder == 1)
  dim(loop)    
    loopn <- subset(x = nsvt, subset = nsvt$Loop.recorder == 1)
    dim(loopn)
    loopno <- subset(x = nonsvt, subset = nonsvt$Loop.recorder == 1)
    dim(loopno)
    loopp <- matrix(c(8 ,7 , 108, 94), byrow = T, ncol = 2, nrow = 2)
    colnames(loopp) <- c("NSVT", "No NSVT")
    rownames(loopp) <- c("loopp", "no loopp")
    loopp
    chisq.test(loopp)
    fisher.test(loopp)
  
    
    
    # Total number of holter monitors
    totalholter <- rep(1,217)
    totalholter
    
    for (i in 1:217) {
      totalholter[i] <- totalholter[i] + sum(!is.na(ansvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholter
    summary(totalholter)
    
    sum(totalholter >1)
    
    totalholternsvt <- rep(1,116)
    totalholternsvt
    
    for (i in 1:116) {
      totalholternsvt[i] <- totalholternsvt[i] + sum(!is.na(nsvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholternsvt
    summary(totalholternsvt)
    
    
    View(nonsvt)
    dim(nonsvt)
    totalholternonsvt <- rep(1,101)
    totalholternonsvt
    
    for (i in 1:101) {
      totalholternonsvt[i] <- totalholternonsvt[i] + sum(!is.na(nonsvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholternonsvt
    summary(totalholternonsvt)
    
    t.test(totalholternsvt, totalholternonsvt)
    
    # VE
    totalholternsvt <- rep(1,111)
    totalholternsvt
    
    for (i in 1:111) {
      totalholternsvt[i] <- totalholternsvt[i] + sum(!is.na(nsvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholternsvt
    summary(totalholternsvt)
    totalholternonsvt <- rep(1,106)
    totalholternonsvt
    
    for (i in 1:106) {
      totalholternonsvt[i] <- totalholternonsvt[i] + sum(!is.na(nonsvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholternonsvt
    summary(totalholternonsvt)
    t.test(totalholternsvt, totalholternonsvt)
    
    # ef
    totalholternsvt <- rep(1,24)
    totalholternsvt
    
    for (i in 1:24) {
      totalholternsvt[i] <- totalholternsvt[i] + sum(!is.na(nsvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholternsvt
    summary(totalholternsvt)
    totalholternonsvt <- rep(1,193)
    totalholternonsvt
    
    for (i in 1:193) {
      totalholternonsvt[i] <- totalholternonsvt[i] + sum(!is.na(nonsvt[i, c(95, 101, 107, 113, 119, 125)]))
    }
    totalholternonsvt
    summary(totalholternonsvt)
    t.test(totalholternsvt, totalholternonsvt)
    

    
newms <- rep(1,217)    
trop <- rep(0,217)
bnp <- rep(0,217)
for (i in 1:217){
  if (ansvt[i,34] > 0.05){
    trop[i] <- 1
  } else {
    bnp[i] <- 0
  }
  if (ansvt[i,35] > 3000){
    bnp[i] <- 1
  }
  bnp[i] <- bnp[i] + trop[i] + bnp[i]
}    
    
    
#####
# Figure 1A
  
  runsper24h <- nsvt$Runs.per.monitoring.period
  runsper24h
  summary(runsper24h)
  
  runsper1h <- nsvt$Runs.normalized
  summary(runsper1h)
  freq(runsper1h)

  # Working plot
  ggplot(data = nsvt, aes(x = Runs.normalized)) +
    geom_histogram() +
    theme_bw() +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25),
      limits = c(0, 1.25)
    ) +
    scale_y_continuous(
      breaks = c(0, 10, 20, 30, 40),
      limits = c(0, 40)) +
    labs(y = "Frequency", x = "NSVT runs per 24h") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18)
    )
  
  # Final Figure 1A
  a <- c(0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 
         0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.80, 0.84, 0.88, 0.92, 0.96, 1.00,
         1.04, 1.08, 1.12, 1.16, 1.20)
  b <- c(4, 39, 20, 12, 9, 2, 4, 3, 2, 0, 1, 0 , 0, 1, 0, 0, 1, 0, 1, 2, 1,
         0, 1, 0, 2, 0, 1, 0, 0, 1)
  y <- b / 116 * 100
  plot1a <- matrix(data = c(a, b, y), nrow = 30)
  colnames(plot1a) <- c("x", "Y", "Percent")
  plot1a <- as.data.frame(plot1a)
  
  p1a <- ggplot(data = plot1a, aes(x = plot1a[,1], y = plot1a[,3])) +
    geom_col(fill = "steelblue") +
    theme_bw() +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
                       limits = c(0, 1.2)) +
    labs(y = "Percentage of population", x = "NSVT runs per hour") +
    scale_y_continuous(breaks = c(0, 10, 20, 30),
                       limits = c(0, 35)) +
    theme(axis.text =  element_text(size = 14),
          axis.title = element_text(size = 18),
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"))
p1a

  
# Figure 1b
nsvt$Longest
freq(nsvt$Longest)
a <- c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
       16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
b <- c(35, 24, 14, 8, 6, 2, 3, 5, 3, 1, 2, 3, 1, 3, 2, 0, 1, 0, 1, 1, 0, 1, 0)
length(a)
length(b)
yb <- b / 116 * 100
plot1b <- matrix(data = c(a, b, yb), nrow = 23)
colnames(plot1b) <- c("x", "Y", "Percent")
plot1b <- as.data.frame(plot1b)
p1b <-  ggplot(data = plot1b, aes(x = plot1b[,1],y = plot1b[,3])) +
  geom_col(fill = "steelblue") +
  theme_bw() +
  scale_x_continuous(
    breaks = c(3, 8, 16, 24),
    limits = c(2, 24)
  ) +
  scale_y_continuous(
    breaks = c(0, 10, 20, 30),
    limits = c(0, 35)) +
  labs(y = "Percentage of population", x = "Longest NSVT run (beats)") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")
  )
p1b


# Figure 1c
freq(ansvt$VE..)

c <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
d <- c(143, 22, 11, 10, 4, 3, 5, 1, 1, 0, 1)
yd <- d / 217 * 100

plot1c <- matrix(data = c(c, d, yd), nrow = 11)
colnames(plot1c) <- c("x", "Y", "Percent")
plot1c <- as.data.frame(plot1c)
p1c <- ggplot(data = plot1c, aes(x = plot1c[,1], y = plot1c[,3])) +
  geom_col(fill = "steelblue",just = 1) +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
                     limits = c(0, 11)) +
  labs(y = "Percentage of population", x = "% Ventricular Ectopy") +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70),
                     limits = c(0, 70)) +
  theme(axis.text =  element_text(size = 14),
        axis.title = element_text(size = 18),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

p1c

x11()
ggarrange(p1a, p1b, p1c, nrow = 1,labels = "AUTO")


# Univariate analysis
  # Age
  agemod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Age)
  summary(agemod)
  cox.zph(agemod)
  
  # Time since diagnosis - not using
  ttd <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Dxtoindex)
  summary(ttd)
  cox.zph(ttd)

  # AF
  af <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Atrial.fibrillation)
  summary(af)
  cox.zph(af)

  # NYHA
  nyha <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ as.factor(ansvt$NYHA))
  summary(nyha)
  cox.zph(nyha)
  
  freq(ansvt$NYHA)
  nyha34 <- ansvt$NYHA
  save <- nyha34
  for (i in 1:217) {
    if (is.na(nyha34[i])) {
      nyha34[i] <- 1
    }
  }
  for (i in 1:217) {
    if (nyha34[i] >= 3) {
      nyha34[i] <- 1
    } else {
      nyha34[i] <- 0
    }
  }
  freq(nyha34)
  
  nyha <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ nyha34)
  summary(nyha)
  
  # Baseline ICD - don't use - cannot use
  icd <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$ICD.baseline)
  summary(icd)
  
  # For all pacemakers - don't use - significant but not clinically relevant
  ppm <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ppmbaseline)
  summary(ppm)
  
  # NSVT
  presnsvt <- nsvtrunsperhour <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$NSVT)
  summary(presnsvt)
  cox.zph(presnsvt)
  
  # Runs per monitoring period
  run24h <- nsvtrunsperhour <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Runs.per.monitoring.period)
  summary(run24h)
  cox.zph(run24h)
  run1h <- nsvtrunsperhour <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Runs.normalized)
  summary(run1h)
  
  # Longest run
  longestrun <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Longest)
  summary(longestrun)
  cox.zph(longestrun)
  
  # Ventricular ectopy
  ve <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE..)
  summary(ve)
  cox.zph(ve)
  
  ve0.5 <- ansvt$VE..
  outputve <- rep(0,217)
  for (i in 1:217) {
    if (ve0.5[i] >= 0.5) {
      outputve[i] <- 1
    } else {
      outputve[i] <- 0
    }
  }
  ve0.5mod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ outputve)
  summary(ve0.5mod)
  cox.zph(ve0.5mod)
  
  # EF
  # EF continuous
  nsvtef <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$EF)
  summary(nsvtef)
  cox.zph(nsvtef)
  
  # EF < 35%
  ef35 <- ansvt$EF
  output <- rep(0,217)
  freq(ansvt$EF) # Zero NAs

  for (i in 1:217) {
    if (ef35[i] < 35) {
      output[i] <- 1
    } else {
      output[i] <- 0
    }
  }
  ef35
  output
  
  nsvtef35 <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ output)
  summary(nsvtef35)
  
  # septal wall thickness
  freq(ansvt$IVST)
  nsvtseptum <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$IVST)
  summary(nsvtseptum)
  cox.zph(nsvtseptum)
  
  # Posterior wall thickness
  nsvtpwt <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$PWT)
  summary(nsvtpwt)
  cox.zph(nsvtpwt)
  
  # LAVI
  nsvtlavi <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$LAVI)
  summary(nsvtlavi)
  cox.zph(nsvtlavi)
  
  # Mayo stage
  ms <- ansvt$Mayo.stage
  for (i in 1:217) {
    if (is.na(ms[i])) {
      ms[i] <- 1
    }
  }
  ms
  
  mayostage  <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Mayo.stage)
  summary(mayostage)
  cox.zph(nsvtstage)
  
  # HS-TNT
  nsvthst <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Hs.cTnT)
  summary(nsvthst)
  
  hstnt <- ansvt$Hs.cTnT
  output <- rep(0,217)
  freq(ansvt$Hs.cTnT) # 1 NA
  # 12 NAs
  for (i in 1:217) {
    if (is.na(hstnt[i])) {
      hstnt[i] <- 0
    }
  }
  hstnt
  freq(hstnt)
  for (i in 1:217) {
    if (hstnt[i] > 0.05) {
      output[i] <- 1
    } else {
      output[i] <- 0
    }
  }
  freq(output)
  
  hscutoff <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ output)
  summary(hscutoff)
  
  
  # BNP
  bnp <- ansvt$NT.proBNP
  output <- rep(0,217)
  freq(ansvt$NT.proBNP) # 1 NA
  
  for (i in 1:217) {
    if (is.na(bnp[i])) {
      bnp[i] <- 0
    }
  }
  bnp
  for (i in 1:217) {
    if (bnp[i] > 3000) {
      output[i] <- 1
    } else {
      output[i] <- 0
    }
  }
  
  bnp3000 <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ output)
  summary(bnp3000)
  
  # Beta blocker
  nsvtbb <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$BB)
  summary(nsvtbb)
  cox.zph(nsvtbb)
  
  # CCB
  nsvtccb <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$CCB)
  summary(nsvtccb)
  cox.zph(nsvtccb)
  
  # Tafamadis at baseline - has alot of NAs
  tafbaseline <- ansvt$Tafamidis.baseline 
  for(i in 1:217) {
    if (is.na(tafbaseline[i])) {
      tafbaseline[i] <- 0
    } 
  }
  tafbaseline
  nsvttafamadis <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ tafbaseline)
  summary(nsvttafamadis)
  
  # Antiarrhythmic
  aa <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$AA...holter.1)
  summary(aa)
  cox.zph(aa)
  
# Multivariate model
# Interesting that total runs is significant, but When you use total runs
  # It makes NSVT better than VE. But should be using a normalized runs
  # Per hour which we have done
    
  # Multivariate
  finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT)
  summary(finmod) # VE and NSVT alone are both sig
  
  
  finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + ansvt$Age + nyha34 + ansvt$EF)
  summary(finmod) # VE and NSVT alone are both sig
  
  # Tentative one to use for 11 events
  finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + nyha34 + ansvt$Longest)
  summary(finmod)
  
  
  
  
  # Models to place just in text and not in the column
  # For NSVT
  finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Age + 
                    ansvt$NSVT + nyha34)
  summary(finmod) 
  
  # For VE and not NSVT
  finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$Age + 
                    ansvt$VE.. + nyha34)
  summary(finmod) 
  
  
  
  # VE and sustained VA
  finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE..)
  summary(finmod) 
  
  # VE and mortality
  finmod <- coxph(Surv(time = ansvt[,75], event = ansvt$Mortality) ~ ansvt$VE..)
  summary(finmod) 
  
  
  
  # The rest we don't use
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + ansvt$Longest)
    summary(finmod) # add longest and just VE
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + ansvt$Longest + nyha34)
    summary(finmod) # add longest and nyha3/4 and now none are significant
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT)
    summary(finmod)
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + ansvt$Runs.normalized)
    summary(finmod)
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + nyha34 + ansvt$PPM.baseline)
    summary(finmod)
    cox.zph(finmod)  
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT + ansvt$Longest)
    summary(finmod)
    cox.zph(finmod)  
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                    ansvt$NSVT)
    summary(finmod)
    
    finmod <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ ansvt$VE.. + 
                       ansvt$EF + ansvt$Age + ansvt$NSVT)
    summary(finmod)
  
  ansvt$VT.or.VF.WS[]  

#####
# Figures
# Figure 2
nsvt <- ansvt$NSVT
km.model1 <- survfit(Surv(time = ansvt[,72] , event = ansvt[,63]) ~nsvt, data = ansvt)
km.model1
summary(km.model1)
p2 <- ggsurvplot(km.model1,
                 surv.scale = "percent",
                 risk.table = TRUE,
                 pval=TRUE,
                 pval.method = TRUE,
                 xlab = c("Months"),
                 ylab = c("Survival free from VT/VF"),
                 font.legend.labs = c(18),
                 legend.labs = c("No NSVT", "NSVT"),
                 legend.title = c(""),
                 xscale = c("d_m"),
                 xlim = c(0,70*30.437),
                 break.x.by = 10* 30.437,
                 risk.table.title = c("At Risk"),
                 risk.table.fontsize = 4,
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE,
                 #tables.height = 0.15,
                 censor = TRUE,
                 font.x = 16,
                 font.y = 16,
                 font.tickslab = 16,
                 palette = c("steelblue", "red"),
)
p2
x11()
p2

# Old Figure 3A - Mortality by NSVT figure
ansvt[,75]
ansvt[,75] / 30.437


km.model3a <- survfit(Surv(time = ansvt[,75] , event = ansvt$Mortality) ~ nsvt, data = ansvt)
km.model3a
summary(km.model3a)
p3a <- ggsurvplot(km.model3a,
                 surv.scale = "percent",
                 risk.table = TRUE,
                 pval=TRUE,
                 pval.method = TRUE,
                 font.legend.labs = c(18),
                 xlab = c("Months"),
                 ylab = c("Survival free from all-cause mortality"),
                 legend.labs = c("No NSVT", "NSVT"),
                 legend.title = c(""),
                 xscale = c("d_m"),
                 xlim = c(0,70 * 30.437),
                 break.x.by = 10 * 30.437,
                 risk.table.title = c("At Risk"),
                 risk.table.fontsize = 4,
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE,
                 #tables.height = 0.15,
                 censor = TRUE,
                 font.x = 16,
                 font.y = 16,
                 font.tickslab = 16,
                 palette = c("steelblue", "red"),
)
p3a

# Composite by NSVT figure
km.model3b <- survfit(Surv(time = ansvt[,77] , event = ansvt[,76]) ~ nsvt, data = ansvt)
km.model3b
summary(km.model3b)
p3b <- ggsurvplot(km.model3b,
                 surv.scale = "percent",
                 risk.table = TRUE,
                 pval=TRUE,
                 pval.method = TRUE,
                 xlab = c("Months"),
                 ylab = c("Survival free from VT/VF or mortality"),
                 legend.labs = c("No NSVT", "NSVT"),
                 font.legend.labs = c(18),
                 legend.title = c(""),
                 font.title = c(20, "bold", "black"),
                 xscale = c("d_m"),
                 xlim = c(0,70 * 30.437),
                 break.x.by = 10 * 30.437,
                 risk.table.title = c("At Risk"),
                 risk.table.fontsize = 4,
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE,
                 tables.height = 0.15,
                 censor = TRUE,
                 font.x = 16,
                 font.y = 16,
                 font.tickslab = 16,
                 palette = c("steelblue", "red"),
)
p3b

survplotlist <- list()
survplotlist[[1]] <- p3a + labs(tag = 'A')
survplotlist[[2]] <- p3b + labs(tag = 'B')

survplotlist
arrange_ggsurvplots(survplotlist, ncol = 2, nrow = 1, labels = "AUTO", 
                    risk.table.height = 0.15)


# Playing around with VE and kaplan meier curves

  # VE cutoff of 0.5
  ve <- ansvt$VE..
  summary(ve)
  
  freq(ansvt$VE..)
  ve <- ansvt$VE..
  output <- rep(0,217)
  for (i in 1:217) {
    if (ve[i] >= 0.5) {
      output[i] <- 1
    } else {
      output[i] <- 0
    }
  }
  freq(output)
  
  km.model4 <- survfit(Surv(time = ansvt[,72] , event = ansvt[,63]) ~output, data = ansvt)
  km.model4
  summary(km.model4)
  
  p4 <- ggsurvplot(km.model4,
                   surv.scale = "percent",
                   risk.table = TRUE,
                   pval=TRUE,
                   pval.method = TRUE,
                   xlab = c("Months"),
                   ylab = c("Survival free from VT/VF"),
                   font.legend.labs = c(14), 
                   legend.labs = c("VE < 0.5%", "VE ≥ 0.5%"),
                   legend.title = c(""),
                   xscale = c("d_m"),
                   xlim = c(0,70*30.437),
                   break.x.by = 10* 30.437,
                   risk.table.title = c("At Risk"),
                   risk.table.fontsize = 4,
                   risk.table.y.text = FALSE,
                   risk.table.y.text.col = TRUE,
                   tables.y.text = FALSE,
                   #tables.height = 0.15,
                   censor = TRUE,
                   font.x = 16,
                   font.y = 16,
                   font.tickslab = 16,
                   palette = c("steelblue", "red"),
  )
  x11()
  p4

  km.model5 <- survfit(Surv(time = ansvt[,72] , event = ansvt[,63]) ~nyha34, data = ansvt)
  km.model5
  summary(km.model5)
  
  p5 <- ggsurvplot(km.model5,
                   surv.scale = "percent",
                   risk.table = TRUE,
                   pval=TRUE,
                   pval.method = TRUE,
                   xlab = c("Months"),
                   ylab = c("Survival free from VT/VF"),
                   font.legend.labs = c(14), 
                   legend.labs = c("No", "NYHA 3/4"),
                   legend.title = c(""),
                   xscale = c("d_m"),
                   xlim = c(0,70*30.437),
                   break.x.by = 10* 30.437,
                   risk.table.title = c("At Risk"),
                   risk.table.fontsize = 4,
                   risk.table.y.text = FALSE,
                   risk.table.y.text.col = TRUE,
                   tables.y.text = FALSE,
                   tables.height = 0.15,
                   censor = TRUE,
                   font.x = 16,
                   font.y = 16,
                   font.tickslab = 16,
                   palette = c("steelblue", "red"),
  )
p5  

# Look at mortality and compositve across VE


km.model6a <- survfit(Surv(time = ansvt[,75] , event = ansvt$Mortality) ~ output, data = ansvt)
km.model6a
summary(km.model6a)
p6a <- ggsurvplot(km.model6a,
                  surv.scale = "percent",
                  risk.table = TRUE,
                  pval=TRUE,
                  pval.method = TRUE,
                  font.legend.labs = c(14),
                  xlab = c("Months"),
                  ylab = c("Survival free from all-cause mortality"),
                  legend.labs = c("VE < 0.5%", "VE ≥ 0.5%"),
                  legend.title = c(""),
                  xscale = c("d_m"),
                  xlim = c(0,70 * 30.437),
                  break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"),
                  risk.table.fontsize = 4,
                  risk.table.y.text = FALSE,
                  risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE,
                  #tables.height = 0.15,
                  censor = TRUE,
                  font.x = 16,
                  font.y = 16,
                  font.tickslab = 16,
                  palette = c("steelblue", "red"),
)
p6a

# Composite by VE figure
km.model6b <- survfit(Surv(time = ansvt[,77] , event = ansvt[,76]) ~ output, data = ansvt)
km.model6b
summary(km.model6b)
p6b <- ggsurvplot(km.model6b,
                  surv.scale = "percent",
                  risk.table = TRUE,
                  pval=TRUE,
                  pval.method = TRUE,
                  xlab = c("Months"),
                  ylab = c("Survival free from VT/VF or mortality"),
                  legend.labs = c("VE < 0.5%", "VE ≥ 0.5%"),
                  font.legend.labs = c(14),
                  legend.title = c(""),
                  font.title = c(20, "bold", "black"),
                  xscale = c("d_m"),
                  xlim = c(0,70 * 30.437),
                  break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"),
                  risk.table.fontsize = 4,
                  risk.table.y.text = FALSE,
                  risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE,
                  tables.height = 0.15,
                  censor = TRUE,
                  font.x = 16,
                  font.y = 16,
                  font.tickslab = 16,
                  palette = c("steelblue", "red"),
)
p6b

survplotlist <- list()
survplotlist[[1]] <- p4 + labs(tag = 'A')
survplotlist[[2]] <- p6a + labs(tag = 'B')
survplotlist[[3]] <- p6b + labs(tag = 'C')

survplotlist
arrange_ggsurvplots(survplotlist, ncol = 3, nrow = 1, labels = "AUTO", 
                    risk.table.height = 0.15)

centfig <- list()
centfig[[1]] <- p2 + labs(tag = 'A')
centfig[[2]] <- p3a + labs(tag = 'B')
centfig[[3]] <- p4 + labs(tag = 'C')
centfig[[4]] <- p6a + labs(tag = 'D')

arrange_ggsurvplots(centfig, arrange.by.row = TRUE , nrow = 2, labels = "AUTO")

?arrange_ggsurvplots


layoutmatrix <- matrix(c(1,2,3,4), nrow = 2, ncol = 2, byrow = TRUE)
arranged_plots <- arrange_ggsurvplots(list(p2, p3a, p4, p6a),layout_matrix = layoutmatrix)
print(arranged_plots)



myplot %>% ggsave(device="png", filename="Figure.png", width = 15, height = 5, units = "in")

plota <- p2$plot
plotb <- p3a$plot
plotc <- p4$plot
plotd <- p6a$plot
x11()
grid.arrange(plota, plotb, plotc, plotd, nrow = 2)


# Reviewer commments
#####


# Indications for ICD
freq(ansvt$Incidcation.for.ICD)

# Antiarrhythmics at first holter
freq(ansvt$Antiarrhythmic....name)

# Morphology of EP and EPS


# LGE analysis
  # Number with CMR
freq(ansvt$Cardiac.MR.obtained.1...yes..0...no)


# LGE in cox analysis

lge <- subset(ansvt,subset = ansvt$Cardiac.MR.obtained.1...yes..0...no == 1)
freq(lge$VT_VF_WS)
freq(lge$Subendo.LGE)
freq(lge$Midwall.LGE)
freq(lge$Subepi.LGE)

        # For the 25
        lge <- subset(ansvt, subset = ansvt$Cardiac.MR.obtained.1...yes..0...no == 1 & ansvt$CMR.to.index > -1)
        dim(lge)
        lge$VT_VF_WS

    # Time Index to CMR
    summary(lge$CMR.to.index) / 365.25
    lge$CMR.to.index / 365.25
    
    sum(lge$CMR.to.index >= 0)
    
    nsvtcmr <- subset(x = lge,subset = lge$NSVT ==1)
    dim(nsvtcmr)
    nonsvtcmr <- subset(x = lge, subset = lge$NSVT == 0)
    dim(nonsvtcmr)    
    
    summary(nsvtcmr$CMR.to.index) / 365.25  
    summary(nonsvt$CMR.to.index) / 365.25
    t.test(x = nsvtcmr$CMR.to.index, y = nonsvt$CMR.to.index)
        
subendo <- ifelse(lge$Subendo.LGE == 1 & lge$Midwall.LGE == 0 & lge$Subepi.LGE == 0,1,0)
subendo  
sum(subendo)

lgeendoonly <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ subendo)
summary(lgeendoonly)
cox.zph(lgeendoonly)

transmural <- rep(0,93) 
transmural <- ifelse(lge$Subendo.LGE == 1 & lge$Midwall.LGE == 1 & lge$Subepi.LGE == 1,1,0)
transmural
sum(transmural)

lgetransmural <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ transmural)
summary(lgetransmural)
cox.zph(lgetransmural)

extensivelge <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ lge$LGE2nd)
summary(extensivelge)
cox.zph(extensivelge)


lgescore <- rep(0,93)
for (i in 1:93) {
  if (subendo[i] == 1){
    lgescore[i] <- 1
  } else if (transmural[i] == 1) {
    lgescore[i] <- 2
  }
}
lgescore

# Based on transmural vs. subepi
tmcompare <- lge

lgetm <- rep(0,93)
for (i in 1:93) {
  if(transmural[i] == 1) {
    lgetm[i] <- 1
  }
}
lgetm
sum(lgetm)

lgeendo <- rep(0,93)
for (i in 1:93) {
  if(subendo[i] == 1) {
    lgeendo[i] <- 1
  }
}
lgeendo
sum(lgeendo)

km.model1 <- survfit(Surv(time = lge[,72] , event = lge[,63]) ~ lgetm, data = lge)
km.model1

km.model1 <- survfit(Surv(time = lge[,72] , event = lge[,63]) ~ lgeendo, data = lge)
km.model1

summary(km.model1)
p2 <- ggsurvplot(km.model1, surv.scale = "percent", risk.table = TRUE,
                 pval=TRUE, pval.method = TRUE, xlab = c("Months"),
                 ylab = c("Survival free from VT/VF"), font.legend.labs = c(18),
                 legend.labs = c("other","transmural"), legend.title = c(""),
                 xscale = c("d_m"), xlim = c(0,70*30.437), break.x.by = 10* 30.437,
                 risk.table.title = c("At Risk"), risk.table.fontsize = 4,
                 risk.table.y.text = FALSE, risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE, #tables.height = 0.15,
                 censor = TRUE, font.x = 16, font.y = 16,
                 font.tickslab = 16, palette = c("steelblue", "red"),
)
p2
km.model3a <- survfit(Surv(time = lge[,75] , event = lge$Mortality) ~ lgetm, data = lge)
km.model3a
km.model3a <- survfit(Surv(time = lge[,75] , event = lge$Mortality) ~ lgeendo, data = lge)
km.model3a
summary(km.model3a)
p3a <- ggsurvplot(km.model3a, surv.scale = "percent", risk.table = TRUE,
                  pval=TRUE, pval.method = TRUE, font.legend.labs = c(18),
                  xlab = c("Months"), ylab = c("Survival free from all-cause mortality"),
                  legend.labs = c("other", "transmural"), legend.title = c(""),
                  xscale = c("d_m"), xlim = c(0,70 * 30.437), break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"), risk.table.fontsize = 4,
                  risk.table.y.text = FALSE, risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE, #tables.height = 0.15,
                  censor = TRUE, font.x = 16, font.y = 16, font.tickslab = 16,
                  palette = c("steelblue", "red"),
)
p3a

# Using second column - extensive or diffuse vs patchy/focal
km.model1 <- survfit(Surv(time = lge[,72] , event = lge[,63]) ~lge$LGE2nd, data = lge)
km.model1
summary(km.model1)

p2 <- ggsurvplot(km.model1,
                 surv.scale = "percent",
                 risk.table = TRUE,
                 pval=TRUE,
                 pval.method = TRUE,
                 xlab = c("Months"),
                 ylab = c("Survival free from VT/VF"),
                 font.legend.labs = c(18),
                 legend.labs = c("focal/none", "extensive"),
                 legend.title = c(""),
                 xscale = c("d_m"),
                 xlim = c(0,70*30.437),
                 break.x.by = 10* 30.437,
                 risk.table.title = c("At Risk"),
                 risk.table.fontsize = 4,
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE,
                 #tables.height = 0.15,
                 censor = TRUE,
                 font.x = 16,
                 font.y = 16,
                 font.tickslab = 16,
                 palette = c("steelblue", "red"),
)
p2

km.model3a <- survfit(Surv(time = lge[,75] , event = lge$Mortality) ~ lge$LGE2nd, data = lge)
km.model3a
summary(km.model3a)
p3a <- ggsurvplot(km.model3a, surv.scale = "percent", risk.table = TRUE,
                  pval=TRUE, pval.method = TRUE, font.legend.labs = c(18),
                  xlab = c("Months"), ylab = c("Survival free from all-cause mortality"),
                  legend.labs = c("subendo", "transmural"), legend.title = c(""),
                  xscale = c("d_m"), xlim = c(0,70 * 30.437), break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"), risk.table.fontsize = 4,
                  risk.table.y.text = FALSE, risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE, #tables.height = 0.15,
                  censor = TRUE, font.x = 16, font.y = 16, font.tickslab = 16,
                  palette = c("steelblue", "red"),
)
p3a



# Broken down by where LGE is
lgeendo <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ lge$Subendo.LGE)
summary(lgeendo)
cox.zph(lgeendo)

lgemid <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ lge$Midwall.LGE)
summary(lgemid)
cox.zph(lgemid)

lgeepi <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ lge$Subepi.LGE)
summary(lgeepi)
cox.zph(lgeepi)

# Cox for lge score
lgescore <- lge$Subendo.LGE + lge$Midwall.LGE + lge$Subepi.LGE
lgescore

length(lgescore)
dim(lge)


lgequant <- coxph(Surv(time = lge[,72], event = lge[,63]) ~ lgescore)
summary(lgequant)
cox.zph(lgequant)

# Kaplan meier for those with a CMR and LGE score
km.model1 <- survfit(Surv(time = lge[,72] , event = lge[,63]) ~lgescore, data = lge)
km.model1
summary(km.model1)
p2 <- ggsurvplot(km.model1,
                 surv.scale = "percent",
                 risk.table = TRUE,
                 pval=TRUE,
                 pval.method = TRUE,
                 xlab = c("Months"),
                 ylab = c("Survival free from VT/VF"),
                 font.legend.labs = c(18),
                 legend.labs = c("0", "1", "2", "3"),
                 legend.title = c(""),
                 xscale = c("d_m"),
                 xlim = c(0,70*30.437),
                 break.x.by = 10* 30.437,
                 risk.table.title = c("At Risk"),
                 risk.table.fontsize = 4,
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE,
                 #tables.height = 0.15,
                 censor = TRUE,
                 font.x = 16,
                 font.y = 16,
                 font.tickslab = 16,
                 palette = c("black", "steelblue", "red", "green"),
)
p2

km.model3a <- survfit(Surv(time = lge[,75] , event = lge$Mortality) ~ lgescore, data = lge)
km.model3a
summary(km.model3a)
p3a <- ggsurvplot(km.model3a,
                  surv.scale = "percent",
                  risk.table = TRUE,
                  pval=TRUE,
                  pval.method = TRUE,
                  font.legend.labs = c(18),
                  xlab = c("Months"),
                  ylab = c("Survival free from all-cause mortality"),
                  legend.labs = c("0", "1", "2", "3"),
                  legend.title = c(""),
                  xscale = c("d_m"),
                  xlim = c(0,70 * 30.437),
                  break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"),
                  risk.table.fontsize = 4,
                  risk.table.y.text = FALSE,
                  risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE,
                  #tables.height = 0.15,
                  censor = TRUE,
                  font.x = 16,
                  font.y = 16,
                  font.tickslab = 16,
                  palette = c("black", "steelblue", "red", "green"),
)
p3a

# Composite by NSVT figure
km.model3b <- survfit(Surv(time = lge[,77] , event = lge[,76]) ~ lgescore, data = lge)
km.model3b
summary(km.model3b)
p3b <- ggsurvplot(km.model3b,
                  surv.scale = "percent",
                  risk.table = TRUE,
                  pval=TRUE,
                  pval.method = TRUE,
                  xlab = c("Months"),
                  ylab = c("Survival free from VT/VF or mortality"),
                  legend.labs = c("0", "1", "2", "3"),
                  font.legend.labs = c(18),
                  legend.title = c(""),
                  font.title = c(20, "bold", "black"),
                  xscale = c("d_m"),
                  xlim = c(0,70 * 30.437),
                  break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"),
                  risk.table.fontsize = 4,
                  risk.table.y.text = FALSE,
                  risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE,
                  tables.height = 0.15,
                  censor = TRUE,
                  font.x = 16,
                  font.y = 16,
                  font.tickslab = 16,
                  palette = c("black", "steelblue", "red", "green"),
)
p3b


# Crossover question
  # Including variable for people that crossed over
nsvtall <- rep(0,217)  
nsvtany <- ifelse(!is.na(ansvt[, 16]) & ansvt[, 16] == 1 |
                    !is.na(ansvt[, 97]) & ansvt[, 97] == 1 |
                    !is.na(ansvt[, 103]) & ansvt[, 103] == 1 |
                    !is.na(ansvt[, 109]) & ansvt[, 109] == 1 |
                    !is.na(ansvt[, 115]) & ansvt[, 115] == 1 |
                    !is.na(ansvt[, 121]) & ansvt[, 121] == 1 |
                    !is.na(ansvt[, 127]) & ansvt[, 127] == 1, 1, 0)
nsvtany
sum(nsvtany)



nsvtallholter <- coxph(Surv(time = ansvt[,72], event = ansvt[,63]) ~ nsvtany)
summary(nsvtallholter)
cox.zph(nsvtallholter)

nsvtyany <- ifelse(!is.na(nsvt[, 16]) & nsvt[, 16] == 1 |
                    !is.na(nsvt[, 97]) & nsvt[, 97] == 1 |
                    !is.na(nsvt[, 103]) & nsvt[, 103] == 1 |
                    !is.na(nsvt[, 109]) & nsvt[, 109] == 1 |
                    !is.na(nsvt[, 115]) & nsvt[, 115] == 1 |
                    !is.na(nsvt[, 121]) & nsvt[, 121] == 1 |
                    !is.na(nsvt[, 127]) & nsvt[, 127] == 1, 1, 0)
nsvtyany
sum(nsvtyany)

nsvtnany <- ifelse(!is.na(nonsvt[, 16]) & nonsvt[, 16] == 1 |
                     !is.na(nonsvt[, 97]) & nonsvt[, 97] == 1 |
                     !is.na(nonsvt[, 103]) & nonsvt[, 103] == 1 |
                     !is.na(nonsvt[, 109]) & nonsvt[, 109] == 1 |
                     !is.na(nonsvt[, 115]) & nonsvt[, 115] == 1 |
                     !is.na(nonsvt[, 121]) & nonsvt[, 121] == 1 |
                     !is.na(nonsvt[, 127]) & nonsvt[, 127] == 1, 1, 0)
nsvtnany
sum(nsvtnany)

# VE - stats for supplement
vensvt <- matrix(c(85, 45, 26, 61), byrow = T, ncol = 2, nrow = 2)
colnames(vensvt) <- c("NSVT", "No NSVT")
rownames(vensvt) <- c("eflow", "no eflow")
vensvt
chisq.test(vensvt)
fisher.test(vensvt)

# EF - stats for supplement
vensvt <- matrix(c(17, 113, 7, 80), byrow = T, ncol = 2, nrow = 2)
colnames(vensvt) <- c("NSVT", "No NSVT")
rownames(vensvt) <- c("eflow", "no eflow")
vensvt
chisq.test(vensvt)
fisher.test(vensvt)



km.model1 <- survfit(Surv(time = ansvt[,72] , event = ansvt[,63]) ~nsvtany, data = ansvt)
km.model1
summary(km.model1)
p2 <- ggsurvplot(km.model1,
                 surv.scale = "percent",
                 risk.table = TRUE,
                 pval=TRUE,
                 pval.method = TRUE,
                 xlab = c("Months"),
                 ylab = c("Survival free from VT/VF"),
                 font.legend.labs = c(18),
                 legend.labs = c("No NSVT", "NSVT"),
                 legend.title = c(""),
                 xscale = c("d_m"),
                 xlim = c(0,70*30.437),
                 break.x.by = 10* 30.437,
                 risk.table.title = c("At Risk"),
                 risk.table.fontsize = 4,
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 tables.y.text = FALSE,
                 #tables.height = 0.15,
                 censor = TRUE,
                 font.x = 16,
                 font.y = 16,
                 font.tickslab = 16,
                 palette = c("steelblue", "red"),
)
p2


# Old Figure 3A - Mortality by NSVT figure
ansvt[,75]
ansvt[,75] / 30.437


km.model3a <- survfit(Surv(time = ansvt[,75] , event = ansvt$Mortality) ~ nsvtany, data = ansvt)
km.model3a
summary(km.model3a)
p3a <- ggsurvplot(km.model3a,
                  surv.scale = "percent",
                  risk.table = TRUE,
                  pval=TRUE,
                  pval.method = TRUE,
                  font.legend.labs = c(18),
                  xlab = c("Months"),
                  ylab = c("Survival free from all-cause mortality"),
                  legend.labs = c("No NSVT", "NSVT"),
                  legend.title = c(""),
                  xscale = c("d_m"),
                  xlim = c(0,70 * 30.437),
                  break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"),
                  risk.table.fontsize = 4,
                  risk.table.y.text = FALSE,
                  risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE,
                  #tables.height = 0.15,
                  censor = TRUE,
                  font.x = 16,
                  font.y = 16,
                  font.tickslab = 16,
                  palette = c("steelblue", "red"),
)
p3a

# Composite by NSVT figure
km.model3b <- survfit(Surv(time = ansvt[,77] , event = ansvt[,76]) ~ nsvtany, data = ansvt)
km.model3b
summary(km.model3b)
p3b <- ggsurvplot(km.model3b,
                  surv.scale = "percent",
                  risk.table = TRUE,
                  pval=TRUE,
                  pval.method = TRUE,
                  xlab = c("Months"),
                  ylab = c("Survival free from VT/VF or mortality"),
                  legend.labs = c("No NSVT", "NSVT"),
                  font.legend.labs = c(18),
                  legend.title = c(""),
                  font.title = c(20, "bold", "black"),
                  xscale = c("d_m"),
                  xlim = c(0,70 * 30.437),
                  break.x.by = 10 * 30.437,
                  risk.table.title = c("At Risk"),
                  risk.table.fontsize = 4,
                  risk.table.y.text = FALSE,
                  risk.table.y.text.col = TRUE,
                  tables.y.text = FALSE,
                  tables.height = 0.15,
                  censor = TRUE,
                  font.x = 16,
                  font.y = 16,
                  font.tickslab = 16,
                  palette = c("steelblue", "red"),
)
p3b

# morphology question
setwd('~/Siontis research/')
vtablation <- read.csv('IVIEW_PROCEDURES_IDENTIFIED (1) - 2024.02.02.csv', na.strings = c("", "NA"))
dim(vtablation)
View(vtablation)

common_patients <- inner_join(vtablation, ansvt,by = c("Current.Clinic" = "Mayo.MRN"))
View(common_patients)

for (i in 1:217) {
  for (j in 1:1406) {
    if (ansvt[i,2] == vtablation[j,1]){
      print(ansvt[i,2])
    }
  }
}

# MFDFA
Application of the MultiFractal Detrended Fluctuation Analysis (MFDFA) on Time Series

Version: 0.1.0

Author: Mohamed Laib, Luciano Telesca and Mikhail Kanevski

Maintainer: Mohamed Laib <Mohamed.Laib@unil.ch>

URL: https://sites.google.com/site/mohamedlaibwebpage/

License: GPL-3

Note: This R code was developed and used for the following paper:
   Long-range fluctuations and multifractality in connectivity density
   time series of a wind speed monitoring network, submitted.

## MFDFA package installation: from github ####

install.packages("devtools")

devtools::install_github("mlaib/MFDFA")

library(MFDFA)

## Example #####

a<-0.9

N<-1024

tsx<-MFsim(N,a)

scale=10:100

q<--10:10

m<-1

mfdfa<-MFDFA(tsx, scale, m, q)


## Results plot ####

dev.new()

par(mai=rep(1, 4))

plot(q, mfdfa$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent",
     ylim=c(min(mfdfa$Hq),max(mfdfa$Hq)))
     
grid(col="midnightblue")

axis(1)

axis(2)


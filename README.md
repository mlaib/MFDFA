# MFDFA: MultiFractal Detrended Fluctuation Analysis for Time Series
Applies the MultiFractal Detrended Fluctuation Analysis (MFDFA) to time series. The package contains some suggestion plot of the MFDFA results.

#### Version 
1.0

#### Authors 
Mohamed Laib, Luciano Telesca and Mikhail Kanevski

#### Maintainer
Mohamed Laib <Mohamed.Laib@unil.ch>

#### URL
https://mlaib.github.io

#### License
GPL-3

#### Note
This R code was developed and used for the following papers:
  
 * M. Laib, L. Telesca, M. Kanevski, Long-range fluctuations and multifractality in connectivity density time series of a wind speed monitoring network, Chaos: An Interdisciplinary Journal of Nonlinear Science 28 (3), 033108.
   
 * M. Laib, J. Golay, L. Telesca, M. Kanevski, Multifractal analysis of the time series of daily means of wind speed in complex regions, Chaos, Solitons & Fractals, 109 (2018) pp. 118-127.

### MFDFA package installation: from github 
```{r}
install.packages("devtools")
devtools::install_github("mlaib/MFDFA")
library(MFDFA)
```

#### Example 
```{r}
a<-0.9
N<-1024
tsx<-MFsim(N,a)
scale=10:100
q<--10:10
m<-1
mfdfa<-MFDFA(tsx, scale, m, q)
```

#### Results plot 
```{r}
dev.new()
par(mai=rep(1, 4))
plot(q, mfdfa$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
     cex.axis=1.8, main="Hurst exponent",
     ylim=c(min(mfdfa$Hq),max(mfdfa$Hq)))
grid(col="midnightblue")
axis(1)
axis(2)
```

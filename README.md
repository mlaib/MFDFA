## MFDFA: MultiFractal Detrended Fluctuation Analysis for Time Series
Applies the MultiFractal Detrended Fluctuation Analysis (MFDFA) to time series. The package contains some suggestion plot of the MFDFA results.

The MFDFA R library is now available on CRAN. Further update will be added soon.

A new file is available [Here](https://gist.github.com/mlaib/bb0c09df9593dad16ae270334ec3e7d7). It proposes the MFDFA with a parallel version (MFDFA2.R). Useful for long time series. It can be used as the first one with same parameters. It uses (N-1) of workers of your computer. ENJOY ...

#### Version 
1.0

#### Authors 
Mohamed Laib, Luciano Telesca and Mikhail Kanevski

#### Maintainer
Mohamed Laib [mohamed.laib (at) unil.ch] or 
             [laib.med (at) gmail.com]

#### URL
[https://cran.r-project.org/package=MFDFA](https://cran.r-project.org/package=MFDFA)

[https://mlaib.github.io/MFDFA/](https://mlaib.github.io/MFDFA/)

[https://mlaib.github.io](https://mlaib.github.io)



#### License
GPL-3

#### Note
This R code was developed and used for the following papers:
  
 * M. Laib, L. Telesca, M. Kanevski, Long-range fluctuations and multifractality in connectivity density time series of a wind speed monitoring network, Chaos: An Interdisciplinary Journal of Nonlinear Science 28 (3), 033108. [Paper](https://www.researchgate.net/publication/319121707_Long-range_fluctuations_and_multifractality_in_connectivity_density_time_series_of_a_wind_speed_monitoring_network)
   
 * M. Laib, J. Golay, L. Telesca, M. Kanevski, Multifractal analysis of the time series of daily means of wind speed in complex regions, Chaos, Solitons & Fractals, 109 (2018) pp. 118-127. [Paper](https://www.researchgate.net/publication/320223480_Multifractal_analysis_of_the_time_series_of_daily_means_of_wind_speed_in_complex_regions)

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

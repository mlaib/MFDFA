#' Multiscale Multifractal Analysis
#'
#' Applies the Multiscale Multifractal Analysis (MMA) on time series.
#' @usage MMA(tsx, scale, qminmax, ovlap=0, m=2)
#' @param tsx  Univariate time series (must be a vector or a ts object).
#' @param scale Vector of scales.
#' @param qminmax Vector of two values min and max of q-order of the moment.
#' @param ovlap Overlapping parameter (By default ovlap=0: no overlapping).
#' @param m Polynomial order for the detrending (by defaults m=2).
#' @return A matrix with three columns (q-order, scale (s), and the scale exponent).
#'
#' @note The original code of this function is in Matlab, you can find it on the
#' following website \href{https://physionet.org/physiotools/mma/}{Physionet}. See
#' references below.
#'
#'
#'
#' @examples
#'
#' \dontrun{
#' library(MFDFA)
#' library(plotly)
#' library(plot3D)
#'
#' a<-0.6
#' N<-800
#' tsx<-MFsim(N,a)
#' scale=10:100
#' res<-MMA(tsx, scale, qminmax=c(-10,10), ovlap=0, m=2)
#'
#' ## Visualisation 1:
#' S_exponent <- matrix(res[,3], nrow=length(unique(res[,1])), ncol=length(min(scale):(max(scale)/5)))
#' m_scale <- unique(res[,2])
#' q <- unique(res[,1])
#' plot_ly() %>% add_surface(x = ~m_scale, y = ~q,
#'                          z = ~S_exponent)
#'
#' ## Visualisation 2:
#' image2D(S_exponent, xlab="q", ylab="scale", axes=F)
#' axis(1, seq(0,1,0.1), round(quantile(q, seq(0, 1, 0.1)), 2))
#' axis(2, seq(0,1,0.1), round(quantile(m_scale, seq(0, 1, 0.1)), 2))
#' }
#'
#' @references
#' J. Feder, Fractals, Plenum Press, New York, NY, USA, 1988.
#'
#' J. Gieraltowski, J. J. Zebrowski, and R. Baranowski,
#' Multiscale multifractal analysis of heart rate variability recordings
#' http://dx.doi.org/10.1103/PhysRevE.85.021915
#'
#' Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
#' Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit,
#' and PhysioNet: Components of a New Research Resource for Complex
#' Physiologic Signals. Circulation 101(23):e215-e220.
#'
#' J. W. Kantelhardt, S. A. Zschiegner, E. Koscielny-Bunde, S. Havlin,
#' A. Bunde, H. Stanley, Multifractal detrended fluctuation analysis of
#' nonstationary time series, Physica A: Statistical Mechanics and its
#' Applications, 316 (1) (2002) 87 – 114.
#'
#' J. Gierałtowski, J. J. Żebrowski, and R. Baranowski, "Multiscale
#' multifractal analysis of heart rate variability recordings with a
#' large number of occurrences of arrhythmia," Phys. Rev. E 85, 021915 (2012)
#'
#'
#' @importFrom numbers mod
#' @export
#'
#'

MMA <- function(tsx, scale, qminmax, ovlap=0, m=2){
  qrange <- seq(qminmax[1], qminmax[2], 0.1)
  qrange[which(qrange==0)]<-0.0001
  prof = cumsum(tsx)
  slength = length(prof)
  fqs<-c()
  for (i in scale){
    if (ovlap==1){
      vec <- 0:(i-1)
      ind <- 1:(slength-i+1)
      A<-matrix(rep(vec, each = length(ind)), length(ind), length(vec))
      coorxy <- apply(A, 2, FUN=function(x) (x+ind))
    } else {
      nd<-(slength-mod(slength,i))
      nc<-(slength-mod(slength,i))/i
      coorxy<-matrix(1:nd, nc, i, byrow = TRUE)
    }
    segments <- apply(coorxy, 2, FUN=function(x)(prof[x]))
    xbs <- 1:i
    f2nis <- c()

    for (j in 1:nrow(segments)){
      seg <- segments[j,]
      ft <- poly_fit.val(xbs, seg, m)
      fit<-ft$polyval
      f2nis <-c(f2nis, mean((seg- fit)^2))

    }

    for (qq in qrange){
      fqs<- c(fqs, qq, i, mean(f2nis^(qq/2))^(1/qq))
    }


  }

  fqs1 <- matrix(fqs, ncol=3, byrow = T)
  fqs1 <- matrix(append(fqs1, log(fqs1[,2])), ncol=4)
  fqs1 <- matrix(append(fqs1, log(fqs1[,3])), ncol=5)


  hqs <- c()

  for (i in min(scale):(max(scale)/5)){
    for (qq in qrange){
      tempfqs <- fqs1[which(fqs1[,1]==qq & fqs1[,2] >= i & fqs1[,2] <= 5*i),]
      hft <- poly_fit.val(tempfqs[,4], tempfqs[,5], 1)
      hqs <- c(hqs, qq, i, as.numeric(hft$polyfit[1]))
    }
  }

  hqs1 <- matrix(hqs, ncol=3, byrow = T)
  colnames(hqs1)<- c("q", "s", "S_Exp")
  return(hqs1)
}

poly_fit.val<-function(x,y,n){
  formule<-lm(as.formula(paste('y~',paste('I(x^',1:n,')', sep='',collapse='+'))))
  res1<-coef(formule)
  poly.res<-res1[length(res1):1]
  suppressWarnings(res2<-predict(formule, interval = "prediction"))
  poly.eva<-res2[,1]
  allres<-list(polyfit=round(poly.res,4), polyval=poly.eva)
  return(allres)}

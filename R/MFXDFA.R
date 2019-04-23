#' Multifractal detrended cross-correlation analysis
#'
#' Applies the MultiFractal Detrended Fluctuation cross-correlation Analysis (MFXDFA) on two time series.
#' @usage MFXDFA(tsx1, tsx2, scale, m=1, q)
#' @param tsx1  Univariate time series (must be a vector or a ts object).
#' @param tsx2  Univariate time series (must be a vector or a ts object).
#' @param scale Vector of scales.
#' @param m Polynomial order for the detrending (by default m=1).
#' @param q q-order of the moment. There is no default value
#' for this parameter, please add values.
#' @return A list of the following elements:
#'  \itemize{
#'   \item \code{Hq} Hurst exponent.
#'   \item \code{h} Holder exponent.
#'   \item \code{Dh} Multifractal spectrum.
#'   \item \code{Fq} Fluctuation function in log.
#'   }
#'
#' @note The original code of this function is in Matlab, you can find it on the
#' following website \href{https://ch.mathworks.com/matlabcentral/fileexchange/38262-multifractal-detrended-fluctuation-analyses?focused=5247306&tab=function}{Mathworks}.
#'
#'
#'
#' @examples
#'
#' library(MFDFA)
#' a<-0.6
#' N<-1024
#'
#' tsx1<-MFsim(N,a)
#' b<-0.8
#' N<-1024
#'
#' tsx2<-MFsim(N,b)
#' scale=10:100
#' q<--10:10
#' m<-1
#'
#' \dontrun{
#' b<-MFXDFA(tsx1, tsx2, scale, m=1, q)
#'
#' ## Supplementary functions: #####
#' reset <- function(){
#' par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
#' plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)}
#'
#' poly_fit<-function(x,y,n){
#'   formule<-lm(as.formula(paste('y~',paste('I(x^',1:n,')', sep='',collapse='+'))))
#'   res1<-coef(formule)
#'   poly.res<-res1[length(res1):1]
#'   allres<-list(polyfit=poly.res, model1=formule)
#'   return(allres)}
#'
#' ## Plot results: #####
#' dev.new()
#' layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(4, 4))
#' ## b : mfdfa output
#' par(mai=rep(0.8, 4))
#'
#' ## 1st plot: Fluctuations function
#' p1<-which(q==2)
#' plot(log(scale),b$Fq[,p1],  pch=16, col=1, axes = FALSE, xlab = "s",
#'      ylab=expression('log'*'(F'[2]*')'), cex=1, cex.lab=1.6, cex.axis=1.6,
#'      main= "Fluctuation function F for q=2",
#'      ylim=c(min(b$Fq[,c(p1)]),max(b$Fq[,c(p1)])))
#' lines(log(scale),b$line[,p1], type="l", col=1, lwd=2)
#' grid(col="midnightblue")
#' axis(2)
#' lbl<-scale[c(1,floor(length(scale)/8),floor(length(scale)/4),
#'              floor(length(scale)/2),length(scale))]
#' att<-log(lbl)
#' axis(1, at=att, labels=lbl)
#'
#' ## 2nd plot: q-order Hurst exponent
#' plot(q, b$Hq, col=1, axes= FALSE, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
#'      cex.axis=1.8, main="Hurst exponent", ylim=c(min(b$Hq),max(b$Hq)))
#' grid(col="midnightblue")
#' axis(1, cex=4)
#' axis(2, cex=4)
#'
#' ## 3rd plot: Spectrum
#' plot(b$h, b$Dh, col=1, axes=FALSE, pch=16, main="Multifractal spectrum",
#'      ylab=bquote("f ("~alpha~")"),cex.lab=1.8, cex.axis=1.8,
#'      xlab=bquote(~alpha))
#' grid(col="midnightblue")
#' axis(1, cex=4)
#' axis(2, cex=4)
#'
#' x1=b$h
#' y1=b$Dh
#' rr<-poly_fit(x1,y1,4)
#' mm1<-rr$model1
#' mm<-rr$polyfit
#' x2<-seq(min(x1),max(x1)+1,0.01)
#' curv<-mm[1]*x2^4+mm[2]*x2^3+mm[3]*x2^2+mm[4]*x2+mm[5]
#' lines(x2,curv, col="red", lwd=2)
#' }
#'
#' @references
#' J. Feder, Fractals, Plenum Press, New York, NY, USA, 1988.
#'
#' Espen A. F. Ihlen, Introduction to multifractal detrended fluctuation analysis
#' in matlab, Frontiers in Physiology: Fractal Physiology, 3 (141),(2012) 1-18.
#'
#' J. W. Kantelhardt, S. A. Zschiegner, E. Koscielny-Bunde, S. Havlin,
#' A. Bunde, H. Stanley, Multifractal detrended fluctuation analysis of
#' nonstationary time series, Physica A: Statistical Mechanics and its
#' Applications, 316 (1) (2002) 87 – 114.
#'
#' Kantelhardt J.W. (2012) Fractal and Multifractal Time Series. In: Meyers R. (eds) Mathematics
#' of Complexity and Dynamical Systems. Springer, New York, NY.
#'
#' @export
#'

MFXDFA <- function(tsx1, tsx2, scale, m=1, q){
  if ((length(tsx1)==length(tsx2))==FALSE){
    stop("check the length of these time series")
  }
  N <- length(tsx1)
  X<-cumsum(tsx1-mean(tsx1))
  Y<-cumsum(tsx2-mean(tsx2))

  q0 <- which(q==0)
  Fq <- matrix(0, nrow=length(scale), ncol=length(q))
  for (ns in 1:length(scale)){
    Ns <- floor(N/scale[ns])
    VR <- matrix(0, nrow=Ns, ncol=length(q))

    for (v in 1:Ns){
      SegInd <- ((((v-1)*scale[ns])+1):(v*scale[ns]))

      # 1st time serie
      SegX <- X[SegInd]
      fitX <- poly_fit.val(SegInd, SegX, m)$polyval

      # 2nd time serie
      SegY <- Y[SegInd]
      fitY <- poly_fit.val(SegInd, SegY, m)$polyval


      for (nq in 1:length(q)){
        VR[v,nq] <- ((sum(((SegX-fitX)^2)*((SegY-fitY)^2)))/scale[ns])^(q[nq]/4)
      }

    }
    for (nq in 1:length(q)){
      if ((q[nq]==0)==FALSE){
        Fq[ns,nq] <- ((sum(VR[,nq])+sum(VR[,nq]))/(2*Ns))^(1/q[nq])
      }
    }
    Fq[ns,q0] <- (Fq[ns,(q0+1)]+Fq[ns,(q0-1)])/2


  }

  LogF <- log2(Fq)

  qRegLine <- list()
  Hq <- c()
  for (nq in 1:length(q)){
    polyft<-poly_fit.val(log2(scale),LogF[,nq], 1)

    Hq[nq] <- polyft$polyfit[1]
    qRegLine[[nq]]<-polyft$polyval
  }
  qRegLine <- as.data.frame(Reduce("cbind", qRegLine))
  tau <- (q*Hq)-1

  H <- diff(tau)/diff(q)
  Dh <- (q[-length(q)]*H)-tau[-length(tau)]
  h <- H-1

return(list(Hq=Hq, h=h, Dh=Dh, Fq=LogF, line=qRegLine))
}






poly_fit.val<-function(x,y,n){
  formule<-lm(as.formula(paste('y~',paste('I(x^',1:n,')', sep='',
                                          collapse='+'))))
  res1<-coef(formule)
  poly.res<-res1[length(res1):1]
  suppressWarnings(res2<-predict(formule, interval = "prediction"))
  poly.eva<-res2[,1]
  allres<-list(polyfit=round(poly.res,4), polyval=poly.eva)
  return(allres)}

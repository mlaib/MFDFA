#' MultiFractal Detrended Fluctuation Analysis
#'
#' Applies the MultiFractal Detrended Fluctuation Analysis (MFDFA) to time series.
#' @usage MFDFA(tsx, scale, m=1, q)
#' @param tsx  Univariate time series (must be a vector).
#' @param scale Vector of scales. There is no default value
#' for this parameter, please add values.
#' @param m An integer of the polynomial order for the detrending (by default m=1).
#' @param q q-order of the moment. There is no default value
#' for this parameter, please add values.
#'
#' @return A list of the following elements:
#'  \itemize{
#'   \item \code{Hq} Hurst exponent.
#'   \item \code{tau_q} Mass exponent.
#'   \item \code{spec} Multifractal spectrum (\eqn{\alpha}{\alpha} and
#'   \eqn{f(\alpha)}{f(\alpha)})
#'   \item \code{Fq} Fluctuation function.
#'   }
#'
#' @details The original code of this function is in Matlab, you can find it on the
#' following website \href{https://ch.mathworks.com/matlabcentral/fileexchange/38262-multifractal-detrended-fluctuation-analyses?focused=5247306&tab=function}{Mathworks}.
#'
#'
#' @examples
#'
#' \dontrun{
#' ## MFDFA package installation: from github ####
#' install.packages("devtools")
#' devtools::install_github("mlaib/MFDFA")
#'
#' ## Get the Parellel version:
#' devtools::source_gist("bb0c09df9593dad16ae270334ec3e7d7", filename = "MFDFA2.r")
#' }
#'
#' library(MFDFA)
#' a<-0.9
#' N<-1024
#' tsx<-MFsim(N,a)
#' scale=10:100
#' q<--10:10
#' m<-1
#' b<-MFDFA(tsx, scale, m, q)
#'
#'\dontrun{
#' ## Results plot ####
#' dev.new()
#' par(mai=rep(1, 4))
#' plot(q, b$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
#'      cex.axis=1.8, main="Hurst exponent",
#'      ylim=c(min(b$Hq),max(b$Hq)))
#' grid(col="midnightblue")
#' axis(1)
#' axis(2)
#'
#' ##################################
#' ## Suggestion of output plot: ####
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
#' ##################################
#' ## Output plots: #################
#' dev.new()
#' layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(4, 4))
#' ## b : mfdfa output
#' par(mai=rep(0.8, 4))
#'
#' ## 1st plot: Scaling function order Fq (q-order RMS)
#' p1<-c(1,which(q==0),which(q==q[length(q)]))
#' plot(log2(scale),log2(b$Fqi[,1]),  pch=16, col=1, axes = F, xlab = "s (days)",
#'      ylab=expression('log'[2]*'(F'[q]*')'), cex=1, cex.lab=1.6, cex.axis=1.6,
#'      main= "Fluctuation function Fq",
#'      ylim=c(min(log2(b$Fqi[,c(p1)])),max(log2(b$Fqi[,c(p1)]))))
#' lines(log2(scale),b$line[,1], type="l", col=1, lwd=2)
#' grid(col="midnightblue")
#' axis(2)
#' lbl<-scale[c(1,floor(length(scale)/8),floor(length(scale)/4),
#'              floor(length(scale)/2),length(scale))]
#' att<-log2(lbl)
#' axis(1, at=att, labels=lbl)
#' for (i in 2:3){
#'   k<-p1[i]
#'   points(log2(scale), log2(b$Fqi[,k]),  col=i,pch=16)
#'   lines(log2(scale),b$line[,k], type="l", col=i, lwd=2)
#' }
#' legend("bottomright", c(paste('q','=',q[p1] , sep=' ' )),cex=2,lwd=c(2,2,2),
#'  bty="n", col=1:3)
#'
#' ## 2nd plot: q-order Hurst exponent
#' plot(q, b$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
#'     cex.axis=1.8, main="Hurst exponent", ylim=c(min(b$Hq),max(b$Hq)))
#' grid(col="midnightblue")
#' axis(1, cex=4)
#' axis(2, cex=4)
#'
#' ## 3rd plot: q-order Mass exponent
#' plot(q, b$tau_q, col=1, axes=F, cex.lab=1.8, cex.axis=1.8,
#'      main="Mass exponent",
#'      pch=16,ylab=expression(tau[q]))
#' grid(col="midnightblue")
#' axis(1, cex=4)
#' axis(2, cex=4)
#'
#' ## 4th plot: Multifractal spectrum
#' plot(b$spec$hq, b$spec$Dq, col=1, axes=F, pch=16, #main="Multifractal spectrum",
#'      ylab=bquote("f ("~alpha~")"),cex.lab=1.8, cex.axis=1.8,
#'      xlab=bquote(~alpha))
#' grid(col="midnightblue")
#' axis(1, cex=4)
#' axis(2, cex=4)
#'
#' x1=b$spec$hq
#' y1=b$spec$Dq
#' rr<-poly_fit(x1,y1,4)
#' mm1<-rr$model1
#' mm<-rr$polyfit
#' x2<-seq(0,max(x1)+1,0.01)
#' curv<-mm[1]*x2^4+mm[2]*x2^3+mm[3]*x2^2+mm[4]*x2+mm[5]
#' lines(x2,curv, col="red", lwd=2)
#' reset()
#' legend("top", legend="MFDFA Plots", bty="n", cex=2)
#' }
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
#' Kantelhardt J.W. (2012) Fractal and Multifractal Time Series. In: Meyers R. (eds)
#' Mathematics of Complexity and Dynamical Systems. Springer, New York, NY.
#'
#' M. Laib, L. Telesca and M. Kanevski, Long-range fluctuations and
#' multifractality in connectivity density time series of a wind speed
#' monitoring network, Chaos: An Interdisciplinary Journal of Nonlinear
#' Science, 28 (2018) p. 033108, \href{https://aip.scitation.org/doi/10.1063/1.5022737}{Paper}.
#'
#' M. Laib, J. Golay, L. Telesca, M. Kanevski, Multifractal
#' analysis of the time series of daily means of wind speed
#' in complex regions, Chaos, Solitons & Fractals, 109 (2018)
#' pp. 118-127, \href{https://www.sciencedirect.com/science/article/pii/S0960077918300699}{Paper}.
#'
#' @importFrom stats var
#'
#' @export

MFDFA<-function(tsx, scale, m=1, q){
  if (!is.numeric(tsx) | var(tsx)==0){
    stop("Check your time series")
  }

  if (m%%1 != 0){
    stop("m must be an integer")
  }

  if (!is.numeric(scale) | length(scale) <= 1 | any(scale%%1 != 0)){
    stop("scale must be a vector containing integers")
  }

  if (!is.numeric(q) | any(q%%1 != 0)){
    stop("q must contain integers")
  }


  X<-cumsum(tsx-mean(tsx))
  seg<-list()
  qRMS<-list()
  Fq<-c()
  Fqi<-list()
  Hq<-c()
  qRegLine<-list()
  RMSvi<-list()
  for (i in 1:length(scale)){
    seg[[i]]<-floor(length(X)/scale[i])
    rmvi<-c()
    for (vi in 1:seg[[i]]){
      Index=((((vi-1)*scale[i])+1):(vi*scale[i]))
      polyft<-poly_fit.val( Index,X[Index], m)
      C<-polyft$polyfit
      fit<-polyft$polyval
      rmvi[vi]<-sqrt(mean((fit-X[Index])^2))
      RMSvi[[i]]<-rmvi
    }


    for (nq in 1:length(q)){
      rod<-RMSvi[[i]]^q[nq]
      qRMS[[nq]]<-rod
      Fq[nq]<-mean(qRMS[[nq]])^(1/q[nq])

    }
    if (any(q==0)){Fq[which(q==0)]<-exp(0.5*mean(log(RMSvi[[i]]^2)))}

    Fqi[[i]]<-Fq

  }
  Fqi<-Reduce("rbind", Fqi)
  for (nq in 1:length(q)){
    polyft<-poly_fit.val( log2(scale),log2(Fqi[,nq]), 1)
    C<-polyft$polyfit
    Hq[nq]<-C[1]
    qRegLine[[nq]]<-polyft$polyval
  }
  qRegLine<-as.data.frame(Reduce("cbind", qRegLine))
  tq<-Hq*q-1
  hq<-diff(tq)/(q[2]-q[1])
  Dq<-(q[1:(length(q)-1)]*hq)-tq[1:(length(tq)-1)]
  return(list(Hq=Hq, tau_q=tq, spec=data.frame(hq=hq, Dq=Dq),
              Fqi=Fqi, line=qRegLine))}


## intern functions: polyfit ####
poly_fit.val<-function(x,y,n){
  formule<-lm(as.formula(paste('y~',paste('I(x^',1:n,')', sep='',collapse='+'))))
  res1<-coef(formule)
  poly.res<-res1[length(res1):1]
  suppressWarnings(res2<-predict(formule, interval = "prediction"))
  poly.eva<-res2[,1]
  allres<-list(polyfit=round(poly.res,4), polyval=poly.eva)
  return(allres)}

reset <- function(){
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

poly_fit<-function(x,y,n){
  formule<-lm(as.formula(paste('y~',paste('I(x^',1:n,')', sep='',collapse='+'))))
  res1<-coef(formule)
  poly.res<-res1[length(res1):1]
  allres<-list(polyfit=poly.res, model1=formule)
  return(allres)}


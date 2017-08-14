#' Simulated multifractal series.
#'
#'
#' Generates series using the binomial multifractal model (see references).
#' @usage MFsim(N,a)
#' @param N The length of the generated multifractal series.
#' @param a The Hurst exponent, which takes values in [0.5, 1].
#'
#' @return A vector containing the multifractal series
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' a<-0.9
#' N<-1024
#' tsx<-MFsim(N,a)
#'
#' scale=10:100
#' q<--10:10
#' m<-1
#' b<-MFDFA(tsx, scale, m, q)
#' dev.new()
#' par(mai=rep(1, 4))
#' plot(q, b$Hq, col=1, axes= F, ylab=expression('h'[q]), pch=16, cex.lab=1.8,
#'      cex.axis=1.8, main="q-order Hurst exponent", ylim=c(min(b$Hq),max(b$Hq)))
#' grid(col="midnightblue")
#' axis(1)
#' axis(2)
#' }
#' @references
#'
#' J. Feder, Fractals, Plenum Press, New York, NY, USA, 1988.
#'
#' E.L. Flores-Márquez, A. Ramírez-Rojas, L. Telesca, Multifractal detrended
#' fluctuation analysis of earthquake magnitude series of Mexican South Pacific
#' Region, Applied Mathematics and Computation, Volume 265, 2015,
#' Pages 1106-1114, ISSN 0096-3003.
#'
#' @importFrom graphics par plot
#' @importFrom stats as.formula coef lm predict
#' @export

MFsim<-function(N,a){
  m<-1:N
  b1<-a^(nbit(m-1))
  b2<-(1-a)^(16-nbit(m-1))
  XM<-b1*b2
  return(XM)
}

# intern function: nbit
nbit<-function(num){
  num<-as.matrix(num)
  s<-list()
  s<-apply(num, 1, intToBits)
  rs<-apply(s, 2, FUN = function(xa) (length(which(xa==1))))
  return(rs)
}

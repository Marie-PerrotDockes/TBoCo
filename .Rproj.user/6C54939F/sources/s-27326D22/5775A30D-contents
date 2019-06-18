
sf.test<-function (x)
{
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  if ((n < 5 || n > 5000))
    stop("sample size must be between 5 and 5000")
  y <- qnorm(ppoints(n, a = 3/8))
  W <- cor(x, y)^2
  u <- log(n)
  v <- log(u)
  mu <- -1.2725 + 1.0521 * (v - u)
  sig <- 1.0308 - 0.26758 * (v + 2/u)
  z <- (log(1 - W) - mu)/sig
  pval <- pnorm(z, lower.tail = FALSE)
  RVAL <- list(statistic = c(W = W), p.value = pval, method = "Shapiro-Francia normality test",
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}



#' Title
#'
#' @param i optional name of the vector which is being transform
#' @param data the vector to transform
#' @param method which method should be used. Default to "all"
#' @param lam the sequence of lambda to try
#' @param plotit logical
#'
#' @return a numeric vector transform using boxcox transformation
#' @export
#' @import tseries
#' @examples
#' x <- rnorm(100)
#' boxcoxnc2(data =exp(x),lam = seq(-2, 2, 0.5))
boxcoxnc2<-function (i = "data to transform",
                     data, lam = seq(-2, 2, 0.01), plotit = TRUE,
                      p.method = "BY")
{
  data <- as.numeric(data)
  if (is.na(min(data)) == TRUE)
    stop("Data include NA")
  if (min(data) <= 0)
    stop("Data must include positive values")

    jb <- NULL
    for (i in 1:length(lam)) {
      if (round(lam[i], 2) != 0)
        jb <- rbind(jb, c(lam[i], jarque.bera.test((data^(lam[i]) -
                                                      1)/(lam[i]))$statistic))
      if (round(lam[i], 2) == 0)
        jb <- rbind(jb, c(lam[i], jarque.bera.test(log(data))$statistic))
    }
    jblam <- jb[which.min(jb[, 2]), 1]
    jbstat <- jb[which.min(jb[, 2]), 2]
    if (plotit == TRUE) {
      plot(jb[, 1], jb[, 2], ylab = "test statistic", xlab = expression(lambda),
           main = "Jarque-Bera")
      abline(v = jblam, lty = 2)
    }
    if (jblam == max(lam))
      print(paste( i ,'not optimal'))
    if (jblam == min(lam))
      print(paste( i ,'not optimal'))
    if (jblam != 0)
      data.jb <- ((data^jblam) - 1)/jblam
    if (jblam == 0)
      data.jb <- log(data)
    sw.pvalue_data.jb <- p.adjust(c(shapiro.test(data.jb)$p.value,
                                    sf.test(data.jb)$p.value, jarque.bera.test(data.jb)$p.value),
                                  p.method)[1]
    sf.pvalue_data.jb <- p.adjust(c(shapiro.test(data.jb)$p.value,
                                    sf.test(data.jb)$p.value, jarque.bera.test(data.jb)$p.value),
                                  p.method)[2]
    jb.pvalue_data.jb <- p.adjust(c(shapiro.test(data.jb)$p.value,
                                    sf.test(data.jb)$p.value, jarque.bera.test(data.jb)$p.value),
                                  p.method)[3]
    result <- matrix(0, 4, 1)
    colnames(result) <- c("jb")
    rownames(result) <- c("lambda.hat", "sw.pvalue", "sf.pvalue",
                          "jb.pvalue")
    result[, 1] <- c(jblam, sw.pvalue_data.jb, sf.pvalue_data.jb,
                     jb.pvalue_data.jb)
    out <- list()
    out$title <- "Implementation of Box-Cox Power Transformation"
    out$method = "Jarque-Bera"
    out$date <- date()
    out$result <- result
    out

}

#' Title
#'
#' @param bas matrix to transform ( each columns is transform using a boxCox transformation)
#' and then it is scaled
#' @return the transformed matrix ( after a scaling step)
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rpois(100,3)
#' z <- rnorm(100,10,2)
#' z[z<0] <- 0
#' bas <- cbind(x = exp(x), y = exp(y), z = z)
#' bxccx_transf(bas)
bxccx_transf<-function(bas){
  bas <-bas+1
 if (is.null(colnames(bas))) colnames(bas) <- paste("V",1:ncol(bas))
  lambda<-function(i){
    tac<-boxcoxnc2(data =bas[,i],lam=seq(-2,2,0.5), plotit=F)
    return(tac$result["lambda.hat",] )
  }
  Lambda<-sapply(colnames(bas),lambda)
  print(Lambda)


  transf<-function(i){
    l<-Lambda[i]
    if(l !=0){
      atransf<-bas[,i]
      return((atransf^l-1)/l)

    }
    else {
      return(log(bas[,i]))

    }
  }

  bxcx<-sapply(colnames(bas),transf)
  return(scale(bxcx))
}

## ---- echo = FALSE-------------------------------------------------------
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(pdSpecEst)
## Generate example time series and periodogram
set.seed(123)
d <- 3
n <- 2^9 
example <- rExamples(2 * n, example = "bumps")
freq <- example$freq
str(example)

## ---- echo = F, out.width='75%', fig.width = 8, fig.height = 3, dpi=100, fig.align = "center"----
## Plot time series observations
colnames(example$ts) <- c("X1", "X2", "X3")
plot.ts(Re(example$ts), ylab = "Time", main = "Time series components")

## ---- echo = F, out.width='75%', fig.width = 8, fig.height = 5, dpi=100, fig.align = "center", fig.cap = "Figure 1: Generating (3 x 3)-dimensional spectral matrix"----
## Plot spectral matrix
plotspec <- function(i, data, data1 = NULL){
  
  ylim <- range(Re(data$f), Im(data$f))
  if(!is.null(data1)){
    ylim <- range(range(Re(data$f), Im(data$f)), range(Re(data1$f), Im(data1$f)))
  }
  if(i[1] == i[2]){
    plot(range(data$freq), range(ylim), type = "n", main = paste0("Auto-spectrum (", i[1], ",", i[1], ")"),
         xlab = "", yaxt = "n", ylab = "", ylim = ylim, mgp = c(3,0.5,0))
    abline(h = 0, lty = 3)
    if(length(dim(data$f)) == 3){
      if(!is.null(data1)){
        lines(data1$freq, Re(data1$f[i[1], i[1],]))
      }
      lines(data$freq, Re(data$f[i[1], i[1],]), lty = ifelse(is.null(data1), 1, 2))
    } else if(length(dim(data$f)) == 4){
      for(i1 in 1:dim(data$f)[4]){
        for(i2 in 1:dim(data$f)[4]){
          if(!is.null(data1)){
            lines(data1$freq, Re(data1$f[i[1], i[1], , i1]), col = i1)
          }
          lines(data$freq, Re(data$f[i[1], i[1], , i1]), col = i1, lty = ifelse(is.null(data1), 1, 2))
        }
      }
    }
    title(xlab = expression(paste("Frequency (", omega, ")")), line = 2)
    
  } else{
    plot(range(data$freq), range(ylim), type = "n", main = paste0("Cross-spectrum (", i[1], ",", i[2], ")"),
         xlab = "", yaxt = "n", xaxt = "n", ylab = "", ylim = ylim, mgp = c(3,0.5,0))
    abline(h = 0, lty = 3)
    title(xlab = expression(paste("Frequency (", omega, ")")), line = 0.5)
    if(!is.null(data1)){
      lines(data1$freq, Re(data1$f[i[1], i[2], ]), col = 1)
      lines(data1$freq, Im(data1$f[i[1], i[2], ]), col = 2)
    }
    lines(data$freq, Im(data$f[i[1], i[2], ]), lty = ifelse(is.null(data1), 1, 2), col = 2)
    lines(data$freq, Re(data$f[i[1], i[2], ]), lty = ifelse(is.null(data1), 1, 2), col = 1)
    legend("topright", legend = c("Real", "Imag."), lty = c(1,1), col = c(1,2), bty = "n", cex = 1, y.intersp = 1)
    
  }
}
par(mfrow=c(d, d), mar = c(3.5,1,1.5,1))
invisible(apply(expand.grid(1:d, 1:d), 1, function(i) plotspec(i, data = example)))
invisible(dev.off())

## ------------------------------------------------------------------------
wt.f <- WavTransf1D(example$f, periodic = T)

## ---- echo = F, out.width='75%', fig.width = 8, fig.height = 3, dpi=100, fig.align = "center"----
## Plot wavelet coefficients
plotCoeff <- function(D, title){
  
   if (!requireNamespace("ggplot2", quietly = T) | 
      !requireNamespace("viridis", quietly = T) | 
      !requireNamespace("ggthemes", quietly = T) | 
      !requireNamespace("reshape2", quietly = T)) {
    cat("Packages 'ggplot2', 'viridis', 'ggthemes' and 'reshape2' needed for this function to work. Please install missing packages.")
  } else{
    
    L_b <- (dim(D[[1]])[3] - 1) / 2
    J <- length(D)
    D <- lapply(1:J, function(j) D[[j]][, , L_b + 1:2^(j - 1), drop = F])
    norms <- lapply(1:J, function(j) apply(D[[j]], 3, function(D) pdSpecEst:::NormF(D)))
    longData <- reshape2::melt(sapply(1:J, function(j) rep(norms[[j]], each = 2^(J - j))))
    
    gg <- ggplot2::ggplot(longData, ggplot2::aes(x = longData$Var1, y = longData$Var2, fill = longData$value)) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis() +
      ggplot2::scale_x_continuous(breaks = NULL, labels = NULL, expand=c(0,0)) +
      ggplot2::scale_y_reverse(breaks = 1:J, labels=0:(J-1), expand = c(0, 0)) +
      ggplot2::ggtitle(title) +
      ggplot2::labs(x = "Location", y = "scale") +
      ggthemes::theme_tufte(base_family="sans") +
      ggplot2::theme(plot.title=ggplot2::element_text(size=14, hjust=0), axis.text = ggplot2::element_text(size=12), axis.title=ggplot2::element_text(size=14), legend.key.width=grid::unit(0.4, "cm"), legend.title=ggplot2::element_blank(), legend.key.height = grid::unit(1.25, "cm"), legend.text = ggplot2::element_text(size=12))
    
    print(gg)
  }
}
invisible(plotCoeff(wt.f$D, title = "Frobenius norm of target wavelet coefficients"))

## ------------------------------------------------------------------------
f.hat <- pdSpecEst1D(example$per)
str(f.hat)

## ---- echo=FALSE, out.width='75%', fig.width = 8, fig.height = 3, dpi=100, fig.align = "center"----
wt.per <- WavTransf1D(example$per, periodic = T)
invisible(plotCoeff(wt.per$D, title = "Frobenius norm of noisy periodogram wavelet coefficients"))
invisible(plotCoeff(f.hat$D, title = "Frobenius norm of denoised wavelet coefficients"))

## ---- echo=FALSE, out.width='75%', fig.width = 8, fig.height = 5, dpi=100, fig.align = "center", fig.cap = "Figure 2: target spectral matrix (dashed lines) and estimated spectral matrix (continuous lines)."----
## Plot estimated spectral matrix
par(mfrow=c(d, d), mar = c(3.5,1,1.5,1))
invisible(apply(expand.grid(1:d, 1:d), 1, function(i) plotspec(i, data = example, 
                                                               data1 = list(freq = example$freq, f = f.hat$f))))
invisible(dev.off())

## ---- eval = F-----------------------------------------------------------
#  ## Not run
#  boot.ci <- pdConfInt1D(f.hat$f, alpha = c(0.1, 0.05, 0.01), ci.region = c(0.5, 1),
#                         boot.samples = 1E3, f.0 = example$f)

## ---- echo = F-----------------------------------------------------------
cat("The depth-based confidence balls:")
pdSpecEst:::vignette1.data$depth.CI

cat("Verification that the target spectral matrix is covered:")
pdSpecEst:::vignette1.data$cover.f

cat("Depth of the target spectrum w.r.t. cloud of bootstrap spectral estimates")
pdSpecEst:::vignette1.data$depth.f


## ------------------------------------------------------------------------
## Generate example pseudo time-varying periodogram observations
set.seed(17)
d <- 2
n <- c(2^7, 2^7) 
example <- rExamples2D(n, d, example = "smiley", snr = 0.5)
tf.grid <- example$tf.grid ## time-frequency grid
str(example)

## ---- echo = F, out.width='75%', fig.width = 8, fig.height = 6, dpi=100, fig.align = "center", fig.cap = "Figure 3: Matrix-log's of target time-varying spectrum"----
## Plot 2D spectral matrices
plotspec2D <- function(P, lim = T, Log = F){
  
  if (!requireNamespace("ggplot2", quietly = T) | 
      !requireNamespace("viridis", quietly = T) | 
      !requireNamespace("ggthemes", quietly = T) | 
      !requireNamespace("reshape2", quietly = T)) {
    cat("Packages 'ggplot2', 'viridis', 'ggthemes' and 'reshape2' needed for this function to work. Please install missing packages.")
  } else{
    
    d <- dim(P)[1]
    x_n <- min(dim(P)[3], 32)
    y_n <- min(dim(P)[4], 32)
    P <- P[,,as.integer(seq(from=1,to=dim(P)[3],len=x_n)),as.integer(seq(from=1,to=dim(P)[4],len=y_n))]
    grid_n <- expand.grid(1:x_n, 1:y_n)
    if(Log){
      P <- array(apply(P, c(3,4), function(P) Logm(diag(d), P)), dim = c(d, d, x_n, y_n))
    }
    ylim <- range(mapply(function(i1, i2) range(Re(P[,,i1,i2]), Im(P[,,i1,i2])), grid_n$Var1, grid_n$Var2))
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2*d, d)))
    define_region <- function(row, col){
      grid::viewport(layout.pos.row = row, layout.pos.col = col)
    }
    
    marg <- 1/(2*d)
    
    for(d1 in 1:d){
      for(d2 in 1:d){
        if(d1 == d2){
          data <- Re(P[d1,d1,,])
          longdata <- reshape2::melt(data)
          longdata$Var1 <- rep(seq(from=1/x_n,to=1,len=x_n), times=y_n)
          longdata$Var2 <- rep(seq(from=0.5/y_n,to=0.5,len=y_n), each=x_n)
          
          gg <- ggplot2::ggplot(longdata, ggplot2::aes(x = longdata$Var1, y = longdata$Var2, fill = longdata$value)) +
            ggplot2::geom_tile() +
            viridis::scale_fill_viridis(name = "", limits = (if(lim) ylim else NULL)) +
            ggplot2::labs(x = "Time (t)", y = expression(paste("Freq. (", omega, ")", sep ="")),
                          title = paste0("Auto-spec. (", d1, ",", d1, ")")) +
            ggthemes::theme_tufte(base_family="sans") +
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0, size = 10), 
                           axis.ticks=ggplot2::element_blank(), axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = -2)), axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = -6)), axis.text = ggplot2::element_blank(), legend.key.width=grid::unit(0.3, "cm"), legend.title=ggplot2::element_blank(), legend.key.height = grid::unit(1.1, "cm"), legend.text = ggplot2::element_text(size=8), legend.position = c(1.05, 0.5), 
                           plot.margin = grid::unit(c(5.5,30.5,5.5,5.5), "points"))
          
          print(gg, vp = define_region((2*d1-1):(2*d1), d1))
        } else{
          if(d2 > d1){
            data1 <- Re(P[d1,d2,,])
            longdata1 <- reshape2::melt(data1)
            longdata1$Var1 <- rep(seq(from=1/x_n,to=1,len=x_n), times=y_n)
            longdata1$Var2 <- rep(seq(from=pi/y_n,to=pi,len=y_n), each=x_n)
            
            gg1 <- ggplot2::ggplot(longdata1, ggplot2::aes(x = longdata1$Var1, y = longdata1$Var2, fill = longdata1$value)) +
              ggplot2::geom_tile() +
              viridis::scale_fill_viridis(name = "", limits = (if(lim) ylim else NULL)) +
              ggplot2::labs(x = "Time (t)", y = expression(paste("Freq. (", omega, ")", sep = "")), title = paste0("Real cross. (", d1, ",", d2, ")")) +
              ggthemes::theme_tufte(base_family="sans") +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0, size = 10), axis.ticks=ggplot2::element_blank(), axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = -2)), legend.position = "none",
                             axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = -6)), axis.text=ggplot2::element_blank())
            
            print(gg1, vp = define_region(2*(d1-1)+1,d2))
            
            data2 <- Im(P[d1,d2,,])
            longdata2 <- reshape2::melt(data2)
            longdata2$Var1 <- rep(seq(from=1/x_n,to=1,len=x_n), times=y_n)
            longdata2$Var2 <- rep(seq(from=pi/y_n,to=pi,len=y_n), each=x_n)
            
            gg2 <- ggplot2::ggplot(longdata2, ggplot2::aes(x = longdata2$Var1, y = longdata2$Var2, fill = longdata2$value)) +
              ggplot2::geom_tile() +
              viridis::scale_fill_viridis(name = "", limits = (if(lim) ylim else NULL)) +
              ggplot2::labs(x = "Time (t)", y = expression(paste("Freq. (", omega, ")", sep = "")), title = paste0("Imag. cross. (", d1, ",", d2, ")")) +
              ggthemes::theme_tufte(base_family="sans") +
              ggplot2::theme(plot.title=ggplot2::element_text(hjust=0, size = 10), legend.position = "none",
                             axis.ticks=ggplot2::element_blank(), axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = -2)),
                             axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = -6)), axis.text=ggplot2::element_blank())
            
            print(gg2, vp = define_region(2*d1,d2))
          }
        }
      }
    }
  }
}

invisible(plotspec2D(example$f, lim = T, Log = T))

## ---- echo = F, out.width='75%', fig.width = 8, fig.height = 6, dpi=100, fig.align = "center", fig.cap = "Figure 4: Matrix-log's of pseudo-periodogram matrix"----
invisible(plotspec2D(example$per, lim = T, Log = T))

## ------------------------------------------------------------------------
f.hat <- pdSpecEst2D(example$per, order = c(1, 1), progress = F)
str(f.hat)

## ---- echo = F, out.width='75%', fig.width = 8, fig.height = 6, dpi=100, fig.align = "center", fig.cap = "Figure 5: Matrix-log's of estimated time-varying spectrum"----
invisible(plotspec2D(f.hat$f, lim = T, Log = T))

## ------------------------------------------------------------------------
## Fix parameters
set.seed(123)
Phi1 <- array(c(0.5, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
Phi2 <- array(c(0.7, 0, 0, 0.4, rep(0, 4)), dim = c(2, 2, 2))
Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)

## Generate periodogram data for 10 subjects
pgram <- function(Phi) pdPgram(rARMA(2^10, 2, Phi, Theta, Sigma)$X)$P
P <- array(c(replicate(5, pgram(Phi1)), replicate(5, pgram(Phi2))), dim=c(2,2,2^9,10))

pdSpecClust1D(P, K = 2, metric = "logEuclidean")$cl.prob


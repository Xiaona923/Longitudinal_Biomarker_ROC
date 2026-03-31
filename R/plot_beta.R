plot_beta <- function(beta.t, beta.t.CI, visit.time, my.title = NULL){
  #the sequence of beta and time.mat need to be consistent
  x <- visit.time
  #prepare data to plot
  min.y <- floor(min(beta.t.CI[, 1])) 
  max.y <- ceiling(max(beta.t.CI[, 2]))
  
  plot.dat <- unique(data.frame(x, beta.t, beta.t.CI))
  reorder <- order(plot.dat$x)
  plot(x = plot.dat$x[reorder], y = plot.dat$beta.t[reorder], 
       ylim = c(min.y, max.y), 
       type = "l",
       lwd = 1,
       main = my.title,
       font.main = 2,
       yaxt='n',xaxt='n',
       xlab = "s",
       ylab = ifelse(grepl("gamma",my.title), expression(paste(gamma[tau], " (s)")), expression(paste(beta[tau], " (s)"))))
  axis(1, cex.axis=1, tck=-0.02, lwd = 0.8)
  axis(2, cex.axis=1, tck=-0.02, lwd = 0.8)
  
  lines(x = plot.dat$x[reorder], plot.dat[,3][reorder], col = "black", lwd = 0.8, lty = 3)
  lines(x = plot.dat$x[reorder], plot.dat[,4][reorder], col = "black", lwd = 0.8, lty = 3)
  abline(h = 0, col = "#4292C6", lty = 2, lwd = 1)
}

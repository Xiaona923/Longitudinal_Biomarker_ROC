plot_ROC <- function(data, my.add, my.col, my.lty, my.main =""){
  f <- approxfun(data$FalsePos, data$new_meas, method = "constant", f = 0)
  par(pty="s")
  curve(f(x), from = 0, to = 1, ylim = c(0, 1), lwd = 1.5, lty = my.lty,
        xlab = "",
        ylab = "", main = my.main,
        col = my.col,
        add = my.add,
        axes=FALSE,
        frame=TRUE)
  title(xlab = "1 - Specificity", line = 3)
  title(ylab = "Sensitivity", line = 3)
  axis(1, cex.axis=1, tck=-0.02, lwd = 0.8)
  axis(2, cex.axis=1, tck=-0.02, lwd = 0.8)
}

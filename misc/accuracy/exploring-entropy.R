library(aqp)

equal.prob.H <- function(n, b) {
  p <- rep(1, times = n) / n
  shannonEntropy(p, b = b)
}


x <- seq(0, 1, by = 0.01)
y <- 1 - x

p <- cbind(x, y)

H <- apply(p, 1, shannonEntropy)

s <- seq_along(x)

par(mar = c(3, 4.5, 1, 1), bg = 'black', fg = 'white')

plot(s, x, type = 'n', las = 1, axes = FALSE)
# grid()
axis(side = 2, col.axis = 'white', las = 1, cex.axis = 1, at = seq(0, 1, by = 0.1))
mtext(side = 2, text = 'Probability | Shannon H', line = 3, font = 2)

lines(s, x, lwd = 2, col = grey(0.33))
lines(s, y, lwd = 2)
lines(s, H, col = 2, lwd = 2)

abline(h = 0)

abline(h = equal.prob.H(n = 2, b = 2), lty = 2, col = 5)
abline(h = equal.prob.H(n = 2, b = 2) * 0.95, lty = 2, col = 4)
abline(h = equal.prob.H(n = 2, b = 2) * 0.90, lty = 2, col = 3)

idx <- which.min(H < equal.prob.H(n = 2, b = 2))
segments(x0 = 0, y0 = y[idx], x1 = s[idx], y1 = y[idx], col = 5)
segments(x0 = s[idx], y0 = y[idx], x1 = s[idx], y1 = H[idx], col = 5)
text(x = 0, y = y[idx], labels = sprintf('100%%: %s', y[idx]), cex = 0.85, col = 5, adj = c(0, -0.8), font = 2)

idx <- which.min(H < equal.prob.H(n = 2, b = 2) * 0.95)
segments(x0 = 0, y0 = y[idx], x1 = s[idx], y1 = y[idx], col = 4)
segments(x0 = s[idx], y0 = y[idx], x1 = s[idx], y1 = H[idx], col = 4)
text(x = 0, y = y[idx], labels = sprintf('95%%: %s', y[idx]), cex = 0.85, col = 4, adj = c(0, -0.8), font = 2)

idx <- which.min(H < equal.prob.H(n = 2, b = 2) * 0.90)
segments(x0 = 0, y0 = y[idx], x1 = s[idx], y1 = y[idx], col = 3)
segments(x0 = s[idx], y0 = y[idx], x1 = s[idx], y1 = H[idx], col = 3)
text(x = 0, y = y[idx], labels = sprintf('90%%: %s', y[idx]), cex = 0.85, col = 3, adj = c(0, -0.8), font = 2)



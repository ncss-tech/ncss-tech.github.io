# normalized Shannon entropy (Kempen et al, 2009)
shannon.H <- function(i, b) {
  res <- -1 * sum(i * log(i, base=b), na.rm=TRUE)
  return(res)
}

# confusion index (Burrough et al., 1997)
confusion.index <- function(i) {
  i <- sort(i, decreasing = TRUE)
  res <- 1 - (i[1] - i[2])
  return(res)
}

# x.i: data.frame, rows are predictions/observations, columns contain classes
# class.labels: vector of class labels, corrosponding to column names in x.i
# actual: name of column containing the observed class
brierScore <- function(x.i, class.labels, actual='actual') {
  # save the gen hz probabilities into new df
  x.pr <- x.i[, class.labels, drop=FALSE]
  # init new matrix to store most-likely gen hz class
  m <- matrix(0, ncol=ncol(x.pr), nrow=nrow(x.pr))
  # same structure as x.pr
  dimnames(m)[[2]] <- names(x.pr)
  # set appropriate genhz to 1
  for(i in 1:nrow(x.i)) {
    ml.hz.i <- x.i[[actual]][i]
    m[i, ml.hz.i] <- 1
  }
  # compute brier score
  bs <- sum((x.pr - m)^2, na.rm=TRUE) / nrow(x.pr)
  return(bs)
}




# simulate predicted class probabilities via draws from dirichlet 
# (n pixel * k classes)
# n: number of simulated "pixels" 
# alpha: dirichlet hyperparameter (class likelihood)
simulatePredictions <- function(n, alpha) {
  # number of classes
  k <- length(alpha)
  class.labels <- toupper(letters[1:k])
  
  # generate simulated probabilities
  x <- sample_dirichlet(n, alpha = alpha)
  x <- t(x)
  
  # add class labels and id
  d <- as.data.frame(x)
  names(d) <- class.labels
  d$id <- rownames(d)
  
  # simulate actual classes using predicted probabilities
  d <- ddply(d, 'id', function(i) {
    i$actual <- sample(class.labels, size = 1, prob = i[, class.labels])
    return(i)
  })
  
  # reshape for plotting
  m <- melt(d, id.vars = 'id', measure.vars = class.labels)
  
  # compute H and CI for each simulated case
  H <- apply(x, 1, shannon.H, b=k)
  CI <- apply(x, 1, confusion.index)
  
  # reshape for plotting
  z <- data.frame(Shannon.H=H, CI=CI)
  z.long <- make.groups(Shannon.H=z$Shannon.H, CI=z$CI)
  
  return(list(predictions=d, predictions.long=m, stats=z, stats.long=z.long, classes=class.labels))
}

# get an example of predictions and associated uncertainty stats
extractExample <- function(x, n=1) {
  p <- x$predictions[1:n, ]
  p$id <- NULL
  stats <- x$stats[1:n, ]
  d <- data.frame(p, stats)
  return(d)
}


# generate some performance metrics
# x: results from simulatePredictions()
performance <- function(x, w=NULL) {
  # predictions
  p <- x$predictions
  
  # default weights matrix
  if(is.null(w)) {
    w <- outer(1:length(x$classes), 1:length(x$classes))
    w[] <- 0
    diag(w) <- 1
    dimnames(w) <- list(x$classes, x$classes)
  }
  
  
  # confusion matrix
  tab <- crossTabProbs(x)
  
  # simple brier score
  bs <- brierScore(p, x$classes)
  
  ## TODO: decide on use of priors
  # weighted tau
  tau.res <-tauW(tab, w)
  
  
  res <- data.frame(brier.score=bs, tau=tau.res$tau, tau.w=tau.res$tau.w, PCC=tau.res$overall.naive)
  return(res)
}

# generate confusion matrix from simulatePredictions()
# matrix is based on this most likely class vs. actual class
# x: output from simulatePredictions()
crossTabProbs <- function(x) {
  p <- x$predictions
  # most likely class lables
  max.pr.class <- x$classes[apply(p[, x$classes], 1, which.max)]
  # upgrade to factors for cross tabulation
  max.pr.class <- factor(max.pr.class, levels = x$classes)
  p$actual <- factor(p$actual, levels = x$classes)
  # confusion matrix, rows are predictions, columns are actual
  tab <- table(predictions=max.pr.class, actual=p$actual)
  
  return(tab)
}

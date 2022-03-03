#' @export

ebcd.update1 <- function(ebcd.fit){

  ebcd.fit <- ebcd.update.L(ebcd.fit)
  ebcd.fit <- ebcd.update.F(ebcd.fit)
  ebcd.fit <- ebcd.update.tau(ebcd.fit)
  ebcd.fit <- ebcd.update.obj(ebcd.fit)

  return(ebcd.fit)
}

# individual update functions --------------------

ebcd.update.L <- function(ef){

  temp.svd <- svd(ef$Y %*% ef$EF[[2]])
  ef$EF[[1]] <- tcrossprod(temp.svd$u, temp.svd$v)
  ef$EF2[[1]] <- ef$EF[[1]]^2
  ef$KL[[1]] <- rep(0, times=ncol(ef$EF[[1]]))

  return(ef)
}

ebcd.update.F <- function(ef){

  for (k in 1:ncol(ef$EF[[2]]) ){

    x <- crossprod(ef$Y, ef$EF[[1]][, k])
    s <- rep(1/sqrt(ef$tau), times=length(x))

    temp.ebnm.fn <- ef$ebnm.fn[[k]][[2]]
    temp.ebnm <- temp.ebnm.fn(x = x, s = s)

    ef$g[[k]][[2]] <- temp.ebnm$fitted_g
    ef$EF[[2]][, k] <- temp.ebnm$posterior$mean
    ef$EF2[[2]][, k] <- temp.ebnm$posterior$sd^2 + temp.ebnm$posterior$mean^2
    ef$KL[[2]][k] <- temp.ebnm$log_likelihood +
      -flashier_normal.means.loglik(x, s, ef$EF[[2]][, k], ef$EF2[[2]][, k])
  }
  return(ef)
}

ebcd.update.tau <- function(ef){

  ef$tau <- ef$n.nonmissing/(ef$Y2 - 2 * sum((ef$Y %*% ef$EF[[2]]) * ef$EF[[1]] ) +
                               sum(ef$EF2[[2]]) )
  ef$est.tau <- ef$tau

  return(ef)
}

ebcd.update.obj <- function(ef){

  ef$obj <- (-ef$n.nonmissing)/2 * (log(2*pi) - log(ef$tau)) +
    -ef$tau/2 * (ef$Y2 - 2 * sum((ef$Y %*% ef$EF[[2]]) * ef$EF[[1]] ) +
                   sum(ef$EF2[[2]]) ) + sum(ef$KL[[2]])
  return(ef)
}


#' @export

ebcd.block <- function(ebcd.obj,
                       maxiter = 5000,
                       tol = flashier_set.default.tol(ebcd.obj$flash.fit)
                       ){

  ebcd.fit <- ebcd.obj$flash.fit

  verbose.lvl <- ebcd.fit$verbose.lvl
  announce.block(verbose.lvl, n.factors = ncol(ebcd.fit$EF[[2]]), tol)

  iter <- 0
  obj.old <- -Inf
  obj.diff <- Inf
  next.tol.target <- NULL
  while (iter < maxiter && abs(obj.diff) > tol) {

    iter <- iter + 1

    ebcd.fit <- ebcd.update1(ebcd.fit)
    obj.diff <- ebcd.fit$obj - obj.old
    obj.old <- ebcd.fit$obj

    # 'next.tol.target' construction/reporting is from flashier
    if (is.null(next.tol.target) && abs(obj.diff) > 0 && abs(obj.diff) < Inf) {
      # Set the first target.
      next.tol.target <- 10^floor(log10(abs(obj.diff)))
    } else if (!is.null(next.tol.target) && abs(obj.diff) < next.tol.target) {
      # Report progress and set the next target.
      flashier_report.backfit.progress(verbose.lvl, next.tol.target)
      next.tol.target <- next.tol.target / 10
    }

  }

  if (iter == maxiter) {
    flashier_report.maxiter.reached(verbose.lvl)
  }
  flashier_report.backfit.complete(verbose.lvl, ebcd.fit$obj)

  flashier_announce.wrapup(verbose.lvl)
  if (!is.null(ebcd.fit$U)) {
    ebcd.fit$Y        <- ebcd.fit$U %*% ebcd.fit$Y
    ebcd.fit$EF[[1]]  <- ebcd.fit$U %*% ebcd.fit$EF[[1]]
    ebcd.fit$EF2[[1]] <- ebcd.fit$EF[[1]]^2
  }
  ebcd.obj <- flashier_wrapup.flash(ebcd.fit, output.lvl = 3L)
  ebcd.obj$L.psd  <- NULL
  ebcd.obj$L.lfsr <- NULL
  ebcd.obj$L.ghat <- NULL

  flashier_report.completion(verbose.lvl)

  return(ebcd.obj)
}


#================

announce.block <- function(verbose.lvl, n.factors, tol) {
  if (verbose.lvl > 0)
    cat(paste0("Block-fitting ", n.factors, " factors (tolerance: ",
               formatC(tol, format = "e", digits = 2), ")...\n"))
}


#' @references https://stat.ethz.ch/pipermail/r-devel/2013-August/067210.html
#' fix to use internal functions of flashier
#'
`%:::%` = function(pkg, fun) get(fun, envir = asNamespace(pkg), inherits = FALSE)

flashier_announce.wrapup         <- 'flashier'%:::%'announce.wrapup'
flashier_normal.means.loglik     <- 'flashier'%:::%'normal.means.loglik'
flashier_report.backfit.complete <- 'flashier'%:::%'report.backfit.complete'
flashier_report.backfit.progress <- 'flashier'%:::%'report.backfit.progress'
flashier_report.completion       <- 'flashier'%:::%'report.completion'
flashier_report.maxiter.reached  <- 'flashier'%:::%'report.maxiter.reached'
flashier_set.default.tol         <- 'flashier'%:::%'set.default.tol'
flashier_wrapup.flash            <- 'flashier'%:::%'wrapup.flash'




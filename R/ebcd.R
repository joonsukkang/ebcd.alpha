#' @export

ebcd <- function(data,
                 init = 'deflation',
                 compact = TRUE,
                 S = NULL,
                 ebnm.fn = ebnm::ebnm_point_laplace,
                 greedy.Kmax = 10L,
                 sampleratio = 1,
                 verbose = 1L){

  if (any(is.na(data))) { stop('missing values not allowed') }
  if ( (compact==TRUE) & (nrow(data) < ncol(data)) ) {
    compact <- FALSE
    if(verbose>=1){
      cat(paste0("compactifying not used (nrow < ncol) \n"))
      }
    }
  greedy.Kmax <- min(greedy.Kmax, dim(data)-1) # to satisfy orthogonality of L

  ebcd.obj <- flashier::flash.init(data, S = S, var.type = 0)
  ebcd.obj$flash.fit$n.nonmissing <- sampleratio * ebcd.obj$flash.fit$n.nonmissing

  if (init == 'deflation'){
    ebcd.obj <- flashier::flash.add.greedy(ebcd.obj,
                                           Kmax = greedy.Kmax,
                                           ebnm.fn = c(ebnm_norm1, ebnm.fn),
                                           verbose = verbose)
  }else if (init == 'svd'){
    ebcd.obj <- flashier::flash.init.factors(ebcd.obj,
                                             init = svd(data,
                                                        nu = greedy.Kmax,
                                                        nv = greedy.Kmax))
  }else if (init == 'varimax'){
    data.svd <- svd(data, nu = greedy.Kmax, nv = greedy.Kmax)
    L <- data.svd$v %*% diag(x = data.svd$d[1:greedy.Kmax], nrow=greedy.Kmax)
    vm <- varimax(L)

    L <- vm$loadings
    Z <- data.svd$u %*% vm$rotmat
    list.ZL <- list(Z, L)

    ebcd.obj <- flashier::flash.init.factors(ebcd.obj,
                                             init = list.ZL)
  }

  ebcd.obj$flash.fit$verbose.lvl <- verbose

  if (compact) {
    svdY <- svd(ebcd.obj$flash.fit$Y) # Y \approx UC
    ebcd.obj$flash.fit$Y <- diag(svdY$d) %*% t(svdY$v)
    ebcd.obj$flash.fit$U <- svdY$u
    ebcd.obj$flash.fit$EF[[1]] <- crossprod(svdY$u, ebcd.obj$flash.fit$EF[[1]])
    ebcd.obj$flash.fit$EF2[[1]] <- ebcd.obj$flash.fit$EF[[1]]^2
  }


  ebcd.obj <- ebcd.scale(ebcd.obj)
  ebcd.obj <- ebcd.block(ebcd.obj)
  #ebcd.obj <- ebcd.nullcheck(ebcd.obj)

  return(ebcd.obj)
}


#=================
ebnm_norm1 <- function(x,
                       s = 1, mode = 0, scale = "estimate", g_init = NULL,
                       fix_g = FALSE, output = ebnm::output_default(), optmethod = NULL,
                       control = NULL
){
  ebnm.res <- ebnm::ebnm_normal(x) # temporary fix to fill out the 'fitted_g' slot
  x <- x/sqrt(sum(x^2))

  if('lfsr' %in% output){
    out.list <- list(posterior = list(lfsr = rep(0, length(x))))
  } else {
    out.list <- list(data = list(x = x),
                     posterior = list(mean = x,
                                      second_moment = x^2),
                     fitted_g = ebnm.res$fitted_g,
                     log_likelihood = flashier_normal.means.loglik(x, s, x, x^2),
                     #KL=0,
                     posterior_sampler = NULL)

  }
  return(out.list)
}

ebcd.scale <- function(ebcd.obj){

  ebcd.fit <- ebcd.obj$flash.fit
  size <- sqrt(colSums(ebcd.fit$EF[[1]]^2))

  for (k in 1:ncol(ebcd.fit$EF[[1]])) {
    ebcd.fit$EF[[1]][,k] <- ebcd.fit$EF[[1]][,k] / size[k]
    ebcd.fit$EF[[2]][,k] <- ebcd.fit$EF[[2]][,k] * size[k]
    ebcd.fit$EF2[[1]][,k] <- ebcd.fit$EF2[[1]][,k] / size[k]^2
    ebcd.fit$EF2[[2]][,k] <- ebcd.fit$EF2[[2]][,k] * size[k]^2
  }

  ebcd.obj <- flashier_wrapup.flash(ebcd.fit, output.lvl = 3L)
  ebcd.obj$L.psd  <- NULL
  ebcd.obj$L.lfsr <- NULL
  ebcd.obj$L.ghat <- NULL

  return(ebcd.obj)
}

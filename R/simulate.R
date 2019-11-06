#' Title
#'
#' @param J
#' @param n_j
#' @param a
#' @param b
#' @param c_p
#' @param theta_ab
#'
#' @return
#' @export
#'
#' @examples
simulate_mediation_data <- function(J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2){
  requireNamespace('lme4', quietly = TRUE)
  re_a = T
  re_b = T
  newdata <- data.frame(x = rnorm(J*n_j),
                        cov.m = rep(rnorm(J), each = n_j),
                        cov.y = rep(rnorm(J), each = n_j),
                        id = rep(factor(1:J), each = n_j))
  if(re_a){
    re_term.m <- '(1 + x | id)'
    if(re_b){
      re_term.y <- '(1 + m + x | id)'
    } else {
      re_term.y <- re_term.m
    }
  } else {
    re_term.m <- '(1 | id)'
    if(re_b){
      re_term.y <- '(1 + m | id)'
    } else {
      re_term.y <- re_term.m
    }
  }
  form.m <- as.formula(paste0('~ 1 + cov.m + x + ', re_term.m))
  params.m <- list(theta = c(`id.(Intercept)` = 1,
                             `id.x.(Intercept)` = .1,
                             `id.x` = .5),
                   beta = c(`(Intercept)` = 0,
                            cov.m = 0,
                            x = a),
                   sigma = 1)
  newdata$m <- simulate(form.m,
                        newdata = newdata,
                        newparams = params.m,
                        family = 'gaussian', re.form = NA)[[1]]
  form.y <- as.formula(paste0('~ 1 + cov.y + m + x + ', re_term.y))
  params.y <- list(theta = c(`id.(Intercept)` = 1,
                             `id.m.(Intercept)` = .1,
                             `id.x.(Intercept)` = .1,
                             `id.m` = .5,
                             `id.x.m` = theta_ab,
                             `id.x` = .5),
                   beta = c(`(Intercept)` = 0,
                            cov.y = 0,
                            x = c_p,
                            m = b),
                   sigma = 1)
  newdata$y <- simulate(form.y,
                        newdata = newdata,
                        newparams = params.y,
                        family = 'gaussian', re.form = NA)[[1]]
  return(newdata)
}
run_permutation_simulation <- function(nrep, mc.cores, J = 100, n_j = 4, a = 0, b = 0, c_p = 0, theta_ab = .2, optimizer = "bobyqa"){

}

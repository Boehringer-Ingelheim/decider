#' @rdname logit
#' @name logit
#' @aliases log-odds inv_logit
#' @title Logit and inverse logit function
#' @description Logit and inverse logit function that work for vectors,
#' matrices, and other types of arrays that consist of numerical entries.
#' @usage logit(mu)
#' @param mu An array with numerical entries between 0 and 1.
#' @param eta An array with arbitrary numerical entries.
#' @return Returns an array of either the logits of the entries of mu (function \code{logit()}) or the
#' inverse logits of the entries of eta (function \code{inv_logit()}).
#' @details The logit and inverse logit function are defined as
#' \deqn{logit(x)=log(x/(1-x)),}
#' \deqn{inv\_logit(y)=exp(y)/(1+exp(y)).}
#' In particular, one has \eqn{logit(inv\_logit(y))=y} and
#' \eqn{inv\_logit(logit(x))=x}.
#'
#' Both functions are extracted from a \code{\link[stats:family]{stats::family}}
#' object of family \code{binomial()}.
#' @examples
#' y <- logit(0.1)
#' inv_logit(y)
#' logit(c(0.1, 0.2, 0.3, 0.4, 0.5))
#' inv_logit(c(-2, -1, 0, 1, 2))
#' @export
logit <- binomial()$linkfun

#' @export
#' @usage inv_logit(eta)
#' @rdname logit
inv_logit <- binomial()$linkinv

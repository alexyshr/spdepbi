#' Title
#'
#' @param pkg package name
#' @param fun function name
#'
#' @return function or variable value
#' @export
#'
#' @examples
#' library(spdep)
#' 'spdep' %:::% '.spdepOptions'
`%:::%` = function(pkg, fun){
  get(fun, envir = asNamespace(pkg), inherits = FALSE)
}

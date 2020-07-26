#' Print message
#'
#' @param dbg boolean whether to print or not
#' @param mess message to print
#' @param ... additional parameters
#'
#' @return message

spit = function (dbg, mess, ...) {
  if (dbg) { print( sprintf(mess,...) ) }
}

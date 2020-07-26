#' Print message and variable
#'
#' @param dbg boolean whether to print or not
#' @param name message to print
#' @param var variable to print
#'
#' @return name and variable
#'
spat = function (dbg, name, var) {
if (dbg) {print(name); print(var)}
}

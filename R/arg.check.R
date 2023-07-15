#' Check for missing arguments in function call
#'
#' @description Check for missing arguments using function call and a provided
#'  vector with argument names to check
#' @param call match.call(). To get function call with all of the specified
#'  arguments and their full names.
#' @param arguments character. Arguments to be checked
#'
#' @return logical
#'
#' @author Neander Marcel Heming
#'
#' @examples
#' geop <- function(x, tree, ...){
#'                 f4 <- arg.check(match.call(),
#'                                 c("LR", "inv.R", "branch.length", "n.descen"))
#'                 f1 <- arg.check(match.call(),
#'                                 c("tree"))
#'                 c(f1, f4)
#'                 }
#' geop(1, 1)
#' geop(1)
#' geop(1, LR=1)
#'
#' @export
arg.check <- function(call,
                      arguments = c("LR", "inv.R", "branch.length", "n.descen", "tree")  ){
  fun <- match.fun(as.list(call)[[1]])
  defined <- methods::formalArgs(args(fun))
  passed <- names(as.list(call)[-1])

  totest <- arguments %in% defined
  absent <- !(arguments %in% passed) & totest

  setNames(sapply(absent, function(x)(ifelse(x, T, F))), arguments) #[totest]
}

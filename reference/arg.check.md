# Check for missing arguments in function call

Check for missing arguments using function call and a provided vector
with argument names to check

## Usage

``` r
arg.check(
  call,
  arguments = c("LR", "inv.R", "branch.length", "n.descen", "tree")
)
```

## Arguments

- call:

  match.call(). To get function call with all of the specified arguments
  and their full names.

- arguments:

  character. Arguments to be checked

## Value

logical

## Author

Neander Marcel Heming

## Examples

``` r
# \donttest{
geop <- function(x, tree, ...){
                f4 <- arg.check(match.call(),
                                c("LR", "inv.R",
                                "branch.length", "n.descen"))
                f1 <- arg.check(match.call(),
                                c("tree"))
                c(f1, f4)
                }
geop(1, 1)
#> Error in get(sub("phyloraster::|phyloraster::", "", as.character(as.list(call)[[1]]))[1],     mode = "function", envir = loadNamespace("phyloraster")): object 'geop' of mode 'function' was not found
geop(1)
#> Error in get(sub("phyloraster::|phyloraster::", "", as.character(as.list(call)[[1]]))[1],     mode = "function", envir = loadNamespace("phyloraster")): object 'geop' of mode 'function' was not found
geop(1, LR=1)
#> Error in get(sub("phyloraster::|phyloraster::", "", as.character(as.list(call)[[1]]))[1],     mode = "function", envir = loadNamespace("phyloraster")): object 'geop' of mode 'function' was not found
# }
```

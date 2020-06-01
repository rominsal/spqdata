#' A funcion for creating symbols
#'
#' This function obtains the m-surroundings by selecting the m-1 nearest neighbors of each observation, allowing for a degree of overlap of s.
#' @param k,m
#' @keywords symbols
#' @export
#' @examples
#' # Obtain symbols for k=2 classes and m-surroundings of size 5
#' symb25 <- cr_symb(2,5)
#' symb25$p_symb #Permutations symbols
#' symb25$c_symb #Combinations-totals symbols

cr_symb <- function(k, m) {
  p_symb <- gtools::permutations(k, m, repeats.allowed = TRUE)  #Symbols by permutation
  c_symb0 <- gtools::combinations(k, m, repeats.allowed = TRUE)  #Symbols by combination
  c_symb <- matrix(rep(0, nrow(c_symb0) * k), ncol = k)
  for (i in 1:nrow(c_symb0)) {
    for (j in 1:k) {
      c_symb[i, j] = sum(c_symb0[i, ] == j)
    }
  }
  symb <- list(p_symb, c_symb)
  names(symb)[1] <- "p_symb"
  names(symb)[2] <- "c_symb"
  return(symb)
}

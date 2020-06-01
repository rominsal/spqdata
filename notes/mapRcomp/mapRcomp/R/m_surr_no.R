#' A funcion for creating m-surroudings
#'
#' This function obtains the m-surroundings by selecting the m-1 nearest neighbors of each observation, allowing for a degree of overlap of s.
#' @param coords
#' @keywords m_surround
#' @export
#' @examples
#' # Load dataset
#' data("Newark")
#'
#' # Use false origin UTM coordinates
#' coords <- cbind(Newark$X_UTM2,Newark$Y_UTM2)
#'
#' # Obtain m-surroundings of size 3 (m=3), with degree of overlap of two (s=2)
#' mh32 <- m_surr_no(coords,3,2)
#'

m_surr_no <- function(coords, m, s) {
  np <- nrow(coords)
  nnlist <- matrix(0, np, m - 1)  #Matrix with list of nearest neighbors
  nn1 <- order(sqrt((coords[1, 1] - coords[, 1])^2 + (coords[1, 2] -
                                                        coords[, 2])^2))
  nnlist[1, ] <- nn1[2:m]
  ns <- trunc((np - m)/(m - s)) + 1
  list <- rep(0, ns)  #zeros(1,ns)
  list[1] = 1
  # nnlist[1,] <- nn[1,1:m-1]
  blacklist <- c(1, nnlist[1, 1:(m - (s + 1))])
  t <- 1
  for (v in 2:ns) {
    list[v] = nnlist[t, m - s]
    h <- list[v]
    a <- order(sqrt((coords[h, 1] - coords[, 1])^2 + (coords[h, 2] -
                                                        coords[, 2])^2))  #a=nn[h,]
    a <- a[2:np]
    k = 0
    for (i in 1:(np - 1)) {
      if (sum(blacklist == a[i]) == 0) {
        k = k + 1
        nnlist[h, k] = a[i]
        t = h
        if (k == m - 1) {
          break
        }
      }
    }
    blacklist <- c(blacklist, h, nnlist[h, 1:(m - (s + 1))])
  }
  nnlist <- cbind(1:np, nnlist)
  nnlist1 <- cbind(nnlist[which(!nnlist[, 2] == 0), ])
  return(nnlist1)
}

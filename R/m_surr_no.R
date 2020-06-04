#' A funcion for creating m-surroudings MEJORADA
#'
#' This function obtains the m-surroundings by selecting the m-1 nearest neighbors of each observation, allowing for a degree of overlap of s.
#' @param coords son coordenadas
#' @keywords m_surround
#' @export
#' @examples
#'
#' # Obtain m-surroundings of size 3 (m=5), with degree of overlap of two (s=2)
#' N <- 100
#' m = 5
#' s = 2
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' mh <- m_surr_no(x, m, s)
#' W <- matrix(0,ncol = N,nrow = N)
#' for (i in 1:dim(mh)[1]){
#'  W[mh[i,1],mh[i,2:4]]<-1
#' }
#' W <- mat2listw(W)
#' plot(W,x)
#'
#'
#'

m_surr_no <- function(x, m, s) {
  coords <- x
  np <- nrow(coords)
  nnlist <- matrix(0, np, m - 1)  #Matrix with list of nearest neighbors
  nn1 <- cbind(1:np,sqrt((coords[1,1] - coords[,1])^2 + (coords[1,2] - coords[,2])^2),coords)
  nn1 <- nn1[order(nn1[,2]),]
  nnlist[1, ] <- nn1[2:m,1]
  ns <- trunc((np - m)/(m - s)) + 1
  list <- rep(0, ns)  #zeros(1,ns)
  list[1] = 1
  blacklist <- c(1, nnlist[1, 1:(m - (s + 1))])
  t <- 1
  for (v in 2:ns) {
    list[v] = nnlist[t, m - s]
    h <- list[v]
    nn1 <- nn1[!nn1[,1] %in% blacklist,]
    nn1[,2] <-  sqrt((nn1[1,3] - nn1[,3])^2 + (nn1[1,4] - nn1[,4])^2)
    nn1 <- nn1[order(nn1[,2]),]
    nnlist[h,] <- nn1[2:(m),1]
    t = h
    blacklist <- c(blacklist, h, nnlist[h, 1:(m - (s + 1))])
  }
  nnlist <- cbind(1:np, nnlist)
  nnlist1 <- cbind(nnlist[which(!nnlist[,2] == 0),])
  return(nnlist1)
}

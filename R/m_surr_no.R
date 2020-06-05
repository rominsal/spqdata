#'
#' @title A funcion for creating m-surroudings (mejoradas en tiempo de cómputo)
#' @usage m_surr_no(x = x, m = m, s = s, control = NULL)
#' @param x son coordenadas
#' @param m amplitud m-historia
#' @param s solapamiento máximo entre cualesquiera dos m-historias
#' @param control Argumento opcional. Elimina aquellas m-historias cuyos elementos
#' no esten contenidos dentro del conjunto de los k=control vecinos más próximos.
#' @description This function obtains the m-surroundings by selecting the m-1 nearest neighbors of each observation, allowing for a degree of overlap of s.
#' @details Aquí Antonio escribe una linda historia
#' @return una matrix cuyas filas corresponden con las m-historias
#' @keywords m_surround
#' @export
#' @examples
#'
#' # Obtain m-surroundings of size 3 (m=5), with degree of overlap of two (s=2)
#' N <- 1000
#' m = 5
#' s = 2
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' mh <- m_surr_no(x = x, m = m, s = s)
#' W <- matrix(0,ncol = N,nrow = N)
#' for (i in 1:dim(mh)[1]){
#'  W[mh[i,1],mh[i,2:m]] <- 1
#' }
#' W <- mat2listw(W)
#' plot(W,x)
#'
#' mh <- m_surr_no(x = x, m = m, s = s, control = 20)
#' W <- matrix(0,ncol = N,nrow = N)
#' for (i in 1:dim(mh)[1]){
#'  W[mh[i,1],mh[i,2:m]] <- 1
#' }
#' W <- mat2listw(W)
#' plot(W,x)
#'

m_surr_no <- function(x = x, m = m, s = s, control = NULL) {
  if (s<1 || s >(m-1))
    stop("mínimo grado de solapamiento es 1 y menor que la amplitud de la m-historia")

  np <- nrow(x)
  nnlist <- matrix(0, np, m - 1)  # Matrix with list of nearest neighbors
  nn1 <- cbind(1:np,sqrt((x[1,1] - x[,1])^2 + (x[1,2] - x[,2])^2),x)
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
  # control debería ser una lista
  # Se eliminan aquellas m-historias que contengan vecinos que no estén
  # dentro de los k vecinos más próximos
  if (is.null(control)==FALSE){
    knn <- cbind(1:N,spdep::knearneigh(x, control)$nn)
    int <- numeric()
    for (i in 1:dim(nnlist1)[1]){
      int[i] <- length(intersect(nnlist1[i,],knn[nnlist1[i,1],]))
    }
    nnlist1 <- nnlist1[int==m,]
  }
  return(nnlist1)
}

#'
#' @title A funcion to avoid overlaps of order greater than s in m-surroundings
#' @usage m_surr_maxs(x = x, s = s, control = NULL)
#' @param x matrix of m-surroundings
#' @param s maximum overlap between any two m-surroundings
#' @param control Argumento opcional. Por definir
#'                VIP: ¿SE PUEDE PONER UNA DISTANCIA MÁXIMA PARA ELIMINAR
#'                SURROUDINGS QUE SUPEREN ESE UMBRAL?.
#' @description This function prunes the m-surroundings .
#' @details Aquí Antonio escribe una linda historia
#' @return matrix of m-surroundings excluding cases with more than s overlaps
#' @keywords m_surround
#' @export
#' @examples
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' lmh3 <- m_surr_cbl(nc, m = 3)
#' mh3s <- m_surr_maxs(lmh3$mh, s = 1)
#' head(lmh3$mh); head(mh3s)
#' dim(lmh3$mh); dim(mh3s)
m_surr_maxs <- function(x, s = 1, control = NULL) {
  mh <- x
  N <- nrow(mh)
  m <- ncol(mh)
  if(is.null(rownames(mh))) rownames(mh) <- 1:N
  mhs <- mh
  Ns <- dim(mhs)[1]
  i <- 1
  while (i < Ns) {
    idxsij <- NULL # Vector de índices con solapamientos para feature i
    for (j in (i+1):Ns) {
      matchij <- match(mhs[i,], mhs[j,])
      if ((any(!is.na(matchij)))) {
        sij <- sum(!is.na(matchij))
      } else sij <- 0
      if (sij > s) idxsij <- c(idxsij, j)
    }
    if (!is.null(idxsij)) {
      mhs <- mhs[-idxsij, ]
      Ns <- dim(mhs)[1]
    }
    i <- i + 1
  }
  return(mhs)
}

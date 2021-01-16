#'
#' @title A funcion to avoid overlaps of order greater than r in m-surroundings
#' @usage m_surr_maxr(x = x, r = r, control = NULL)
#' @param x matrix of m-surroundings
#' @param r maximum overlap between any two m-surroundings
#' @param control Argumento opcional. Por definir
#'                VIP: ¿SE PUEDE PONER UNA DISTANCIA MÁXIMA PARA ELIMINAR
#'                SURROUDINGS QUE SUPEREN ESE UMBRAL?.
#' @description This function prunes the m-surroundings .
#' @details Aquí Antonio escribe una linda historia
#' @return matrix of m-surroundings excluding cases with more than r overlaps
#' @keywords m_surround
#' @export
#' @examples
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' lmh3 <- m_surr_cbl(nc, m = 3)
#' mh3s <- m_surr_maxr(lmh3$mh, r = 1)
#' head(lmh3$mh); head(mh3s)
#' dim(lmh3$mh); dim(mh3s)
m_surr_maxr <- function(x, r = 1, control = NULL) {
  mh <- x
  N <- nrow(mh)
  m <- ncol(mh)
  if(is.null(rownames(mh))) rownames(mh) <- 1:N
  mhr <- mh
  Ns <- dim(mhr)[1]
  i <- 1
  while (i < Ns) {
    idxrij <- NULL # Vector de índices con solapamientos para feature i
    for (j in (i+1):Ns) {
      matchij <- match(mhr[i,], mhr[j,])
      if ((any(!is.na(matchij)))) {
        rij <- sum(!is.na(matchij))
      } else rij <- 0
      if (rij > r) idxrij <- c(idxrij, j)
    }
    if (!is.null(idxrij)) {
      mhr <- mhr[-idxrij, ]
      Ns <- dim(mhr)[1]
    }
    i <- i + 1
  }
  return(mhr)
}

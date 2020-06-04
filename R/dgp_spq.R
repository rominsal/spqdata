#' @title A funcion to generate qualitative process with spatial structure
#' @description This function genera un proceso espacial de datos cualitativos
#' @param x matrix of point coordinate
#' @param listw A \code{listw} object created for example by
#'   \code{\link[spdep]{nb2listw}} from \pkg{spatialreg} package; if
#'   \code{\link[spdep]{nb2listw}} not given, set to
#'   the same spatial weights as the \code{listw} argument. It can
#'   also be a spatial weighting matrix of order \emph{(NxN)} instead of
#'   a \code{listw} object. Default = \code{NULL}.
#' @param rho nivel de dependencia espacial (un valor entre 0 y 1)
#' @param p un vector indicando el porcentaje de cada clase. Su longitud debe coincidir con el número de clases. Su suma debe ser 1.
#' @return Devuelve un FACTOR codificado con los primeros números naturales
#' @details Aquí Antonio escribe una linda historia
#' @seealso
#' \code{\link{q_symb}}, \code{\link{m_surr_no}}
#' @keywords m_surround
#' @export
#' @examples
#' #
#' N <- 1000
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- spdep::nb2listw(knn2nb(knearneigh(cbind(cx,cy), k=4)))
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' QY <- dgp_spq(x = x, p = p, listw = listw, rho = rho)

dgp_spq <- function(x,p,listw = NULL, rho) {

  if (sum(p)!=1)  stop("ojo, las probabilidades deben sumar 1")

listw <- listw2mat(listw)
N <- dim(listw)[1]
cx <- x[,1]
cy <- x[,2]
y <- Matrix::solve(diag(N)-rho*listw)%*%rnorm(N,1)
Y <- cut(y,quantile(y,c(0,cumsum(p))),include.lowest=TRUE)
levels(Y) <- as.character(1:length(p))
return(Y)
}

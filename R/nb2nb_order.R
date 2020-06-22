#'
#' @title Una función para ordenar elementos de las m_i-historias
#'
#' @description Es una función auxiliar. En el caso de obtener la lista de vecinos de clase nb con \code{\link{poly2nb}}
#' es necesario reordenar los elementos en función de la distancia y/o del angulo
#' @param listw una lista de vecinos tipo nb.
#' @param sf el objeto sf utilizado para .
#' @usage nb2nb_order <- function(listw = listw, sf = sf)
#' @keywords spatial association, qualitative variable, runs test
#' @details Ordena los elementos de una lista nb. En primer lugar por distancia y en segundo lugar (a igualdad de distancias) por angulo.
#' @return devuelve un objeto de la clase nb (elementos ordenados).
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paez@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     A test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}
#' @export
#' @examples
#'
#' # With a sf object (irregular lattice)
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' listw <- spdep::poly2nb(as(nc,"Spatial"), queen = FALSE)
#' listw_order <- nb2nb_order(listw = listw, sf = nc)
#'
#' # With a sf object (regular lattice: hexagons)
#' sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
#' hexs <- st_make_grid(sfc, cellsize = 0.1, square = FALSE)
#' hexs.sf <- st_sf(hexs)
#' listw  <- poly2nb(as(hexs.sf, "Spatial"), queen = FALSE)
#' listw_order <- nb2nb_order(listw = listw, sf = hexs.sf)
#'
nb2nb_order <- function(listw = listw, sf = sf){

# n <- length(listw)
# co <- as.data.frame(sf::st_coordinates(st_centroid(sf)))
# co <- st_as_sf(co,coords=c("X","Y"), crs = 32630)
# dis <- 1/matrix(as.numeric(st_distance(co,co)),ncol=n,nrow=n)
# diag(dis) <- 0
# matw <- nb2mat(listw,style = 'B', zero.policy = TRUE)
# m <- rowSums(matw)
# matwdis <- dis*matw
# NB <- list()
# for (i in 1:n){
#   or <-  order(matwdis[i,], decreasing = TRUE)
#   NB[[i]]<- or[1:m[i]]
#   if (m[i]==0){NB[[i]]<- as.integer(0)}
# }
# class(NB)<- 'nb'
# }
  n <- length(listw)
  matw <- nb2mat(listw,style = 'B', zero.policy = TRUE)
  m <- rowSums(matw)
  co <- sf::st_coordinates(st_centroid(sf))
  NB <- list()
for (i in 1:dim(co)[1]){
  a <- co[listw[[i]],] - t(matrix(co[i,],nrow = 2,ncol = length(listw[[i]])))
  a <- round(cbind(sqrt(a[,1]^2+a[,2]^2),atan2(a[,1],a[,2])*180/pi),digits = 6)
  a <- order(a[,1],a[,2])
  # a <- order(a,decreasing = FALSE) # sort(a,index.return=TRUE)$ix
  NB[[i]] <- listw[[i]][a]
  if (m[i]==0){NB[[i]]<- as.integer(0)}
}
  class(NB)<- 'nb'
return <- NB
}


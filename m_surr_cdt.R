#'
#' @title A funcion for creating m-surroudings minimizing the distances
#' @usage m_surr_cdt(x = x, m = m, r = r, control = NULL)
#' @param x input sf object with points/multipolygons geometry or matrix
#'          of spatial coordinates
#' @param m length of m-surrounding
#' @param r maximum overlap between any two m-surroundings
#' @param control Argumento opcional. Por definir
#' @description This function obtains the m-surroundings by
#'              choosing the m-1 nearest centroids. If there are less
#'              than m-1 neighbors the corresponding m-surrounding
#'              is filled with NA's
#' @details Aquí Antonio escribe una linda historia
#' @return REPASAR A list with SOME matrices: an m-surrounding matrix, named ms,
#'        whose rows correspond to the m-surroundings and a distance matrix,
#'        named mdtms, whose rows correspond to distance between the elements
#'        of m-surroundings.....
#' @keywords m_surrounding
#' @export
#' @examples
#' N <- 1000
#' m = 5
#' r = 2
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' msurr_points <- m_surr_cdt(x = x, m = m, r = r,
#'                            control = list(dtmaxabs = 0.05))
#' ms <- msurr_points$ms
#' mdtms <- msurr_points$mdtms
#' dim(ms); dim(mdtms)
#' ms[1:10,]
#' mdtms[1:10,]
#' # Examples with multipolygons
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' plot(sf::st_geometry(nc))
#' lms3 <- m_surr_cdt(nc, m = 3, r = 2,
#'                    control = list(dtmaxpc = 0.1))
#' lms5 <- m_surr_cdt(nc, m = 5, r = 1,
#'                    control = list(dtmaxpc = 0.1))
#' # m-surroundings or m-histories
#' lms3$ms[1:3, ]; lms5$ms[1:3, ]
#' # Distances in m-surroundings
#' lms3$mdtms[1:3,]; lms5$mdtms[1:3,]
m_surr_cdt <- function(x, m, r = 1, control = list()) {
  con <- list(dtmaxabs = 0, dtmaxpc = 0)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  # Transform matrix coordinates into SpatialPoints class
  if (inherits(x, "matrix")) x <- sp::SpatialPoints(x)
  # Transform Spatial classes into sf class
  if (inherits(x, "Spatial")) x <- as(x, "sf")
  # Compute centroids/coordinates/distances from sf objects
  # VIP: NO ME QUEDA CLARO SI HACE FALTA ESTABLECER CONTROLES PARA IMPEDIR
  #      CENTROIDES DE GEOMETRÍAS. SUPONGO QUE NO...
  if (inherits(x, "sf")) {
    xct <- suppressWarnings(sf::st_centroid(sf::st_geometry(x)))
  } else stop("object must be either sf, sp or matrix class")
  N <- length(xct) # Number of features
  # Compute nearest m-1 neighbords
  lknn <- spdep::knearneigh(xct, k = m - 1)
  ms <- cbind(c(1:N), lknn$nn)
  rownames(ms) <- c(1:N)
  colnames(ms) <- c(1:m)
  # Compute distances in m-surroundings
  mdtfull <- sf::st_distance(xct)
  mdtms <- matrix(0, nrow = N, ncol = m)
  rownames(mdtms) <- c(1:N)
  colnames(mdtms) <- c(1:m)
  for (i in 1:N) {
    for (j in 2:ncol(ms)) {
      mdtms[i,j] <- mdtfull[i,ms[i,j]]
    }
  }
  # Prune the r overlappings in m-surroundings
  msr <- m_surr_maxr(ms, r = r)
  # Build the distance for msr matrix
  mdtmsr <- mdtms[rownames(msr), ]
  ms <- msr
  mdtms <- mdtmsr
  rownames(ms) <- rownames(mdtms) <- NULL
  if (is.null(con$dtmaxabs)) dtmaxabs <- 0 else dtmaxabs <- con$dtmaxabs
  if (is.null(con$dtmaxpc)) dtmaxpc <- 0 else dtmaxpc <- con$dtmaxpc
  lms <- prunemsdthr(dtmaxabs = dtmaxabs, dtmaxpc = dtmaxpc,
                     mdtfull = mdtfull, ms = ms,
                     mdtms = mdtms)
  lms$N <- N
  lms$m <- m
  lms$r <- r
  return(lms)
}

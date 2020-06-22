#'
#' @title A funcion for creating m-surroudings minimizing the distances
#' @usage m_surr_cdt(x = x, m = m, control = NULL)
#' @param x input sf object with points/multipolygons geometry or matrix
#'          of spatial coordinates
#' @param m length of m-surrounding
#' @param s maximum overlap between any two m-surroundings
#' @param control Argumento opcional. Por definir
#' @description This function obtains the m-surroundings by
#'              choosing the m-1 nearest centroids. If there are less
#'              than m-1 neighbors the corresponding m-surrounding
#'              is filled with NA's
#' @details Aquí Antonio escribe una linda historia
#' @return REPASAR A list with SOME matrices: an m-surrounding matrix, named mh,
#'        whose rows correspond to the m-surroundings and a distance matrix,
#'        named dtmh, whose rows correspond to distance between the elements
#'        of m-surroundings.....
#' @keywords m_surrounding
#' @export
#' @examples
#' N <- 1000
#' m = 5
#' s = 2
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
# msurr_points <- m_surr_cdt(x = x, m = m, s = s,
#                            control = list(dthr = 0.05))
#' mh <- msurr_points$mh
#' dtmh <- msurr_points$dtmh
#' dim(mh); dim(dtmh)
#' mh[1:10,]
#' dtmh[1:10,]
#' # Examples with multipolygons
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' plot(sf::st_geometry(nc))
#' lmh3 <- m_surr_cdt(nc, m = 3, s = 2,
#'                    control = list(dthrpc = 0.1))
#' lmh5 <- m_surr_cdt(nc, m = 5, s = 1,
#'                    control = list(dthrpc = 0.1))
#' # m-surroundings or m-histories
#' lmh3$mh[1:3, ]; lmh5$mh[1:3, ]
#' # Distances in m-surroundings
#' lmh3$dtmh[1:3,]; lmh5$dtmh[1:3,]
m_surr_cdt <- function(x, m, s = 1, control = list()) {
  con <- list(dthr = 0, dthrpc = 0)
  #OJO: TAMBIÉN PODEMOS METER EN EL CONTROL EL TIPO DE DISTANCIA, ¿NO?.
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
    xct <- suppressWarnings(sf::st_centroid(sf::st_geometry(x),
                           of_largest_polygon = TRUE))
  } else stop("object must be either sf, sp or matrix class")
  N <- length(xct) # Number of features
  # Compute nearest m-1 neighbords
  lknn <- spdep::knearneigh(xct, k = m - 1)
  mh <- cbind(c(1:N), lknn$nn)
  rownames(mh) <- c(1:N)
  colnames(mh) <- c(1:m)
  # Compute distances in m-surroundings
  dtfull <- sf::st_distance(xct)
  dtmh <- matrix(0, nrow = N, ncol = m)
  rownames(dtmh) <- c(1:N)
  colnames(dtmh) <- c(1:m)
  for (i in 1:N) {
    for (j in 2:ncol(mh)) {
      dtmh[i,j] <- dtfull[i,mh[i,j]]
    }
  }
  # Prune the s overlappings in m-surroundings
  mhs <- m_surr_maxs(mh, s = s)
  # Build the distance for mhs matrix
  dtmhs <- dtmh[rownames(mhs), ]
  mh <- mhs
  dtmh <- dtmhs
  rownames(mh) <- rownames(dtmh) <- NULL
  # Check for threshold in distances dthr = 0, dthrpc = 0
  if (con$dthr > 0 || con$dthrpc > 0) {
    if (con$dthrpc <= 1)
      dthrpc <- as.numeric(con$dthrpc*max(dtfull)) else
        stop("The value of dthrpc must be in [0,1]")
    if ((dthrpc > 0)  && (dthrpc < con$dthr))
      dthr <- dthrpc else if (con$dthr > 0) dthr <- con$dthr
      else dthr <- dthrpc
  }
  if (dthr > 0) {
    excl <- dtmh > dthr
    indexcl <- which(excl==TRUE, arr.ind = TRUE)
    rowexcl <- unique(indexcl[,"row"])
    cat("\n Threshold distance: ", dthr)
    if (length(rowexcl) > 0) {
      cat("\n Number of m-surroundings excluded for exceeding
        the threshold distance: ",length(rowexcl),"\n")
      mh <- mh[-rowexcl,]
      dtmh <- dtmh[-rowexcl,]
    } else cat("\n None m-surrounding excluded for exceeding
        the threshold distance","\n")
  }
  lmh <- list(mh = mh, dtmh = dtmh)
  return(lmh)
}

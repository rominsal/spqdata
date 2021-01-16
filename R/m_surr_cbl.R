#'
#' @title A funcion for creating m-surroudings from polygonal
#'        features maximizing the common borderline
#' @usage m_surr_cbl(x = x, m = m, control = NULL)
#' @param x input sf object with polygons or multipolygons features
#' @param m length of m-surrounding
#' @param r maximum overlap between any two m-surroundings
#' @param control Argumento opcional. Por definir
#' @description This function obtains the m-surroundings by
#'              choosing the m-1 neighboring polygons. If there are more than
#'              (m-1) first order neighbors the algorithm chooses the polygons with the
#'              longest borderlines.
#'              **VIP** If there are less than (m-1) first order neighbors the
#'              corresponding m-surrounding is filled with second order
#'              neighbors in increasing order between centroid distances.
#' @details Aquí Antonio escribe una linda historia
#' @return REPASAR A list with four matrices: an m-surrounding matrix, named ms,
#'         whose rows correspond to the m-surroundings; a matrix, named msr,
#'         containing the m-surroundings with maximum overlapping of r; a matrix including
#'         the length of borderlines between m-surroundings and a distance matrix,
#'         named mdtms, whose rows correspond to distance between the elements
#'         of m-surroundings.
#' @keywords m_surrounding
#' @export
#' @examples
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' plot(sf::st_geometry(nc))
#' # VIP: THE THRESHOLD OF THE DISTANCE NEED TO BE VERY SMALL TO
#' # FIND M-SURROUNDINGS EXCEDING THE THRESHOLD.
#' lms3 <- m_surr_cbl(nc, m = 3, r = 1, control = list(dtmaxpc = 0.05))
#' lms5 <- m_surr_cbl(nc, m = 5, r = 2, control = list(dtmaxpc = 0.1))
#' dim(lms3$ms); dim(lms5$ms)
#' # m-surroundings or m-histories
#'  lms3$ms[1:5, ]; lms5$ms[1:5, ]
#' # length between borderlines in m-surroundings
#' lms3$mscbl[1:5, ]; lms5$mscbl[1:3, ]
#' # Distances in m-surroundings
#' lms3$mdtms[1:5,]; lms5$mdtms[1:5,]
#'
m_surr_cbl <- function(x, m, r = 1, control = list()) {
  con <- list(dtmaxabs = 0, dtmaxpc = 0)
  #OJO: TAMBIÉN PODEMOS METER EN EL CONTROL EL TIPO DE DISTANCIA, ¿NO?.
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  # Transform objects of class Spatial (sp package) into class sf
  if (inherits(x, "Spatial")) x <- as(x, "sf")
  # Check if the object is sf
  if (!inherits(x, "sf")) stop("Object must be of class sf")
  # Check if the geometry is Polygon/Multipolygon
  typegeo <- sf::st_geometry_type(x)
  if (any(!grepl("POLYGON", typegeo)))
    stop("All features must correspond to Polygon/Multipolygon geometries")
  nbx <- spdep::poly2nb(x) # Create nb object
  N <- dim(x)[1] # Number of features
  #W matrix of first order neighb.
  W <- spdep::nb2mat(nbx, style = "W", zero.policy = TRUE)
  W2 <- W %*% W
  # Compute centroids and distances between centroids
  xct <- suppressWarnings(sf::st_centroid(sf::st_geometry(x)))
  mdtfull <- sf::st_distance(xct)
  # Extract geometry
  xgeom <- sf::st_geometry(x)
  ms <- matrix(NA, N, m)
  mscbl <- matrix(NA, N, m)
  rownames(ms) <- rownames(mscbl) <- c(1:N)
  colnames(ms) <- colnames(mscbl) <- c(1:m)
  # m-surroundings sharing borderlines in order of length
  for (i in 1:N) {
    nbi <- c(i, nbx[[i]])
    mi <- length(nbi)
    nbicbl <- rep(0, mi)
    if (mi > 1) {
      for (j in 1:mi) {
        # Se incluye asimismo pero siempre estará en primer lugar
        intnbij <- sf::st_intersection(xgeom[[nbi[1]]],
                                       xgeom[[nbi[j]]])
        if (inherits(intnbij, c("MULTILINESTRING", "POLYGON"))) {
          nbicbl[j] <- sf::st_length(sf::st_multilinestring(intnbij))
        } else if (inherits(intnbij, "LINESTRING")) {
          nbicbl[j] <- sf::st_length(intnbij)
        } else {
          # cat("\n Intersection between regions: ",nbi[1],
          #     " and ",nbi[j]," not MULTILINESTRING/LINESTRING/POLYGON \n")
          # cat("Corresponding length of borderline set to 0 \n")
          nbicbl[j] <- 0
        }
      }
    }
    # Check the length of the first region...
    if (nbicbl[1] < max(nbicbl)) nbicbl[1] <- max(nbicbl) + 0.01
    ord_nbicbl <- order(nbicbl, decreasing = TRUE)
    nbiord <- nbi[ord_nbicbl]
    mscbli <- nbicbl[ord_nbicbl]
    # Add second order neighbors when there is no first order neighb.
    # Criterio: Add the nearest neighboors between centroids
    if (length(nbiord) < m) {
      W2i <- W2[i, ]
      idxordW2i <- order(W2i, decreasing = TRUE)
      idxordW2i <- idxordW2i[W2i[idxordW2i] > 1.0E-9]
      # Excluye los vecinos de segundo orden que ya están en lista.
      idxordW2i <- idxordW2i[!(idxordW2i %in% nbiord)]
      # Ordena los vecinos de segundo orden en orden creciente
      # a partir de distancia entre centroides
      dtcti <- mdtfull[i, idxordW2i]
      idxdtcti <- order(dtcti, decreasing = FALSE)
      idxordW2i <- idxordW2i[idxdtcti]
      nbiord <- c(nbiord, idxordW2i[1:(m-length(nbiord))])
      mscbli <- c(mscbli, rep(0, m-length(mscbli)))
    }
    ms[i, ] <- nbiord[1:m]
    mscbl[i, ] <- mscbli[1:m]
  }
  # Compute distances in m-surroundings
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
  # Build the distance and perimeter matrix for msr matrix
  msrcbl <- mscbl[rownames(msr),]
  dtmsr <- mdtms[rownames(msr), ]
  ms <- msr
  mdtms <- dtmsr
  mscbl <- msrcbl
  rownames(ms) <- rownames(mdtms) <- rownames(mscbl) <- NULL
  if (is.null(con$dtmaxabs)) dtmaxabs <- 0 else dtmaxabs <- con$dtmaxabs
  if (is.null(con$dtmaxpc)) dtmaxpc <- 0 else dtmaxpc <- con$dtmaxpc
  lms <- prunemsdthr(dtmaxabs = dtmaxabs, dtmaxpc = dtmaxpc,
                     mdtfull = mdtfull, ms = ms,
                     mdtms = mdtms, mscbl = mscbl)
  lms$N <- N
  lms$m <- m
  lms$r <- r
  return(lms)
}

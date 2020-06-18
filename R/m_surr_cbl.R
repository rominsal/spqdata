#'
#' @title A funcion for creating m-surroudings from polygonal
#'        features maximizing the common borderline
#' @usage m_surr_cbl(x = x, m = m, control = NULL)
#' @param x input sf object with polygons or multipolygons features
#' @param m length of m-surrounding
#' @param s maximum overlap between any two m-surroundings
#' @param control Argumento opcional. Por definir
#' @description This function obtains the m-surroundings by
#'              choosing the m-1 neighboring polygons. If there are more than
#'              (m-1) first order neighbors the algorithm chooses the polygons with the
#'              longest borderlines.
#'              **VIP** If there are less than (m-1) first order neighbors the
#'              corresponding m-surrounding is filled with second order
#'              neighbors in increasing order between centroid distances.
#' @details Aquí Antonio escribe una linda historia
#' @return REPASAR A list with four matrices: an m-surrounding matrix, named mh,
#'         whose rows correspond to the m-surroundings; a matrix, named mhs,
#'         containing the m-surroundings with maximum overlapping of s; a matrix including
#'         the length of borderlines between m-surroundings and a distance matrix,
#'         named dtmh, whose rows correspond to distance between the elements
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
#' lmh3 <- m_surr_cbl(nc, m = 3, s = 1, control = list(dthrpc = 0.05))
#' lmh5 <- m_surr_cbl(nc, m = 5, s = 2, control = list(dthrpc = 0.1))
#' dim(lmh3$mh); dim(lmh5$mh)
#' # m-surroundings or m-histories
#'  lmh3$mh[1:5, ]; lmh5$mh[1:5, ]
#' # length between borderlines in m-surroundings
#' lmh3$mhcbl[1:5, ]; lmh5$mhcbl[1:3, ]
#' # Distances in m-surroundings
#' lmh3$dtmh[1:5,]; lmh5$dtmh[1:5,]
#'
m_surr_cbl <- function(x, m, s = 1, control = list()) {
  con <- list(dthr = 0, dthrpc = 0)
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
  xct <- suppressWarnings(sf::st_centroid(sf::st_geometry(x),
                         of_largest_polygon = TRUE))
  dtfull <- sf::st_distance(xct)
  # Extract geometry
  xgeom <- sf::st_geometry(x)
  mh <- matrix(NA, N, m)
  mhcbl <- matrix(NA, N, m)
  rownames(mh) <- rownames(mhcbl) <- c(1:N)
  colnames(mh) <- colnames(mhcbl) <- c(1:m)
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
    mhcbli <- nbicbl[ord_nbicbl]
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
      dtcti <- dtfull[i, idxordW2i]
      idxdtcti <- order(dtcti, decreasing = FALSE)
      idxordW2i <- idxordW2i[idxdtcti]
      nbiord <- c(nbiord, idxordW2i[1:(m-length(nbiord))])
      mhcbli <- c(mhcbli, rep(0, m-length(mhcbli)))
    }
    mh[i, ] <- nbiord[1:m]
    mhcbl[i, ] <- mhcbli[1:m]
  }
  # Compute distances in m-surroundings
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
  # Build the distance and perimeter matrix for mhs matrix
  mhscbl <- mhcbl[rownames(mhs),]
  dtmhs <- dtmh[rownames(mhs), ]
  mh <- mhs
  dtmh <- dtmhs
  mhcbl <- mhscbl
  rownames(mh) <- rownames(dtmh) <- rownames(mhscbl) <- NULL
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
      mhcbl <- mhcbl[-rowexcl,]
    } else cat("\n None m-surrounding excluded for exceeding
        the threshold distance","\n")
  }

  lmh <- list(mh = mh, dtmh = dtmh, mhcbl = mhcbl)
  return(lmh)
}

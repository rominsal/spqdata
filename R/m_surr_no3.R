#'
#' @title A funcion for creating m-surroudings (mejoradas en tiempo de cómputo)
#' @usage m_surr_no_new(x = x, m = m, s = s, control = NULL)
#' @param x input sf object with points/multipolygons geometry or matrix
#'          of spatial coordinates
#' @param m amplitud m-historia
#' @param s solapamiento máximo entre cualesquiera dos m-surroundings
#' @param control Argumento opcional. Elimina aquellas m-surroundings cuyos elementos
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
#' system.time( msurr_points <- m_surr_no3(x = x, m = m, s = s,
#'                                         control = list(niter = 5,
#'                                         seedinit = 123456,
#'                                         dthr = 0.3)) )
#' mh <- msurr_points$mh
#' dtmh <- msurr_points$dtmh
#' dim(mh); dim(dtmh)
#' mh[1:10,]
#' dtmh[1:10,]
#'
#'
#' # Examples with multipolygons
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' plot(sf::st_geometry(nc))
#' system.time( msurr_polygonsf <- m_surr_no3(x = nc, m = 5, s = 2,
#'                      control = list(niter = 50, seedinit = 123456,
#'                      dthrpc = 0.20)) )
#' mh <- msurr_polygonsf$mh
#' dtmh <- msurr_polygonsf$dtmh
#' dim(mh); dim(dtmh)
#' mh[1:10,]
#' dtmh[1:10,]

m_surr_no3 <- function(x, m, s = 1, control = list()) {
  con <- list(niter = 20, seedinit = 1111,
              dthr = 0, dthrpc = 0)
  #OJO: TAMBIÉN PODEMOS METER EN EL CONTROL EL TIPO DE DISTANCIA, ¿NO?.
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  if (s < 1 || s > (m-1))
    stop("mínimo grado de solapamiento es 1 y menor que la amplitud de la m-historia")
  # Transform matrix coordinates into SpatialPoints class
  if (is.matrix(x)) x <- sp::SpatialPoints(x)
  # Transform Spatial classes into sf class
  if (inherits(x, "Spatial")) x <- as(x, "sf")
  # Compute centroids/coordinates/distances from sf objects
  # VIP: NO ME QUEDA CLARO SI HACE FALTA ESTABLECER CONTROLES PARA IMPEDIR
  #      CENTROIDES DE GEOMETRÍAS. SUPONGO QUE NO...
  if (inherits(x, "sf")) {
    xct <- suppressWarnings(sf::st_centroid(sf::st_geometry(x),
                           of_largest_polygon = TRUE))
  } else stop("object must be either sf, sp or matrix class")
  N <- length(xct)
  Ns <- trunc((N - m)/(m - s)) + 1
  mcoor <- sf::st_coordinates(xct)
  dtfull <- sf::st_distance(xct)
  rownames(dtfull) <- colnames(dtfull) <- as.character(1:N)
  set.seed(con$seedinit)
  niter <- con$niter
  lmh <- vector(mode = "list", length = niter)
  ldtmh <- vector(mode = "list", length = niter)
  for (k in 1:niter) {
    mh <- matrix(NA, Ns, m)
    dt <- matrix(NA, Ns, m)
    Si <- as.character(1:N)
    s0i <- as.character(sample(1:N, 1)) # s0
    dti <- dtfull
    for (i in 1:Ns) {
      dsts0i <- dti[s0i,] # Vector of distances to s0i
      idxnbs0i <- sort(dsts0i, decreasing = FALSE, index.return = TRUE)
      # Distances to s0i ordered
      if (is.list(idxnbs0i)) nbs0i <- names(idxnbs0i$x[2:m])
      if (inherits(idxnbs0i, "units")) nbs0i <- names(idxnbs0i[2:m])
      # Set of (m-1) neighbors to s0i (excluding itself)
      mh[i,] <- c(s0i, nbs0i)
      dt[i,] <- dsts0i[c(s0i,nbs0i)]
      Ai <- c(s0i,nbs0i[1:(m-s-2)])
      Si <- Si[!(Si %in% Ai)]
      dti <- dti[Si, Si]
      s0i <- nbs0i[m-s-1]
      #mhs <- m_surr_maxs(mh, s = s) # CREO QUE NO HACE FALTA
      #dtmhs <- dt[rownames(mhs), ] # CREO QUE NO HACE FALTA
    }
    lmh[[k]] <- mh
    #lmhs[[k]] <- mhs
    ldtmh[[k]] <- dt
    #ldtmhs[[k]] <- dtmhs
  }
  mh <- NULL
  #mhs <- NULL
  dtmh <- NULL
  #dtmhs <- NULL
  for (k in 1:niter) {
    mh <- rbind(mh, lmh[[k]])
    #mhs <- rbind(mhs, lmhs[[k]])
    dtmh <- rbind(dtmh, ldtmh[[k]])
    #dtmhs <- rbind(dtmhs, ldtmhs[[k]])
  }
  # Check for threshold in distances dthr = 0, dthrpc = 0
  if (con$dthr > 0 || con$dthrpc > 0) {
    if (con$dthrpc <= 1)
      dthrpc <- as.numeric(con$dthrpc*max(dtfull)) else
        stop("The value of dthrpc must be in [0,1]")
    if ((dthrpc > 0)  && (dthrpc < con$dthr))
         dthr <- dthrpc else if (con$dthr > 0) dthr <- con$dthr
         else dthr <- dthrpc
  } else dthr <- 0
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
  lmh <- list(mh = matrix(as.integer(mh), nrow=nrow(mh), ncol = ncol(mh)),
              dtmh = dtmh)
  return(lmh)
}
  # nnlist <- matrix(0, N, m - 1)  # Matrix with list of nearest neighbors
  # nn1 <- cbind(1:N, dtfull[1,], mcoor)
  # #nn1 <- cbind(1:N, sqrt((mcoor[1,1] - mcoor[,1])^2 +
  # #                          (mcoor[1,2] - mcoor[,2])^2), mcoor)
  # nn1 <- nn1[order(nn1[,2]), ]
  # nnlist[1, ] <- nn1[2:m, 1]
  # ns <- trunc((N - m)/(m - s)) + 1
  # list <- rep(0, ns)  #zeros(1,ns)
  # list[1] = 1
  # blacklist <- c(1, nnlist[1, 1:(m - (s + 1))])
  # t <- 1
  # for (v in 2:ns) {
  #   list[v] = nnlist[t, m - s]
  #   h <- list[v]
  #   nn1 <- nn1[!nn1[,1] %in% blacklist,]
  #   #**VIP**: SALE DISTINTO PERO DEBERÍA DAR IGUAL... CONSULTAR CON
  #   # FERNANDO
  #   #nn1[,2] <- dtfull[nn1[1,1],nn1[,1]]
  #   nn1[,2] <-  sqrt((nn1[1,3] - nn1[,3])^2 + (nn1[1,4] - nn1[,4])^2)
  #   nn1 <- nn1[order(nn1[,2]),]
  #   nnlist[h,] <- nn1[2:(m),1]
  #   t = h
  #   blacklist <- c(blacklist, h, nnlist[h, 1:(m - (s + 1))])
  # }
  # nnlist <- cbind(1:N, nnlist)
  # nnlist1 <- cbind(nnlist[which(!nnlist[,2] == 0),])
  # # control debería ser una lista
  # # Se eliminan aquellas m-historias que contengan vecinos que no estén
  # # dentro de los k vecinos más próximos
  # if (!is.null(control)){
  #   knn <- cbind(1:N, spdep::knearneigh(mcoor, control)$nn)
  #   int <- numeric()
  #   for (i in 1:dim(nnlist1)[1]){
  #     int[i] <- length(intersect(nnlist1[i,],knn[nnlist1[i,1],]))
  #   }
  #   nnlist1 <- nnlist1[int==m,]
  # }
  #return(nnlist1)
#}

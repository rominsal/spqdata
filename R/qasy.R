#'
#' @title A funcion to compute q test for spatial qualitative data using
#'        asymptotic distribution
#' @usage qasy(formula = NULL, data = NULL, na.action,
#'             m = 3, s = 1, type = "p", control = NULL,
#'             xf = NULL, coord = NULL)
#' @param xf input sf object with points/multipolygons geometry or matrix
#'          of spatial coordinates
#' @param m length of m-surrounding
#' @param s solapamiento máximo entre cualesquiera dos m-surroundings
#' @param type "p" or "c"
#' @param typems type of m-surrounding ("cbl","cdt","no"). Default = "no"
#' @param control Argumento opcional. Por definir
#'
#' @description This function ....
#' @details Aquí Antonio escribe una linda historia
#' @return Por definir, probablemente un objeto tipo spq que herede también
#'         de h-test, ANOVA test o similar...
#' @keywords ...
#' @export
#' @examples
#'
#' # Examples with multipolygons
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' qb79 <- quantile(nc$BIR79)
#' nc$QBIR79 <- (nc$BIR79 > qb79[2]) + (nc$BIR79 > qb79[3]) +
#' (nc$BIR79 >= qb79[4]) + 1
#' nc$QBIR79 <- as.factor(nc$QSID79)
#' plot(nc["QBIR79"], pal = c("#FFFEDE","#FFDFA2", "#FFA93F", "#D5610D"),
#'      main = "BIR79 (Quartiles)")
#' sid79 <- quantile(nc$SID79)
#' nc$QSID79 <- (nc$SID79 > sid79[2]) + (nc$SID79 > sid79[3]) +
#' (nc$SID79 >= sid79[4]) + 1
#' nc$QSID79 <- as.factor(nc$QSID79)
#' plot(nc["QSID79"], pal = c("#FFFEDE","#FFDFA2", "#FFA93F", "#D5610D"),
#'      main = "SID79 (Quartiles)")
#' f1 <- ~ QSID79 + QBIR79
#' lq <- qasy(formula = f1, data = nc, m = 5, s = 2,
#'            typems = "no",
#'            control = list(niter = 10, seedinit = 1111,
#'                           dthrpc = 0.5) )
#' lq$QSID79; lq$QBIR79

#' # Examples with points and matrix of variables
#'
#' xf <- matrix(c(nc$QBIR79, nc$QSID79), ncol = 2, byrow = TRUE)
#' mctr <- suppressWarnings(sf::st_centroid(sf::st_geometry(nc),
#'                                          of_largest_polygon = TRUE))
#' mcoor <- sf::st_coordinates(mctr)[,c("X","Y")]
#' lq <- qasy(xf = xf, mcoor = mcoor, m = 5, s = 2, typems = "no",
#'            control = list(niter = 10, seedinit = 1111,
#'                           dthrpc = 0.5))
#' lq
qasy <- function(formula = NULL, data = NULL, na.action,
                 m = 3, s = 1, typems = "no", control = list(),
                 xf = NULL, mcoord = NULL) {
  con <- list(niter = 20, seedinit = 1111, dthr = 0, dthrpc = 0)
  #OJO: TAMBIÉN PODEMOS METER EN EL CONTROL EL TIPO DE DISTANCIA, ¿NO?.
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

  cl <- match.call()
  # Lectura Datos

  if (!is.null(formula) && !is.null(data)) {
    if (inherits(data, "Spatial")) data <- as(data, "sf")
    mxf <- get_all_vars(formula, data)
    # for (i in 1:ncol(mxf)) {
    #   if (!is.factor(mxf[[i]])) mxf[[i]] <- as.factor(mxf[[i]])
    # }
    # xsf <- dplyr::left_join(data, mxf)
  #   mt <- terms(formula, data = data)
  #   mf <- lm(formula, data = data, na.action = na.action,
  #            method = "model.frame")
  #   mf$drop.unused.levels <- TRUE
  #   na.act <- attr(mf, "na.action")
  } else if (!is.null(xf) && !is.null(mcoord)) {
    mxf <- xf
    if (!is.matrix(mxf)) mxf <- as.matrix(mxf, ncol = 1)
    mxf <- as.data.frame(mxf)
    for (i in 1:ncol(mxf)) {
      if (!is.factor(mxf[,i])) mxf[,i] <- as.factor(mxf[,i])
    }
    mxf <- as.matrix(mxf)
    if (is.matrix(mcoord)) mcoord <- sp::SpatialPoints(mcoord)
    if (inherits(mcoord, "Spatial")) mcoord <- as(mcoord, "sf")
    data <- mcoord #sf object
  } else stop("Data wrong")

  # CÁLCULO M-SURROUNDINGS...
  if (typems == "no") {
    lmh <- m_surr_no3(x = data, m = m, s = s,
                      control = list(niter = con$niter,
                                     seedinit = con$seedinit,
                                     dthr = con$dthr,
                                     dthrpc = con$dthrpc))
  } else if (typems == "cdt") {
    lmh <- m_surr_cdt(x = data, m = m, s = s,
                      control = list(dthr = con$dthr,
                                     dthrpc = con$dthrpc))
  } else if (typems == "cbl") {
    lmh <- m_surr_cbl(x = data, m = m, s = s,
                      control = list(dthr = con$dthr,
                                     dthrpc = con$dthrpc))
  } else stop("Value for argument 'typems' not valid.")
  mh <- lmh$mh
  dtmh <- lmh$dtmh

  lres <- vector(mode = "list", length = ncol(mxf))
  names(lres) <- colnames(mxf)
  for (i in 1:ncol(mxf)) {
    xfi <- mxf[, i]
    if (!is.factor(xfi)) xfi <- as.factor(xfi)
    nlevelsi <- length(levels(xfi))
    symbi <- cr_symb(nlevelsi, m)
    qi <- q_symb_A(xfi, mh, symbi)
    statistic <-  c(qi$qp, qi$qc)
    names(statistic) <- c("Qq", "Qc")
    method <- c("Qp for combinations-totals symbols",
                "Qc for combinations-totals symbols")
    if (!is.null(colnames(mxf)[i]))
      data.name <- rep(colnames(mxf)[i], 2) else data.name <- c(NULL, NULL)
    parameter <- c(qi$qp_df, qi$qc_df)
    names(parameter) <- rep("df", 2)
    p.value <- rep(0, length(parameter))
    for (k in 1:length(p.value)) {
      p.value[k] <- pchisq(statistic[k], df = parameter[k],
                           lower.tail = FALSE)
    }
    lres[[i]] <- vector(mode = "list", length = length(statistic))
    for (j in 1:length(statistic)) {
      lres[[i]][[j]] <- list(statistic = statistic[j],
                             parameter = parameter[j],
                             p.value = p.value[j],
                             method = method[j],
                             data.name = data.name[j])
      class(lres[[i]][[j]]) <- c("spqtest", "htest")
    }
  }
  return(lres)
}

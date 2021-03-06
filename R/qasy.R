#'
#' @title A funcion to compute q test for spatial qualitative data using
#'        asymptotic distribution
#' @usage qasy(formula = NULL, data = NULL, na.action,
#'             m = 3, r = 1, type = "p", control = NULL,
#'             xf = NULL, coord = NULL)
#' @param xf input sf object with points/multipolygons geometry or matrix
#'          of spatial coordinates
#' @param m length of m-surrounding
#' @param r solapamiento máximo entre cualesquiera dos m-surroundings
#' @param distr distribución entre "asymptotic" (default) o "bootstrap"
#' @param type "p" or "c"
#' @param typems type of m-surrounding ("cbl","cdt","no"). Default = "no"
#' @param control Argumento opcional.
#' @param distance Dentro de control: character dando tipo de
#'                distancia a calcular. Default = "Euclidean"
#'                for Cartesian coordinates only: one of Euclidean,
#'                 Hausdorff or Frechet; for geodetic coordinates,
#'                 great circle distances are computed (see sf::st_distance())
#' @param nsim Dentro de control: número de simulaciones
#'             para distribución bootstrap. Default = 999
#' @param seedinit Dentro de control: semilla inicial para selección
#'                 observación inicial en algoritmo para obtener m-surroundings
#'                 o para iniciar proceso bootstrap
#' @description This function ....
#' @details Aquí Antonio escribe una linda historia
#' @return Por definir, probablemente un objeto tipo spq que herede también
#'         de h-test, ANOVA test o similar...
#' @keywords ...
#' @export
#' @examples
#'
#' # Ejemplos Fernando
#' # Ejemplo 1
#' data("FastFood")
#' # Add other dependent variable to test...
#' FastFood.sf$TypeII <- sample(FastFood.sf$Type,
#'   length(FastFood.sf$Type), replace = TRUE)
#' f1 <- ~ Type + TypeII
#' vm <- c(3, 4)
#' vr <- c(1, 2, 3)
#' lq <- qasy(formula = f1, data = FastFood.sf, m = vm, r = vr,
#'            distr = "asymptotic",
#'            #typems = "no",
#'            control = list(distance = "Euclidean",
#'            seedinit = 1111))
#' summary(lq)
#' summary.spqtest(lq)
#' lq
#' lqboot <- qasy(formula = f1, data = FastFood.sf, m = vm, r = vr,
#'            distr = "bootstrap",
#'            #typems = "no",#'
#'            control = list(distance = "Euclidean",
#'            seedinit = 1111, nsim = 50))


#'  # Ejemplo 2
#'
#' data("Spain")
#' f1 <- ~ COSTA
#' lq <- qasy(formula = f1, data = spain.sf, m = 4, r = 3,
#'            typems = "no", control = list(seedinit = 1111))
#' summary(lq)
#' lq


#' # Examples with multipolygons
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' qb79 <- quantile(nc$BIR79)
#' nc$QBIR79 <- (nc$BIR79 > qb79[2]) + (nc$BIR79 > qb79[3]) +
#' (nc$BIR79 >= qb79[4]) + 1
#' nc$QBIR79 <- as.factor(nc$QBIR79)
#' plot(nc["QBIR79"], pal = c("#FFFEDE","#FFDFA2", "#FFA93F", "#D5610D"),
#'      main = "BIR79 (Quartiles)")
#' sid79 <- quantile(nc$SID79)
#' nc$QSID79 <- (nc$SID79 > sid79[2]) + (nc$SID79 > sid79[3]) +
#' (nc$SID79 >= sid79[4]) + 1
#' nc$QSID79 <- as.factor(nc$QSID79)
#' plot(nc["QSID79"], pal = c("#FFFEDE","#FFDFA2", "#FFA93F", "#D5610D"),
#'      main = "SID79 (Quartiles)")
#' f1 <- ~ QSID79 + QBIR79
#' lq1nc <- qasy(formula = f1, data = nc, m = 5, r = 2,
#'            typems = "no", control = list(seedinit = 1111,
#'                                          dtmaxpc = 0.5) )
#' lq1nc$QSID79; lq1nc$QBIR79
#'
#' lq2nc <- qasy(formula = f1, data = nc, m = 5, r = 2,
#'               typems = "cbl", control = list(dtmaxpc = 0.2) )
#' lq2nc$QSID79; lq2nc$QBIR79
#'
#' lq3nc <- qasy(formula = f1, data = nc, m = 5, r = 2,
#'               typems = "cdt", control = list(dtmaxpc = 0.2) )
#' lq3nc$QSID79; lq3nc$QBIR79
#' # Examples with points and matrix of variables
#'
#' xf <- matrix(c(nc$QBIR79, nc$QSID79), ncol = 2, byrow = TRUE)
#' mctr <- suppressWarnings(sf::st_centroid(sf::st_geometry(nc)))
#' mcoor <- sf::st_coordinates(mctr)[,c("X","Y")]
#' lq <- qasy(xf = xf, mcoor = mcoor, m = 5, r = 2, typems = "no",
#'            control = list(seedinit = 1111,
#'                           dtmaxpc = 0.5))
#' lq
#'
#'
qasy <- function(formula = NULL, data = NULL, na.action,
                 m = 3, r = 1, distr = "asymptotic",
                 #typems = "no",
                 control = list(),
                 xf = NULL, mcoord = NULL) {
  con <- list(distance = "Euclidean",
              nsim = 999, seedinit = NULL,
              dtmaxabs = 0, dtmaxpc = 0)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  distance <- con$distance
  nsim <- con$nsim
  seedinit <- con$seedinit
  dtmaxabs <- con$dtmaxabs
  dtmaxpc <- con$dtmaxpc
  cl <- match.call()
   # Lectura Datos.
  if (!is.null(formula) && !is.null(data)) {
    if (inherits(data, "Spatial")) data <- as(data, "sf")
    mxf <- model.frame(formula, data, na.action = na.action)
    #mxf <- get_all_vars(formula, data)
  } else if (!is.null(xf) && !is.null(mcoord)) {
    mxf <- xf
    if (!is.matrix(mxf)) mxf <- as.matrix(mxf, ncol = 1)
    mxf <- as.data.frame(mxf)
    if (is.matrix(mcoord)) mcoord <- sp::SpatialPoints(mcoord)
    if (inherits(mcoord, "Spatial")) mcoord <- as(mcoord, "sf")
    data <- mcoord #sf object
  } else stop("input data wrong")
  for (i in 1:ncol(mxf)) {
    if (!is.factor(mxf[,i])) mxf[,i] <- as.factor(mxf[,i])
  }
  mxf <- as.matrix(mxf)
  lres <- NULL # Initialize the list of results
  # Same random initial observation for every m-sorrounding
  if (!is.null(seedinit)) set.seed(seedinit)
  initobs <- sample(1:nrow(mxf), 1)
  for (i in 1:ncol(mxf)) {
    xfi <- mxf[, i]
    if (!is.factor(xfi)) xfi <- as.factor(xfi)
    N <- length(xfi)
    ki <- length(levels(xfi))
    lms <- vector(mode = "list", length = length(m)) #List of m-sorroundings
    for (j in seq_len(length(m))) {
      mj <- m[j]
      symbi <- cr_symb(ki, mj)
      if (distr == "asymptotic") {
        rcut <- r[!(r > (mj - 1))]
        lms[[j]] <- vector(mode = "list",
                           length = length(rcut))
        for (h in seq_len(length(rcut))) {
          rh <- rcut[h]
           ### M-SURROUNDINGS...
          lms[[j]][[h]] <- m_surr_no3(x = data, m = mj, r = rh,
                                      distance = distance,
                                      control = list(initobs = initobs,
                                                     dtmaxabs = dtmaxabs,
                                                     dtmaxpc = dtmaxpc))
          # if (typems == "no") {
          R <- nrow(lms[[j]][[h]]$ms)
          #   lms[[j]][[h]] <- m_surr_no3(x = data, m = mj, r = rh,
          #                        control = list(seedinit = seedinit,
          #                                       dtmaxabs = dtmaxabs,
          #                                       dtmaxpc = dtmaxpc))
          #   #lms__jkpr <- m_surr_no(x = sf::st_coordinates(data), m = mj, s = rh)
          # } else if (typems == "cdt") {
          #   lms[[j]][[h]] <- m_surr_cdt(x = data, m = mj, r = rh,
          #                        control = list(dtmaxabs = dtmaxabs,
          #                                       dtmaxpc = dtmaxpc))
          # } else if (typems == "cbl") {
          #   lms[[j]][[h]] <- m_surr_cbl(x = data, m = mj, r = rh,
          #                        control = list(dtmaxabs = dtmaxabs,
          #                                       dtmaxpc = dtmaxpc))
          # } else stop("Value for argument 'typems' not valid.")
          qi <- q_symb_A2(xfi, lms[[j]][[h]]$ms, symbi)
          qi$qp_df <- nrow(symbi$p_symb) - 1
          qi$qc_df <- nrow(symbi$c_symb) - 1
          statistic <-  c(qi$qp, qi$qc)
          names(statistic) <- c("Qp", "Qc")
          method <- c("Qp (asymptotic distrib.) for symbolization based on permutations",
                      "Qc (asymptotic distrib.) for symbolization based on combinations")
          if (!is.null(colnames(mxf)[i])) {
            data.name <- rep(colnames(mxf)[i], 2)
          } else data.name <- c(NULL, NULL)
          parameter <- c(qi$qp_df, qi$qc_df)
          names(parameter) <- rep("df", 2)
           p.value <- rep(0, length(parameter))
          for (l in 1:length(p.value)) {
            p.value[l] <- pchisq(statistic[l], df = parameter[l],
                                 lower.tail = FALSE)
          }
          lres_i <- vector(mode = "list", length = length(statistic))
          for (t in 1:length(statistic)) {
            lres_i[[t]]<- list(statistic = statistic[t],
                               parameter = parameter[t],
                               p.value = p.value[t],
                               method = method[t],
                               data.name = paste("Variable",
                                                 data.name[t],
                                                 " m = ",mj," r = ", rh),
                               var.name = data.name[t],
                               N = N,
                               R = R,
                               m = mj, r = rh,
                               k = ki, df = parameter[t])
            lres_i[[t]]$distr <- distr
            lres_i[[t]]$ms <- lms[[j]][[h]]$ms
            lres_i[[t]]$mdtms <- lms[[j]][[h]]$mdtms
            lres_i[[t]]$initobs <- lms[[j]][[h]]$initobs
            lres_i[[t]]$distance <- lms[[j]][[h]]$distance
            if(names(statistic)[t] == "Qp") {
              lres_i[[t]]$type <- "permutations"
              lres_i[[t]]$symb <- symbi$p_symb
              lres_i[[t]]$efp_symb <- qi$efp_symb
              lres_i[[t]]$qp_symb <- qi$qp_symb
              lres_i[[t]]$PSymb <- qi$PSymb
            } else if(names(statistic)[t] == "Qc") {
              lres_i[[t]]$type <- "combinations"
              lres_i[[t]]$symb <- symbi$c_symb
              lres_i[[t]]$efc_symb <- qi$efc_symb
              lres_i[[t]]$qc_symb <- qi$qc_symb
              lres_i[[t]]$CSymb <- qi$CSymb
            } else { stop("statistic neither Qp or Qc") }
            lres_i[[t]]$n <- nrow(lres_i[[t]]$symb)
            class(lres_i[[t]]) <- c("htest")
          }
          lres <- c(lres, lres_i)
        } # end for (h in seq_len(length(rcut)))
      } else if (distr == "bootstrap") {
        qi <- q_boots(fx = xfi, x = data,
                      m = mj, nsim = nsim,
                      distance = distance,
                      seedinit = seedinit)
        statistic <-  c(qi$qp, qi$qc)
        names(statistic) <- c("Qp", "Qc")
        method <- c("Qp (bootstrap distrib.) for symbolization based on permutations",
                    "Qc (bootstrap distrib.) for symbolization based on combinations")
        if (!is.null(colnames(mxf)[i])) {
          data.name <- rep(colnames(mxf)[i], 2)
        } else data.name <- c(NULL, NULL)
        parameter <- c("NA", "NA")
        names(parameter) <- rep("df", 2)
        p.value <- c(qi$pvalueboot_qp,
                     qi$pvalueboot_qc)
        lres_i <- vector(mode = "list",
                         length = length(statistic))
        for (t in 1:length(statistic)) {
          lres_i[[t]]<- list(statistic = statistic[t],
                             parameter = parameter[t],
                             p.value = p.value[t],
                             method = method[t],
                             data.name = paste("Variable",
                                               data.name[t],
                                               " m = ",mj,
                                               " r = ", mj-1),
                             var.name = data.name[t],
                             N = N,
                             R = nrow(qi$ms),
                             m = mj, r = mj-1,
                             k = ki, df = parameter[t])
          lres_i[[t]]$distr <- distr
          lres_i[[t]]$ms <- qi$ms
          lres_i[[t]]$mdtms <- qi$mdtms
          lres_i[[t]]$distance <- qi$distance
          if(names(statistic)[t] == "Qp") {
            lres_i[[t]]$type <- "permutations"
            lres_i[[t]]$symb <- qi$symb$p_symb
            lres_i[[t]]$efp_symb <- qi$efp_symb_qp
            lres_i[[t]]$qp_boot <- qi$qpboot
            lres_i[[t]]$efp_symb_boot <- qi$efp_symb_boot
            lres_i[[t]]$qp_symb <- qi$qp_symb
          } else if(names(statistic)[t] == "Qc") {
            lres_i[[t]]$type <- "combinations"
            lres_i[[t]]$symb <- qi$symb$c_symb
            lres_i[[t]]$efc_symb <- qi$efc_symb_qc
            lres_i[[t]]$qc_boot <- qi$qcboot
            lres_i[[t]]$efc_symb_boot <- qi$efc_symb_boot
          } else { stop("statistic neither Qp or Qc") }
          lres_i[[t]]$n <- nrow(lres_i[[t]]$symb)
          class(lres_i[[t]]) <- c("htest")
        }
        lres <- c(lres, lres_i)
      }
    } # end for (j in seq_len(length(m)))
  } # end for (i in 1:ncol(mxf))
  class(lres) <- c("spqtest", class(lres))
  return(lres)
}

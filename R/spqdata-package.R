#' @docType package
#' @name spqdata-package
#' @rdname spqdata-package
#'
#' @title Spatial Qualitative Data.
#'
#' @description
#'  \pkg{spqdata} A nice r package to analyze the spatial independence of cathegorical variables.
#'
#' @details
#'  Some functionalities that have been included in \pkg{spqdata} package are:
#'
#' @section 1. Test Q
#' @section 2. Spatial run test
#' @section Datasets: (example Spain, FastFood)
#'   \pkg{spqdata} includes two different datasets: spain and Fastfood. These
#'   sets are used to illustrate the capabilities of different functions.
#'   Briefly, their
#'   main characteristics are the following \cr
#'   \itemize{
#'     \item The \emph{FastFood} A sf object with points of the localization of Fastfood restaurants in Toronto.
#'      \item The \emph{spain.sp} An object of the class sf with the geometry of spanish provinces and several data sets.
#'    }
#'
#' @references
#'   \itemize{
#'     \item Breusch T, Pagan A (1980). The Lagrange multiplier test
#'       and its applications to model specification in econometrics.
#'       \emph{Review of Economic Studies} 47: 239-254.
#'
#'      \item LeSage, J., and Pace, R. K. (2009). \emph{Introduction to
#'        spatial econometrics}. Chapman and Hall/CRC.
#'
#'      \item Lopez, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science},
#'        53(1), 197-220.
#'
#'      \item López, F.A., Martínez-Ortiz, P.J., and Cegarra-Navarro, J.G.
#'        (2017). Spatial spillovers in public expenditure on a municipal
#'        level in Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#'
#'      \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'        effects in seemingly unrelated regressions. \emph{Spatial Economic
#'        Analysis}, 5(4), 399-440.
#'   }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @importFrom dplyr group_by rename select
#' @importFrom gtools combinations permutations
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar
#' @importFrom ggplot2 aes element_text labs theme
#' @importFrom magrittr %>%
#' @importFrom Matrix solve
#' @importFrom methods as
#' @importFrom rsample bootstraps analysis
#' @importFrom sf st_coordinates st_distance
#' @importFrom spdep knearneigh knn2nb poly2nb
NULL

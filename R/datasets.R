#' España por provincias
#'
#' Conjunto de datos de prueba. PARA PROBAR ELEMENTOS SIN VECINOS
#' @docType data
#'
#' @usage data(Spain)
#'
#' @format a sf object with 50 rows and 7 variables:
#' \describe{
#'   \item{Province}{Province name.}
#'   \item{Population}{Population.}
#'   \item{Older65}{Factor 3 categorias en funcion del % mayores.}
#'   \item{AverageAge}{Edad Media.}
#'   \item{MenWoman}{Factor 2 categorias "men" si % hombres > muejres}
#'   \item{MassTransitSystems}{Tiene metro.}
#'   \item{Coast}{Tiene costa}
#' }
#'
#' @source \url{https://onlinelibrary.wiley.com/doi/full/10.1111/gean.12241}
#'
#' @references
#'   \itemize{
#'     \item Paez, A., Lopez, F. A., Menezes, T., Cavalcanti, R., & Pitta, M. (2020). \emph{A Spatio‐Temporal Analysis of
#'      the Environmental Correlates of COVID‐19 Incidence in Spain.}. Geographical Analysis.
#'   }
"spain.sf"

#' Distribution Fast-Food restaurants in Toronto
#'
#' @docType data
#'
#' @usage data(FastFood)
#'
#' @format A sf object with 877 rows and 4 variables:
#'
#' \describe{
#'   \item{ID}{Identification}
#'   \item{Lat}{Latitude.}
#'   \item{Lon}{Longitude.}
#'   \item{Type}{Factor with 3 types of restaurant, "H" Hamburguer "P" Pizza "S" Shanwich.}
#' }
#'
#' @source El papel con A y M en JGS
#'
#' @references
#'   \itemize{
#'     \item PRuiz M, López FA, A Páez. (2010). \emph{Testing for spatial association of qualitative
#'     data using symbolic dynamics}. Journal of Geographical Systems. 12 (3) 281-309
#'   }
"FastFood.sf"

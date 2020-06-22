#' España por provincias 
#'
#' Conjunto de datos de prueba. PARA PROBAR ELEMENTOS SIN VECINOS
#'
#' @format A sf object with 50 rows and X variables:
#' 
#' \describe{
#'   \item{PROVINCIA}{Province name.}
#'   \item{POPULATION}{Population.}
#'   \item{MAY_65}{Factor 3 categorias en funcion del % mayores.}
#'   \item{EDAD_MED}{Edad Media.}
#'   \item{HOM_MUJ}{Factor 2 categorias "men" si % hombres > muejres}
#'   \item{METRO}{Tiene metro.}
#'   \item{COSTA}{Tiene costa}
#'   \item{AREA}{Area.}
#'   \item{ACRES}{Acres}
#'   \item{PERIMETER}{Perimetro}
#' }
#'
#' @source El papel del COVID-19 con A
#'
#' @references
#'   \itemize{
#'     \item Paez, A., Lopez, F. A., Menezes, T., Cavalcanti, R., & Pitta, M. (2020). \emph{A Spatio‐Temporal Analysis of
#'      the Environmental Correlates of COVID‐19 Incidence in Spain.}. Geographical Analysis.
#'   }
"spain.sf"

#' Distribution Fast-Food restaurant in Toronto 
#'
#' Conjunto de datos de prueba. PARA PROBAR ELEMENTOS SIN VECINOS
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
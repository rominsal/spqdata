#' @title A funcion to calculate Q boots
#'
#' @description This function calculates Qboot, a measure of spatial association based on symbolic entropy.
#' @param Y a factor of the same length as the coordinates x
#' @param x coordenadas asociadas a la observación de Y
#' @param m amplitud de las m-historias
#' @param s grado de solapamiento
#' @param nsim number of permutations
#' @usage q_boots(Y, x, m, s, nsim=999)
#' @keywords spatial association, qualitative variable, symbolic entropy, symbols
#' @details Aquí Antonio escribe una linda historia
#' @return decir que cosas son las que devuelve
#'   \tabular{ll}{
#'     \code{Qboots_p} \tab value of the statistic of the observed distribution. Symbols p\cr
#'     \code{p.value.p} \tab  the pseudo p-value of the Qboots_p test \cr
#'     \code{efp_symb} \tab frecuencia empírica de los simbolos "p" de cada permutación\cr
#'     \code{Qboots_c} \tab value of the statistic of the observed distribution. Symbols c\cr
#'     \code{p.value.c} \tab  the pseudo p-value of the Qboots_c test \cr
#'     \code{efc_symb} \tab frecuencia empírica de los simbolos "c" de cada permutación\cr
#'     \code{nsim} \tab  nsim simulated values of statistic \cr
#'     }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz  \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     Testing for spatial association of qualitative data using symbolic dynamics.
#'       \emph{Journal of Geographical Systems}, 12(3), 281-309.
#'     \item López, F., and Páez, A. (2012).
#'     Distribution-free inference for Q(m) based on permutational bootstrapping: an application
#'     to the spatial co-location pattern of firms in Madrid.
#'       \emph{Estadística Española}, 177, 135-156.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}, \code{\link{m_surr_no}},\code{\link{q_symb}}
#' @export
#' @examples
#' #
#' N <- 1000
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- spdep::nb2listw(knn2nb(knearneigh(cbind(cx,cy), k=4)))
#' p <- c(1/6,3/6,2/6)
#' rho = 0.0
#' QY <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
#' Qboots <- q_boots(QY,x,m=3,s=1,nsim=99)


q_boots <- function(fx, x, m, s, nsim=999) {

  if (is.factor(fx)){
    levels(fx) <- as.character(1:length(levels(fx)))
    Y <- as.numeric(fx)
  }
  if (length(Y) != dim(x)[1])
    stop("La longitud e Y no coincide con la dimensión de las coordenadas")
  if (s<1 || s >(m-1))
    stop("mínimo grado de solapamiento es 1 y menor que la amplitud de la m-historia")

 num_clases <- length(table(Y))
 mh <- m_surr_no(x,m,s)
 symb <- cr_symb(num_clases, m)
 Q0 <- q_symb(Y,mh,symb)
 Q0r <- c(Q0$qp,Q0$qc)
 set.seed(123)
 qlist <- list()
 Qpr <- as.numeric()
 efp_symb <- matrix(0,ncol =dim(symb$p_symb)[1],nrow = nsim)
 efc_symb <- matrix(0,ncol =dim(symb$c_symb)[1],nrow = nsim)
 for (i in 1:nsim){
   Yp <- Y[sample(length(Y))]
   Qp <- q_symb(Yp,mh,symb)
   Qpr <- rbind(Qpr,c(Qp$qp,Qp$qc))
   efp_symb[i,] <- Qp$efp_symb
   efc_symb[i,] <- Qp$efc_symb
 }
 pseudovalor_p=sum(Qpr[,1]>Q0r[1])/(nsim+1)
 pseudovalor_c=sum(Qpr[,2]>Q0r[2])/(nsim+1)
 results <- list(Q0r[1], pseudovalor_p, Q0r[2], pseudovalor_c,efp_symb,efc_symb,nsim) #, p_symb_plot, c_symb_plot)
 names(results) <- c("Q_p","p.value.p","Q_c","p.value.c","efp_symb","efc_symb","nsim")
 return(results)
}

#' @title A funcion to calculate QMap
#'
#' @description This function COMPARA MAPAS based on symbolic entropy.
#' @param df Es un data frame con tantas columnas como mapas a comparar
#' @param m longitud m-historia
#' @param r solapamiento
#' @param type Type of symbols: "standard" or "compact". Default "standard"
#' @usage QMap(df, mh, symb, type)
#' @keywords spatial association, qualitative variable, symbolic entropy, symbols
#' @details Aquí Antonio escribe una linda historia
#' @return decir que cosas son las que devuelve
#'   \tabular{ll}{
#'     \code{QMap} \tab valor del estadístico \cr
#'     \code{pmap_pval} \tab  p-value for QMap \cr
#'     \code{gl} \tab grados de libertad \cr
#'     \code{Tm} \tab  numero de mapas \cr
#'     \code{R} \tab  tamaño de cada mapa \cr
#'     }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio  \tab \email{paez@@gmail.com} \cr
#'   Manolo  \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz M, López FA and A Páez (2011).
#'     Comparison of Thematic Maps Using Symbolic Entropy.
#'       \emph{International Journal of Geographical Information Science},  26, 413-439.
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     Testing for spatial association of qualitative data using symbolic dynamics.
#'       \emph{Journal of Geographical Systems}, 12(3), 281-309.0.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}, \code{\link{m_surr_no}}
#' @export
#' @examples
#' # Load dataset
#' N <- 1000
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- spdep::nb2listw(spdep::knn2nb(knearneigh(cbind(cx,cy), k=4)))
#' p <- c(2/6,2/6,2/6)
#' rho = 0.5
#' QY1 <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
#' rho = 0.7
#' QY2 <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
#' df = data.frame(cbind(QY1,QY2))
#' df$QY1 <- as.factor(df$QY1)
#' df$QY2 <- as.factor(df$QY2)
#' m=3
#' s=1
#' myqmap <- QMap(df, m=m, s=s, type ="standard")
#' myqmap$qmap
#' myqmap$qmap_pval
#' myqmap <- QMap(df, m=m, s=s, type ="compact")
#' myqmap$qmap
#' myqmap$qmap_pval


QMap <- function(df, m=m, s=s, type=type) {

  R <- dim(df)[1]
  Tm <- dim(df)[2]
  # Transformar factores en numeros
  for ( t in 1: Tm){
    levels(df[,1]) <- as.character(1:length(levels(df[,t])))
    df[,t] <- as.numeric(df[,t])
  }

  # Calcular simbolos permutation
  symb <- cr_symb(m=m,max(df[,1]))$p_symb
  mh <- m_surr_no(x=x,m=m,s=s)

  Z <- list()
  for (t in 1:Tm){
    Z[[t]]<- matrix(df[,t][mh], ncol = m)
  }

  nusi <- dim(symb)[1]
  nsk <- matrix(0,ncol = nusi,nrow = Tm)
  for (t in 1:Tm){
    for (i in 1:nusi){
      nsk[t,i] <- sum(rowSums(Z[[t]]==symb[i,])==m)
    }
  }
  ##################################
  # Calcular simbolos combinaciones
  # Simbolizo la serie a partir de Z
  Z2 <- list()
  # for (t in 1:Tm){
  # kk <- numeric()
  # for (i in 1:3){
  #   kk <- cbind(kk,rowSums(Z[[t]]==i))
  # }
  # Z2[[t]] <- kk
  # }

  for (t in 1:Tm){
    Z2[[t]]<- t(apply(Z[[t]], 1, function(x) table(factor(x, levels = unique(sort(c(Z[[t]])))))))
  }


  symb2 <- cr_symb(m=m,max(df[,1]))$c_symb
  nusi2 <- dim(symb2)[1]
  nsk2 <- matrix(0,ncol = nusi2,nrow = Tm)
  for (t in 1:Tm){
    for (i in 1:nusi2){
      nsk2[t,i] <- sum(rowSums(Z2[[t]]==symb2[i,])==3)
    }
  }
  #####
  # AUNQUE CALCULA LAS ENTROPIAS CON LOS DOS TIPOS DE SIMBOLOS
  # SOLO DA EL OUTPUT DE UNO DE ELLOS
  # ESTO SE PUEDE CAMBIAR CLARO

  if (type=="compact"){
    nsk <- nsk2
    nusi <- nusi2
    symb <- symb2
  }
  #####
  fnk <- nsk/(R*Tm)
  fnT=colMeans(nsk)/(R*Tm)

  lnsk <- matrix(0,ncol = nusi,nrow = Tm)
  for (t in 1:Tm){
    lnsk[t,] <-  log(fnT/fnk[t,])
  }
  lnsk[lnsk==Inf] <- 0
  lnsk[is.na(lnsk)] <- 0
  QMap=-2*sum(nsk*lnsk)
  gl <- (Tm-1)*nusi
  qmap_pval = pchisq(QMap, df = gl, lower.tail = FALSE)
  results <- list(QMap,qmap_pval,gl,Tm,R,symb)
  names(results) <- c("qmap","qmap_pval","gl","Tm","R","symb")
  #####

  return(results)
}

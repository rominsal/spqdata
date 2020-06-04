#' A funcion for creating m-surroudings ALTERNATIVA
#'
#' This function obtains the m-surroundings by selecting the m-1 nearest neighbors of each observation, allowing for a degree of overlap of s.
#' @param x son coordenadas
#' @param m amplitud m-historia
#' @param s solapamiento máximo entre cualesquiera dos m-historias
#' @description Esta función calcula m-historias. El algoritmo basado en m_surr_no2 obtiene m-historias
#' cuyos elementos pueden estar muy distantes (esto distorsiona el resultado del contraste). Este NUEVO algoritmo está basado en los knn.
#' Parte de una m-historia semilla (la correspondiente a la observación que se encuentra en primer lugar). La siguiente m-historia se añade siempre
#' que no solape en más de "s" elementos a todas las anteriores. El proceso se repite, incorporando nuevas m-historias siempre que no solapen
#' en más de "s" elementos a todas las anteriores.
#' SE OBTIENE MENOR NÚMERO DE M-HISTORIAS QUE CON LA VERSION CLÁSICA
#' PERO SON MÁS COMPACTAS.
#' NO HAY UN NÚMERO PREFIJADO DE M-HISTORIAS. EL ALGORITMO DEPENDE DE LA SEMILLA Y DEL ORDEN EN QUE SE INCLUYAN LAS M-HISTORIAS.
#' @keywords m_surround
#' @export
#' @examples
#'
#' # Obtain m-surroundings of size 3 (m=5), with degree of overlap of two (s=2)
#' N <- 100
#' m = 5
#' s = 2
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' mh <- m_surr_no2(x, m, s)
#' W <- matrix(0,ncol = N,nrow = N)
#' for (i in 1:dim(mh)[1]){
#'  W[mh[i,1],mh[i,2:m]]<-1
#' }
#' W <- mat2listw(W)
#' plot(W,x)
#'
#'
#'

m_surr_no2 <- function(x, m, s) {

N <- dim(x)[1]
knn <- cbind(1:N,spdep::knearneigh(x, k=(m-1))$nn)
Knn <- knn[1:2,]
for (f in 3:N){
mis <- 0
for (j in 1:m){
  mis <- mis + rowSums(Knn==knn[f,j])
}
if (sum(mis<(s+1))==dim(Knn)[1]){
  Knn <- rbind(Knn,knn[f,])
}
}

# # Calculando el índice de solapamiento REAL
# moverlap <- matrix(0,ncol = dim(Knn)[1],nrow = dim(Knn)[1])
# overlap <- 0
# for (f in 1:(dim(Knn)[1]-2)){
# a <- Knn[f,]
# b <- Knn[-(1:f),]
# for (g in 1:dim(b)[1]){
#   overlap  <- overlap  + length(intersect(b[g,],a))
#   moverlap[f,(g+1)] <- length(intersect(b[g,],a))
# }
# f=dim(Knn)[1]-1
# a <- Knn[f,]
# b <- Knn[-(1:f),]
# overlap  <- overlap  + length(intersect(b,a))
# }
#
# results <- list(Knn,overlap)
# names(results)[1] <- "mh"
# names(results)[1] <- "overlap"
return(Knn)
}

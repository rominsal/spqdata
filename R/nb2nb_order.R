# Esta es una funcion auxiliar para creation_nvar_SR
#  Toma un objeto sf con estructura poligonal
#  calcula los vecinos con frontera comun creando el objeto nb
#  los ordena en función de la DISTANCIA a centroides
#  devuelve una matriz con los elementos ordenados
#  rellena con -99 para usarse en creation_nvar_SR

nb2nb_order <- function(listw = listw, sf = sf){
n <- length(listw)
co <- as.data.frame(sf::st_coordinates(st_centroid(sf)))
co <- st_as_sf(co,coords=c("X","Y"), crs = 32630)
dis <- 1/matrix(as.numeric(st_distance(co,co)),ncol=n,nrow=n)
diag(dis) <- 0
matw <- nb2mat(listw,style = 'B', zero.policy = TRUE)
m <- rowSums(matw)
matwdis <- dis*matw
NB <- list()
for (i in 1:n){
  or <-  order(matwdis[i,], decreasing = TRUE)
  NB[[i]]<- or[1:m[i]]
  if (m[i]==0){NB[[i]]<- as.integer(0)}
}
class(NB)<- 'nb'
return <- NB
 }


#'
#' @title A function to calculate the global spatial runs test (SRQ)
#'
#' @description This function calculates the global spatial runs test for spatial independence of a categorical spatial data set.
#' @param data an (optional) data frame or a sf object containing the variable to testing for.
#' @param listw una lista de vecinos (tipo knn o nb) o una matrix W que indique el orden de cada $m_i-entorno$
#' (por ejemplo de inversa distancia). Para calcular el numero de rachas en cada m_i-entorno debe establecerse
#' un orden, por ejemplo del vecino más próximo al más lejano.
#' @param xf un factor. vector de observaciones de la clase factor.
#' @param alternative a character string specifying the alternative hypothesis, must be one
#' of "two.sided", "greater" or "less" (default).
#' @param nsim Number of permutations to optain confidence intervals (CI).
#' Default value is NULL to don`t get CI of number of runs.
#' @usage srq_test(formula = NULL, data = NULL, na.action, xf = NULL,
#' listw = listw, alternative = "greater" , nsim = NULL, control = list())
#' @keywords spatial association, qualitative variable, runs test
#' @details El orden de las vecindades (m_i-entornos) es crítico.
#' El objeto listw puede ser de tres tipos:
#'   \tabular{ll}{
#'     \code{knn} \tab Un objeto tipo knn obtenida utilizando el criterio del vecino más próximo, obtenida por ejemplo con \code{\link{knearneigh}}\cr
#'     \code{nb} \tab Un objeto tipo nb obtenido con \code{\link{poly2nb}}\cr
#'     \code{matrix} \tab Una matriz indicando el orden de vecindad de cada m_i-entorno.
#'     Por ejemplo basada en inversa de la distancia con cierto punto de corte.\cr
#'     }
#' @return decir que devuelve
#'   \tabular{ll}{
#'     \code{SR} \tab numero total de rachas \cr
#'     \code{dnr} \tab distribución empírica del numero de rachas \cr
#'     \code{SRQ} \tab Test de homogeneidad. Negative sign indicates global homogeneity \cr
#'     \code{p.valueSRQ} \tab p-value of the SRQ \cr
#'     \code{p.valueSRQB} \tab the pseudo p-value of the SRQ test\cr
#'     \code{MeanNeig} \tab Número medio de vecinos\cr
#'     \code{MaxNeig} \tab Número máximo de vecinos\cr
#'     \code{listw} \tab la lista de vecinos\cr
#'     \code{nsim} \tab numero de repeticiones\cr
#'     \code{SRGP} \tab vector con los 'nsim' valores de SRQ por remuestreo permutacional\cr
#'     \code{SRLP} \tab matrix con el número de rachas en cada localización para cada permutación\cr
#'     }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paez@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     A test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}
#' @export
#' @examples
#'
#' # SRQ test based on knn
#' rm(list = ls())
#' N <- 1000
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- knearneigh(cbind(cx,cy), k=3)
#' p <- c(1/6,3/6,2/6)
#' rho <- 0.5
#' xf <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
#' srq <- srq_test(xf = xf, listw = listw)
#' srq$SRQ
#' plot.srq(srq)
#'
#' # Fastfood example. sf (points)
#' data("FastFood")
#' x <- cbind(FastFood.sf$Lon,FastFood.sf$Lat)
#' listw <- spdep::knearneigh(x, k=4)
#' formula <- ~ Type
#' srq <- srq_test(formula = formula, data = FastFood.sf, listw = listw)
#' srq$SRQ
#' plot.srq(srq)
#'
#' # With a sf object (poligons)
#' library(sf)
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' listw <- spdep::poly2nb(as(nc,"Spatial"), queen = FALSE)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' co <- sf::st_coordinates(st_centroid(nc))
#' nc$xf <- dgp_spq(x = co, p = p, listw = listw, rho = rho)
#' plot(nc["xf"])
#' formula <- ~ xf
#' srq <- srq_test(formula = formula, data = nc, listw = listw, nsim=399)
#' c(srq$SRQ,srq$p.valueSRQ,srq$p.valueSRQB)
#' plot.srq(srq)
#'
#'# With a sf object (poligons) WITH ZEROS
#' data(Spain)
#' listw <- spdep::poly2nb(spain.sf, queen = FALSE)
#' listw <- nb2nb_order(listw,spain.sf)
#'
#' # SRQ test based on inverse distance
#' rm(list = ls())
#' N <- 500
#' cx <- runif(N)
#' cy <- runif(N)
#' x <-cbind(cx,cy)
#' xx <- as.data.frame(cbind(cx,cy))
#' xx <- st_as_sf(xx,coords = c("cx","cy"))
#' n = dim(xx)[1]
#' dis <- 1/matrix(as.numeric(st_distance(xx,xx)),ncol=n,nrow=n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.05))*dis
#' rowSums(dis>0)

#' p <- c(1/6,3/6,2/6)
#' rho <- 0.5
#' xf <- dgp_spq(x = x, p = p, listw = dis, rho = rho)
#' srq <- srq_test(xf = xf, listw = dis)
#' srq$SRQ
#' plot.srq(srq)
#'
#' data("FastFood")
#' n = dim(FastFood.sf)[1]
#' dis <- 1000000/matrix(as.numeric(st_distance(FastFood.sf,FastFood.sf)),ncol=n,nrow=n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.005))*dis
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' co <- sf::st_coordinates(st_centroid(FastFood.sf))
#' FastFood.sf$xf <- dgp_spq(x = co, p = p, listw = dis, rho = rho)
#' plot(FastFood.sf["xf"])
#' formula <- ~ xf
#' srq <- srq_test(formula = formula, data = FastFood.sf, listw = dis)
#' srq$SRQ
#' plot.srq(srq)
#'
#'
#' # SRQ test based on inverse distance
#' data("FastFood")
#' n = dim(FastFood.sf)[1]
#' dis <- 1000000/matrix(as.numeric(st_distance(FastFood.sf,FastFood.sf)),ncol=n,nrow=n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.005))*dis
#' formula <- ~ Type
#' srq <- srq_test(formula = formula, data = FastFood.sf, listw = dis)
#' srq$SRQ
#' plot.srq(srq)
#'
#'

srq_test <-  function(formula = NULL, data = NULL, na.action, xf = NULL,
                       listw = listw, alternative = "two-sided" , regular = FALSE, nsim = NULL, control = list()){

  # Solo admite matrices knn o nb
  if (class(listw)[1] != "knn"){
    if (class(listw)[1] != "nb"){
      if (class(listw)[1] != "matrix"){
      stop ("tiene que ser de la clase knn, nb o matrix")
      }
    }
  }
################################################
## Tratamiento de la matrix
################################################
# Si se trata de un objeto sf y la matrix es tipo 'nb' / "matrix" hay que ordenar los m_i-entornos
if (sum(class(data)[1]=="sf")==1){
if (class(listw)[1]=='nb'){ # hay que ordenar los elementos
  listw <- nb2nb_order(listw=listw, sf = data)
}
if (class(listw)[1]=='matrix'){ # hay que ordenar los elementos
      listw <- mat2listw(listw)$neighbours
      class(listw) <- "nb"
      listw2 <- nb2nb_order(listw=listw, sf = data)
}
}
################################################
## Tratamiento del input de los datos
################################################
# Selecciona los argumentos. Bien con (formula + data) o bien incluye la variable (xf)
  if (!is.null(formula) && !is.null(data)) {
    if (inherits(data, "Spatial")) data <- as(data, "sf")
    mxf <- get_all_vars(formula, data)
  } else if (!is.null(xf)) {
    mxf <- xf
    # if (!is.matrix(mxf)) mxf <- as.matrix(mxf, ncol = 1)
    mxf <- as.data.frame(mxf)
    for (i in 1:ncol(mxf)) {
      if (!is.factor(mxf[,i])) mxf[,i] <- as.factor(mxf[,i])
    }
  } else stop("data wrong")

  # xf debe ser un factor. Lo transformo en var numerica para calcular
  if (is.factor(mxf[,1])){
    levels(mxf[,1]) <- as.character(1:length(levels(mxf[,1])))
    y <- as.numeric(mxf[,1])
  }
  if (is.character(mxf[,1])){
    y <- as.factor(mxf[,1])
    levels(fx[,1]) <- as.character(1:length(levels(mxf[,1])))
    y <- as.numeric(mxf[,1])
  }
  if (is.numeric(mxf[,1])){
    stop("Only factors are admitted")
  }

################################################
## Empezamos las cuentas del test
################################################

# Calculo valores previos para obtener media y varianza estadístico
nv <- creation_nvar_SR(listw = listw)
q <- max(y)
n <- length(y)
# Cont is a binary variable that takes on the value of 1 if data are
# continuous and 0 if data are categorical.

m <- numeric() # matrix(0,nrow=q,ncol=1)
pprod <- numeric() # matrix(0,nrow=q,ncol=1)
for (i in 1:q){
  m[i] <- sum(y==i)
  pprod[i]<- m[i]*(n-m[i])
}
if (class(listw)[1]=="knn"){
lnnb <- matrix(dim(listw$nn)[2],ncol = 1,nrow = dim(listw$nn)[1])}
if (class(listw)[1]=="nb"){
lnnb <- rowSums(nb2mat(listw,style='B',zero.policy = TRUE))
}
MaxNeig <- max(lnnb)+1 # El elemento en cuestion + sus vecinos

# here we categorize the original data set y into the q categories
# compute the m_k needed for the computation of mean and variance
# pprod is needed for the computation of p
p <- sum(pprod)/(n*(n-1))


##### COMPUTING THE VARIANCE #####
## case 1 ##
aux1 <- numeric()
aux31 <- numeric()
aux3 <- numeric()

t1=0;
for (k in 1:q){
  for (c in 1:q){
    t1=t1+1
    aux1[t1]=m[k]*m[c]*(n-m[c]-1)
    aux31[t1]=m[k]*m[c]*((m[k]-1)*(n-m[k]-1)+(m[c]-1)*(n-m[c]-1))
    if(k==c){
      aux1[t1]=0
      aux31[t1]=0
    }
  }
}

t3=0
aux3<-numeric()
for (k in 1:q){
  for (c in 1:q){
    for (d in 1:q){
      t3=t3+1
      aux3[t3]=m[k]*m[c]*m[d]*(n-m[d]-2)
      if (c==k){aux3[t3]=0}
      if (d==k){aux3[t3]=0}
      if (d==c){aux3[t3]=0}
    }
  }
}

var1=1/(n*(n-1)*(n-2)*(n-3))*(sum(aux3)+sum(aux31));
var2=1/(n*(n-1)*(n-2))*sum(aux1);
var3=p;

varSR=p*(1-p)*sum(lnnb)+nv[1]*var1+nv[2]*var2+nv[3]*var3-(nv[1]+nv[2]+nv[3])*p^2

# Here we compute the runs starting at each location and it sum is the total number of runs
nruns <- matrix(0,ncol = 1,nrow = n)
for (i in 1:n){
  if (lnnb[i]!= 0){ # Solo calcula los test locales si el elemento tiene vecinos
  if (class(listw)[1]=="knn"){
    runs <- y[c(i,listw$nn[i,])]}
  if (class(listw)[1]=="nb"){
    runs <- y[c(i,listw[[i]])]}
nruns[i] <- 1 + sum(abs(diff(runs))>0)
  }
}

# La distribución del numero de rachas
dnr <- table(factor(nruns, levels = c(1:MaxNeig))) # dnr <- table(SRQlocal[,1],exclude = 0) # Excluimos localizaciones con 0 vecinos

SR <- sum(nruns)

# The mean of the statistic
meanSR <- n+p*sum(lnnb)

# The SRQ global test statistic which is N(0,1) distributed
SRQ <- (SR-meanSR)/sqrt(varSR)
if (alternative =="two-sided"){
p.valueSRQ <- 2*(1-pnorm(abs(SRQ), mean = 0, sd = 1))
} else if (alternative =="less"){
p.valueSRQ <- pnorm(SRQ, mean = 0, sd = 1)
} else if (alternative =="greater"){
p.valueSRQ <- 1-pnorm(SRQ, mean = 0, sd = 1)
}
############################################################################
# Para la obtención de los intervalos de confianza por boots permutacional
############################################################################
if (is.null(nsim) == FALSE){
  set.seed(123)
SRGP <- matrix(0,ncol = 1,nrow = nsim)
SRLP <- matrix(0,ncol = nsim, nrow = n)
    for (i in 1:nsim){
    yp <- y[sample(1:n)]
    srqp <- SR_test_boots(xf = yp, listw = listw, nv = nv)
    SRGP[i] <- srqp$SRglobal
    SRLP[,i] <- srqp$nruns
    }
if (alternative =="greater"){
p.valueSRQB <- sum(SRGP > SRQ)/(nsim+1)
}
else if (alternative =="less"){
  p.valueSRQB <- sum(SRGP < SRQ)/(nsim+1)
}
else if (alternative =="two-sided"){
  p.valueSRQB <- (sum(SRGP < SRQ) + sum(SRGP > SRQ))/(nsim+1)
} else stop("alternative wrong")
}

# El número total de rachas es la suma de las rachas en cada m_i-entorno
# colSums(SRLP) da el número de rachas de las nsim repeticiones
# SRGP es un vector con los valores de SRGlobal bajo aleatoriedad
############################################################
# Salida en función si se piden intervalos IC o no
if (is.null(nsim) == FALSE){
return <- list(SR = SR, dnr = dnr, SRQ  = SRQ, p.valueSRQ = p.valueSRQ, alternative = alternative,
               SRGP = SRGP, p.valueSRQB = p.valueSRQB,
               nsim = nsim, SRLP = SRLP, MeanNeig=sum(lnnb)/n,
               listw = listw, MaxNeig = MaxNeig)
}
else
{
return <- list(SR = SR, dnr = dnr, SRQ = SRQ, p.valueSRQ = p.valueSRQ, alternative = alternative,
               MeanNeig=sum(lnnb)/n,listw = listw, MaxNeig = MaxNeig)
}
}

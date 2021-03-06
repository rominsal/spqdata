#'
#' @title A function to calculate the local spatial runs tests.
#'
#' @description This function calculates the local spatial runs tests for each localization.
#' @param formula An (optional) formula with the factor included in \code{data}
#' @param data An (optional) data frame or a sf object containing the variable to testing for.
#' @param xf An (optional) factor of observations with the same length as the neighbors list in \code{listw}
#' @param listw A neighbours list of the class knn or nb. Alternatively a matrix W que indique el orden de cada $m_i-entorno$
#' (por ejemplo de inversa distancia). Para calcular el numero de rachas en cada $m_i-entorno$ debe establecerse
#' un orden, por ejemplo del vecino más próximo al más lejano.
#' @param alternative a character string specifying the alternative hypothesis, must be one
#' of "two.sided" (default), "greater" or "less".
#' @param distr a character string specifying the distribution "asymptotic" (default) or "bootstrap"
#' @param nsim Default value is NULL to obtain the asymptotic version of the local test.
#' To obtain the boots version the number of permutations to obtain the statistics.
#' @usage local.sp.runs.test(formula = NULL, data = NULL, na.action, xf = NULL,
#' distr = "asymptotic",listw = listw, alternative = "two.sided" , nsim = NULL,
#' control = list())
#' @keywords local spatial association, qualitative variable, spatial runs test
#' @details The object \code{listw} can be the class:
#'   \tabular{ll}{
#'     \code{knn} \tab Un objeto tipo knn obtenida utilizando el criterio del vecino más próximo,
#'      obtenida por ejemplo con \code{\link{knearneigh}} \cr
#'     \code{nb} \tab Un objeto tipo nb obtenido con \code{\link{spdep::poly2nb}}.\cr
#'     \code{matrix} \tab Una matriz indicando el orden de vecindad de cada $m_i-entorno$.
#'      Por ejemplo basada en inversa de la distancia con cierto punto de corte. \cr
#'     }
#' @return The output is an object of the class localsrq \cr
#' \cr
#'   \code{local.SRQ} A matrix with \cr
#'   \tabular{ll}{
#'     \code{runs.i} \tab número de rachas en la localización i. \cr
#'     \code{E.i} \tab expectation of local local runs statistic. \cr
#'     \code{Sd.i} \tab standard deviate of local runs statistic. \cr
#'     \code{z.value} \tab standard value of local runs statistic (only for asymptotic version). \cr
#'     \code{p.value} \tab p-value of local local runs statistic (only for asymptotic version). \cr
#'     \code{zseudo.value} \tab standard value of local runs statistic (only for boots version). \cr
#'     \code{pseudo.value} \tab p-value of local local runs statistic (only for boots version). \cr
#'     }
#'     \code{MeanNeig} Mean of run.i \cr
#'     \code{MaxNeig} Max of run.i  \cr
#'     \code{listw} the object \code{listw} \cr
#'     \code{alternative} a character string describing the alternative hypothesis \cr
#'
#' @section Control arguments:
#'   \tabular{ll}{
#'     \code{seedinit} \tab Numerical value for the seed in boot version. Default value seedinit = 123 \cr
#'       }
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paezha@@gmail.com} \cr
#'   Manuel Ruiz \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'   @references
#'   \itemize{
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     A test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#' @seealso
#' \code{\link{sp.runs.test}}, \code{\link{dgp.spq}}
#' @export
#' @examples
#'
#' # Case 1: Local SRQ test based on knn
#' rm(list = ls())
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- knearneigh(cbind(cx,cy), k = 10)
#' p <- c(1/6,3/6,2/6)
#' rho <- 0.5
#' xf <- dgp.spq(p = p, listw = listw, rho = rho)
#' # Asymtotic version
#' lsrq <- local.sp.runs.test(xf = xf, listw = listw, alternative = "less")
#' print(lsrq)
#' plot(lsrq, sig = 0.05)
#' # Asymtotic version
#' lsrq <- local.sp.runs.test(xf = xf, listw = listw, alternative = "two.sided", distr ="bootstrap", nsim = 399)
#' print(lsrq)
#' plot(lsrq, sig = 0.1)
#'
#' # Case 2:Fastfood example. sf (points)
#' rm(list = ls())
#' data("FastFood")
#' x <- cbind(FastFood.sf$Lat,FastFood.sf$Lon)
#' listw <- spdep::knearneigh(x, k = 10)
#' formula <- ~ Type
#' lsrq <- local.sp.runs.test(formula = formula, data = FastFood.sf, listw = listw)
#' print(lsrq)
#' plot(lsrq, sf = FastFood.sf, sig = 0.05)
#'
#' # Case 3: With a sf object (poligons)
#' rm(list = ls())
#' fname <- system.file("shape/nc.shp", package="sf")
#' nc <- st_read(fname)
#' listw <- spdep::poly2nb(as(nc,"Spatial"), queen = FALSE)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' nc$xf <- dgp.spq(p = p, listw = listw, rho = rho)
#' plot(nc["xf"])
#' formula <- ~ xf
#' lsrq <- local.sp.runs.test(formula = formula, data = nc, listw = listw)
#' print(lsrq)
#' plot(lsrq, sf = nc)
#' # Version boot
#' lsrq <- local.sp.runs.test(formula = formula, data = nc, listw = listw, distr ="bootstrap", nsim = 399)
#' print(lsrq)
#' plot(lsrq, sf = nc)
#'
#' # Case 4: With isolated areas
#' rm(list = ls())
#' data(Spain)
#' listw <- spdep::poly2nb(as(spain.sf,"Spatial"), queen = FALSE)
#' plot(spain.sf["MenWoman"])
#' formula <- ~ MenWoman
#' lsrq <- local.sp.runs.test(formula = formula, data = spain.sf, listw = listw)
#' print(lsrq)
#' plot(lsrq, sf = spain.sf, sig = 0.1)
#' # Boots Version
#' lsrq <- local.sp.runs.test(formula = formula, data = spain.sf, listw = listw, distr ="bootstrap", nsim = 199)
#' print(lsrq)
#' plot(lsrq, sf = spain.sf, sig = 0.10)
#'
#' # Case 5: SRQ test based on a distance matrix (inverse distance)
#' rm(list = ls())
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' coor <- as.data.frame(cbind(cx,cy))
#' coor <- st_as_sf(coor,coords = c("cx","cy"))
#' n = dim(coor)[1]
#' dis <- 1/matrix(as.numeric(st_distance(coor,coor)), ncol = n, nrow = n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.10))*dis
#' p <- c(1/6,3/6,2/6)
#' rho <- 0.5
#' xf <- dgp.spq(p = p, listw = dis, rho = rho)
#' lsrq <- local.sp.runs.test(xf = xf, listw = dis)
#' print(lsrq)
#' plot(lsrq, coor = cbind(cx,cy), sig = 0.05)
#' lsrq <- local.sp.runs.test(xf = xf, listw = dis, data = )
#' print(lsrq)
#' plot(lsrq, sf = coor)
#' # Version boots
#' lsrq <- local.sp.runs.test(xf = xf, listw = dis, data = coor, distr ="bootstrap", nsim = 299)
#' print(lsrq)
#' plot(lsrq, sf = coor)
#'
#' # SRQ test based on inverse distance
#' rm(list = ls())
#' data("FastFood")
#' n = dim(FastFood.sf)[1]
#' dis <- 1000000/matrix(as.numeric(st_distance(FastFood.sf,FastFood.sf)),ncol=n,nrow=n)
#' diag(dis) <- 0
#' dis <- (dis < quantile(dis,.01))*dis
#' formula <- ~ Type
#' lsrq <- local.sp.runs.test(formula = formula, data = FastFood.sf, listw = dis)
#' print(lsrq)
#' plot(lsrq, sf = FastFood.sf)
#'


local.sp.runs.test <-  function(formula = NULL, data = NULL, na.action, xf = NULL, distr = "asymptotic",
                       listw = listw, alternative = "two.sided" , nsim = NULL, control = list()){

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  distr <- match.arg(distr, c("asymptotic", "bootstrap"))

  ################################################
  # Solo admite matrices knn, nb o matrix
  if (class(listw)[1] != "knn"){
    if (class(listw)[1] != "nb"){
      if (class(listw)[1] != "matrix"){
        stop ("listw must be is a knn, nb o matrix class element")
      }
    }
  }
  ################################################
  # Controls
  con <- list(seedinit = 123)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  seedinit <- con$seedinit
  ################################################

  ## LOS M_I-ENTORNOS DEBEN ESTAR ORDENADOS POR ORDEN DE PROXIMIDAD
# Si se trata de un objeto sf y la matrix es tipo 'nb' hay que ordenar los m_i-entornos

  if (sum(class(data)[1]=="sf")==1){
    if (class(listw)[1]=='nb'){ # hay que ordenar los elementos
      listw <- nb2nb_order(listw=listw, sf = data)
    }
    if (class(listw)[1]=='matrix'){ # hay que ordenar los elementos
      listw <- mat2listw(listw)$neighbours
      class(listw) <- "nb"
      listw <- nb2nb_order(listw=listw, sf = data)
    }
  }

  if (sum(class(data)[1]=="sf")==0){
    if (class(listw)[1]=='matrix'){ # hay que ordenar los elementos
      listw <- nb2nb_order(listw=listw)
    }
  }

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

#### EMPEZAMOS LAS CUENTAS

# Calculo valores previos para obtener media y varianza estadístico
nv <- creation_nvar_SR(listw = listw)
q <- max(y)
n <- length(y)
# Cont is a binary variable that takes on the value of 1 if data are
# continuous and 0 if data are categorical.

m <- numeric()
pprod <- numeric()
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

############################################################################
# Here we compute the runs starting at each location and it sum is the total number of runs
############################################################################
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

MeanR <- 1 + lnnb*p
StdR <- sqrt(lnnb*p*(1-p)+2*(lnnb-1)*(var2-p^2)+(lnnb-1)*(lnnb-2)*(var1-p^2))
ZZ <- (nruns-MeanR)/StdR
if (alternative =="two.sided"){
pZ <- 2*(1-pnorm(abs(ZZ), mean = 0, sd = 1))
} else if (alternative =="less"){
pZ <- pnorm(ZZ, mean = 0, sd = 1)
} else if (alternative =="greater"){
pZ <- 1-pnorm(ZZ, mean = 0, sd = 1)
}

local.SRQ <- cbind(nruns,MeanR,StdR,ZZ,pZ)

local.SRQ  <- as.data.frame(local.SRQ)
names(local.SRQ) <- c("runs.i","E.i","Std.i","z.value","p.value")

# # El test de rachas da NaN en caso de una sola racha. Pongo Z=99
# # OJO VER QUE PASA CON RACHAS CORTAS EN HEXAGONOS
# local.SRQ[is.na(local.SRQ[,4]),4] <- 99

############################################################################
# Para la obtención de los intervalos de confianza por boots permutacional
############################################################################
if (is.null(nsim) == FALSE && (distr != "asymptotic")){
  if (!is.null(seedinit)) set.seed(seedinit)
  LSRQP <- matrix(0,ncol = nsim, nrow = n)
    for (i in 1:nsim){
    yp <- y[sample(1:n)]
    lsrqp <- local.sp.runs.test.boots(xf = yp, listw = listw, nv = nv)
    LSRQP[,i] <- lsrqp$nruns
    }

  local.SRQP <- as.data.frame(local.SRQ[,1])
  names(local.SRQP) <- "SRQ"
  local.SRQP$EP.i <- rowMeans(LSRQP)
  local.SRQP$SdP.i <- apply(LSRQP,1, sd, na.rm = TRUE)

  local.SRQP$zseudo.value <- (local.SRQP$SRQ - local.SRQP$EP.i)/local.SRQP$SdP.i
  if (alternative =="two.sided"){
    pZ <- 2*(1-pnorm(abs(local.SRQP$zseudo.value), mean = 0, sd = 1))
  } else if (alternative =="less"){
    pZ <- pnorm(local.SRQP$zseudo.value, mean = 0, sd = 1)
  } else if (alternative =="greater"){
    pZ <- 1-pnorm(local.SRQP$zseudo.value, mean = 0, sd = 1)
  }
  local.SRQP$pseudo.value <- pZ
}

############################################################
# Salida en función si se piden intervalos IC o no
############################################################
if (is.null(nsim) == FALSE){
local <- list(local.SRQP = local.SRQP, LSRQP = LSRQP, nsim = nsim , MeanNeig=sum(lnnb)/n,
               listw = listw, MaxNeig = MaxNeig, alternative = alternative)
}
else
{
local <- list(local.SRQ = local.SRQ, MeanNeig = sum(lnnb)/n, listw = listw,
              MaxNeig = MaxNeig, alternative = alternative)
}
class(local) <- "localsrq"
return <- local
}


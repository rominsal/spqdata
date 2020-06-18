#' @title A funcion to calculate runs test
#'
#' @description This function calculates the runs test for spatial independence
#' @param fx vector de observaciones Factor
#' @param W matrix W
#' @usage SR_test(fx=fx,W=W)
#' @keywords spatial association, qualitative variable, run test
#' @details Aquí Antonio escribe una linda historia
#' @return decir que cosas son las que devuelve
#'   \tabular{ll}{
#'     \code{pc} \tab Qp for combinations-totals symbols \cr
#'     \code{pc_pval} \tab  p-value for Qp \cr
#'     \code{qc} \tab Qc for combinations-totals symbols \cr
#'     \code{qc_pval} \tab  p-value for Qc \cr
#'     \code{p_symb} \tab  Matriz que lista los símbolos sin compactar \cr
#'     \code{efp_symb} \tab  Frecuancia absoluta de cada símbolo (no compactados)\cr
#'     \code{c_symb} \tab  Matriz que lista los símbolos COMPACTOS \cr
#'     \code{efc_symb} \tab  Frecuancia absoluta de cada símbolo (Compactado)\cr
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
#'     \item Ruiz, M., López, F., and Páez, A. (2010).
#'     TA test for global and local homogeneity of categorical data based on spatial runs.
#'       \emph{Geographical Analysis}.
#'   }
#' @seealso
#' \code{\link{dgp_spq}}, \code{\link{m_surr_no}}
#' @export
#' @examples
#'
#' N <- 100
#' cx <- runif(N)
#' cy <- runif(N)
#' x <- cbind(cx,cy)
#' listw <- spdep::nb2listw(knn2nb(knearneigh(cbind(cx,cy), k=4)))
#' listw <- knearneigh(cbind(cx,cy), k=4)
#' p <- c(1/6,3/6,2/6)
#' rho = 0.5
#' fx <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
#' SR_test(fx=fx,listw = listw)

SR_test <- function(fx=fx,listw=listw){
####

  nv <- creation_nvar_SR(listw=listw)

  if (is.factor(fx)){
    levels(fx) <- as.character(1:length(levels(fx)))
    y <- as.numeric(fx)
  }
  if (is.character(fx)){
    y <- as.factor(fx)
    levels(fx) <- as.character(1:length(levels(fx)))
    y <- as.numeric(fx)
  }
  if (is.numeric(fx)){
    y <- fx
  }

q <- max(y)
n <- length(y)
# Cont is a binary variable that takes on the value of 1 if data are
# continuous and 0 if data are categorical.

m <- numeric() # matrix(0,nrow=q,ncol=1)
pprod <- numeric() # matrix(0,nrow=q,ncol=1)
for (i in 1:q){
  m[i]<- sum(y==i)
  pprod[i]<- m[i]*(n-m[i])
}

if (class(listw)=='knn'){
  lnnb <- matrix(dim(listw$nn)[2],ncol = 1,nrow = dim(listw$nn)[1])
  NNB <- cbind(1:R,listw$nn)
}

# if (class(W)=='matrix'){
# mlnnb <- max(rowSums(W>0))+1
# NNB <- matrix(0,ncol = (mlnnb+1),nrow = n)
# lnnb <- matrix(0,ncol = 1,nrow = n)
# for (i in 1:n) {
#   tmp=sort(W[i,],index.return=TRUE)
#   a <- (tmp$x>0)*(tmp$ix)
#   a <- a[a>0]
#   lnnb[i] <- length(a)
#   NNB[i,] <- c(i,a,repmat(-99,1,mlnnb-lnnb[i]))
# }
# }

# here we categorize the original data set y into the q categories
# compute the m_k needed for the computation of mean and variance
# pprod is needed for the computation of p

p=sum(pprod)/(n*(n-1))

# Here we compute the runs starting at each location and it sum is the total number of runs
nruns <- numeric()
SRlocal <- matrix(0,ncol = 4,nrow = n)
for (i in 1:n){
#   if (class(W)=='matrix'){
# tmp=sort(W[i,],index.return=TRUE)
# a <- (tmp$x>0)*(tmp$ix)
# a <- rev(a[a>0])
# runs <- y[c(i,a)]
# }
if (class(W)=='knn'){
    runs <- y[NNB[i,]]
}
# if (class(W)=='nb'){
#     runs <- y[W[[i]]]
# }
rrun <- abs(diff(runs))
nruns[i] <- 1+sum(abs(rrun)>0)
qqsi <- length(unique(runs))
bb <- runs_test0(runs,qqsi)
SRlocal[i,1] <- bb$R
SRlocal[i,2] <- bb$meanR
SRlocal[i,3] <- bb$stdR
SRlocal[i,4] <- bb$Z
}
# El test de rachas da NaN en caso de una sola racha. Pongo Z=99

# OJO VER QUE PASA CON RACHAS CORTAS EN HEXAGONOS
SRlocal[is.na(SRlocal[,4]),4] <- 99

SR=sum(nruns)
#The mean of the statistic
meanSR=n+p*sum(lnnb)

##### COMPUTING THE VARIANCE #####
## case 1 ##
aux2=m*(n-m)*(n-m-1)
aux31=m*(m-1)*(n-m)*(n-m-1)

sum1=sum(pprod)
sum2=sum(aux2)
sum3=2*sum(aux31)

t3=0
aux32<-numeric()
for (k in 1:q){
for (c in 1:q){
for (d in 1:q){
t3=t3+1
aux32[t3]=m[k]*m[c]*m[d]*(n-m[d]-2)
if (c==k){aux32[t3]=0}
if (d==k){aux32[t3]=0}
if (d==c){aux32[t3]=0}
}
}
}

sum32=sum(aux32)
var1=1/(n*(n-1)*(n-2)*(n-3))*(sum3+sum32)
var2=1/(n*(n-1)*(n-2))*(sum2)
var3=1/(n*(n-1))*(sum1)
# nv <- creation_nvar_SR(W)
varSR=p*(1-p)*sum(lnnb)+2*sum(lnnb-1)*(var2-p^2)+sum((lnnb-1)*(lnnb-2))*(var1-p^2)+
  nv[1]*var1+nv[2]*var2+nv[3]*var3-(nv[1]+nv[2]+nv[3])*p^2

# The SR global test statistic which is N(0,1) distributed
SRglobal=(SR-meanSR)/sqrt(varSR)
return <- list(SR=SR,SRglobal=SRglobal,SRlocal=SRlocal,MeanNeig=sum(lnnb)/n,W=W)
}

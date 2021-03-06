% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local.sp.runs.test.R
\name{local.sp.runs.test}
\alias{local.sp.runs.test}
\title{A function to calculate the local spatial runs tests.}
\usage{
local.sp.runs.test(formula = NULL, data = NULL, na.action, xf = NULL,
distr = "asymptotic",listw = listw, alternative = "two.sided" , nsim = NULL,
control = list())
}
\arguments{
\item{formula}{An (optional) formula with the factor included in \code{data}}

\item{data}{An (optional) data frame or a sf object containing the variable to testing for.}

\item{xf}{An (optional) factor of observations with the same length as the neighbors list in \code{listw}}

\item{distr}{a character string specifying the distribution "asymptotic" (default) or "bootstrap"}

\item{listw}{A neighbours list of the class knn or nb. Alternatively a matrix W que indique el orden de cada $m_i-entorno$
(por ejemplo de inversa distancia). Para calcular el numero de rachas en cada $m_i-entorno$ debe establecerse
un orden, por ejemplo del vecino más próximo al más lejano.}

\item{alternative}{a character string specifying the alternative hypothesis, must be one
of "two.sided" (default), "greater" or "less".}

\item{nsim}{Default value is NULL to obtain the asymptotic version of the local test.
To obtain the boots version the number of permutations to obtain the statistics.}
}
\value{
The output is an object of the class localsrq \cr
\cr
  \code{local.SRQ} A matrix with \cr
  \tabular{ll}{
    \code{runs.i} \tab número de rachas en la localización i. \cr
    \code{E.i} \tab expectation of local local runs statistic. \cr
    \code{Sd.i} \tab standard deviate of local runs statistic. \cr
    \code{z.value} \tab standard value of local runs statistic (only for asymptotic version). \cr
    \code{p.value} \tab p-value of local local runs statistic (only for asymptotic version). \cr
    \code{zseudo.value} \tab standard value of local runs statistic (only for boots version). \cr
    \code{pseudo.value} \tab p-value of local local runs statistic (only for boots version). \cr
    }
    \code{MeanNeig} Mean of run.i \cr
    \code{MaxNeig} Max of run.i  \cr
    \code{listw} the object \code{listw} \cr
    \code{alternative} a character string describing the alternative hypothesis \cr
}
\description{
This function calculates the local spatial runs tests for each localization.
}
\details{
The object \code{listw} can be the class:
  \tabular{ll}{
    \code{knn} \tab Un objeto tipo knn obtenida utilizando el criterio del vecino más próximo,
     obtenida por ejemplo con \code{\link{knearneigh}} \cr
    \code{nb} \tab Un objeto tipo nb obtenido con \code{\link{spdep::poly2nb}}.\cr
    \code{matrix} \tab Una matriz indicando el orden de vecindad de cada $m_i-entorno$.
     Por ejemplo basada en inversa de la distancia con cierto punto de corte. \cr
    }
}
\section{Control arguments}{

  \tabular{ll}{
    \code{seedinit} \tab Numerical value for the seed in boot version. Default value seedinit = 123 \cr
      }
}

\examples{

# Case 1: Local SRQ test based on knn
rm(list = ls())
N <- 100
cx <- runif(N)
cy <- runif(N)
x <- cbind(cx,cy)
listw <- knearneigh(cbind(cx,cy), k = 10)
p <- c(1/6,3/6,2/6)
rho <- 0.5
xf <- dgp.spq(p = p, listw = listw, rho = rho)
# Asymtotic version
lsrq <- local.sp.runs.test(xf = xf, listw = listw, alternative = "less")
print(lsrq)
plot(lsrq, sig = 0.05)
# Asymtotic version
lsrq <- local.sp.runs.test(xf = xf, listw = listw, alternative = "two.sided", distr ="bootstrap", nsim = 399)
print(lsrq)
plot(lsrq, sig = 0.1)

# Case 2:Fastfood example. sf (points)
rm(list = ls())
data("FastFood")
x <- cbind(FastFood.sf$Lat,FastFood.sf$Lon)
listw <- spdep::knearneigh(x, k = 10)
formula <- ~ Type
lsrq <- local.sp.runs.test(formula = formula, data = FastFood.sf, listw = listw)
print(lsrq)
plot(lsrq, sf = FastFood.sf, sig = 0.05)

# Case 3: With a sf object (poligons)
rm(list = ls())
fname <- system.file("shape/nc.shp", package="sf")
nc <- st_read(fname)
listw <- spdep::poly2nb(as(nc,"Spatial"), queen = FALSE)
p <- c(1/6,3/6,2/6)
rho = 0.5
nc$xf <- dgp.spq(p = p, listw = listw, rho = rho)
plot(nc["xf"])
formula <- ~ xf
lsrq <- local.sp.runs.test(formula = formula, data = nc, listw = listw)
print(lsrq)
plot(lsrq, sf = nc)
# Version boot
lsrq <- local.sp.runs.test(formula = formula, data = nc, listw = listw, distr ="bootstrap", nsim = 399)
print(lsrq)
plot(lsrq, sf = nc)

# Case 4: With isolated areas
rm(list = ls())
data(Spain)
listw <- spdep::poly2nb(as(spain.sf,"Spatial"), queen = FALSE)
plot(spain.sf["MenWoman"])
formula <- ~ MenWoman
lsrq <- local.sp.runs.test(formula = formula, data = spain.sf, listw = listw)
print(lsrq)
plot(lsrq, sf = spain.sf, sig = 0.1)
# Boots Version
lsrq <- local.sp.runs.test(formula = formula, data = spain.sf, listw = listw, distr ="bootstrap", nsim = 199)
print(lsrq)
plot(lsrq, sf = spain.sf, sig = 0.10)

# Case 5: SRQ test based on a distance matrix (inverse distance)
rm(list = ls())
N <- 100
cx <- runif(N)
cy <- runif(N)
coor <- as.data.frame(cbind(cx,cy))
coor <- st_as_sf(coor,coords = c("cx","cy"))
n = dim(coor)[1]
dis <- 1/matrix(as.numeric(st_distance(coor,coor)), ncol = n, nrow = n)
diag(dis) <- 0
dis <- (dis < quantile(dis,.10))*dis
p <- c(1/6,3/6,2/6)
rho <- 0.5
xf <- dgp.spq(p = p, listw = dis, rho = rho)
lsrq <- local.sp.runs.test(xf = xf, listw = dis)
print(lsrq)
plot(lsrq, coor = cbind(cx,cy), sig = 0.05)
lsrq <- local.sp.runs.test(xf = xf, listw = dis, data = )
print(lsrq)
plot(lsrq, sf = coor)
# Version boots
lsrq <- local.sp.runs.test(xf = xf, listw = dis, data = coor, distr ="bootstrap", nsim = 299)
print(lsrq)
plot(lsrq, sf = coor)

# SRQ test based on inverse distance
rm(list = ls())
data("FastFood")
n = dim(FastFood.sf)[1]
dis <- 1000000/matrix(as.numeric(st_distance(FastFood.sf,FastFood.sf)),ncol=n,nrow=n)
diag(dis) <- 0
dis <- (dis < quantile(dis,.01))*dis
formula <- ~ Type
lsrq <- local.sp.runs.test(formula = formula, data = FastFood.sf, listw = dis)
print(lsrq)
plot(lsrq, sf = FastFood.sf)

}
\seealso{
\code{\link{sp.runs.test}}, \code{\link{dgp.spq}}
}
\author{
\tabular{ll}{
  Fernando López  \tab \email{fernando.lopez@upct.es} \cr
  Román Mínguez  \tab \email{roman.minguez@uclm.es} \cr
  Antonio Páez \tab \email{paezha@gmail.com} \cr
  Manuel Ruiz \tab \email{manuel.ruiz@upct.es} \cr
  }
  @references
  \itemize{
    \item Ruiz, M., López, F., and Páez, A. (2010).
    A test for global and local homogeneity of categorical data based on spatial runs.
      \emph{Geographical Analysis}.
  }
}
\keyword{association,}
\keyword{local}
\keyword{qualitative}
\keyword{runs}
\keyword{spatial}
\keyword{test}
\keyword{variable,}

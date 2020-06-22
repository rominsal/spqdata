#' @name plot.lsrq
#' @rdname plot.lsrq
#'
#' @title Plot the empirical distribution of runs.
#'
#' @param x object of class \emph{spq}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Antonio Páez \tab \email{paez@@gmail.com} \cr
#'   Manuel Ruiz  \tab \email{manuel.ruiz@@upct.es} \cr
#'   }
#'
#' @seealso
#' \code{\link{srq_test}}.
#'
#'
#' @examples
#' # Fastfood example. sf (points)
#' data("FastFood")
#' x <- cbind(FastFood.sf$Lon,FastFood.sf$Lat)
#' listw <- spdep::knearneigh(x, k = 2)
#' formula <- ~ Type
#' lsrq <- local_srq_test(formula = formula, data = FastFood.sf, listw = listw)
#' plot.lsrq(lsrq,sf=FastFood.sf)
#'
#' # Spain example (poligons with 0 neinghbourhood)
#' data("Spain")
#' listw <- spdep::poly2nb(as(spain.sf,"Spatial"), queen = FALSE)
#' formula <- ~ HOM_MUJ
#' plot(spain.sf["HOM_MUJ"])
#' lsrq <- local_srq_test(formula = formula, data = spain.sf, listw = listw)
#' plot.lsrq(lsrq = lsrq, sf = spain.sf)
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
#' lsrq <- local_srq_test(formula = formula, data = nc, listw = listw)
#' lsrq$SRQlocal
#' plot.lsrq(lsrq,sf=nc)

#' @export
#'

plot.lsrq <- function(lsrq = lsrq, sf = NULL, coor = NULL,  sig = 0.05 ){

#####################
### Plot Q Local
#####################

  if (!is.null(sf)){
    if (is.null(lsrq$nsim)){
# if (sum(is.na(lsrq$SRQlocal$`z-value`))==0){
a <- as.factor((lsrq$SRQlocal$p.menor < sig)*1+(lsrq$SRQlocal$p.mayor < sig)*2)
sf$plot <- addNA(a)
levels(sf$plot) <- c("no-sig","sig +","sig -","NA")
ggplot(sf) +
  geom_sf(aes(fill = plot), color = "black")+
  theme_bw() +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank()) +
  xlab(paste0("Significance p-value = ", sig)) +
  scale_fill_manual(values = c("white", "red","blue","gray"))

# }
    }
  }

  # if (!is.null(sf)){
  #   if (is.null(lsrq$nsim)){
  #     # if (sum(is.na(lsrq$SRQlocal$`z-value`))==0){
  #     a <- as.factor((lsrq$SRQlocal$p.menor < sig)*1+(lsrq$SRQlocal$p.mayor < sig)*2)
  #     sf$Signif <- addNA(a)
  #     mylevel <- levels(sf$Signif)
  #     mylevel[mylevel==0]="no-sig"
  #     mylevel[mylevel==1]="sig +"
  #     mylevel[mylevel==2]="sig -"
  #     mylevel[mylevel=="NA"]="NA"
  #     mycolor[a==0]="gray"
  #     mycolor[a==1]="red"
  #     mycolor[a==2]="blue"
  #     sf$Signif <- mycolor
  #     levels(sf$Signif) <- mylevel
  #     ggplot(sf) +
  #     geom_sf(aes(fill = Signif, color = Signif), size=3,color = mycolor) +
  #     theme_bw() +
  #     theme(axis.text.x=element_blank(),axis.text.y=element_blank()) +
  #     xlab(paste0("Significance p-value = ", sig)) +
  #     scale_fill_manual(values = c("blue","gray"))
  #
  #     # }
  #   }}
}






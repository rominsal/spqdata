B <- read.table("/Users/fernandoair/Downloads/20200601_maestra_1_mitma_municipio.txt",header = TRUE,sep="|")
BC <- B[B$origen=="01001_AM",]






N <- 1000
cx <- runif(N)
cy <- runif(N)
x <- cbind(cx,cy)
listw <- knearneigh(cbind(cx,cy), k=3)
p <- c(1/6,3/6,2/6)
rho <- -0.9
xf <- dgp_spq(x = x, p = p, listw = listw, rho = rho)
srq <- srq_test(xf = xf, listw = listw)

W <- nb2mat(knn2nb(listw))
g1 = igraph::graph.adjacency(W)

plot(g1,layout=layout.norm(x),edge.arrow.mode=0,edge.width=0,vertex.size=2,vertex.label="",vertex.color="red",vertex.label.font=0)


sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
hexs <- st_make_grid(sfc, cellsize = .1, square = FALSE)
hexs.sf <- st_sf(hexs)
plot(hexs.sf)
listw <- spdep::poly2nb(as(hexs.sf,"Spatial"), queen = FALSE)
# Antes
co <- sf::st_coordinates(st_centroid(hexs.sf))
W <- poly2nb(as(hexs.sf, "Spatial"), queen = FALSE)
W <- hex_ord(co=co,W=W)
# ahora


listw2 <- nb2nb_order(listw,sf=hexs.sf)
co <- sf::st_coordinates(st_centroid(hexs.sf))
listw3 <- hex_ord(co=co,W=listw)
mm <- matrix(0,ncol=1,nrow=137)
mm[listw3[[32]]] <- c(1:6)
hexs.sf$mm <- mm
plot(hexs.sf["mm"])


sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))))
hexs <- st_make_grid(sfc, cellsize = RR[r], square = FALSE)
hexs.sf <- st_sf(hexs)
W <- poly2nb(as(hexs.sf, "Spatial"), queen = FALSE)
co <- sf::st_coordinates(st_centroid(hexs.sf))
R <- length(W)
W <- hex_ord(co=co,W=W)
nv <- creation_nvar_SR_hex(W)

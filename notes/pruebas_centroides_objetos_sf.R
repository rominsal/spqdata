# sfc POINT
point1 = st_point(c(5, 2))
point2 = st_point(c(1, 3))
points_sfc = st_sfc(point1, point2)
points_sfc
plot(points_sfc)
coorpoints_sfc <- sf::st_coordinates(points_sfc)
coorpoints_sfc
ctpoints_sfc <- sf::st_centroid(points_sfc, of_largest_polygon = TRUE)
ctpoints_sfc
plot(ctpoints_sfc)



















multilinestring_list = list(rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2)),
                            rbind(c(1, 2), c(2, 4)))
pr1 <- sf::st_multilinestring((multilinestring_list))
plot(pr1)
coorpr1 <- sf::st_coordinates(pr1)
ctpr1 <- sf::st_centroid(pr1)
plot(ctpr1)

multipoint_matrix = rbind(c(5, 2), c(1, 3), c(3, 4), c(3, 2))
pr2 <- sf::st_multipoint(multipoint_matrix)
plot(pr2)
coorpr2 <- sf::st_coordinates(pr3)
coorpr2
ctpr2 <- sf::st_centroid(pr2)
plot(ctpr2)

linestring_matrix = rbind(c(1, 5), c(4, 4), c(4, 1), c(2, 2), c(3, 2))
pr3 <- st_linestring(linestring_matrix)
plot(pr3)
coorpr3 <- sf::st_coordinates(pr3)
ctpr3 <- sf::st_centroid(pr3)
plot(ctpr3)


polygon_list = list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5)))
pr4 <- st_polygon(polygon_list)
plot(pr4)
coorpr4 <- sf::st_coordinates(pr4)
coorpr4
ctpr4 <- sf::st_centroid(pr4)
plot(ctpr4)
sf::st_coordinates(ctpr4)

multipolygon_list = list(list(rbind(c(1, 5), c(2, 2), c(4, 1), c(4, 4), c(1, 5))),
                         list(rbind(c(0, 2), c(1, 2), c(1, 3), c(0, 3), c(0, 2))))
pr5 <- st_multipolygon(multipolygon_list)
plot(pr5)
coorpr5 <- sf::st_coordinates(pr5)
coorpr5
ctpr5 <- sf::st_centroid(pr5)
plot(ctpr5)
sf::st_coordinates(ctpr5)

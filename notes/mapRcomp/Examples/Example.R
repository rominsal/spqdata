Newark <- Newark
library(ggplot2)

coords <- cbind(Newark$X_UTM2,Newark$Y_UTM2)
Y1 <- 1*Newark$NW + 2*Newark$IRISH + 3*Newark$GERMAN;#Ethnicity only
ETHNICITY <- factor(Y1,labels = c("Yankee","Irish","German"))

Newark <- cbind(Newark,ETHNICITY)

ggplot(Newark, aes(x = X_UTM2, y = Y_UTM2)) + geom_point(aes(color=ETHNICITY))

mh51 <- m_surr_no(coords,5,1)
symb35 <- cr_symb(3,5)
results.35 <- q_symb(Y1,mh51,symb35)
results.35$c_plot


rm(list = ls())

library(readxl)
Synthetic_Maps <- read_excel("D:/Antonio/Cases 2017/mapRcomp/Examples/Synthetic Maps.xlsx")
View(Synthetic_Maps)

coords <- cbind(Synthetic_Maps$XCOORD,Synthetic_Maps$YCOORD)

ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y1))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y2))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y3))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y4))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y6))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y7))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y8))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y9))) + geom_point(shape = 15, size = 3.5)
ggplot(Synthetic_Maps,aes(x = XCOORD, y = YCOORD, color = as.factor(Y10))) + geom_point(shape = 15, size = 3.5)


mh41 <- m_surr_no(coords,4,1)
symb34 <- cr_symb(3,4)

map1.35 <- q_symb(Synthetic_Maps$Y1,mh41,symb34)
map1.35$c_plot

map2.35 <- q_symb(Synthetic_Maps$Y2,mh41,symb34)
map2.35$c_plot

map3.35 <- q_symb(Synthetic_Maps$Y3,mh41,symb34)
map3.35$c_plot

map4.35 <- q_symb(Synthetic_Maps$Y4,mh41,symb34)
map4.35$c_plot

map6.35 <- q_symb(Synthetic_Maps$Y6,mh41,symb34)
map6.35$c_plot

map7.35 <- q_symb(Synthetic_Maps$Y7,mh41,symb34)
map7.35$c_plot

map8.35 <- q_symb(Synthetic_Maps$Y8,mh41,symb34)
map8.35$c_plot

map9.35 <- q_symb(Synthetic_Maps$Y9,mh41,symb34)
map9.35$c_plot

map10.35 <- q_symb(Synthetic_Maps$Y10,mh41,symb34)
map10.35$c_plot

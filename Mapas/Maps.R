
## Configuração Inicial
require(tikzDevice)
require(RColorBrewer) # Para Paletas de cores
require(maptools)
library(rgdal)
library(raster)
require(SDMTools) # para legend.gradient
<<<<<<< HEAD
library(extrafont)
loadfonts()
=======
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

getwd()
setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Mapas")
## rm(list = ls()) ## Para remover todos os objetos

options(OutDec=",",digits=10,scipen=5)
<<<<<<< HEAD
par(mar = c(4,4,0.1,0.1))
plot(1)
=======
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a
##-----------------------------------------------------------------------------##
## Funções para colocar Norte e legenda nos mapas

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
# checking arguments
if(missing(loc)) stop("loc is missing")
if(missing(size)) stop("size is missing")
# default colors are white and black
if(missing(cols)) cols <- rep(c("white","black"),8)
# calculating coordinates of polygons
radii <- rep(size/c(1,4,2,4),4)
x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
# drawing polygons
for (i in 1:15) {
x1 <- c(x[i],x[i+1],loc[1])
y1 <- c(y[i],y[i+1],loc[2])
polygon(x1,y1,col=cols[i])
}
# drawing the last polygon
polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
# drawing letters
b <- c("E","N","W","S")
for (i in 0:3) text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
(size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
cex=cex)
}


scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
if(missing(loc)) stop("loc is missing")
if(missing(length)) stop("length is missing")
x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
y <- c(0,length/(10*3:1))+loc[2]
cols <- rep(c("black","white"),2)
for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
for (i in 1:5) segments(x[i],y[2],x[i],y[3])
labels <- x[c(1,3)]-loc[1]
labels <- append(labels,paste(x[5]-loc[1],unit))
text(x[c(1,3,5)],y[4],labels=labels,adj=.5,cex=division.cex)
}

##-----------------------------------------------------------------------------##

## Importando dados
dir()

## Projeção
proj <- CRS("+proj=utm +zone=22 +south +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 +units=m +no_defs")
## lendo dados tipo shapefile
limite <- readShapeSpatial("Shapes/SC_limite.shp",proj4string = proj)
FluviEST <- readShapeSpatial("Shapes/Fluvi_Estat.shp",proj4string = proj)
PluviEST <- readShapeSpatial("Shapes/Pluvi_Estat.shp",proj4string = proj)
Hidro <- readShapeSpatial("Shapes/Hidro.shp",proj4string = proj)
Bacias <- readShapeSpatial("Shapes/Bacias/merged/merged.shp",proj4string = proj)
Centroids <- readShapeSpatial("Shapes/Centroids.shp",proj4string = proj)
MBacias <- readShapeSpatial("Shapes/Macro_Bacias.shp",proj4string = proj)

## lendo dados tipo Raster
Cota <- raster("Rasters/cotas.tif")

PluviANO <- mask(raster("Rasters/PluviANO.ascii",crs=proj),limite)
PluviDJF <- mask(raster("Rasters/PluviDJF.ascii",crs=proj),limite)
PluviMAM <- mask(raster("Rasters/PluviMAM.ascii",crs=proj),limite)
PluviJJA <- mask(raster("Rasters/PluviJJA.ascii",crs=proj),limite)
PluviSON <- mask(raster("Rasters/PluviSON.ascii",crs=proj),limite)

IESANO <- mask(raster("Rasters/IESANO.ascii",crs=proj),limite)
IEBANO <- mask(raster("Rasters/IEBANO.ascii",crs=proj),limite)
IEBDJF <- mask(raster("Rasters/IEBDJF.ascii",crs=proj),limite)
IEBMAM <- mask(raster("Rasters/IEBMAM.ascii",crs=proj),limite)
IEBJJA <- mask(raster("Rasters/IEBJJA.ascii",crs=proj),limite)
IEBSON <- mask(raster("Rasters/IEBSON.ascii",crs=proj),limite)

MuANO <- mask(raster("Rasters/MuANO.ascii",crs=proj),limite)
MuDJF <- mask(raster("Rasters/MuDJF.ascii",crs=proj),limite)
MuMAM <- mask(raster("Rasters/MuMAM.ascii",crs=proj),limite)
MuJJA <- mask(raster("Rasters/MuJJA.ascii",crs=proj),limite)
MuSON <- mask(raster("Rasters/MuSON.ascii",crs=proj),limite)

SigmaANO <- mask(raster("Rasters/SigmaANO.ascii",crs=proj),limite)
SigmaDJF <- mask(raster("Rasters/SigmaDJF.ascii",crs=proj),limite)
SigmaMAM <- mask(raster("Rasters/SigmaMAM.ascii",crs=proj),limite)
SigmaJJA <- mask(raster("Rasters/SigmaJJA.ascii",crs=proj),limite)
SigmaSON <- mask(raster("Rasters/SigmaSON.ascii",crs=proj),limite)

Q98ANO <- mask(raster("Rasters/Q98ANO.ascii",crs=proj),limite)
Q98DJF <- mask(raster("Rasters/Q98DJF.ascii",crs=proj),limite)
Q98MAM <- mask(raster("Rasters/Q98MAM.ascii",crs=proj),limite)
Q98JJA <- mask(raster("Rasters/Q98JJA.ascii",crs=proj),limite)
Q98SON <- mask(raster("Rasters/Q98SON.ascii",crs=proj),limite)

QmANO <- mask(raster("Rasters/QmANO.ascii",crs=proj),limite)
VR <- mask(raster("Rasters/VrQfmax.ascii",crs=proj),limite)

##-----------------------------------------------------------------------------##


## Paleta exportada do GRASS
color_cota<-colorRampPalette(c(rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(1,0,1),rgb(1,0,0)))(255)
color_hidro <- brewer.pal(9,"Blues")

<<<<<<< HEAD
#tikz("Figuras/DEM.tex",width=15/2.54, height=12/2.54)

pdf("Figuras/DEM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
par(mar = c(4,4,0.1,1.5))
=======
tikz("Figuras/DEM.tex",width=18/2.54, height=18/2.54)
#postscript("Figuras/DEM.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")

>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a
image(Cota,col=color_cota,asp = 1,xlab = "",ylab = "")
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
<<<<<<< HEAD
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

legend(2e5,6.93e6, legend = c(0,366.2,732.4,1098.6,1464.8,1831), fill =
    c(rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(1,0,1),rgb(1,0,0)),border = "white", cex =0.8,
=======
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 1,border = "white")

legend(2e5,6.93e6, legend = c(0,366.2,732.4,1098.6,1464.8,1831), fill =
    c(rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(1,0,1),rgb(1,0,0)),border = "white", cex =1,
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a
bg = "white",title="Cotas (m)",bty = "n")

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

<<<<<<< HEAD
text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()
embed_fonts("Figuras/DEM.pdf",outfile = "Figuras/DEM.pdf")


pdf("Figuras/QmANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family= "CM Roman")
=======
text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()

postscript("Figuras/QmANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(QmANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(QmANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(QmANO),2),round(maxValue(QmANO),2)),title =
<<<<<<< HEAD
    expression(bar(Q)~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=0.8)
=======
    expression(bar(Q)~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()
<<<<<<< HEAD
embed_fonts("Figuras/QmANO.pdf",outfile = "Figuras/QmANO.pdf")


pdf("Figuras/VR.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======

postscript("Figuras/VR.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(VR,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(VR,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(VR),2),round(maxValue(VR),2)),title =
<<<<<<< HEAD
    expression(V[r]~(m^{3}*km^{-2}*10^{5})),cex=0.8)
=======
    expression(V[r]~(m^{3}*km^{-2}*10^{5})),cex=1)
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()
<<<<<<< HEAD
embed_fonts("Figuras/VR.pdf",outfile = "Figuras/VR.pdf")


pdf("Figuras/Q98ANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======


postscript("Figuras/Q98ANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(Q98ANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98ANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98ANO),2),round(maxValue(Q98ANO),2)),title =
    expression(Q[98]~Anual~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()
<<<<<<< HEAD
embed_fonts("Figuras/Q98ANO.pdf",outfile = "Figuras/Q98ANO.pdf")


pdf("Figuras/Q98DJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======

postscript("Figuras/Q98DJF.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(Q98DJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98DJF,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98DJF),2),round(maxValue(Q98DJF),2)),title =
    expression(Q[98]~Verão~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/Q98MAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/Q98MAM.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(Q98MAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98MAM,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98MAM),2),round(maxValue(Q98MAM),2)),title =
    expression(Q[98]~Outono~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/Q98JJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/Q98JJA.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(Q98JJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98JJA,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98JJA),2),round(maxValue(Q98JJA),2)),title =
    expression(Q[98]~Inverno~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/Q98SON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/Q98SON.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(Q98SON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98SON,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98SON),2),round(maxValue(Q98SON),2)),title =
    expression(Q[98]~Primavera~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()

 
<<<<<<< HEAD
pdf("Figuras/MuANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/MuANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(MuANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuANO),2),round(maxValue(MuANO),2)),title =
    expression(mu~Anual),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()

<<<<<<< HEAD
pdf("Figuras/MuDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/MuDJF.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(MuDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuDJF,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuDJF),2),round(maxValue(MuDJF),2)),title =
    expression(mu~Verão),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/MuMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/MuMAM.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(MuMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuMAM,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuMAM),2),round(maxValue(MuMAM),2)),title =
    expression(mu~Outono),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/MuJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/MuJJA.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(MuJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuJJA,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuJJA),2),round(maxValue(MuJJA),2)),title =
    expression(mu~Inverno),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/MuSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/MuSON.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(MuSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuSON,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuSON),2),round(maxValue(MuSON),2)),title =
    expression(mu~Primavera),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()



<<<<<<< HEAD
pdf("Figuras/SigmaANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/SigmaANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(SigmaANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaANO),2),round(maxValue(SigmaANO),2)),title =
    expression(sigma~Anual),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()

<<<<<<< HEAD
pdf("Figuras/SigmaDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/SigmaDJF.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(SigmaDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaDJF,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaDJF),2),round(maxValue(SigmaDJF),2)),title =
    expression(sigma~Verão),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/SigmaMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/SigmaMAM.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(SigmaMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaMAM,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaMAM),2),round(maxValue(SigmaMAM),2)),title =
    expression(sigma~Outono),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/SigmaJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/SigmaJJA.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(SigmaJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaJJA,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaJJA),2),round(maxValue(SigmaJJA),2)),title =
    expression(sigma~Inverno),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/SigmaSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/SigmaSON.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(SigmaSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaSON,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaSON),2),round(maxValue(SigmaSON),2)),title =
    expression(sigma~Primavera),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
#tikz("Figuras/PluviANO.tex",width=15/2.54, height=12/2.54)
pdf("Figuras/PluviANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======

postscript("Figuras/PluviANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(PluviANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviANO),2),round(maxValue(PluviANO),2)),title =
<<<<<<< HEAD
    "Precipitação Anual (mm)",cex=0.8)
=======
    "Precipitação Anual (mm)",cex=1)
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()
<<<<<<< HEAD
embed_fonts("Figuras/PluviANO.pdf",outfile = "Figuras/PluviANO")


pdf("Figuras/PluviDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======



postscript("Figuras/PluviDJF.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(PluviDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviDJF,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviDJF),2),round(maxValue(PluviDJF),2)),title =
    "Precipitação Verão (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/PluviMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/PluviMAM.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(PluviMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviMAM,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviMAM),2),round(maxValue(PluviMAM),2)),title =
    "Precipitação Outono (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/PluviJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/PluviJJA.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(PluviJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviJJA,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviJJA),2),round(maxValue(PluviJJA),2)),title =
    "Precipitação Inverno (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/PluviSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/PluviSON.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(PluviSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviSON,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviSON),2),round(maxValue(PluviSON),2)),title =
    "Precipitação Primavera (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()

<<<<<<< HEAD
pdf("Figuras/IEBANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/IEBANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(IEBANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IEBANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IEBANO),2),round(maxValue(IEBANO),2)),title =
    "IEB",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()

<<<<<<< HEAD
pdf("Figuras/IESANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/IESANO.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

image(IESANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IESANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IESANO),2),round(maxValue(IESANO),2)),title =
    "IES",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)

dev.off()



##Amarzenando as bordas do Estado em um objeto
##bordas<-rbind(limite@polygons[[1]]@Polygons[[1]]@coords,limite@polygons[[2]]@Polygons[[1]]@coords)


##tikz("Estat.tex", width=12/2.54, height=12/2.54)

<<<<<<< HEAD
pdf("Figuras/BHs.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/BHs.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

plot(Hidro,col="lightblue",axes = T,asp=1)
plot(limite,add=T,cex=2)
plot(Bacias,lty=2,add=T,border="red")
points(FluviEST,col=1,pch=19)
points(Centroids,col=2,pch=19)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)


legend(2e5,6.815e6, legend=c("Estado de Santa Catarina","Hidrografia","Bacias Hidrográficas",
"Exutório",
"Centroide","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", inset=0.05,
    lty=c(1,1,2,NA,NA,NA,NA), pch=c(NA,NA,NA,19, 19,NA,NA),
    col=c(1,"lightblue", "red", 1,"red",NA,NA),cex=0.8)

dev.off()

<<<<<<< HEAD
pdf("Figuras/MacroBHs.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/MacroBHs.eps",onefile = T,horizontal = F, width=20/2.54, height=20/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a

plot(Hidro,col="lightblue",axes = T,asp=1)
plot(MBacias,lty=2,add=T,border="red")
plot(limite,add=T,cex=2)
text(getSpPPolygonsLabptSlots(MBacias), labels=as.character(MBacias$bacia), cex=0.7,font=2)

## colocar labels 
SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)


legend(2e5,6.785e6, legend=c("Estado de Santa Catarina","Hidrografia","Bacias Hidrográficas","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", inset=0.05,
    lty=c(1,1,2,NA,NA),
    col=c(1,"lightblue", "red",NA,NA),cex=0.8)

dev.off()


<<<<<<< HEAD
pdf("Figuras/Estat1.pdf",onefile = T, width=10/2.54, height=10/2.54,paper = "special",family
    = "CM Roman")
=======
postscript("Figuras/Estat1.eps",onefile = T,horizontal = F, width=10/2.54, height=10/2.54,paper = "special")
>>>>>>> 5080a23dbc1bfc104475844919d3270c46243e1a


plot(Hidro,col="lightblue",axes = T,asp=1)
plot(limite,add=T,cex=2)
points(FluviEST,col=1,pch=19)
points(PluviEST,col="darkblue",pch=19)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)


legend(2e5,6.8e6, legend=c("Estado de Santa Catarina","Hidrografia",
"Estações Fluviométricas",
"Estações Pluviométricas","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", inset=0.05,
    lty=c(1,1,NA,NA,NA,NA), pch=c(NA,NA,19, 19,NA,NA),
    col=c(1,"lightblue", 1, "darkblue",NA,NA),cex=0.8)

dev.off()










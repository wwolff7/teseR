## Configuração Inicial
require(tikzDevice)
require(RColorBrewer) # Para Paletas de cores
require(maptools)
library(sp)
library(rgdal)
library(raster)
require(SDMTools) # para legend.gradient
library(extrafont)
loadfonts()

getwd()
setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Mapas")
#rm(list = ls()) ## Para remover todos os objetos
citation("maptools")
options(OutDec=".",digits=10,scipen=5)
#par(mar = c(4,4,0.1,0.1))

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

map<- get_map(location = 'US', zoom = 4)

map <- get_map('Santa Catarina',zoom = 6,mapty="hybrid")
ggmap(map)
qmap('Santa Catarina',zoom = 6,mapty="hybrid")

## Projeção
proj <- CRS("+proj=utm +zone=22 +south +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 +units=m +no_defs")
## lendo dados tipo shapefile
Sul <- readShapeSpatial("Shapes/Sul_Brasil.shp",proj4string = proj)
limite <- readShapeSpatial("Shapes/SC_limite.shp",proj4string = proj)
FluviEST <- readShapeSpatial("Shapes/Fluvi_Estat.shp",proj4string = proj)
PluviEST <- readShapeSpatial("Shapes/Pluvi_Estat.shp",proj4string = proj)
Hidro <- readShapeSpatial("Shapes/Hidro.shp",proj4string = proj)
Bacias <- readShapeSpatial("Shapes/Bacias/merged/merged.shp",proj4string = proj)
Centroids <- readShapeSpatial("Shapes/Centroids.shp",proj4string = proj)
MBacias <- readShapeSpatial("Shapes/Macro_Bacias.shp",proj4string = proj)


## lendo dados tipo Raster
Cota <- raster("Rasters/cotas.tif")
Decliv <- terrain(Cota, opt = "slope", unit = "degrees", df=F)

coords <- locator(type="l") 

MBacias_sub <- MBacias[1:3, ]

image(Cota)
plot(MBacias_sub, add=T)

extract(Cota, MBacias_sub, fun = mean, na.rm = T, small = T, df = T)


PluviANO <- mask(raster("Rasters/PluviANO.tiff",crs=proj),limite)
PluviDJF <- mask(raster("Rasters/PluviDJF.tiff",crs=proj),limite)
PluviMAM <- mask(raster("Rasters/PluviMAM.tiff",crs=proj),limite)
PluviJJA <- mask(raster("Rasters/PluviJJA.tiff",crs=proj),limite)
PluviSON <- mask(raster("Rasters/PluviSON.tiff",crs=proj),limite)

IESANO <- mask(raster("Rasters/IESANO.tiff",crs=proj),limite)
IEBANO <- mask(raster("Rasters/IEBANO.tiff",crs=proj),limite)
IEBDJF <- mask(raster("Rasters/IEBDJF.tiff",crs=proj),limite)
IEBMAM <- mask(raster("Rasters/IEBMAM.tiff",crs=proj),limite)
IEBJJA <- mask(raster("Rasters/IEBJJA.tiff",crs=proj),limite)
IEBSON <- mask(raster("Rasters/IEBSON.tiff",crs=proj),limite)

MuANO <- mask(raster("Rasters/MuANO.tiff",crs=proj),limite)
MuDJF <- mask(raster("Rasters/MuDJF.tiff",crs=proj),limite)
MuMAM <- mask(raster("Rasters/MuMAM.tiff",crs=proj),limite)
MuJJA <- mask(raster("Rasters/MuJJA.tiff",crs=proj),limite)
MuSON <- mask(raster("Rasters/MuSON.tiff",crs=proj),limite)

SigmaANO <- mask(raster("Rasters/SigmaANO.tiff",crs=proj),limite)
SigmaDJF <- mask(raster("Rasters/SigmaDJF.tiff",crs=proj),limite)
SigmaMAM <- mask(raster("Rasters/SigmaMAM.tiff",crs=proj),limite)
SigmaJJA <- mask(raster("Rasters/SigmaJJA.tiff",crs=proj),limite)
SigmaSON <- mask(raster("Rasters/SigmaSON.tiff",crs=proj),limite)

## Mapas de vazões
## funções
Qp <- function(mu,sigma,area,p){
    Qp<-(exp(mu+(sigma*qnorm(1-p))))*area
}


Q98ANO <- Qp(MuANO,SigmaANO,1,0.98)*1000
Q98DJF <- Qp(MuDJF,SigmaDJF,1,0.98)*1000
Q98MAM <- Qp(MuMAM,SigmaMAM,1,0.98)*1000
Q98JJA <- Qp(MuJJA,SigmaJJA,1,0.98)*1000
Q98SON <- Qp(MuSON,SigmaSON,1,0.98)*1000

QmANO <- mask(raster("Rasters/QmANO.tiff",crs=proj),limite)
VR <- mask(raster("Rasters/Vr-Qfmax.tiff",crs=proj),limite)


##-----------------------------------------------------------------------------##


## Paleta exportada do GRASS
color_cota<-colorRampPalette(c(rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(1,0,1),rgb(1,0,0)))(255)
color_hidro <- brewer.pal(9,"Blues")

#tikz("Figuras/DEM.tex",width=15/2.54, height=12/2.54)

pdf("Figuras/Estat-FluviDEM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(Cota,col=color_cota,asp = 1,xlab = "",ylab = "")
plot(Sul,add=T,cex=2)
#points(PluviEST,col=1,pch=19)
points(FluviEST,col=1,pch=19)

## legend("bottomleft",
##        legend=c("States","Rain gauges","Coordinate reference system - UTM ","Datum: SAD/69 Zone 22S"),
##     title="Legend", bg="white", lty=c(1,NA,NA,NA), pch=c(NA,19,NA,NA),
##     col=c(1,1,NA,NA), cex = 0.8,border = "white")

legend("bottomleft",
       legend=c("Estados","Estações Fluviométricas","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zona 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA,NA), pch=c(NA,19,NA,NA),
    col=c(1,1,NA,NA), cex = 0.8,border = "white")

legend(2e5,6.945e6, legend = c(0,366.2,732.4,1098.6,1464.8,1831), fill =
    c(rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(1,0,1),rgb(1,0,0)),border = "white", cex =1,
bg = "white",title="Cotas (m)",bty = "n")

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

text(5.0e5, 7.2e6, "Paraná", cex= 1,font=2)
text(5.0e5, 7.015e6, "Santa Catarina", cex= 1,font=2)
text(4.2e5, 6.85e6, "Rio Grande do Sul", cex= 1,font=2)
text(7.18e5, 6.77e6, "Oceano atlântico", cex= 1,font=2)
#text(7.18e5, 6.77e6, "Atlantic ocean", cex= 1,font=2)

dev.off()
embed_fonts("Figuras/Estat-FluviDEM.pdf",outfile = "Figuras/Estat-FluviDEM.pdf")


pdf("Figuras/MacroBHs.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

#plot(Hidro,col="lightblue",axes = T,asp=1)
plot(MBacias,lty=2,axes = T,asp=1,border="black")
plot(limite,add=T,cex=2)
points(FluviEST,col="grey",pch=19)
text(coordinates(MBacias)[-c(16,20),], labels=as.character(MBacias$bacia)[-c(16,20)], cex=0.7,font=2,col="black")
text(coordinates(MBacias)[c(16,20),]*0.998, labels=as.character(MBacias$bacia)[c(16,20)], cex=0.7,font=2,col="black")
##text(coordinates(FluviEST), labels=as.character(FluviEST$Estação), cex=0.7,font=2,col="grey")

## colocar labels 
SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 0.8)
text(7.5e5, 6.72e6, "200 km", cex= 0.8)


legend(2e5,6.785e6, legend=c("Estado de Santa Catarina","Macrobacias","Estações Fluviométricas","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", inset=0.05,
    lty=c(1,2,NA,NA,NA),
    col=c(1,1, "grey",NA,NA),cex=0.8,
    pch=c(NA,NA,19,NA,NA))

dev.off()
embed_fonts("Figuras/MacroBHs.pdf",outfile = "Figuras/MacroBHs.pdf")



pdf("Figuras/QmANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family= "CM Roman")

image(QmANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(QmANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(QmANO),2),round(maxValue(QmANO),2)),title =
    expression(bar(Q)~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/QmANO.pdf",outfile = "Figuras/QmANO.pdf")


pdf("Figuras/VR.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(VR,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(VR,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(VR),2),round(maxValue(VR),2)),title =
    expression(V[r]~(m^{3}*km^{-2}*10^{6})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/VR.pdf",outfile = "Figuras/VR.pdf")


pdf("Figuras/Q98ANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(Q98ANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98ANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98ANO),2),round(maxValue(Q98ANO),2)),title =
    expression(Q[98]~Anual~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/Q98ANO.pdf",outfile = "Figuras/Q98ANO.pdf")

pdf("Figuras/Q98ALL.pdf",onefile = T, width=21/2.54, height=29/2.54,paper = "special",family
    = "CM Roman")

par(mfrow = c(3,2),mar = c(2,2.5,2.5,1.5),mgp = c(1,1,0))


##pdf("Figuras/Q98DJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(Q98DJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98DJF,add = T)
plot(limite,add=T,cex=2)
mtext("(a)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98DJF),2),round(maxValue(Q98DJF),2)),title =
    expression(Q[98]~Verão~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/Q98DJF.pdf",outfile = "Figuras/Q98DJF.pdf")

##pdf("Figuras/Q98MAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(Q98MAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98MAM,add = T)
plot(limite,add=T,cex=2)
mtext("(b)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98MAM),2),round(maxValue(Q98MAM),2)),title =
    expression(Q[98]~Outono~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/Q98MAM.pdf",outfile = "Figuras/Q98MAM.pdf")

##pdf("Figuras/Q98JJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(Q98JJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98JJA,add = T)
plot(limite,add=T,cex=2)
mtext("(c)",adj=0,line = 0.5, cex=1.5)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98JJA),2),round(maxValue(Q98JJA),2)),title =
    expression(Q[98]~Inverno~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/Q98JJA.pdf",outfile = "Figuras/Q98JJA.pdf")
## 
##pdf("Figuras/Q98SON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(Q98SON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(Q98SON,add = T)
plot(limite,add=T,cex=2)
mtext("(d)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.88e6,6.88e6,6.78e6,6.78e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(Q98SON),2),round(maxValue(Q98SON),2)),title =
    expression(Q[98]~Primavera~(m^{3}*s^{-1}*km^{-2}*10^{-3})),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/Q98SON.pdf",outfile = "Figuras/Q98SON.pdf")

plot(1,axes = F,ann = F,type = "n")
legend("topleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 1)

dev.off()
embed_fonts("Figuras/Q98ALL.pdf",outfile = "Figuras/Q98ALL.pdf")

pdf("Figuras/MuANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(MuANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuANO),2),round(maxValue(MuANO),2)),title =
    expression(mu~Anual),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/MuANO.pdf",outfile = "Figuras/MuANO.pdf")

pdf("Figuras/MuALL.pdf",onefile = T, width=21/2.54, height=29/2.54,paper = "special",family
    = "CM Roman")

par(mfrow = c(3,2),mar = c(2,2.5,2.5,1.5),mgp = c(1,1,0))


##pdf("Figuras/MuDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(MuDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuDJF,add = T)
plot(limite,add=T,cex=2)
mtext("(a)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuDJF),2),round(maxValue(MuDJF),2)),title =
    expression(mu~Verão),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/MuDJF.pdf",outfile = "Figuras/MuDJF.pdf")
## 
##pdf("Figuras/MuMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(MuMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuMAM,add = T)
plot(limite,add=T,cex=2)
mtext("(b)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuMAM),2),round(maxValue(MuMAM),2)),title =
    expression(mu~Outono),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/MuMAM.pdf",outfile = "Figuras/MuMAM.pdf")
## 
##pdf("Figuras/MuJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(MuJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuJJA,add = T)
plot(limite,add=T,cex=2)
mtext("(c)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuJJA),2),round(maxValue(MuJJA),2)),title =
    expression(mu~Inverno),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/MuJJA.pdf",outfile = "Figuras/MuJJA.pdf")
## 
##pdf("Figuras/MuSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(MuSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(MuSON,add = T)
plot(limite,add=T,cex=2)
mtext("(d)",adj=0,line = 0.5, cex=1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(MuSON),2),round(maxValue(MuSON),2)),title =
    expression(mu~Primavera),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/MuSON.pdf",outfile = "Figuras/MuSON.pdf")

plot(1,axes = F,ann = F,type = "n")
legend("topleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 1)

dev.off()
embed_fonts("Figuras/MuALL.pdf",outfile = "Figuras/MuALL.pdf")



pdf("Figuras/SigmaANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(SigmaANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaANO),2),round(maxValue(SigmaANO),2)),title =
    expression(sigma~Anual),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/SigmaANO.pdf",outfile = "Figuras/SigmaANO.pdf")

pdf("Figuras/SigmaALL.pdf",onefile = T, width=21/2.54, height=29/2.54,paper = "special",family
    = "CM Roman")

par(mfrow = c(3,2),mar = c(2,2.5,2.5,1.5),mgp = c(1,1,0))


##pdf("Figuras/SigmaDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(SigmaDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaDJF,add = T)
plot(limite,add=T,cex=2)
mtext("(a)",line = 0.5, adj=0,cex = 1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaDJF),2),round(maxValue(SigmaDJF),2)),title =
    expression(sigma~Verão),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/SigmaDJF.pdf",outfile = "Figuras/SigmaDJF.pdf")
## 
##pdf("Figuras/SigmaMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(SigmaMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaMAM,add = T)
plot(limite,add=T,cex=2)
mtext("(b)",line = 0.5, adj=0,cex = 1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaMAM),2),round(maxValue(SigmaMAM),2)),title =
    expression(sigma~Outono),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/SigmaMAM.pdf",outfile = "Figuras/SigmaMAM.pdf")
## 
##pdf("Figuras/SigmaJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(SigmaJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaJJA,add = T)
plot(limite,add=T,cex=2)
mtext("(c)",line = 0.5, adj=0,cex = 1.5)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaJJA),2),round(maxValue(SigmaJJA),2)),title =
    expression(sigma~Inverno),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/SigmaJJA.pdf",outfile = "Figuras/SigmaJJA.pdf")

##pdf("Figuras/SigmaSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(SigmaSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(SigmaSON,add = T)
plot(limite,add=T,cex=2)
mtext("(d)",line = 0.5, adj=0,cex = 1.5)

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(SigmaSON),2),round(maxValue(SigmaSON),2)),title =
    expression(sigma~Primavera),cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/SigmaSON.pdf",outfile = "Figuras/SigmaSON.pdf")

plot(1,axes = F,ann = F,type = "n")
legend("topleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 1)

dev.off()
embed_fonts("Figuras/SigmaALL.pdf",outfile = "Figuras/SigmaALL.pdf")
plot(1)

#tikz("Figuras/PluviANO.tex",width=18/2.54, height=18/2.54)
pdf("Figuras/PluviANO_ing.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(PluviANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Santa Catarina state","Coordinate reference system - UTM ","Datum: SAD/69 Zone 22S"),
    title="Legend", bg="white", lty=c(1,NA,NA), 
    col=c(1,NA,NA), cex = 0.8,border = "white",bty = "n")

## legend("bottomleft",
##        legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##     title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8 ,bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviANO),2),round(maxValue(PluviANO),2)),title =
    "Annual precipitation (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/PluviANO_ing.pdf",outfile = "Figuras/PluviANO_ing.pdf")

pdf("Figuras/PluviALL_ing.pdf",onefile = T, width=16/2.54, height=21/2.54,paper = "special",family
    = "CM Roman")

par(mfrow = c(3,2),mar = c(2,2.5,2.5,1),mgp = c(1,1,0))


##pdf("Figuras/PluviDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(PluviDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviDJF,add = T)
plot(limite,add=T,cex=2)
mtext("(a)",line = 0.5,adj=0,cex=1)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8, bty="n")

pnts <- cbind(x=c(2.0e5,2.2e5,2.2e5,2.0e5),y=c(6.85e6,6.85e6,6.75e6,6.75e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviDJF),2),round(maxValue(PluviDJF),2)),title =
    "Summer precipitation (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/PluviDJF.pdf",outfile = "Figuras/PluviDJF.pdf")
## 
##pdf("Figuras/PluviMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(PluviMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviMAM,add = T)
plot(limite,add=T,cex=2)
mtext("(b)",line = 0.5,adj=0,cex=1)
##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8, bty="n")

pnts <- cbind(x=c(2.0e5,2.2e5,2.2e5,2.0e5),y=c(6.85e6,6.85e6,6.75e6,6.75e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviMAM),2),round(maxValue(PluviMAM),2)),title =
    "Autumn precipitation (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/PluviMAM.pdf",outfile = "Figuras/PluviMAM.pdf")
## 
##pdf("Figuras/PluviJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(PluviJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviJJA,add = T)
plot(limite,add=T,cex=2)
mtext("(c)",line = 0.5,adj=0,cex=1)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8, bty="n")

pnts <- cbind(x=c(2.0e5,2.2e5,2.2e5,2.0e5),y=c(6.85e6,6.85e6,6.75e6,6.75e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviJJA),2),round(maxValue(PluviJJA),2)),title =
    "Winter precipitation (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/PluviJJA.pdf",outfile = "Figuras/PluviJJA.pdf")

##pdf("Figuras/PluviSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(PluviSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(PluviSON,add = T)
plot(limite,add=T,cex=2)
mtext("(d)",line = 0.5,adj=0,cex=1)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,bty="n")

pnts <- cbind(x=c(2.0e5,2.2e5,2.2e5,2.0e5),y=c(6.85e6,6.85e6,6.75e6,6.75e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(PluviSON),2),round(maxValue(PluviSON),2)),title =
    "Spring precipitation (mm)",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/PluviSON.pdf",outfile = "Figuras/PluviSON.pdf")

plot(1,axes = F,ann = F,type = "n")

legend("topleft",
       legend=c("Santa Catarina state","Coordinate reference system - UTM ","Datum: SAD/69 Zone 22S"),
    title="Legend", bg="white", lty=c(1,NA,NA), 
    col=c(1,NA,NA),bty = "n", cex = 1)

## legend("topleft",
##        legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##     title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 1)

dev.off()
embed_fonts("Figuras/PluviALL_ing.pdf",outfile = "Figuras/PluviALL_ing.pdf")


pdf("Figuras/IEBANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(IEBANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IEBANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IEBANO),2),round(maxValue(IEBANO),2)),title =
    "IEB Anual",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/IEBANO.pdf",outfile = "Figuras/IEBANO.pdf")


pdf("Figuras/IEBALL.pdf",onefile = T, width=21/2.54, height=29/2.54,paper = "special",family
    = "CM Roman")

par(mfrow = c(3,2),mar = c(2,2.5,2.5,1.5),mgp = c(1,1,0))


##pdf("Figuras/IEBDJF.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(IEBDJF,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IEBDJF,add = T)
plot(limite,add=T,cex=2)
mtext("(a)",line = 0.5,adj=0,cex=1.5)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IEBDJF),2),round(maxValue(IEBDJF),2)),title =
    "IEB Verão",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/IEBDJF.pdf",outfile = "Figuras/IEBDJF.pdf")
## 
## 
##pdf("Figuras/IEBMAM.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(IEBMAM,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IEBMAM,add = T)
plot(limite,add=T,cex=2)
mtext("(b)",line = 0.5,adj=0,cex=1.5)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IEBMAM),2),round(maxValue(IEBMAM),2)),title =
    "IEB Outono",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/IEBMAM.pdf",outfile = "Figuras/IEBMAM.pdf")
## 
## 
##pdf("Figuras/IEBJJA.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(IEBJJA,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IEBJJA,add = T)
plot(limite,add=T,cex=2)
mtext("(c)",line = 0.5,adj=0,cex=1.5)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IEBJJA),2),round(maxValue(IEBJJA),2)),title =
    "IEB Inverno",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/IEBJJA.pdf",outfile = "Figuras/IEBJJA.pdf")
## 
## 
##pdf("Figuras/IEBSON.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
##    = "CM Roman")

image(IEBSON,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IEBSON,add = T)
plot(limite,add=T,cex=2)
mtext("(d)",line = 0.5,adj=0,cex=1.5)

##legend("bottomleft",
##       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
##    title="Legenda", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IEBSON),2),round(maxValue(IEBSON),2)),title =
    "IEB Primavera",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

##dev.off()
##embed_fonts("Figuras/IEBSON.pdf",outfile = "Figuras/IEBSON.pdf")

plot(1,axes = F,ann = F,type = "n")
legend("topleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty = "n", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 1)

dev.off()
embed_fonts("Figuras/IEBALL.pdf",outfile = "Figuras/IEBALL.pdf")

pdf("Figuras/IESANO.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

image(IESANO,col=color_hidro,asp = 1,xlab = "",ylab = "")
contour(IESANO,add = T)
plot(limite,add=T,cex=2)

legend("bottomleft",
       legend=c("Estado de Santa Catarina","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bg="white", lty=c(1,NA,NA), col=c(1,NA,NA),cex = 0.8,border = "white",bty="n")

pnts <- cbind(x=c(2.3e5,2.5e5,2.5e5,2.3e5),y=c(6.9e6,6.9e6,6.8e6,6.8e6))
legend.gradient(pnts,cols=color_hidro,limits = c(round(minValue(IESANO),2),round(maxValue(IESANO),2)),title =
    "IES Anual",cex=1)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("white", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)

dev.off()
embed_fonts("Figuras/IESANO.pdf",outfile = "Figuras/IESANO.pdf")


##Amarzenando as bordas do Estado em um objeto
##bordas<-rbind(limite@polygons[[1]]@Polygons[[1]]@coords,limite@polygons[[2]]@Polygons[[1]]@coords)


##tikz("Estat.tex", width=12/2.54, height=12/2.54)

pdf("Figuras/BHs.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")

plot(Hidro,col="lightblue",axes = T,asp=1)
plot(limite,add=T,cex=2)
plot(Bacias,lty=2,add=T,border="red")
points(FluviEST,col=1,pch=19)
points(Centroids,col=2,pch=19)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)


legend(2e5,6.815e6, legend=c("Estado de Santa Catarina","Hidrografia","Bacias Hidrográficas",
"Exutório",
"Centroide","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n",
    lty=c(1,1,2,NA,NA,NA,NA), pch=c(NA,NA,NA,19, 19,NA,NA),
    col=c(1,"lightblue", "red", 1,"red",NA,NA),cex=0.8)

dev.off()
embed_fonts("Figuras/BHs.pdf",outfile = "Figuras/BHs.pdf")

pdf("Figuras/Estat-Pluvi.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")


plot(Hidro,col="lightblue",axes = T,asp=1)
plot(limite,add=T,cex=2)
points(PluviEST,col=1,pch=19)


SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)


legend(2e5,6.8e6, legend=c("Estado de Santa Catarina","Hidrografia",
"Estações Pluviométricas","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda",bty = "n",
    lty=c(1,1,NA,NA,NA), pch=c(NA,NA,19,NA,NA),
    col=c(1,"lightblue", 1, NA,NA),cex=0.8)

dev.off()
embed_fonts("Figuras/Estat-Pluvi.pdf",outfile = "Figuras/Estat-Pluvi.pdf")


pdf("Figuras/Estat-Fluvi.pdf",onefile = T, width=18/2.54, height=18/2.54,paper = "special",family
    = "CM Roman")


plot(Hidro,col="lightblue",axes = T,asp=1)
plot(limite,add=T,cex=2)
points(FluviEST,col=1,pch=19)

SpatialPolygonsRescale(layout.north.arrow(1), offset= c(2.5e5,7.13e6), scale = 5e4, plot.grid=F)
SpatialPolygonsRescale(layout.scale.bar(), offset= c(5.5e5,6.7e6), scale= 2e5, fill= c("transparent", "black"), plot.grid= F)

text(5.5e5, 6.72e6, "0", cex= 1)
text(7.5e5, 6.72e6, "200 km", cex= 1)


legend(2e5,6.8e6, legend=c("Estado de Santa Catarina","Hidrografia",
"Estações Fluviométricas","Sistemas de Coordenadas Projetadas UTM","Datum: SAD/69 Zonas 22S"),
    title="Legenda", bty="n",
    lty=c(1,1,NA,NA,NA), pch=c(NA,NA,19,NA,NA),
    col=c(1,"lightblue", 1, NA,NA),cex=0.8)

dev.off()
embed_fonts("Figuras/Estat-Fluvi.pdf",outfile = "Figuras/Estat-Fluvi.pdf")


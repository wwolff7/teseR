
## Configuração Inicial

require(RColorBrewer) # Para Paletas de cores
require(hydroTSM) 
colour<-colorRampPalette(c("white","blue","dark Blue"))
require(MASS)
require(rgdal)
require(maptools)
require(geoR)
require(xtable)
##rm(list = ls())

setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Pluvigeo")

options(OutDec=",",digits=10,scipen=5)
##-----------------------------------------------------------------------------##

## Funções Norte e Legenda

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

dados <- read.csv("Pluvi.csv", head=T,dec=",")[,c(20,19,14,15,16,17,18,21)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC_conti.csv",h=T);limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
floripa <- read.csv("Limite_Floripa.csv",h=T);floripa$Longitude<-floripa$Longitude/1000;floripa$Latitude<-floripa$Latitude/1000


head(limite)
head(floripa)
dim(limite)
## lendo dados tipo shapefile
limiteSHP <- readShapeSpatial("SC.shp",proj4string = CRS("+proj=utm +zone=22 +south +ellps=aust_SA +towgs84=-57,1,-41,0,0,0,0 +units=m +no_defs"))
## Amarzenando as bordas do Estado em um objeto
## bordas<-rbind(limite@polygons[[1]]@Polygons[[1]]@coords,limite@polygons[[2]]@Polygons[[1]]@coords)

plot(limiteSHP,axes = T,)

##-----------------------------------------------------------------------------##

## verificando argumentos e documentação da função utilizada
args(read.table)
help(read.table)

## visualizando dados importados
View(dados) 
dados
##-----------------------------------------------------------------------------##

## verificando parte dos dados e caracteristicas do conjunto de dados
head(dados)
class(dados)
dim(dados) 
names(dados)
smry(dados)

mean(dados[,3])
mean(dados$ANUAL)

##-----------------------------------------------------------------------------##
## Curvas de permanência das Estações Fluviométricas
## Probabilidades testadas, lembrar que refere-se a 1-p !! 


# Adding the sp.layout parameter shows the locations of the measurements

postscript("Figuras/Estacoes.eps",onefile = T, horizontal = F, width=15/2.54, height=15/2.54)
plot(limite,type="l",asp=1,xlim = c(min(limite[,1]),max(limite[,1])),ylim = c(min(limite[,2]),max(limite[,2])),ylab = "Latitude (km)",xlab = "Longitude (km)")
lines(floripa)
points(dados$Longitude,dados$Latitude,pch=19,col=1)
northarrow(loc=c(300,7000),size = 30,cex =0.8)
scalebar(loc=c(300,7100),length=400)


dev.off()


##-----------------------------------------------------------------------------##
## Análise exploratória dos dados 

##-----------------------------------------------------------------------------##
## Pluvi Anual

## Análise exploratória geoestatística
names(dados)
Pluvi.geo<-as.geodata(dados,coords.col = 1:2, data.col = 3,covar.col = 8,)
Pluvi.geo$borders <- limite

plot(Pluvi.geo,low=T)
plot(Pluvi.geo,low=T,trend="1st")

postscript("Figuras/BoxCox_PluviANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Pluvi.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(Pluvi.geo$data~1,lambda = seq(-3,3,l=20),ylab = "Log-Verossimilhança")
dev.off()

## Não precisa transformar
## Pluvi.geo$borders<-limite
##hist(Pluvi.geo$data)
##boxcox(Pluvi.geo$data~1) 
##boxcox(Pluvi.geo$data~1,lambda=seq(-1,0,l=20))
## 
##trans<-boxcox(Pluvi.geo$data~1,lambda=seq(-0.6,-0.4,l=20));trans # dando zoom para verqual lambda usar
##lambdaANO<-with(trans, x[which.max(y)]);lambdaANO # Valor máximo de Lambda
## 
##pluvi.n<-(((dados$ANUAL^(lambdaANO)) - 1)/lambdaANO);pluvi.n #normalizando - (X^lambda)-1/lambda
##pluvi.inv<-(((lambdaANO*pluvi.n)+1)^(1/lambdaANO));pluvi.inv #inverter e diminuir com 6 para encontrar a varivel original
##Pluvi.geo$data <- pluvi.n
## 
##smry(dados$ANUAL)

hist(Pluvi.geo$data)
hist(dados$ANUAL)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normal e sem tendência
##-----------------------------------------------------------------------------##
## Conclusão: há dependência espacial
##-----------------------------------------------------------------------------##
##-----------------------------------------------------------------------------##

## Interpolação espacial (krigagem)

## Modelagem Precipitação média anual

ef.pluvi<-eyefit(v.pluvi)

Pluvi <- list()
Pluvi$lf0 <- likfit(Pluvi.geo, ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf1 <- likfit(Pluvi.geo, trend=~Cota,  ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf2 <- likfit(Pluvi.geo, trend=~coords[,1],ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf3 <- likfit(Pluvi.geo, trend=~coords[,2],ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf4 <- likfit(Pluvi.geo, trend="1st",  ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf5 <- likfit(Pluvi.geo, trend="2nd",  ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf6 <- likfit(Pluvi.geo, trend=~coords[,1]+Cota,  ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf7 <- likfit(Pluvi.geo, trend=~coords[,2]+Cota,  ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf8 <- likfit(Pluvi.geo, trend=~coords[,1]+poly(Cota,2), ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf9 <- likfit(Pluvi.geo, trend=~coords[,2]+poly(Cota,2), ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf10 <- likfit(Pluvi.geo, trend=~poly(coords[,1],2)+poly(Cota,2), ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf11 <- likfit(Pluvi.geo, trend=~poly(coords[,2],2)+poly(Cota,2), ini=c(2.80e-6,24.32), nug=3.5e-7)
Pluvi$lf12 <- likfit(Pluvi.geo, trend=~poly(coords[,1],2)+poly(coords[,2],2)+poly(Cota,2),ini=c(2.80e-6,24.32), nug=3.5e-7)


sort(sapply(Pluvi,AIC))
sort(sapply(Pluvi,function(x) x$BIC))
sapply(Pluvi,function(x) x$value)
plot(Pluvi.geo,low=T,trend="1st")
plot(Pluvi.geo,low=T,trend="2nd")
plot(Pluvi.geo,low=T,trend=~coords[,1])


##-----------------------------------------------------------------------------##


Pluvilf4<-list()
Pluvilf4$exp <- likfit(Pluvi.geo, ini=c(2.80e-6,24.32), nug=3.5e-7,trend="1st")
Pluvilf4$gau <- likfit(Pluvi.geo, cov.model = "gau", ini=c(2.80e-6,36.49), nug=3.5e-7,trend="1st")  
Pluvilf4$sph <- likfit(Pluvi.geo, cov.model = "sph", ini=c(2.80e-6,72.97), nug=3.5e-7,trend="1st") 
Pluvilf4$cir <- likfit(Pluvi.geo, cov.model = "cir", ini=c(2.80e-6,72.97), nug=3.5e-7,trend="1st")
Pluvilf4$kappa1.5 <- likfit(Pluvi.geo, cov.model = "mat", kappa= 1.5, ini=c(2.80e-6,12.16), nug=3.5e-7,trend="1st")  
Pluvilf4$kappa2.5 <- likfit(Pluvi.geo, cov.model = "mat", kappa= 2.5, ini=c(2.80e-6,12.16), nug=3.5e-7,trend="1st") 

##-----------------------------------------------------------------------------

v.pluviANO<-variog(Pluvi.geo,max.dist=200,uvec=seq(0, 200, by=10),trend="1st",option = "cloud")
postscript("Figuras/semivarioPaANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
plot(v.pluviANO,xlab="Distância (km)",ylab=expression(gamma(u)))#, ylim = c(0,2e4),mgp=c(2.5,1,0))
lines(lowess(v.pluviANO$u,v.pluviANO$v),col = 2,cex = 2)
dev.off()


sort(sapply(Pluvilf4, AIC))

plot(v.pluviANO,xlab="Distance (km)")#,ylim =c(0,4e-6),ylab=expression(gamma(u)))
lines(Pluvilf4$exp,col=1)
lines(Pluvilf4$gau,col=2)
lines(Pluvilf4$sph,col=3)
lines(Pluvilf4$cir,col=1)
lines(Pluvilf4$kappa1.5,col=5)
lines(Pluvilf4$kappa2.5,col=6)

## modelo escolhido, pela MV o melhor foi o circular
## escolhendo o modelo circular
summary(Pluvilf4$cir)
Pluvilf4$cir$sigmasq
Pluvilf4$cir$tausq
Pluvilf4$cir$beta
Pluvilf4$cir$phi
Pluvilf4$cir$parameters.summary


##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Predição na área

##rm(list = ls())

grid <- pred_grid(c(200,800),c(6700,7200), by=1)
points(Pluvi.geo)
points(grid, pch=19, cex=0.25, col=2)
rm(gr)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
kc.pluvi <- krige.conv(Pluvi.geo, loc=grid, krige=krige.control(obj=Pluvilf4$cir))
attributes(kc.pluvi)
kc.pluvi$predict
##kc.pluvi$predict<-(backtransform.moments(lambda=lambdaANO,mean=kc.pluvi$predict,variance=kc.pluvi$krige.var)$mean)

Pluvi_krige <- cbind(grid*1000,kc.pluvi$predict)
write.table(Pluvi_krige,"PluviANO.ascii",col.names = F,row.names=F,quote=F)

##-----------------------------------------------------------------------------
## Validação cruzada
xv.pluvi<-xvalid(Pluvi.geo,model=Pluvilf4$cir)
names(xv.pluvi)
xv.pluvi$data
dados$ANUAL

##xv.pluvi$data <- dados$ANUAL
##xv.pluvi$predicted <- (backtransform.moments(lambda=lambdaANO,mean=xv.pluvi$predicted,variance=xv.pluvi$krige.var)$mean)

plot(xv.pluvi$predicted,xv.pluvi$data) #,xlim=c(0,1.2),ylim=c(0,1.2))
abline(0,1)
abline(lm(xv.pluvi$data~xv.pluvi$predicted)$coef[1],lm(xv.pluvi$data~xv.pluvi$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.pluvi)


smry(Pluvi.geo$data)


postscript("Figuras/xv_pluviANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
plot(xv.pluvi$predicted,xv.pluvi$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)
hist(xv.pluvi$std.error,breaks = 20,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)
dev.off()


postscript("Figuras/kr_pluvi.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.pluvi,val=kc.pluvi$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.pluvi, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)


dev.off()

summary(limite)

pdf("Figuras/sd.pluvi.pdf",onefile = T, width=20/2.54, height=20/2.54)

image(kc.pluvi, val=kc.pluvi$krige.var,col=grey(seq(1,0.2,len=21)),x.leg=c(250,450),y.leg=c(6800,6810),xlab = "Longitude (km)",ylab = "Latitude (km)")

dev.off()


##-----------------------------------------------------------------------------              

## Pluvi Verão

## Análise exploratória geoestatística
head(dados)
names(dados)
PluviDJF.geo<-as.geodata(dados,coords.col = 1:2, data.col = 4,covar.col = 8,)

plot(PluviDJF.geo,low=T)
## PluviDJF.geo$borders<-limite

postscript("Figuras/BoxCox_PluviDJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper="special")
par(mfrow = c(1,2))
hist(PluviDJF.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(PluviDJF.geo$data~1,lambda = seq(-10,10),ylab = "Log-Verossimilhança")
dev.off()

names(PluviDJF.geo)
summary(PluviDJF.geo$data)
plot(PluviDJF.geo,low=T)
boxcox(PluviDJF.geo$data~1,lambda=seq(-10,10,l=20))
boxcox(PluviDJF.geo$data~1,lambda=seq(-4,-2,l=20))

trans<-boxcox(PluviDJF.geo$data~1,lambda=seq(-3,-2.5,l=20));trans # dando zoom para ver qual lambda usar
lambdaDJF<-with(trans, x[which.max(y)]);lambdaDJF # Valor máximo de Lambda
dados$DJF
pluviDJF.n<-(((dados$DJF^(lambdaDJF)) - 1)/lambdaDJF);pluviDJF.n #normalizando - (X^lambda)-1/lambda
pluviDJF.inv<-(((lambdaDJF*pluviDJF.n)+1)^(1/lambdaDJF));pluviDJF.inv #inverter e diminuir com 6 para encontrar a varivel original
PluviDJF.geo$data <- pluviDJF.n

boxcox(pluviDJF.n~1,lambda=seq(,l=20)) 
smry(pluviDJF.n)
hist(pluviDJF.n)
qqnorm(PluviDJF.geo$data)
qqline(PluviDJF.geo$data)
qqnorm(pluviDJF.n)
qqline(pluviDJF.n)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normal e sem tendência
##-----------------------------------------------------------------------------##
## Conclusão: há dependência espacial
##-----------------------------------------------------------------------------##
##-----------------------------------------------------------------------------##

## Interpolação espacial (krigagem)

## Modelagem Precipitação média anual
v.pluviDJF<-variog(PluviDJF.geo,max.dist=200,uvec=seq(0, 200, by=20),trend=~coords[,2]+Cota)
plot(v.pluviDJF)

ef.pluviDJF<-eyefit(v.pluviDJF)

PluviDJF <- list()
PluviDJF$lf0 <- likfit(PluviDJF.geo, ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf1 <- likfit(PluviDJF.geo, trend=~Cota,  ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf2 <- likfit(PluviDJF.geo, trend=~coords[,1],ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf3 <- likfit(PluviDJF.geo, trend=~coords[,2],ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf4 <- likfit(PluviDJF.geo, trend="1st",  ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf5 <- likfit(PluviDJF.geo, trend="2nd",  ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf6 <- likfit(PluviDJF.geo, trend=~coords[,1]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf7 <- likfit(PluviDJF.geo, trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf8 <- likfit(PluviDJF.geo, trend=~coords[,1]+poly(Cota,2),ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf9 <- likfit(PluviDJF.geo, trend=~coords[,2]+poly(Cota,2),ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf10 <- likfit(PluviDJF.geo, trend=~poly(coords[,1],2)+poly(Cota,2),ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf11 <- likfit(PluviDJF.geo, trend=~poly(coords[,2],2)+poly(Cota,2), ini=c(3950,22), nug=952,lambda=)
PluviDJF$lf12 <- likfit(PluviDJF.geo, trend=~poly(coords[,1],2)+poly(coords[,2],2)+poly(Cota,2),ini=c(3950,22), nug=952,lambda=)

sort(sapply(PluviDJF,AIC))
sort(sapply(PluviDJF,function(x) x$BIC))
sapply(PluviDJF,function(x) x$value)
plot(PluviDJF.geo,low=T,lambda=lambdaDJF)
plot(PluviDJF.geo,low=T,trend="2nd")
plot(PluviDJF.geo,low=T,trend=~coords[,2]+Cota)

##-----------------------------------------------------------------------------##
## Melhor modelo trend="coords[,1]"
summary(PluviDJF$lf7)

## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

PluviDJFlf7<-list()
PluviDJFlf7$exp <- likfit(PluviDJF.geo, trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJFlf7$gau <- likfit(PluviDJF.geo, cov.model = "gau", trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJFlf7$sph <- likfit(PluviDJF.geo, cov.model = "sph", trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJFlf7$cir <- likfit(PluviDJF.geo, cov.model = "cir", trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJFlf7$kappa1.5 <- likfit(PluviDJF.geo, cov.model = "mat", kappa= 1.5,trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)
PluviDJFlf7$kappa2.5 <- likfit(PluviDJF.geo, cov.model = "mat", kappa= 2.5,trend=~coords[,2]+Cota,ini=c(3950,22), nug=952,lambda=)

##-----------------------------------------------------------------------------

v.pluviDJF<-variog(PluviDJF.geo,max.dist=200,uvec=seq(0, 200, by=10),trend=~coords[,2]+Cota,lambda=)
postscript("Figuras/semivarioPaDJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
plot(v.pluviDJF,xlab="Distância (km)",ylab=expression(gamma(u))) #, ylim = c(0,4e-6),mgp=c(2.5,1,0))
lines(lowess(v.pluviDJF$u,v.pluviDJF$v))
dev.off()

sort(sapply(PluviDJFlf7, AIC))

plot(v.pluviDJF,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(PluviDJFlf7$exp,col=1)
lines(PluviDJFlf7$gau,col=2)
lines(PluviDJFlf7$sph,col=3)
lines(PluviDJFlf7$cir,col=1)
lines(PluviDJFlf7$kappa1.5,col=5)
lines(PluviDJFlf7$kappa2.5,col=6)


##-----------------------------------------------------------------------------
## modelo escolhido, pela MV o melhor foi o circular
## escolhendo o modelo circular
summary(PluviDJFlf7$sph)
names(PluviDJFlf7$sph)
PluviDJFlf7$sph$sigmasq
PluviDJFlf7$sph$nugget
PluviDJFlf7$sph$phi
PluviDJFlf7$sph$parameters.summary

##-----------------------------------------------------------------------------
## Predição na área

##rm(list = ls())

grid <- pred_grid(c(200,800),c(6700,7200), by=1)
points(PluviDJF.geo)
points(grid, pch=19, cex=0.25, col=2)
##rm(gr)
gr <- locations.inside(grid, rbind(limite,floripa))
points(gr, pch=19, cex=0.24,col=4)
kc.pluviDJF <- krige.conv(PluviDJF.geo, loc=grid, krige=krige.control(obj=PluviDJFlf7$sph))
attributes(kc.pluviDJF)
##kc.pluviDJF$predict<-(backtransform.moments(lambda=lambdaDJF,mean=kc.pluviDJF$predict,variance=kc.pluviDJF$krige.var)$mean)

PluviDJF_krige <- cbind(grid*1000,kc.pluviDJF$predict)
write.table(PluviDJF_krige,"PluviDJF.ascii",col.names = F,row.names=F,quote=F)
rm(PluviDJF.geo)
PluviDJF.geo$borders <- limite
##-----------------------------------------------------------------------------
## Validação cruzada
xv.pluviDJF<-xvalid(PluviDJF.geo,model=PluviDJFlf7$sph)
names(xv.pluviDJF)
xv.pluviDJF$data
dados$DJF

##xv.pluviDJF$data <- dados$DJF
##xv.pluviDJF$predicted <- (backtransform.moments(lambda=lambdaDJF,mean=xv.pluviDJF$predicted,variance=xv.pluviDJF$krige.var)$mean)

plot(xv.pluviDJF$predicted,xv.pluviDJF$data) #,xlim=c(0,1.2),ylim=c(0,1.2))
abline(0,1)
abline(lm(xv.pluviDJF$data~xv.pluviDJF$predicted)$coef[1],lm(xv.pluviDJF$data~xv.pluviDJF$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.pluviDJF)


smry(PluviDJF.geo$data)


postscript("Figuras/xv_pluviDJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.pluviDJF$predicted,xv.pluviDJF$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.pluviDJF$std.error,breaks = 20,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()


postscript("Figuras/kr_pluviDJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.pluviDJF,val=kc.pluviDJF$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.pluviDJF, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)


dev.off()

summary(limite)

pdf("Figuras/sd.pluviDJF.pdf",onefile = T, width=20/2.54, height=20/2.54)

image(kc.pluviDJF, val=kc.pluviDJF$krige.var,col=grey(seq(1,0.2,len=21)),x.leg=c(250,450),y.leg=c(6800,6810),xlab = "Longitude (km)",ylab = "Latitude (km)")

dev.off()


##-----------------------------------------------------------------------------              

## Pluvi Outono

## Análise exploratória geoestatística
head(dados)
names(dados)
PluviMAM.geo<-as.geodata(dados,coords.col = 1:2, data.col = 5,covar.col = 8,)

plot(PluviMAM.geo,low=T)
## PluviMAM.geo$borders<-limite

postscript("Figuras/BoxCox_PluviMAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper="special")
par(mfrow = c(1,2))
hist(PluviMAM.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(PluviMAM.geo$data~1,lambda = seq(-3,3),ylab = "Log-Verossimilhança")
dev.off()

names(PluviMAM.geo)
shapiro.test(PluviMAM.geo$data)
plot(PluviMAM.geo,low=T)
boxcox(PluviMAM.geo$data~1,lambda=seq(-5,5,l=20)) 
boxcox(PluviMAM.geo$data~1,lambda=seq(-2,2,l=20))

trans<-boxcox(PluviMAM.geo$data~1,lambda=seq(-1,0.5,l=20));trans # dando zoom para ver qual lambda usar
lambdaMAM<-with(trans, x[which.max(y)]);lambdaMAM # Valor máximo de Lambda

pluviMAM.n<-(((dados$ANUAL^(lambdaMAM)) - 1)/lambdaMAM);pluviMAM.n #normalizando - (X^lambda)-1/lambda
pluviMAM.inv<-(((lambdaMAM*pluviMAM.n)+1)^(1/lambdaMAM));pluviMAM.inv #inverter e diminuir com 6 para encontrar a varivel original
PluviMAM.geo$data
PluviMAM.geo$data <- pluviMAM.n

boxcox(pluviMAM.n~1,lambda=seq(-5,5,l=20)) 
smry(pluviMAM.n)
hist(PluviMAM.geo$data)
qqnorm(PluviMAM.geo$data)
qqline(PluviMAM.geo$data)
qqnorm(pluviMAM.inv)
qqline(pluviMAM.inv)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normal e sem tendência
##-----------------------------------------------------------------------------##
## Conclusão: há dependência espacial
##-----------------------------------------------------------------------------##
##-----------------------------------------------------------------------------##

## Interpolação espacial (krigagem)

## Modelagem Precipitação média anual
v.pluviMAM<-variog(PluviMAM.geo,max.dist=500,uvec=seq(0, 500, by=10),trend="2nd")

ef.pluviMAM<-eyefit(v.pluviMAM)

PluviMAM <- list()
PluviMAM$lf0 <- likfit(PluviMAM.geo, ini=c(1900,54), nug=292)
PluviMAM$lf1 <- likfit(PluviMAM.geo, trend=~Cota,  ini=c(1900,54), nug=292)
PluviMAM$lf2 <- likfit(PluviMAM.geo, trend=~coords[,1],ini=c(1900,54), nug=292)
PluviMAM$lf3 <- likfit(PluviMAM.geo, trend=~coords[,2],ini=c(1900,54), nug=292)
PluviMAM$lf4 <- likfit(PluviMAM.geo, trend="1st",  ini=c(1900,54), nug=292)
PluviMAM$lf5 <- likfit(PluviMAM.geo, trend="2nd",  ini=c(1900,54), nug=292)
PluviMAM$lf6 <- likfit(PluviMAM.geo, trend=~coords[,1]+Cota,ini=c(1900,54), nug=292)
PluviMAM$lf7 <- likfit(PluviMAM.geo, trend=~coords[,2]+Cota,ini=c(1900,54), nug=292)
PluviMAM$lf8 <- likfit(PluviMAM.geo, trend=~coords[,1]+poly(Cota,2),ini=c(1900,54), nug=292)
PluviMAM$lf9 <- likfit(PluviMAM.geo, trend=~coords[,2]+poly(Cota,2),ini=c(1900,54), nug=292)
PluviMAM$lf10 <- likfit(PluviMAM.geo, trend=~poly(coords[,1],2)+poly(Cota,2),ini=c(1900,54), nug=292)
PluviMAM$lf11 <- likfit(PluviMAM.geo, trend=~poly(coords[,2],2)+poly(Cota,2), ini=c(1900,54), nug=292)
PluviMAM$lf12 <- likfit(PluviMAM.geo, trend=~poly(coords[,1],2)+poly(coords[,2],2)+poly(Cota,2),ini=c(1900,54), nug=292)

sort(sapply(PluviMAM,AIC))
sort(sapply(PluviMAM,function(x) x$BIC))
sapply(PluviMAM,function(x) x$value)
plot(PluviMAM.geo,low=T,trend="1st")
plot(PluviMAM.geo,low=T,trend=)
plot(PluviMAM.geo,low=T,trend="2nd")
plot(PluviMAM.geo,low=T,trend=~coords[,2]+Cota)

##-----------------------------------------------------------------------------##
## Melhor modelo trend="coords[,1]"
summary(PluviMAM$lf5)

## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

PluviMAMlf5<-list()
PluviMAMlf5$exp <- likfit(PluviMAM.geo, trend="2nd",  ini=c(1900,54), nug=292)
PluviMAMlf5$gau <- likfit(PluviMAM.geo, cov.model = "gau", trend="2nd",  ini=c(1900,54), nug=292)
PluviMAMlf5$sph <- likfit(PluviMAM.geo, cov.model = "sph", trend="2nd",  ini=c(1900,54), nug=292)
PluviMAMlf5$cir <- likfit(PluviMAM.geo, cov.model = "cir", trend="2nd",  ini=c(1900,54), nug=292)
PluviMAMlf5$kappa1.5 <- likfit(PluviMAM.geo, cov.model = "mat", kappa= 1.5,trend=~poly(coords[,1],2)+poly(coords[,2],2),  ini=c(1900,54), nug=292)
PluviMAMlf5$kappa2.5 <- likfit(PluviMAM.geo, cov.model = "mat", kappa= 2.5,trend="2nd",  ini=c(1900,54), nug=292)

##-----------------------------------------------------------------------------

v.pluviMAM<-variog(PluviMAM.geo,max.dist=200,uvec=seq(0, 500, by=10),trend="2nd")
postscript("Figuras/semivarioPaMAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
plot(v.pluviMAM,xlab="Distância (km)",ylab=expression(gamma(u)))#, ylim = c(0,4e-6),mgp=c(2.5,1,0))
lines(lowess(v.pluviMAM$u,v.pluviMAM$v))
dev.off()

sort(sapply(PluviMAMlf5, AIC))

plot(v.pluviMAM,xlab="Distance (km)",ylab=expression(gamma(u)),ylim =c(0,5000))
lines(PluviMAMlf5$exp,col=1)
lines(PluviMAMlf5$gau,col=2)
lines(PluviMAMlf5$sph,col=3)
lines(PluviMAMlf5$cir,col=1)
lines(PluviMAMlf5$kappa1.5,col=5)
lines(PluviMAMlf5$kappa2.5,col=6)


##-----------------------------------------------------------------------------
## modelo escolhido, pela MV o melhor foi o circular
## escolhendo o modelo circular
summary(PluviMAMlf5$kappa1.5)
names(PluviMAMlf7$kappa1.5)
PluviMAMlf5$kappa1.5$sigmasq
PluviMAMlf5$kappa1.5$nugget
PluviMAMlf5$kappa1.5$beta
PluviMAMlf5$kappa1.5$phi
PluviMAMlf5$kappa1.5$parameters.summary


##-----------------------------------------------------------------------------
## Predição na área

##rm(list = ls())

grid <- pred_grid(c(200,800),c(6700,7200), by=1)
points(PluviMAM.geo)
points(grid, pch=19, cex=0.25, col=2)
##rm(gr)
gr <- locations.inside(grid, rbind(limite,floripa))
points(gr, pch=19, cex=0.24,col=4)
kc.pluviMAM <- krige.conv(PluviMAM.geo, loc=grid, krige=krige.control(obj=PluviMAMlf5$kappa1.5))
attributes(kc.pluviMAM)
##kc.pluviMAM$predict<-(backtransform.moments(lambda=lambdaMAM,mean=kc.pluviMAM$predict,variance=kc.pluviMAM$krige.var)$mean)
summary(kc.pluviMAM$predict)

PluviMAM_krige <- cbind(grid*1000,kc.pluviMAM$predict)
write.table(PluviMAM_krige,"PluviMAM.ascii",col.names = F,row.names=F,quote=F)
rm(PluviMAM.geo)
PluviMAM.geo$borders <- limite
##-----------------------------------------------------------------------------
## Validação cruzada
xv.pluviMAM<-xvalid(PluviMAM.geo,model=PluviMAMlf5$kappa1.5)
names(xv.pluviMAM)
xv.pluviMAM$data
dados$MAM
xv.pluviMAM$predicted
##xv.pluviMAM$data <- dados$MAM
##xv.pluviMAM$predicted <- (backtransform.moments(lambda=lambdaMAM,mean=xv.pluviMAM$predicted,variance=xv.pluviMAM$krige.var)$mean)


plot(xv.pluviMAM$predicted,xv.pluviMAM$data) #,xlim=c(0,1.2),ylim=c(0,1.2))
abline(0,1)
abline(lm(xv.pluviMAM$data~xv.pluviMAM$predicted)$coef[1],lm(xv.pluviMAM$data~xv.pluviMAM$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.pluviMAM)


smry(PluviMAM.geo$data)


postscript("Figuras/xv_pluviMAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.pluviMAM$predicted,xv.pluviMAM$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.pluviMAM$std.error,breaks = 20,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()


postscript("Figuras/kr_pluviMAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.pluviMAM,val=kc.pluviMAM$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.pluviMAM, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)


dev.off()

summary(limite)

pdf("Figuras/sd.pluviMAM.pdf",onefile = T, width=20/2.54, height=20/2.54)

image(kc.pluviMAM, val=kc.pluviMAM$krige.var,col=grey(seq(1,0.2,len=21)),x.leg=c(250,450),y.leg=c(6800,6810),xlab = "Longitude (km)",ylab = "Latitude (km)")

dev.off()


##-----------------------------------------------------------------------------              

## Pluvi Inverno

## Análise exploratória geoestatística
head(dados)
names(dados)
PluviJJA.geo<-as.geodata(dados,coords.col = 1:2, data.col = 6,covar.col = 8,)

plot(PluviJJA.geo,low=T)
## PluviJJA.geo$borders<-limite


postscript("Figuras/BoxCox_PluviJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper="special")
par(mfrow = c(1,2))
hist(PluviJJA.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(PluviJJA.geo$data~1,lambda = seq(-3,3),ylab = "Log-Verossimilhança")
dev.off()


plot(PluviJJA.geo,low=T)
boxcox(PluviJJA.geo$data~1,lambda=seq(-5,5,l=20)) 
boxcox(PluviJJA.geo$data~1,lambda=seq(-2,-0,l=20))

trans<-boxcox(PluviJJA.geo$data~1,lambda=seq(-1,0.5,l=20));trans # dando zoom para ver qual lambda usar
lambdaJJA<-with(trans, x[which.max(y)]);lambdaJJA # Valor máximo de Lambda

pluviJJA.n<-(((dados$ANUAL^(lambdaJJA)) - 1)/lambdaJJA);pluviJJA.n #normalizando - (X^lambda)-1/lambda
pluviJJA.inv<-(((lambdaJJA*pluviJJA.n)+1)^(1/lambdaJJA));pluviJJA.inv #inerter e diminuir com 6 para encontrar a varivel original
PluviJJA.geo$data
PluviJJA.geo$data <- pluviJJA.n

boxcox(pluviJJA.n~1,lambda=seq(-5,5,l=20)) 
smry(pluviJJA.n)
hist(PluviJJA.geo$data)
qqnorm(PluviJJA.geo$data)
qqline(PluviJJA.geo$data)
qqnorm(pluviJJA.inv)
qqline(pluviJJA.inv)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normal e sem tendência
##-----------------------------------------------------------------------------##
## Conclusão: há dependência espacial
##-----------------------------------------------------------------------------##
##-----------------------------------------------------------------------------##

## Interpolação espacial (krigagem)

## Modelagem Precipitação média anual
v.pluviJJA<-variog(PluviJJA.geo,max.dist=500,uvec=seq(0, 500, by=10),trend="1st")
ef.pluviJJA<-eyefit(v.pluviJJA)

PluviJJA <- list()
PluviJJA$lf0 <- likfit(PluviJJA.geo, ini=c(1900,54), nug=292)
PluviJJA$lf1 <- likfit(PluviJJA.geo, trend=~Cota,  ini=c(1900,54), nug=292)
PluviJJA$lf2 <- likfit(PluviJJA.geo, trend=~coords[,1],ini=c(1900,54), nug=292)
PluviJJA$lf3 <- likfit(PluviJJA.geo, trend=~coords[,2],ini=c(1900,54), nug=292)
PluviJJA$lf4 <- likfit(PluviJJA.geo, trend="1st",  ini=c(1900,54), nug=292)
PluviJJA$lf5 <- likfit(PluviJJA.geo, trend="2nd",  ini=c(1900,54), nug=292)
PluviJJA$lf6 <- likfit(PluviJJA.geo, trend=~coords[,1]+Cota,ini=c(1900,54), nug=292)
PluviJJA$lf7 <- likfit(PluviJJA.geo, trend=~coords[,2]+Cota,ini=c(1900,54), nug=292)
PluviJJA$lf8 <- likfit(PluviJJA.geo, trend=~coords[,1]+poly(Cota,2),ini=c(1900,54), nug=292)
PluviJJA$lf9 <- likfit(PluviJJA.geo, trend=~coords[,2]+poly(Cota,2),ini=c(1900,54), nug=292)
PluviJJA$lf10 <- likfit(PluviJJA.geo, trend=~poly(coords[,1],2)+poly(Cota,2),ini=c(1900,54), nug=292)
PluviJJA$lf11 <- likfit(PluviJJA.geo, trend=~poly(coords[,2],2)+poly(Cota,2), ini=c(1900,54), nug=292)
PluviJJA$lf12 <- likfit(PluviJJA.geo, trend=~poly(coords[,1],2)+poly(coords[,2],2)+poly(Cota,2),ini=c(1900,54), nug=292)

sort(sapply(PluviJJA,AIC))
sort(sapply(PluviJJA,function(x) x$BIC))
sapply(PluviJJA,function(x) x$value)
plot(PluviJJA.geo,low=T,trend="1st")
plot(PluviJJA.geo,low=T,trend=)
plot(PluviJJA.geo,low=T,trend="2nd")
plot(PluviJJA.geo,low=T,trend=~coords[,1]+Cota)

##-----------------------------------------------------------------------------##
## Melhor modelo trend="coords[,1]"
summary(PluviJJA$lf6)

## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

PluviJJAlf6<-list()
PluviJJAlf6$exp <- likfit(PluviJJA.geo, trend=~coords[,1]+Cota,  ini=c(1900,54), nug=292)
PluviJJAlf6$gau <- likfit(PluviJJA.geo, cov.model = "gau", trend=~coords[,1]+Cota,  ini=c(1900,54), nug=292)
PluviJJAlf6$sph <- likfit(PluviJJA.geo, cov.model = "sph", trend=~coords[,1]+Cota,  ini=c(1900,54), nug=292)
PluviJJAlf6$cir <- likfit(PluviJJA.geo, cov.model = "cir", trend=~coords[,1]+Cota,  ini=c(1900,54), nug=292)
PluviJJAlf6$kappa1.5 <- likfit(PluviJJA.geo, cov.model = "mat", kappa= 1.5,trend=~coords[,1]+Cota,  ini=c(1900,54), nug=292)
PluviJJAlf6$kappa2.5 <- likfit(PluviJJA.geo, cov.model = "mat", kappa= 2.5,trend=~coords[,1]+Cota,  ini=c(1900,54), nug=292)

##-----------------------------------------------------------------------------

v.pluviJJA<-variog(PluviJJA.geo,max.dist=300,uvec=seq(0, 300, by=10),trend=~coords[,1]+Cota)
postscript("Figuras/semivarioPaJJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
plot(v.pluviJJA,xlab="Distância (km)",ylab=expression(gamma(u)))#, ylim = c(0,4e-6),mgp=c(2.5,1,0))
lines(lowess(v.pluviJJA$u,v.pluviJJA$v))
dev.off()

sort(sapply(PluviJJAlf6, AIC))

plot(v.pluviJJA,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,5000))
lines(PluviJJAlf6$exp,col=1)
lines(PluviJJAlf6$gau,col=2)
lines(PluviJJAlf6$sph,col=3)
lines(PluviJJAlf6$cir,col=1)
lines(PluviJJAlf6$kappa1.5,col=5)
lines(PluviJJAlf6$kappa2.5,col=6)


##-----------------------------------------------------------------------------
## modelo escolhido, pela MV o melhor foi o circular
## escolhendo o modelo circular
summary(PluviJJAlf6$exp)
PluviJJAlf6$exp$sigmasq
PluviJJAlf6$exp$nugget
PluviJJAlf6$exp$beta
PluviJJAlf6$exp$phi
PluviJJAlf6$exp$parameters.summary


##-----------------------------------------------------------------------------
## Predição na área

##rm(list = ls())

grid <- pred_grid(c(200,800),c(6700,7200), by=1)
points(PluviJJA.geo)
points(grid, pch=19, cex=0.25, col=2)
##rm(gr)
gr <- locations.inside(grid, rbind(limite,floripa))
points(gr, pch=19, cex=0.24,col=4)
kc.pluviJJA <- krige.conv(PluviJJA.geo, loc=grid, krige=krige.control(obj=PluviJJAlf6$exp))
attributes(kc.pluviJJA)
##kc.pluviJJA$predict<-(backtransform.moments(lambda=lambdaJJA,mean=kc.pluviJJA$predict,variance=kc.pluviJJA$krige.var)$mean)
summary(kc.pluviJJA$predict)

PluviJJA_krige <- cbind(grid*1000,kc.pluviJJA$predict)
write.table(PluviJJA_krige,"PluviJJA.ascii",col.names = F,row.names=F,quote=F)
PluviJJA.geo$borders
PluviJJA.geo$borders <- limite
##-----------------------------------------------------------------------------
## Validação cruzada
xv.pluviJJA<-xvalid(PluviJJA.geo,model=PluviJJAlf6$exp)
names(xv.pluviJJA)
xv.pluviJJA$data
dados$JJA
xv.pluviJJA$predicted
##xv.pluviJJA$data <- dados$JJA
##xv.pluviJJA$predicted <- (backtransform.moments(lambda=lambdaJJA,mean=xv.pluviJJA$predicted,variance=xv.pluviJJA$krige.var)$mean)


plot(xv.pluviJJA$predicted,xv.pluviJJA$data) #,xlim=c(0,1.2),ylim=c(0,1.2))
abline(0,1)
abline(lm(xv.pluviJJA$data~xv.pluviJJA$predicted)$coef[1],lm(xv.pluviJJA$data~xv.pluviJJA$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.pluviJJA)


smry(PluviJJA.geo$data)


postscript("Figuras/xv_pluviJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.pluviJJA$predicted,xv.pluviJJA$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.pluviJJA$std.error,breaks = 20,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()


postscript("Figuras/kr_pluviJJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.pluviJJA,val=kc.pluviJJA$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.pluviJJA, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)


dev.off()

summary(limite)

pdf("Figuras/sd.pluviJJA.pdf",onefile = T, width=20/2.54, height=20/2.54)

image(kc.pluviJJA, val=kc.pluviJJA$krige.var,col=grey(seq(1,0.2,len=21)),x.leg=c(250,450),y.leg=c(6800,6810),xlab = "Longitude (km)",ylab = "Latitude (km)")

dev.off()

##-----------------------------------------------------------------------------##


## Pluvi Primavera

## Análise exploratória geoestatística
head(dados)
names(dados)
PluviJJA.geo<-as.geodata(dados,coords.col = 1:2, data.col = 7,covar.col = 8,)

plot(PluviJJA.geo,low=T)
## PluviJJA.geo$borders<-limite


postscript("Figuras/BoxCox_PluviJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper="special")
par(mfrow = c(1,2))
hist(PluviJJA.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(PluviJJA.geo$data~1,lambda = seq(-3,3),ylab = "Log-Verossimilhança")
dev.off()



plot(PluviJJA.geo,low=T)
boxcox(PluviJJA.geo$data~1,lambda=seq(-5,5,l=20)) 
boxcox(PluviJJA.geo$data~1,lambda=seq(-2,-0,l=20))

trans<-boxcox(PluviJJA.geo$data~1,lambda=seq(-1,0.5,l=20));trans # dando zoom para ver qual lambda usar
lambdaJJA<-with(trans, x[which.max(y)]);lambdaJJA # Valor máximo de Lambda

pluviJJA.n<-(((dados$ANUAL^(lambdaJJA)) - 1)/lambdaJJA);pluviJJA.n #normalizando - (X^lambda)-1/lambda
pluviJJA.inv<-(((lambdaJJA*pluviJJA.n)+1)^(1/lambdaJJA));pluviJJA.inv #inverter e diminuir com 6 para encontrar a varivel original
PluviJJA.geo$data
PluviJJA.geo$data <- pluviJJA.n

boxcox(pluviJJA.n~1,lambda=seq(-5,5,l=20)) 
smry(pluviJJA.n)
hist(PluviJJA.geo$data)
qqnorm(PluviJJA.geo$data)
qqline(PluviJJA.geo$data)
qqnorm(pluviJJA.inv)
qqline(pluviJJA.inv)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normal e sem tendência
##-----------------------------------------------------------------------------##
## Conclusão: há dependência espacial
##-----------------------------------------------------------------------------##
##-----------------------------------------------------------------------------##

## Interpolação espacial (krigagem)

## Modelagem Precipitação média anual
v.pluviJJA<-variog(PluviJJA.geo,max.dist=500,uvec=seq(0, 500, by=10),trend="1st")
ef.pluviJJA<-eyefit(v.pluviJJA)

PluviJJA <- list()
PluviJJA$lf0 <- likfit(PluviJJA.geo, ini=c(1900,54), nug=292)
PluviJJA$lf1 <- likfit(PluviJJA.geo, trend=~Cota,  ini=c(1900,54), nug=292)
PluviJJA$lf2 <- likfit(PluviJJA.geo, trend=~coords[,1],ini=c(1900,54), nug=292)
PluviJJA$lf3 <- likfit(PluviJJA.geo, trend=~coords[,2],ini=c(1900,54), nug=292)
PluviJJA$lf4 <- likfit(PluviJJA.geo, trend="1st",  ini=c(1900,54), nug=292)
PluviJJA$lf5 <- likfit(PluviJJA.geo, trend="2nd",  ini=c(1900,54), nug=292)
PluviJJA$lf6 <- likfit(PluviJJA.geo, trend=~coords[,1]+Cota,ini=c(1900,54), nug=292)
PluviJJA$lf7 <- likfit(PluviJJA.geo, trend=~coords[,2]+Cota,ini=c(1900,54), nug=292)
PluviJJA$lf8 <- likfit(PluviJJA.geo, trend=~coords[,1]+poly(Cota,2),ini=c(1900,54), nug=292)
PluviJJA$lf9 <- likfit(PluviJJA.geo, trend=~coords[,2]+poly(Cota,2),ini=c(1900,54), nug=292)
PluviJJA$lf10 <- likfit(PluviJJA.geo, trend=~poly(coords[,1],2)+poly(Cota,2),ini=c(1900,54), nug=292)
PluviJJA$lf11 <- likfit(PluviJJA.geo, trend=~poly(coords[,2],2)+poly(Cota,2), ini=c(1900,54), nug=292)
PluviJJA$lf12 <- likfit(PluviJJA.geo, trend=~poly(coords[,1],2)+poly(coords[,2],2)+poly(Cota,2),ini=c(1900,54), nug=292)

sort(sapply(PluviJJA,AIC))
sort(sapply(PluviJJA,function(x) x$BIC))
sapply(PluviJJA,function(x) x$value)
plot(PluviJJA.geo,low=T,trend="1st")
plot(PluviJJA.geo,low=T,trend=)
plot(PluviJJA.geo,low=T,trend="2nd")
plot(PluviJJA.geo,low=T,trend=~coords[,1]+Cota)

##-----------------------------------------------------------------------------##
## Melhor modelo trend="coords[,1]"
summary(PluviJJA$lf4)

## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

PluviJJAlf4<-list()
PluviJJAlf4$exp <- likfit(PluviJJA.geo, trend="1st",  ini=c(1900,54), nug=292)
PluviJJAlf4$gau <- likfit(PluviJJA.geo, cov.model = "gau", trend="1st",  ini=c(1900,54), nug=292)
PluviJJAlf4$sph <- likfit(PluviJJA.geo, cov.model = "sph", trend="1st",  ini=c(1900,54), nug=292)
PluviJJAlf4$cir <- likfit(PluviJJA.geo, cov.model = "cir", trend="1st",  ini=c(1900,54), nug=292)
PluviJJAlf4$kappa1.5 <- likfit(PluviJJA.geo, cov.model = "mat", kappa= 1.5,trend="1st",  ini=c(1900,54), nug=292)
PluviJJAlf4$kappa2.5 <- likfit(PluviJJA.geo, cov.model = "mat", kappa= 2.5,trend="1st",  ini=c(1900,54), nug=292)

##-----------------------------------------------------------------------------

v.pluviJJA<-variog(PluviJJA.geo,max.dist=200,uvec=seq(0, 200, by=10),trend="1st")
postscript("Figuras/semivarioPaJJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
plot(v.pluviJJA,xlab="Distância (km)",ylab=expression(gamma(u)))#, ylim = c(0,4e-6),mgp=c(2.5,1,0))
lines(lowess(v.pluviJJA$u,v.pluviJJA$v))
dev.off()

sort(sapply(PluviJJAlf4, AIC))

plot(v.pluviJJA,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,5000))
lines(PluviJJAlf4$exp,col=1)
lines(PluviJJAlf4$gau,col=2)
lines(PluviJJAlf4$sph,col=3)
lines(PluviJJAlf4$cir,col=1)
lines(PluviJJAlf4$kappa1.5,col=5)
lines(PluviJJAlf4$kappa2.5,col=6)


##-----------------------------------------------------------------------------
## modelo escolhido, pela MV o melhor foi o circular
## escolhendo o modelo circular
summary(PluviJJAlf4$sph)
names(PluviJJAlf4$sph)
PluviJJAlf4$sph$sigmasq
PluviJJAlf4$sph$nugget
PluviJJAlf4$sph$beta
PluviJJAlf4$sph$phi
PluviJJAlf4$sph$parameters.summary


##-----------------------------------------------------------------------------
## Predição na área

##rm(list = ls())

grid <- pred_grid(c(200,800),c(6700,7200), by=1)
points(PluviJJA.geo)
points(grid, pch=19, cex=0.25, col=2)
##rm(gr)
gr <- locations.inside(grid, rbind(limite,floripa))
points(gr, pch=19, cex=0.24,col=4)
kc.pluviJJA <- krige.conv(PluviJJA.geo, loc=grid, krige=krige.control(obj=PluviJJAlf4$sph))
attributes(kc.pluviJJA)
##kc.pluviJJA$predict<-(backtransform.moments(lambda=lambdaJJA,mean=kc.pluviJJA$predict,variance=kc.pluviJJA$krige.var)$mean)
summary(kc.pluviJJA$predict)

PluviJJA_krige <- cbind(grid*1000,kc.pluviJJA$predict)
write.table(PluviJJA_krige,"PluviJJA.ascii",col.names = F,row.names=F,quote=F)
PluviJJA.geo$borders
PluviJJA.geo$borders <- limite
##-----------------------------------------------------------------------------
## Validação cruzada
xv.pluviJJA<-xvalid(PluviJJA.geo,model=PluviJJAlf4$sph)
names(xv.pluviJJA)
xv.pluviJJA$data
dados$JJA
xv.pluviJJA$predicted
##xv.pluviJJA$data <- dados$JJA
##xv.pluviJJA$predicted <- (backtransform.moments(lambda=lambdaJJA,mean=xv.pluviJJA$predicted,variance=xv.pluviJJA$krige.var)$mean)


plot(xv.pluviJJA$predicted,xv.pluviJJA$data) #,xlim=c(0,1.2),ylim=c(0,1.2))
abline(0,1)
abline(lm(xv.pluviJJA$data~xv.pluviJJA$predicted)$coef[1],lm(xv.pluviJJA$data~xv.pluviJJA$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.pluviJJA)


smry(PluviJJA.geo$data)


postscript("Figuras/xv_pluviJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.pluviJJA$predicted,xv.pluviJJA$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.pluviJJA$std.error,breaks = 20,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()


postscript("Figuras/kr_pluviJJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.pluviJJA,val=kc.pluviJJA$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.pluviJJA, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)


dev.off()

summary(limite)

pdf("Figuras/sd.pluviJJA.pdf",onefile = T, width=20/2.54, height=20/2.54)

image(kc.pluviJJA, val=kc.pluviJJA$krige.var,col=grey(seq(1,0.2,len=21)),x.leg=c(250,450),y.leg=c(6800,6810),xlab = "Longitude (km)",ylab = "Latitude (km)")

dev.off()
##-----------------------------------------------------------------------------##

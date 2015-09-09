
## Configuração Inicial

require(RColorBrewer) # Para Paletas de cores
require(hydroTSM) 
colour<-colorRampPalette(c("white","blue","dark Blue"))
require(MASS)
require(rgdal)
require(geoR)
require(tikzDevice)
require(extrafont)
fonts()
loadfonts()
##font_install("fontcm")

getwd()
setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/IEBgeo")
##rm(list = ls()) ## Para remover todos os objetos

##Ver documentação
cov.spatial()
trend.spatial()
options(OutDec=",",digits=10,tikzSanitizeCharacters =c('$'))


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
#[-22,]
dir()
##dados <- read.csv("Param_geo.csv", head=T,dec=",",sep= ";")[-22,];dados$Centroid.Y<-dados$Centroid.Y/1000;dados$Centroid.X<-dados$Centroid.X/1000

dados <- read.csv("ParamANO.csv", head=T,dec=",",sep= "")[,-c(1,4,5,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC.csv",h=T,sep=",",dec = ".");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000

summary(limite)
##write.table(dados,"Resultados.txt")

names(dados)
dim(dados)
## lendo dados tipo shapefile
## limite <- readShapePoly("SC.shp")
## Amarzenando as bordas do Estado em um objeto
## bordas<-rbind(limite@polygons[[1]]@Polygons[[1]]@coords,limite@polygons[[2]]@Polygons[[1]]@coords)


##-----------------------------------------------------------------------------##

## verificando argumentos e documentação da função utilizada
args(read.table)
help(read.table)

## visualizando dados importados
View(dados) 
dados

##-----------------------------------------------------------------------------##
## Parametros ANO mu

## Análise exploratória geoestatística

names(dados)

dim(dados)
IEB.ANO.geo<-as.geodata(dados,coords.col = 1:2, data.col = 5,borders=T)
IEB.ANO.geo$borders<-limite

plot(IEB.ANO.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

boxcox(IEB.ANO.geo$data~1) 

## 
##trans<-boxcox(ANO.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.ANO<-with(trans, x[which.max(y)]);lambda.mu.ANO # Valor máximo de Lambda
## 
##Mu.ANO.geo$data<-(((ANO.mu^(lambda.mu.ANO)) - 1)/lambda.mu.ANO);Mu.ANO.geo$data #normalizando - (X^lambda.mu.ANO)-1/lambda.mu.ANO
##ANO.mu.inv<-(((lambda.mu.ANO*Mu.ANO.geo$data)+1)^(1/lambda.mu.ANO))-5;ANO.mu.inv #inverter e diminuir com 5 para encontrar a varivel original

bcIEBANO <- boxcox(IEB.ANO.geo$data~1,plotit = T,lambda = seq(-3,3),0.1)

tikz("Figuras/plot.tex",onefile = T, width=18/2.54, height=12/2.54)

##pdf("Figuras/BoxCox_IEBANO.pdf",onefile = T,family="CM Roman Greek", width=25/2.54, height=15/2.54,paper = "special")

par(mfrow = c(1,2))
hist(IEB.ANO.geo$data,main="", xlab="",ylab = "Frequ\\^{e}ncia")
plot(bcIEBANO,ylab = "Log-Verossimilhan\\c{c}a",xlab = "$\\lambda$",type="l")
abline(v=range(bcIEBANO$x[bcIEBANO$y > max(bcIEBANO$y)-qchisq(0.95,1)/2]),lty=2)
legend('topleft',leg = 'IC 95 \\%',lty = 2,bty = 'n')

dev.off()



#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.IEB.ANO<-variog(IEB.ANO.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.IEB.ANO,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
names(dados)
(trend_IEB.ANO <- dados[,c(3,4,6:23)])
names(trend_IEB.ANO)
ncol(trend_IEB.ANO)

## Tendência linear
IEB.ANO <- list()
for(i in 1:ncol(trend_IEB.ANO)){
IEB.ANO[[i]] <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend_IEB.ANO[,i],ini=c(0.04,23.35), nug=0.02,)
IEB.ANO[[ncol(trend_IEB.ANO)+1]] <- likfit(IEB.ANO.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}

sapply(IEB.ANO, AIC)

names(trend_IEB.ANO)
summary(IEB.ANO[[order(sapply(IEB.ANO, AIC))[1]]])
trend1_IEB.ANO<-trend_IEB.ANO[,order(sapply(IEB.ANO, AIC))[1]];trend_IEB.ANO <- trend_IEB.ANO[,-order(sapply(IEB.ANO, AIC))[1]]

IEB.ANO1 <- list()
for(i in 1:ncol(trend_IEB.ANO)){
IEB.ANO1[[i]] <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend1_IEB.ANO+trend_IEB.ANO[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.ANO)
sapply(IEB.ANO1, AIC)
summary(IEB.ANO1[[order(sapply(IEB.ANO1, AIC))[1]]])

## Tendência quadrática
(trend_IEB.ANO <- dados[,c(3,4,6:23)])
names(trend_IEB.ANO)

IEB.ANO2 <- list()
for(i in 1:ncol(trend_IEB.ANO)){
IEB.ANO2[[i]] <- likfit(IEB.ANO.geo, cov.model="gau", trend =~poly(trend_IEB.ANO[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(IEB.ANO2, AIC)
summary(IEB.ANO2[[order(sapply(IEB.ANO2, AIC))[1]]])
trend3_IEB.ANO<-trend_IEB.ANO[,order(sapply(IEB.ANO2, AIC))[1]];trend_IEB.ANO <- trend_IEB.ANO[,-order(sapply(IEB.ANO2, AIC))[1]]

IEB.ANO3 <- list()
for(i in 1:ncol(trend_IEB.ANO)){
IEB.ANO3[[i]] <- likfit(IEB.ANO.geo, cov.model="gau", trend =~poly(trend3_IEB.ANO,2)+poly(trend_IEB.ANO[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.ANO)
sapply(IEB.ANO3, AIC)
summary(IEB.ANO3[[order(sapply(IEB.ANO3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_IEB.ANO <- dados[,c(3,4,6:23)]
names(trend_IEB.ANO)

IEB.ANOfim <- list()
IEB.ANOfim$lf1 <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend_IEB.ANO[,11]+trend_IEB.ANO[,18],ini=c(0.04,23.35), nug=0.02,)
IEB.ANOfim$lf2 <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend_IEB.ANO[,11]+poly(trend_IEB.ANO[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.ANOfim$lf3 <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.ANOfim$lf4 <- likfit(IEB.ANO.geo, cov.model="gau", trend =~poly(trend_IEB.ANO[,18],2)+poly(trend_IEB.ANO[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.ANOfim$lf5 <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend_IEB.ANO[,11]+poly(trend_IEB.ANO[,18],2),ini=c(0.04,23.35), nug=0.02,)
IEB.ANOfim$lf6 <- likfit(IEB.ANO.geo, cov.model="gau", trend =~trend_IEB.ANO[,6]+poly(trend_IEB.ANO[,18],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(IEB.ANOfim, AIC))
plot(IEB.ANO.geo,low=T,trend=~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2))

tikz()
#postscript("Figuras/Vario-IEB_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.IEB.ANO<-variog(IEB.ANO.geo,max.dist=200,uvec=seq(0, 200, by=15),trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2))
plot(v.IEB.ANO,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.IEB.ANO$u,v.IEB.ANO$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_ANO.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

IEB.ANOlf3<-list()
IEB.ANOlf3$exp <- likfit(IEB.ANO.geo, ini=c(0.04,23.35), nug=0.02,trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2))
IEB.ANOlf3$gau <- likfit(IEB.ANO.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2))  
IEB.ANOlf3$sph <- likfit(IEB.ANO.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2)) 
IEB.ANOlf3$cir <- likfit(IEB.ANO.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2))
IEB.ANOlf3$kappa1.5 <- likfit(IEB.ANO.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2))  
IEB.ANOlf3$kappa2.5 <- likfit(IEB.ANO.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~trend_IEB.ANO[,18]+poly(trend_IEB.ANO[,6],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(IEB.ANOlf2, logLik),"AIC"=sapply(IEB.ANOlf2, AIC))
sort(sapply(IEB.ANOlf3, logLik),decreasing = T)
sort(sapply(IEB.ANOlf3, AIC))

summary(IEB.ANOlf3$exp)
IEB.ANOlf3$exp$parameters.summary
IEB.ANOlf3$exp$phi
IEB.ANOlf3$exp$sigmasq
IEB.ANOlf3$exp$tausq
  
plot(v.IEB.ANO,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(IEB.ANOlf3$exp,col=1)
lines(IEB.ANOlf3$gau,col=2)
lines(IEB.ANOlf3$sph,col=3)
lines(IEB.ANOlf3$cir,col=1)
lines(IEB.ANOlf3$kappa1.5,col=5)
lines(IEB.ANOlf3$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(IEB.ANO.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
IEB.ANOlf3$exp
kc.IEB.ANO <- krige.conv(IEB.ANO.geo, locations = grid, krige=krige.control(obj=IEB.ANOlf3$exp))
attributes(kc.IEB.ANO)
##kc.IEB.ANO$predict<-(backtransform.moments(lambda=lambda.IEB.ANO,mean=kc.IEB.ANO$predict,variance=kc.IEB.ANO$krige.var)$mean)-5
(1-kc.IEB.ANO$predict)-kc.IEB.ANO$IESpredict == 0
kc.IEB.ANO$IESpredict <- 1-kc.IEB.ANO$predict

IEBANO_krige <- cbind(grid*1000,kc.IEB.ANO$predict)
write.table(IEBANO_krige,"IEBANO.ascii",col.names = F,row.names=F,quote=F)

IESANO_krige <- cbind(grid*1000,kc.IEB.ANO$IESpredict)
write.table(IESANO_krige,"IESANO.ascii",col.names = F,row.names=F,quote=F)

summary(kc.IEB.ANO$predict)

postscript("Figuras/IEB_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.ANO,val=kc.IEB.ANO$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.ANO, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/IES_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.ANO,val=kc.IEB.ANO$IESpredict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.ANO,val=kc.IEB.ANO$IESpredict, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.IEB.ANO<-xvalid(IEB.ANO.geo,model=IEB.ANOlf3$exp)
names(xv.IEB.ANO)

plot(xv.IEB.ANO$predicted,xv.IEB.ANO$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.IEB.ANO$data~xv.IEB.ANO$predicted)$coef[1],lm(xv.IEB.ANO$data~xv.IEB.ANO$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.IEB.ANO)

smry(IEB.ANO.geo$data)


postscript("Figuras/xv_IEB-ANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.IEB.ANO$predicted,xv.IEB.ANO$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.IEB.ANO$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Verão
dir()
dados <- read.csv("ParamDJF.csv", head=T,dec=",",sep= ";")[,-c(1,4,5,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
names(dados)

##-----------------------------------------------------------------------------##
## Parametros IEB verão

## Análise exploratória geoestatística

dim(dados)
IEB.DJF.geo<-as.geodata(dados,coords.col = 1:2, data.col = 5,borders=T)
IEB.DJF.geo$borders<-limite

plot(IEB.DJF.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

boxcox(IEB.DJF.geo$data~1) 

## 
##trans<-boxcox(DJF.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.DJF<-with(trans, x[which.max(y)]);lambda.mu.DJF # Valor máximo de Lambda
## 
##Mu.DJF.geo$data<-(((DJF.mu^(lambda.mu.DJF)) - 1)/lambda.mu.DJF);Mu.DJF.geo$data #normalizando - (X^lambda.mu.DJF)-1/lambda.mu.DJF
##DJF.mu.inv<-(((lambda.mu.DJF*Mu.DJF.geo$data)+1)^(1/lambda.mu.DJF))-5;DJF.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_IEBDJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(IEB.DJF.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(IEB.DJF.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.IEB.DJF<-variog(IEB.DJF.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.IEB.DJF,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
names(dados)
(trend_IEB.DJF <- dados[,c(3,4,6:23)])
names(trend_IEB.DJF)
ncol(trend_IEB.DJF)

## Tendência linear
IEB.DJF <- list()
for(i in 1:ncol(trend_IEB.DJF)){
IEB.DJF[[i]] <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend_IEB.DJF[,i],ini=c(0.04,23.35), nug=0.02,)
IEB.DJF[[ncol(trend_IEB.DJF)+1]] <- likfit(IEB.DJF.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}

sapply(IEB.DJF, AIC)

names(trend_IEB.DJF)
summary(IEB.DJF[[order(sapply(IEB.DJF, AIC))[1]]])
trend1_IEB.DJF<-trend_IEB.DJF[,order(sapply(IEB.DJF, AIC))[1]];trend_IEB.DJF <- trend_IEB.DJF[,-order(sapply(IEB.DJF, AIC))[1]]

IEB.DJF1 <- list()
for(i in 1:ncol(trend_IEB.DJF)){
IEB.DJF1[[i]] <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend1_IEB.DJF+trend_IEB.DJF[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.DJF)
sapply(IEB.DJF1, AIC)
summary(IEB.DJF1[[order(sapply(IEB.DJF1, AIC))[1]]])

## Tendência quadrática
(trend_IEB.DJF <- dados[,c(3,4,6:23)])
names(trend_IEB.DJF)

IEB.DJF2 <- list()
for(i in 1:ncol(trend_IEB.DJF)){
IEB.DJF2[[i]] <- likfit(IEB.DJF.geo, cov.model="gau", trend =~poly(trend_IEB.DJF[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(IEB.DJF2, AIC)
summary(IEB.DJF2[[order(sapply(IEB.DJF2, AIC))[1]]])
trend3_IEB.DJF<-trend_IEB.DJF[,order(sapply(IEB.DJF2, AIC))[1]];trend_IEB.DJF <- trend_IEB.DJF[,-order(sapply(IEB.DJF2, AIC))[1]]

IEB.DJF3 <- list()
for(i in 1:ncol(trend_IEB.DJF)){
IEB.DJF3[[i]] <- likfit(IEB.DJF.geo, cov.model="gau", trend =~poly(trend3_IEB.DJF,2)+poly(trend_IEB.DJF[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.DJF)
sapply(IEB.DJF3, AIC)
summary(IEB.DJF3[[order(sapply(IEB.DJF3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_IEB.DJF <- dados[,c(3,4,6:23)]
names(trend_IEB.DJF)

IEB.DJF <- list()
IEB.DJF$lf1 <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend_IEB.DJF[,6]+trend_IEB.DJF[,18],ini=c(0.04,23.35), nug=0.02,)
IEB.DJF$lf2 <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend_IEB.DJF[,6]+poly(trend_IEB.DJF[,18],2),ini=c(0.04,23.35), nug=0.02,)
IEB.DJF$lf3 <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.DJF$lf4 <- likfit(IEB.DJF.geo, cov.model="gau", trend =~poly(trend_IEB.DJF[,18],2)+poly(trend_IEB.DJF[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.DJF$lf5 <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend_IEB.DJF[,11]+poly(trend_IEB.DJF[,18],2),ini=c(0.04,23.35), nug=0.02,)
IEB.DJF$lf6 <- likfit(IEB.DJF.geo, cov.model="gau", trend =~trend_IEB.DJF[,6]+poly(trend_IEB.DJF[,18],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(IEB.DJF, AIC))
plot(IEB.DJF.geo,low=T,trend=~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2))

postscript("Figuras/Vario-IEB_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.IEB.DJF<-variog(IEB.DJF.geo,max.dist=200,uvec=seq(0, 200, by=15),trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2))
plot(v.IEB.DJF,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.IEB.DJF$u,v.IEB.DJF$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_DJF.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

IEB.DJFlf3<-list()
IEB.DJFlf3$exp <- likfit(IEB.DJF.geo, ini=c(0.04,23.35), nug=0.02,trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2))
IEB.DJFlf3$gau <- likfit(IEB.DJF.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2))  
IEB.DJFlf3$sph <- likfit(IEB.DJF.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2)) 
IEB.DJFlf3$cir <- likfit(IEB.DJF.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2))
IEB.DJFlf3$kappa1.5 <- likfit(IEB.DJF.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2))  
IEB.DJFlf3$kappa2.5 <- likfit(IEB.DJF.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~trend_IEB.DJF[,18]+poly(trend_IEB.DJF[,6],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(IEB.DJFlf2, logLik),"AIC"=sapply(IEB.DJFlf2, AIC))
sort(sapply(IEB.DJFlf3, logLik),decreasing = T)
sort(sapply(IEB.DJFlf3, AIC))

summary(IEB.DJFlf3$exp)
IEB.DJFlf3$exp$parameters.summary
IEB.DJFlf3$exp$phi
IEB.DJFlf3$exp$sigmasq
IEB.DJFlf3$exp$tausq
  
plot(v.IEB.DJF,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(IEB.DJFlf3$exp,col=1)
lines(IEB.DJFlf3$gau,col=2)
lines(IEB.DJFlf3$sph,col=3)
lines(IEB.DJFlf3$cir,col=1)
lines(IEB.DJFlf3$kappa1.5,col=5)
lines(IEB.DJFlf3$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(IEB.DJF.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
IEB.DJFlf3$exp
kc.IEB.DJF <- krige.conv(IEB.DJF.geo, locations = grid, krige=krige.control(obj=IEB.DJFlf3$exp))
attributes(kc.IEB.DJF)
##kc.IEB.DJF$predict<-(backtransform.moments(lambda=lambda.IEB.DJF,mean=kc.IEB.DJF$predict,variance=kc.IEB.DJF$krige.var)$mean)-5
kc.IEB.DJF$predict

IEBDJF_krige <- cbind(grid*1000,kc.IEB.DJF$predict)
write.table(IEBDJF_krige,"IEBDJF.ascii",col.names = F,row.names=F,quote=F)


summary(kc.IEB.DJF$predict)

postscript("Figuras/IEB_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.DJF,val=kc.IEB.DJF$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.DJF, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/IES_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.DJF,val=kc.IEB.DJF$IESpredict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.DJF,val=kc.IEB.DJF$IESpredict, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.IEB.DJF<-xvalid(IEB.DJF.geo,model=IEB.DJFlf3$exp)
names(xv.IEB.DJF)

plot(xv.IEB.DJF$predicted,xv.IEB.DJF$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.IEB.DJF$data~xv.IEB.DJF$predicted)$coef[1],lm(xv.IEB.DJF$data~xv.IEB.DJF$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.IEB.DJF)

smry(IEB.DJF.geo$data)


postscript("Figuras/xv_IEB-DJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.IEB.DJF$predicted,xv.IEB.DJF$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.IEB.DJF$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Outono
dir()
dados <- read.csv("ParamMAM.csv", head=T,dec=",",sep= ";")[,-c(1,4,5,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
names(dados)

##-----------------------------------------------------------------------------##
## Parametros IEB Outono

## Análise exploratória geoestatística

dim(dados)
IEB.MAM.geo<-as.geodata(dados,coords.col = 1:2, data.col = 5,borders=T)
IEB.MAM.geo$borders<-limite

plot(IEB.MAM.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

boxcox(IEB.MAM.geo$data~1) 

## 
##trans<-boxcox(MAM.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.MAM<-with(trans, x[which.max(y)]);lambda.mu.MAM # Valor máximo de Lambda
## 
##Mu.MAM.geo$data<-(((MAM.mu^(lambda.mu.MAM)) - 1)/lambda.mu.MAM);Mu.MAM.geo$data #normalizando - (X^lambda.mu.MAM)-1/lambda.mu.MAM
##MAM.mu.inv<-(((lambda.mu.MAM*Mu.MAM.geo$data)+1)^(1/lambda.mu.MAM))-5;MAM.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_IEBMAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(IEB.MAM.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(IEB.MAM.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.IEB.MAM<-variog(IEB.MAM.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.IEB.MAM,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
names(dados)
(trend_IEB.MAM <- dados[,c(3,4,6:23)])
names(trend_IEB.MAM)
ncol(trend_IEB.MAM)

## Tendência linear
IEB.MAM <- list()
for(i in 1:ncol(trend_IEB.MAM)){
IEB.MAM[[i]] <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend_IEB.MAM[,i],ini=c(0.04,23.35), nug=0.02,)
IEB.MAM[[ncol(trend_IEB.MAM)+1]] <- likfit(IEB.MAM.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}

sapply(IEB.MAM, AIC)

names(trend_IEB.MAM)
summary(IEB.MAM[[order(sapply(IEB.MAM, AIC))[1]]])
trend1_IEB.MAM<-trend_IEB.MAM[,order(sapply(IEB.MAM, AIC))[1]];trend_IEB.MAM <- trend_IEB.MAM[,-order(sapply(IEB.MAM, AIC))[1]]

IEB.MAM1 <- list()
for(i in 1:ncol(trend_IEB.MAM)){
IEB.MAM1[[i]] <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend1_IEB.MAM+trend_IEB.MAM[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.MAM)
sapply(IEB.MAM1, AIC)
summary(IEB.MAM1[[order(sapply(IEB.MAM1, AIC))[1]]])

## Tendência quadrática
(trend_IEB.MAM <- dados[,c(3,4,6:23)])
names(trend_IEB.MAM)

IEB.MAM2 <- list()
for(i in 1:ncol(trend_IEB.MAM)){
IEB.MAM2[[i]] <- likfit(IEB.MAM.geo, cov.model="gau", trend =~poly(trend_IEB.MAM[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(IEB.MAM2, AIC)
summary(IEB.MAM2[[order(sapply(IEB.MAM2, AIC))[1]]])
trend3_IEB.MAM<-trend_IEB.MAM[,order(sapply(IEB.MAM2, AIC))[1]];trend_IEB.MAM <- trend_IEB.MAM[,-order(sapply(IEB.MAM2, AIC))[1]]

IEB.MAM3 <- list()
for(i in 1:ncol(trend_IEB.MAM)){
IEB.MAM3[[i]] <- likfit(IEB.MAM.geo, cov.model="gau", trend =~poly(trend3_IEB.MAM,2)+poly(trend_IEB.MAM[,i],2),ini=c(0.015,10.35), nug=0.001,)

}

names(trend_IEB.MAM)
sapply(IEB.MAM3, AIC)
summary(IEB.MAM3[[order(sapply(IEB.MAM3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_IEB.MAM <- dados[,c(3,4,6:23)]
names(trend_IEB.MAM)

IEB.MAM <- list()
IEB.MAM$lf1 <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend_IEB.MAM[,11]+trend_IEB.MAM[,18],ini=c(0.04,23.35), nug=0.02,)
IEB.MAM$lf2 <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend_IEB.MAM[,11]+poly(trend_IEB.MAM[,16],2),ini=c(0.04,23.35), nug=0.02,)
IEB.MAM$lf3 <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend_IEB.MAM[,18]+poly(trend_IEB.MAM[,16],2),ini=c(0.04,23.35), nug=0.02,)
IEB.MAM$lf4 <- likfit(IEB.MAM.geo, cov.model="gau", trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.MAM$lf5 <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend_IEB.MAM[,11]+poly(trend_IEB.MAM[,18],2),ini=c(0.04,23.35), nug=0.02,)
IEB.MAM$lf6 <- likfit(IEB.MAM.geo, cov.model="gau", trend =~trend_IEB.MAM[,16]+poly(trend_IEB.MAM[,18],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(IEB.MAM, AIC))
plot(IEB.MAM.geo,low=T,trend=~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2))

postscript("Figuras/Vario-IEB_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.IEB.MAM<-variog(IEB.MAM.geo,max.dist=200,uvec=seq(0, 200, by=15),trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2))
plot(v.IEB.MAM,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.IEB.MAM$u,v.IEB.MAM$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_MAM.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

IEB.MAMlf4<-list()
IEB.MAMlf4$exp <- likfit(IEB.MAM.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2))
IEB.MAMlf4$gau <- likfit(IEB.MAM.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2))  
IEB.MAMlf4$sph <- likfit(IEB.MAM.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2)) 
IEB.MAMlf4$cir <- likfit(IEB.MAM.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2))
IEB.MAMlf4$kappa1.5 <- likfit(IEB.MAM.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2))  
IEB.MAMlf4$kappa2.5 <- likfit(IEB.MAM.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_IEB.MAM[,18],2)+poly(trend_IEB.MAM[,6],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(IEB.MAMlf2, logLik),"AIC"=sapply(IEB.MAMlf2, AIC))
sort(sapply(IEB.MAMlf4, logLik),decreasing = T)
sort(sapply(IEB.MAMlf4, AIC))

summary(IEB.MAMlf4$exp)
IEB.MAMlf4$exp$parameters.summary
IEB.MAMlf4$exp$phi
IEB.MAMlf4$exp$sigmasq
IEB.MAMlf4$exp$tausq
  
plot(v.IEB.MAM,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(IEB.MAMlf4$exp,col=1)
lines(IEB.MAMlf4$gau,col=2)
lines(IEB.MAMlf4$sph,col=3)
lines(IEB.MAMlf4$cir,col=1)
lines(IEB.MAMlf4$kappa1.5,col=5)
lines(IEB.MAMlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(IEB.MAM.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
IEB.MAMlf4$exp
kc.IEB.MAM <- krige.conv(IEB.MAM.geo, locations = grid, krige=krige.control(obj=IEB.MAMlf4$exp))
attributes(kc.IEB.MAM)
##kc.IEB.MAM$predict<-(backtransform.moments(lambda=lambda.IEB.MAM,mean=kc.IEB.MAM$predict,variance=kc.IEB.MAM$krige.var)$mean)-5
kc.IEB.MAM$predict

IEBMAM_krige <- cbind(grid*1000,kc.IEB.MAM$predict)
write.table(IEBMAM_krige,"IEBMAM.ascii",col.names = F,row.names=F,quote=F)



summary(kc.IEB.MAM$predict)

postscript("Figuras/IEB_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.MAM,val=kc.IEB.MAM$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.MAM, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/IES_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.MAM,val=kc.IEB.MAM$IESpredict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.MAM,val=kc.IEB.MAM$IESpredict, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.IEB.MAM<-xvalid(IEB.MAM.geo,model=IEB.MAMlf4$exp)
names(xv.IEB.MAM)

plot(xv.IEB.MAM$predicted,xv.IEB.MAM$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.IEB.MAM$data~xv.IEB.MAM$predicted)$coef[1],lm(xv.IEB.MAM$data~xv.IEB.MAM$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.IEB.MAM)

smry(IEB.MAM.geo$data)


postscript("Figuras/xv_IEB-MAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.IEB.MAM$predicted,xv.IEB.MAM$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.IEB.MAM$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Inverno
dir()
dados <- read.csv("ParamJJA.csv", head=T,dec=",",sep= ";")[,-c(1,4,5,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
names(dados)

##-----------------------------------------------------------------------------##
## Parametros IEB Inverno

## Análise exploratória geoestatística

dim(dados)
IEB.JJA.geo<-as.geodata(dados,coords.col = 1:2, data.col = 5,borders=T)
IEB.JJA.geo$borders<-limite

plot(IEB.JJA.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

boxcox(IEB.JJA.geo$data~1) 

## 
##trans<-boxcox(JJA.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.JJA<-with(trans, x[which.max(y)]);lambda.mu.JJA # Valor máximo de Lambda
## 
##Mu.JJA.geo$data<-(((JJA.mu^(lambda.mu.JJA)) - 1)/lambda.mu.JJA);Mu.JJA.geo$data #normalizando - (X^lambda.mu.JJA)-1/lambda.mu.JJA
##JJA.mu.inv<-(((lambda.mu.JJA*Mu.JJA.geo$data)+1)^(1/lambda.mu.JJA))-5;JJA.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_IEBJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(IEB.JJA.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(IEB.JJA.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.IEB.JJA<-variog(IEB.JJA.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.IEB.JJA,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
names(dados)
(trend_IEB.JJA <- dados[,c(3,4,6:23)])
names(trend_IEB.JJA)
ncol(trend_IEB.JJA)

## Tendência linear
IEB.JJA <- list()
for(i in 1:ncol(trend_IEB.JJA)){
IEB.JJA[[i]] <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend_IEB.JJA[,i],ini=c(0.014,23.35), nug=0.002,)
IEB.JJA[[ncol(trend_IEB.JJA)+1]] <- likfit(IEB.JJA.geo, cov.model="gau",ini=c(0.014,30.35), nug=0.002,)
}

sapply(IEB.JJA, AIC)

names(trend_IEB.JJA)
summary(IEB.JJA[[order(sapply(IEB.JJA, AIC))[1]]])
trend1_IEB.JJA<-trend_IEB.JJA[,order(sapply(IEB.JJA, AIC))[1]];trend_IEB.JJA <- trend_IEB.JJA[,-order(sapply(IEB.JJA, AIC))[1]]

IEB.JJA1 <- list()
for(i in 1:ncol(trend_IEB.JJA)){
IEB.JJA1[[i]] <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend1_IEB.JJA+trend_IEB.JJA[,i],ini=c(0.014,23.35), nug=0.002,)

}

names(trend_IEB.JJA)
sapply(IEB.JJA1, AIC)
summary(IEB.JJA1[[order(sapply(IEB.JJA1, AIC))[1]]])

## Tendência quadrática
(trend_IEB.JJA <- dados[,c(3,4,6:23)])
names(trend_IEB.JJA)

IEB.JJA2 <- list()
for(i in 1:ncol(trend_IEB.JJA)){
IEB.JJA2[[i]] <- likfit(IEB.JJA.geo, cov.model="gau", trend =~poly(trend_IEB.JJA[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(IEB.JJA2, AIC)
summary(IEB.JJA2[[order(sapply(IEB.JJA2, AIC))[1]]])
trend3_IEB.JJA<-trend_IEB.JJA[,order(sapply(IEB.JJA2, AIC))[1]];trend_IEB.JJA <- trend_IEB.JJA[,-order(sapply(IEB.JJA2, AIC))[1]]

IEB.JJA3 <- list()
for(i in 1:ncol(trend_IEB.JJA)){
IEB.JJA3[[i]] <- likfit(IEB.JJA.geo, cov.model="gau", trend =~poly(trend3_IEB.JJA,2)+poly(trend_IEB.JJA[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.JJA)
sapply(IEB.JJA3, AIC)
summary(IEB.JJA3[[order(sapply(IEB.JJA3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_IEB.JJA <- dados[,c(3,4,6:23)]
names(trend_IEB.JJA)

IEB.JJA <- list()
IEB.JJA$lf1 <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend_IEB.JJA[,11]+trend_IEB.JJA[,20],ini=c(0.04,23.35), nug=0.02,)
IEB.JJA$lf2 <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend_IEB.JJA[,11]+poly(trend_IEB.JJA[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.JJA$lf3 <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend_IEB.JJA[,20]+poly(trend_IEB.JJA[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.JJA$lf4 <- likfit(IEB.JJA.geo, cov.model="gau", trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.JJA$lf5 <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend_IEB.JJA[,20]+poly(trend_IEB.JJA[,18],2),ini=c(0.04,23.35), nug=0.02,)
IEB.JJA$lf6 <- likfit(IEB.JJA.geo, cov.model="gau", trend =~trend_IEB.JJA[,6]+poly(trend_IEB.JJA[,18],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(IEB.JJA, AIC))
plot(IEB.JJA.geo,low=T,trend=~trend_IEB.JJA[,18]+poly(trend_IEB.JJA[,6],2))

postscript("Figuras/Vario-IEB_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.IEB.JJA<-variog(IEB.JJA.geo,max.dist=200,uvec=seq(0, 200, by=15),trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2))
plot(v.IEB.JJA,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.IEB.JJA$u,v.IEB.JJA$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_JJA.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

IEB.JJAlf4<-list()
IEB.JJAlf4$exp <- likfit(IEB.JJA.geo, ini=c(0.04,24.45), nug=0.02,trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2))
IEB.JJAlf4$gau <- likfit(IEB.JJA.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2))  
IEB.JJAlf4$sph <- likfit(IEB.JJA.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2)) 
IEB.JJAlf4$cir <- likfit(IEB.JJA.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2))
IEB.JJAlf4$kappa1.5 <- likfit(IEB.JJA.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2))  
IEB.JJAlf4$kappa2.5 <- likfit(IEB.JJA.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_IEB.JJA[,18],2)+poly(trend_IEB.JJA[,6],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(IEB.JJAlf2, logLik),"AIC"=sapply(IEB.JJAlf2, AIC))
sort(sapply(IEB.JJAlf4, logLik),decreasing = T)
sort(sapply(IEB.JJAlf4, AIC))

summary(IEB.JJAlf4$exp)
IEB.JJAlf4$exp$parameters.summary
IEB.JJAlf4$exp$phi
IEB.JJAlf4$exp$sigmasq
IEB.JJAlf4$exp$tausq
  
plot(v.IEB.JJA,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(IEB.JJAlf4$exp,col=1)
lines(IEB.JJAlf4$gau,col=2)
lines(IEB.JJAlf4$sph,col=4)
lines(IEB.JJAlf4$cir,col=1)
lines(IEB.JJAlf4$kappa1.5,col=5)
lines(IEB.JJAlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(IEB.JJA.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
IEB.JJAlf4$exp
kc.IEB.JJA <- krige.conv(IEB.JJA.geo, locations = grid, krige=krige.control(obj=IEB.JJAlf4$exp))
attributes(kc.IEB.JJA)
##kc.IEB.JJA$predict<-(backtransform.moments(lambda=lambda.IEB.JJA,mean=kc.IEB.JJA$predict,variance=kc.IEB.JJA$krige.var)$mean)-5
kc.IEB.JJA$predict

IEBJJA_krige <- cbind(grid*1000,kc.IEB.JJA$predict)
write.table(IEBJJA_krige,"IEBJJA.ascii",col.names = F,row.names=F,quote=F)

summary(kc.IEB.JJA$predict)

postscript("Figuras/IEB_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.JJA,val=kc.IEB.JJA$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.JJA, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/IES_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.JJA,val=kc.IEB.JJA$IESpredict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.JJA,val=kc.IEB.JJA$IESpredict, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.IEB.JJA<-xvalid(IEB.JJA.geo,model=IEB.JJAlf4$exp)
names(xv.IEB.JJA)

plot(xv.IEB.JJA$predicted,xv.IEB.JJA$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.IEB.JJA$data~xv.IEB.JJA$predicted)$coef[1],lm(xv.IEB.JJA$data~xv.IEB.JJA$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.IEB.JJA)

smry(IEB.JJA.geo$data)


postscript("Figuras/xv_IEB-JJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.IEB.JJA$predicted,xv.IEB.JJA$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.IEB.JJA$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Primavera
dir()
dados <- read.csv("ParamSON.csv", head=T,dec=",",sep= ";")[,-c(1,4,5,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
names(dados)

##-----------------------------------------------------------------------------##
## Parametros IEB Primavera

## Análise exploratória geoestatística

dim(dados)
IEB.SON.geo<-as.geodata(dados,coords.col = 1:2, data.col = 5,borders=T)
IEB.SON.geo$borders<-limite

plot(IEB.SON.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

boxcox(IEB.SON.geo$data~1) 

## 
##trans<-boxcox(SON.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.SON<-with(trans, x[which.max(y)]);lambda.mu.SON # Valor máximo de Lambda
## 
##Mu.SON.geo$data<-(((SON.mu^(lambda.mu.SON)) - 1)/lambda.mu.SON);Mu.SON.geo$data #normalizando - (X^lambda.mu.SON)-1/lambda.mu.SON
##SON.mu.inv<-(((lambda.mu.SON*Mu.SON.geo$data)+1)^(1/lambda.mu.SON))-5;SON.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_IEBSON.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(IEB.SON.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(IEB.SON.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.IEB.SON<-variog(IEB.SON.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.IEB.SON,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
names(dados)
(trend_IEB.SON <- dados[,c(3,4,6:23)])
names(trend_IEB.SON)
ncol(trend_IEB.SON)

## Tendência linear
IEB.SON <- list()
for(i in 1:ncol(trend_IEB.SON)){
IEB.SON[[i]] <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend_IEB.SON[,i],ini=c(0.04,23.35), nug=0.02,)
IEB.SON[[ncol(trend_IEB.SON)+1]] <- likfit(IEB.SON.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}

sapply(IEB.SON, AIC)

names(trend_IEB.SON)
summary(IEB.SON[[order(sapply(IEB.SON, AIC))[1]]])
trend1_IEB.SON<-trend_IEB.SON[,order(sapply(IEB.SON, AIC))[1]];trend_IEB.SON <- trend_IEB.SON[,-order(sapply(IEB.SON, AIC))[1]]

IEB.SON1 <- list()
for(i in 1:ncol(trend_IEB.SON)){
IEB.SON1[[i]] <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend1_IEB.SON+trend_IEB.SON[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.SON)
sapply(IEB.SON1, AIC)
summary(IEB.SON1[[order(sapply(IEB.SON1, AIC))[1]]])

## Tendência quadrática
(trend_IEB.SON <- dados[,c(3,4,6:23)])
names(trend_IEB.SON)

IEB.SON2 <- list()
for(i in 1:ncol(trend_IEB.SON)){
IEB.SON2[[i]] <- likfit(IEB.SON.geo, cov.model="gau", trend =~poly(trend_IEB.SON[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(IEB.SON2, AIC)
names(trend_IEB.SON)
summary(IEB.SON2[[order(sapply(IEB.SON2, AIC))[1]]])
trend3_IEB.SON<-trend_IEB.SON[,order(sapply(IEB.SON2, AIC))[1]];trend_IEB.SON <- trend_IEB.SON[,-order(sapply(IEB.SON2, AIC))[1]]

IEB.SON3 <- list()
for(i in 1:ncol(trend_IEB.SON)){
IEB.SON3[[i]] <- likfit(IEB.SON.geo, cov.model="gau", trend =~poly(trend3_IEB.SON,2)+poly(trend_IEB.SON[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_IEB.SON)
sapply(IEB.SON3, AIC)
summary(IEB.SON3[[order(sapply(IEB.SON3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_IEB.SON <- dados[,c(3,4,6:23)]
names(trend_IEB.SON)

IEB.SON <- list()
IEB.SON$lf1 <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend_IEB.SON[,11]+trend_IEB.SON[,18],ini=c(0.04,23.35), nug=0.02,)
IEB.SON$lf2 <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend_IEB.SON[,11]+poly(trend_IEB.SON[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.SON$lf3 <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend_IEB.SON[,18]+poly(trend_IEB.SON[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.SON$lf4 <- likfit(IEB.SON.geo, cov.model="gau", trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2),ini=c(0.04,23.35), nug=0.02,)
IEB.SON$lf5 <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend_IEB.SON[,11]+poly(trend_IEB.SON[,18],2),ini=c(0.04,23.35), nug=0.02,)
IEB.SON$lf6 <- likfit(IEB.SON.geo, cov.model="gau", trend =~trend_IEB.SON[,6]+poly(trend_IEB.SON[,18],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(IEB.SON, AIC))
plot(IEB.SON.geo,low=T,trend=~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2))

postscript("Figuras/Vario-IEB_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.IEB.SON<-variog(IEB.SON.geo,max.dist=200,uvec=seq(0, 200, by=15),trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2))
plot(v.IEB.SON,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.IEB.SON$u,v.IEB.SON$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_SON.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

IEB.SONlf4<-list()
IEB.SONlf4$exp <- likfit(IEB.SON.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2))
IEB.SONlf4$gau <- likfit(IEB.SON.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2))  
IEB.SONlf4$sph <- likfit(IEB.SON.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2)) 
IEB.SONlf4$cir <- likfit(IEB.SON.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2))
IEB.SONlf4$kappa1.5 <- likfit(IEB.SON.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2))  
IEB.SONlf4$kappa2.5 <- likfit(IEB.SON.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_IEB.SON[,18],2)+poly(trend_IEB.SON[,6],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(IEB.SONlf2, logLik),"AIC"=sapply(IEB.SONlf2, AIC))
sort(sapply(IEB.SONlf4, logLik),decreasing = T)
sort(sapply(IEB.SONlf4, AIC))

summary(IEB.SONlf4$exp)
IEB.SONlf4$exp$parameters.summary
IEB.SONlf4$exp$phi
IEB.SONlf4$exp$sigmasq
IEB.SONlf4$exp$tausq
  
plot(v.IEB.SON,xlab="Disiance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(IEB.SONlf4$exp,col=1)
lines(IEB.SONlf4$gau,col=2)
lines(IEB.SONlf4$sph,col=4)
lines(IEB.SONlf4$cir,col=1)
lines(IEB.SONlf4$kappa1.5,col=5)
lines(IEB.SONlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(IEB.SON.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
IEB.SONlf4$exp
kc.IEB.SON <- krige.conv(IEB.SON.geo, locations = grid, krige=krige.control(obj=IEB.SONlf4$exp))
attributes(kc.IEB.SON)
##kc.IEB.SON$predict<-(backtransform.moments(lambda=lambda.IEB.SON,mean=kc.IEB.SON$predict,variance=kc.IEB.SON$krige.var)$mean)-5
kc.IEB.SON$predict

IEBSON_krige <- cbind(grid*1000,kc.IEB.SON$predict)
write.table(IEBSON_krige,"IEBSON.ascii",col.names = F,row.names=F,quote=F)


summary(kc.IEB.SON$predict)

postscript("Figuras/IEB_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.SON,val=kc.IEB.SON$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.SON, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/IES_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.IEB.SON,val=kc.IEB.SON$IESpredict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.IEB.SON,val=kc.IEB.SON$IESpredict, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.IEB.SON<-xvalid(IEB.SON.geo,model=IEB.SONlf4$exp)
names(xv.IEB.SON)

plot(xv.IEB.SON$predicted,xv.IEB.SON$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.IEB.SON$data~xv.IEB.SON$predicted)$coef[1],lm(xv.IEB.SON$data~xv.IEB.SON$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.IEB.SON)

smry(IEB.SON.geo$data)


postscript("Figuras/xv_IEB-SON.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.IEB.SON$predicted,xv.IEB.SON$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.IEB.SON$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

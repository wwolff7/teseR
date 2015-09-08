
## Configuração Inicial

require(RColorBrewer) # Para Paletas de cores
require(hydroTSM) 
colour<-colorRampPalette(c("white","blue","dark Blue"))
require(MASS)
require(rgdal)
require(maptools)
require(geoR)
require(xtable)

getwd()
setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Paramgeo")
## rm(list = ls()) ## Para remover todos os objetos
cite("geoR")
citation("geoR")

##Ver documentação
cov.spatial()
trend.spatial()
options(OutDec=",",digits=10)
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

dados <- read.csv("ParamANO.csv", head=T,dec=",",sep= "")[,-c(1,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC.csv",h=T,sep=",",dec = ".");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
flori <- read.csv("Limite-Floripa.csv",h=T,sep=",",dec = ".")[,-3];flori$Longitude<-flori$Longitude/1000;flori$Latitude<-flori$Latitude/1000

summary(limite)
##write.table(dados,"Resultados.txt")

head(dados)
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

##-----------------------------------------------------------------------------##
## Curvas de permanência das Estações Fluviométricas
## Probabilidades testadas, lembrar que refere-se a 1-p !! 

postscript("Figuras/Estacoes.eps",onefile = T,horizontal = F, width=18/2.54, height=18/2.54,paper="special")
plot(limite$Longitude,limite$Latitude,type="l",asp=1,xlim = c(min(dados$Longitude),max(dados$Longitude)),ylim = c(min(dados$Latitude),max(dados$Latitude)),ylab = "Latitude (km)",xlab = "Longitude (km)")
lines(flori)
points(dados$Longitude,dados$Latitude,pch=19,col=1)


northarrow(loc=c(425,7025),size = 10,cex =0.8)
scalebar(loc=c(600,6780),length=100)

dev.off()

names(dados)


## Gráfico exemplo para volume reservatório
## Vazão média estação exemplo = E70100000
p <- seq(0,0.99,l=1000)
Qm <- rep(35.04155,times = length(p))

q <- qnorm(1-seq(0,0.9999,l=1000))
D(expression(exp(mu+(sigma*q))),"mu")
D(D(expression(exp(mu+(sigma*q))),"sigma"),"sigma")


Qp.int <- function(mu,sigma,area,p){
 Qp.int<-(exp(mu+(sigma*qnorm(1-p))))*area
     
}


summary(dados[,c(7,8)])
mu.param <- seq(-4.50,-3.10,l=20)
sigma.param <- seq(0.60,1.20,l=20)




for (i in 1:length(mu.param)){

x <- 5

plot(p*100,exp(mu.param[i]+(mean(sigma.param)*qnorm(1-p))),type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0)) 
mtext(bquote(.(x)~Variando~mu~para~sigma~fixo))##~~~Q[m]==.(round(Qm.est[i]),4)),cex=2)

Sys.sleep(0.9)
}


Qm.est <- c()
for (i in 1:length(sigma.param)){

Qm.est[i] <- integrate(Qp.int, lower=0, upper=1, mu=mean(mu.param), sigma=sigma.param[i],area=1)
plot(p*100,exp(mean(mu.param)+(sigma.param[i]*qnorm(1-p))),type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0)) 
mtext(bquote(Variando~sigma~para~mu~fixo),cex=2)
abline(h=Qm.est[i],col="red",lwd=2)
points(x=50,y=Qm.est[i],pch=19)


Sys.sleep(0.9)

}

## Ponto de inflexão
plot(p*100,exp(mean(mu.param)+(sigma.param[1]*qnorm(1-p))),type="n",xlab = bquote(p(theta)),ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0)) 

for (i in 1:length(sigma.param)){

lines(p*100,exp(mean(mu.param)+(sigma.param[i]*qnorm(1-p))),lty = i)
mtext(bquote(Variando~sigma~para~mu~fixo),cex=2)
Sys.sleep(0.9)

}


## qlnorm(1-(seq(0,1,by=0.01)))
## plot(qlnorm(1-(seq(0,1,by=0.01)),0,1),seq(0,1,by=0.01),ty="n",xlim = c(0,3))
## for(i in 1:length(sigma.param)){
## lines(qlnorm(1-seq(0,1,by=0.01),0,sigma.param[i]),seq(0,1,by=0.01),col=i)
## Sys.sleep(0.9)
## }
## abline(h=0.5)

## plot(qlnorm(1-(seq(0,1,by=0.01)),1,0.5),seq(0,1,by=0.01),ty="l",xlim = c(0,3))
## for(i in 1:length(mu.param)){
## lines(qlnorm(1-seq(0,1,by=0.01),mu.param[i],0.5),seq(0,1,by=0.01),col=i)
## Sys.sleep(0.9)
## }

postscript("Vr9.eps",onefile = T,horizontal = F, width=18/2.54, height=18/2.54)

plot(p*100,Qp1,type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
lines(p*100,Qm,col=2,lty=2,lwd = 2)

pqf <-1-pnorm((log(Qm[1]/dados$Área_Bacia[4])-dados$Param_ANO.mu[4])/dados$Param_ANO.sigma[4])

T <- 5
pqfseq <- seq(pqf*100,(1-(1/T))*100,l=1000)

Qpseq <- (exp(dados$Param_ANO.mu[4]+(dados$Param_ANO.sigma[4]*qnorm(1-(pqfseq/100)))))*dados$Área_Bacia[4]

coords.x <- c(pqf*100,pqfseq,(1-(1/T))*100)
coords.y <- c(Qm[1],Qpseq,Qm[1])

polygon(coords.x,coords.y,col=1)
text(0,Qm[1]*1.20,expression(Q[f]),cex=1.5)
text(mean(pqfseq),Qm[1]*1.20,bquote(V[r]~para~T==.(T)~anos),cex = 1.5)
text(98,mean(seq(0,Qm[1])),"T",cex = 1.5)
arrows(97,mean(seq(0,Qm[1])),max(pqfseq),mean(seq(0,Qm[1])),lwd = 2)

dev.off()

## Integral vazão média

Qp.int <- function(mu,sigma,area,p){
 Qp.int<-(exp(mu+(sigma*qnorm(1-p))))*area
     
}

## Feito no R
(Qm.est <- integrate(Qp.int, lower=0, upper=1, mu=dados$Param_ANO.mu[4], sigma=dados$Param_ANO.sigma[4],area=dados$Área_Bacia[4]))
Qm[1]
(Qp.int(mu=dados$Param_ANO.mu[4], sigma=dados$Param_ANO.sigma[4],area=dados$Área_Bacia[4],p=0.3))

## Com a função simpson
simpson <- function(fun, a, b, n=100) {
	
	h <- (b-a)/n
	x <- seq(a, b, by=h)
	if (n == 2) {
		s <- fun(x[1]) + 4*fun(x[2]) +fun(x[3])
	} else {
		s <- fun(x[1]) + fun(x[n+1]) + 2*sum(fun(x[seq(2,n,by=2)])) + 4 *sum(fun(x[seq(3,n-1, by=2)]))
	}
	s <- s*h/3
	return(s)
}

Qp.int2 <- function(p){
    Qp.int(dados$Param_ANO.mu[4],dados$Param_ANO.sigma[4],dados$Área_Bacia[4],p)}

simpson(fun=function(p) Qp.int(dados$Param_ANO.mu[4],dados$Param_ANO.sigma[4],dados$Área_Bacia[4],p),0.0000001,0.9999999,n=100000)
Qm[1]
Qm.est

## Integral Volume do Reservatório

Vr <- function(mu, sigma, area, T,Qf){

    pqfi <-1-pnorm((log(Qf/area)-mu)/sigma)
    pqff <- 1-(1/T) 

Vr <- ((((pqff-pqfi)*Qf)-integrate(Qp.int, lower=pqfi, upper=pqff, mu=mu, sigma=sigma,area=area)$value)*(60*60*24*365))/10^6

return(Vr)    
}

T <- (1:100)
Vr.T <- Vr(mu=dados$Param_ANO.mu[4],sigma = dados$Param_ANO.sigma[4],area=22,T=T,Qf = 0.28*0.8)
Qf <- seq(0.28*0.05,0.28*0.80,l=100)
Vr.Qf <- Vr(mu=dados$Param_ANO.mu[4],sigma = dados$Param_ANO.sigma[4],area=22,T=30,Qf = Qf)

par(mfrow = c(2,2))
plot(T,Vr.T,type = "l",ylab = expression(V[r]~~(m^{3}~10^6)),xlab="T (anos)",mgp=c(2.5,1,0))
plot(Vr.Qf,Qf,type = "l",ylab = expression(Q[f]~~(m^3~s^{-1})),xlab=expression(V[r]~~(m^3~10^6)),mgp=c(2.5,1,0))

## Reservatório tamanho físico m³
vol <- 381807*0.95

fun.spline <- splinefun(Vr.Qf, Qf)
Qf.possível <- fun.spline(vol/10^6)

plot(Vr.Qf,Qf,type = "l",ylab = expression(Q[f]~~(m^3~s^{-1})),xlab=expression(V[r]~~(m^3~10^6)),mgp=c(2.5,1,0))
lines(Vr.Qf,fun.spline(Vr.Qf),col=2,)
points(vol/10^6,Qf.possível,pch=19)

Qp_ANO <- Qp(dados[,c(7,8)],22)
Qf.possível/(0.28*0.8)


Qp <- function(Param,area){
p <- seq(0.01,0.99,l=1000)
    Qp<-matrix(NA,length(p),nrow(Param))
    for(i in 1:nrow(Param)){
        Qp[,i]<-(exp(Param[i,1]+(Param[i,2]*qnorm(1-p))))*area[i]
        
    }
colnames(Qp) <- dados$Código
return(Qp)
}
summary(dados[c(7,8)])

Qp_ANO <- Qp(dados[,c(7,8)],dados$Área_Bacia)
Qp_DJF <- Qp(dados[,c(9,10)],dados$Área_Bacia)
Qp_MAM <- Qp(dados[,c(11,12)],dados$Área_Bacia)
Qp_JJA <- Qp(dados[,c(13,14)],dados$Área_Bacia)
Qp_SON <- Qp(dados[,c(15,16)],dados$Área_Bacia)

dim(Qp_SON)
dim(Qp_MAM)

plot.perm <- function(Param,Qp){
for(i in 1:ncol(Qp)){
       p <- seq(0.01,0.99,l=1000)
       par(mar=c(4,5,3,1.5)) ##‘c(bottom, left, top, right)’
       par(mfrow=c(1,2))

       plot(p*100,Qp[,i],type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
       mtext(bquote(Q[p]==group("{",e^group("{",.(round(Param[i,1],4)) +
       group("[",.(round(Param[i,2],4))%.%phi^{-1}%.%(1-p),"]"),"}"),"}")%.%.(dados$Área_Bacia[i])),cex=1.5)

       plot(limite,type="l",asp=1,xlab = "Longitude (km)",ylab = "Latitude (km)",
            xlim = c(min(dados$Centroid.X),max(dados$Centroid.X)),ylim = c(min(dados$Centroid.Y),max(dados$Centroid.Y)),mgp=c(2.5,1,0))
       points(dados$Centroid.X[i],dados$Centroid.Y[i],pch=19,col="red")
       Sys.sleep(0.9)
   }
return(plot.perm)
}

setwd("/home/wagner/MEGA/Doutorado/Apresentação_Quali")
dir.create("frames0")
setwd("frames0")

postscript("curv_mapa%03d.eps",onefile = F,horizontal = T, width=30/2.54, height=20/2.54,)

plot.perm(dados[,c(7,8)],Qp_ANO)

dev.off()

plot.perm(dados[,c(9,10)],Qp_DJF)
plot.perm(dados[,c(11,12)],Qp_MAM)
plot.perm(dados[,c(13,14)],Qp_JJA)
plot.perm(dados[,c(15,16)],Qp_SON)


dir.create("frames1")
setwd("frames1")

postscript("curv_perm%03d.eps",onefile = F,horizontal = T, width=30/2.54, height=20/2.54,)


for(i in 1:ncol(Qp_DJF)){
       p <- seq(0.01,0.99,l=1000)
       par(mar=c(4,5,3,1.5)) ##‘c(bottom, left, top, right)’
       par(mfrow=c(2,2))
       plot(p*100,Qp_DJF[,i],type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
       mtext(bquote(Q[p]==group("{",e^group("{",.(round(dados$Param_DJF.mu[i],4)) +
       group("[",.(round(dados$Param_DJF.sigma[i],4))%.%phi^{-1}%.%(1-p),"]"),"}"),"}")%.%.(dados$Área_Bacia[i])),cex=1.5)   
       legend("topright","Verão",bty = "n")

       plot(p*100,Qp_MAM[,i],type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
       mtext(bquote(Q[p]==group("{",e^group("{",.(round(dados$Param_MAM.mu[i],4)) +
       group("[",.(round(dados$Param_MAM.sigma[i],4))%.%phi^{-1}%.%(1-p),"]"),"}"),"}")%.%.(dados$Área_Bacia[i])),cex=1.5)   
       legend("topright","Outono",bty = "n")

       plot(p*100,Qp_JJA[,i],type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
       mtext(bquote(Q[p]==group("{",e^group("{",.(round(dados$Param_JJA.mu[i],4)) +
       group("[",.(round(dados$Param_JJA.sigma[i],4))%.%phi^{-1}%.%(1-p),"]"),"}"),"}")%.%.(dados$Área_Bacia[i])),cex=1.5)   
       legend("topright","Inverno",bty = "n")

       plot(p*100,Qp_SON[,i],type="l",xlab = "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
       mtext(bquote(Q[p]==group("{",e^group("{",.(round(dados$Param_SON.mu[i],4)) +
       group("[",.(round(dados$Param_SON.sigma[i],4))%.%phi^{-1}%.%(1-p),"]"),"}"),"}")%.%.(dados$Área_Bacia[i])),cex=1.5)   
       legend("topright","Prmavera",bty = "n")
       Sys.sleep(0.8)
   }
dev.off()

setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Paramgeo")
## converte os pngs para um gif usando ImageMagick
## system("convert -delay 10 curv_perm*.png curv_perm.gif")
 
## remove os arquivos png
## file.remove(list.files(pattern=".png"))

## setwd("/home/wagner/Gdrive/Doutorado/Disciplinas/Geoestatística/Geoestatística-2014-PJ/Artigo")
##-----------------------------------------------------------------------------##


par(new=T)
plot(dados[,2:1],asp=1)
plot(dados[c(1, 3, 4, 10, 12, 34, 50),c(1,2)], asp=1)

##-----------------------------------------------------------------------------##
## resumo dos dados
xtable(summary(dados[,3:5]))
xtable(smry(dados[,3:5]))
smry(dados[,7:16])

##-----------------------------------------------------------------------------##

##-----------------------------------------------------------------------------##
## Análise exploratória dos dados 

apply(dados,2,shapiro.test)
apply(dados[,c(5,7:16)],2,shapiro.test)
names(dados[,c(5,7:16)])


##-----------------------------------------------------------------------------##
## Parametros ANO mu

## Análise exploratória geoestatística

names(dados)

dim(dados)
Mu.ANO.geo<-as.geodata(dados,coords.col = 1:2, data.col = 3,borders=T)
Mu.ANO.geo$borders<-limite



plot(Mu.ANO.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

ANO.mu <- Mu.ANO.geo$data+5

boxcox(ANO.mu~1) 
boxcox(ANO.mu~1,lambda=seq(0,1,l=20))
## 
##trans<-boxcox(ANO.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.ANO<-with(trans, x[which.max(y)]);lambda.mu.ANO # Valor máximo de Lambda
## 
##Mu.ANO.geo$data<-(((ANO.mu^(lambda.mu.ANO)) - 1)/lambda.mu.ANO);Mu.ANO.geo$data #normalizando - (X^lambda.mu.ANO)-1/lambda.mu.ANO
##ANO.mu.inv<-(((lambda.mu.ANO*Mu.ANO.geo$data)+1)^(1/lambda.mu.ANO))-5;ANO.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_muANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Mu.ANO.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(ANO.mu~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.mu.ANO<-variog(Mu.ANO.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.mu.ANO,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
(trend_Mu.ANO <- dados[,c(5:25,1,2)])
names(trend_Mu.ANO)
ncol(trend_Mu.ANO)

## Tendência linear
Mu.ANO <- list()
for(i in 1:ncol(trend_Mu.ANO)){
Mu.ANO[[i]] <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend_Mu.ANO[,i],ini=c(0.04,23.35), nug=0.02,)
Mu.ANO[[ncol(trend_Mu.ANO)+1]] <- likfit(Mu.ANO.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Mu.ANO, AIC)

names(trend_Mu.ANO)
summary(Mu.ANO[[order(sapply(Mu.ANO, AIC))[1]]])
trend1_Mu.ANO<-trend_Mu.ANO[,order(sapply(Mu.ANO, AIC))[1]];trend_Mu.ANO <- trend_Mu.ANO[,-order(sapply(Mu.ANO, AIC))[1]]

Mu.ANO1 <- list()
for(i in 1:ncol(trend_Mu.ANO)){
Mu.ANO1[[i]] <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend1_Mu.ANO+trend_Mu.ANO[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.ANO)
sapply(Mu.ANO1, AIC)
summary(Mu.ANO1[[order(sapply(Mu.ANO1, AIC))[1]]])

## Tendência quadrática
trend_Mu.ANO <- dados[,c(5:25,1,2)]
names(trend_Mu.ANO)

Mu.ANO2 <- list()
for(i in 1:ncol(trend_Mu.ANO)){
Mu.ANO2[[i]] <- likfit(Mu.ANO.geo, cov.model="gau", trend =~poly(trend_Mu.ANO[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Mu.ANO2, AIC)
summary(Mu.ANO2[[order(sapply(Mu.ANO2, AIC))[1]]])
trend3_Mu.ANO<-trend_Mu.ANO[,order(sapply(Mu.ANO2, AIC))[1]];trend_Mu.ANO <- trend_Mu.ANO[,-order(sapply(Mu.ANO2, AIC))[1]]

Mu.ANO3 <- list()
for(i in 1:ncol(trend_Mu.ANO)){
Mu.ANO3[[i]] <- likfit(Mu.ANO.geo, cov.model="gau", trend =~poly(trend3_Mu.ANO,2)+poly(trend_Mu.ANO[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.ANO)
sapply(Mu.ANO3, AIC)
summary(Mu.ANO3[[order(sapply(Mu.ANO3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Mu.ANO <- dados[,c(5:25,1,2)]
names(trend_Mu.ANO)

Mu.ANOfim <- list()
Mu.ANOfim$lf1 <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend_Mu.ANO[,3]+trend_Mu.ANO[,9],ini=c(0.04,23.35), nug=0.02,)
Mu.ANOfim$lf2 <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend_Mu.ANO[,9]+poly(trend_Mu.ANO[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.ANOfim$lf3 <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend_Mu.ANO[,3]+poly(trend_Mu.ANO[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.ANOfim$lf4 <- likfit(Mu.ANO.geo, cov.model="gau", trend =~poly(trend_Mu.ANO[,3],2)+poly(trend_Mu.ANO[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.ANOfim$lf5 <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend_Mu.ANO[,9]+poly(trend_Mu.ANO[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.ANOfim$lf6 <- likfit(Mu.ANO.geo, cov.model="gau", trend =~trend_Mu.ANO[,9]+poly(trend_Mu.ANO[,3],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Mu.ANOfim, AIC))
plot(Mu.ANO.geo,low=T,trend=~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2))

postscript("Figuras/Vario-mu_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.mu.ANO<-variog(Mu.ANO.geo,max.dist=180,uvec=seq(0, 180, by=10),trend =~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2))
plot(v.mu.ANO,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.mu.ANO$u,v.mu.ANO$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_ANO.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Mu.ANOlf4<-list()
Mu.ANOlf4$exp <- likfit(Mu.ANO.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2))
Mu.ANOlf4$gau <- likfit(Mu.ANO.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2))  
Mu.ANOlf4$sph <- likfit(Mu.ANO.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2)) 
Mu.ANOlf4$cir <- likfit(Mu.ANO.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2))
Mu.ANOlf4$kappa1.5 <- likfit(Mu.ANO.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2))  
Mu.ANOlf4$kappa2.5 <- likfit(Mu.ANO.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Mu.ANO[,19],2)+poly(trend_Mu.ANO[,3],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(Mu.ANOlf2, logLik),"AIC"=sapply(Mu.ANOlf2, AIC))
sort(sapply(Mu.ANOlf4, logLik),decreasing = T)
sort(sapply(Mu.ANOlf4, AIC))

summary(Mu.ANOlf4$gau)
Mu.ANOlf4$gau$parameters.summary
Mu.ANOlf4$gau$phi
Mu.ANOlf4$gau$sigmasq
Mu.ANOlf4$gau$tausq
  
plot(v.mu.ANO,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Mu.ANOlf4$exp,col=1)
lines(Mu.ANOlf4$gau,col=2)
lines(Mu.ANOlf4$sph,col=3)
lines(Mu.ANOlf4$cir,col=1)
lines(Mu.ANOlf4$kappa1.5,col=5)
lines(Mu.ANOlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Mu.ANO.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Mu.ANOlf4$gau
kc.mu.ANO <- krige.conv(Mu.ANO.geo, locations = grid, krige=krige.control(obj=Mu.ANOlf4$gau))
attributes(kc.mu.ANO)
##kc.mu.ANO$predict<-(backtransform.moments(lambda=lambda.mu.ANO,mean=kc.mu.ANO$predict,variance=kc.mu.ANO$krige.var)$mean)-5
kc.mu.ANO$predict

Mu_ANO_krige <- cbind(grid*1000,kc.mu.ANO$predict)
write.table(Mu_ANO_krige,"MuANO.ascii",col.names = F,row.names=F,quote=F)

summary(kc.mu.ANO$predict)

postscript("Figuras/mu_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.mu.ANO,val=kc.mu.ANO$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.mu.ANO, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.mu.ANO<-xvalid(Mu.ANO.geo,model=Mu.ANOfim$lf4)
names(xv.mu.ANO)

plot(xv.mu.ANO$predicted,xv.mu.ANO$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.mu.ANO$data~xv.mu.ANO$predicted)$coef[1],lm(xv.mu.ANO$data~xv.mu.ANO$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.mu.ANO)

smry(Mu.ANO.geo$data)


postscript("Figuras/xv_mu-ANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.mu.ANO$predicted,xv.mu.ANO$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.mu.ANO$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##








## Inferência Bayesiana Mu ANO
trend=~Prec_media+Param_ANO.sigma
MC <- model.control(trend.d = "1st",cov.model = "gau")
args(prior.control)

OC <- output.control(n.post=100)

PC <- prior.control(beta.prior = "normal", beta = c(0,0,0),
                    beta.var = cbind(c(2,1.5,0),c(1.5,1.8,.2),c(0,0.2,1.5)),
                    phi.discrete=seq(0,40,l=41), 
                    phi.prior="rec",
                    tausq.rel.discrete=seq(0,0.01,l=20),
                    tausq.rel.prior="rec")

prior.control(beta.prior = "normal", beta = c(0,0,0),
                           beta.var = cbind(c(2,1.5,0),c(1.5,1.8,.2),c(0,0.2,1.5)),
                           phi.prior = "exponential", phi = 2.5, phi.discrete = c(2.5,3),
                           sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = 0.5)



kb <- krige.bayes(Mu.ANO.geo, model=MC, prior=PC, output=OC)
par(mfrow=c(1,2))
names(kb)

plot(kb)
kb$predictive

grid <- pred_grid(range(dados$Centroid.X),range(dados$Centroid.Y), by=10)
plot(grid)

kb <- krige.bayes(Mu.ANO.geo, loc=grid, model=model.control(trend.d = "1st",trend.l = "1st",cov.model = "gau"), prior=PC, output=OC)
plot(limite)
names(kb)

names(kb$post)

(kb$pred)

image(kb)

p50bayes <- apply(kb$pred$simul, 1, function(x) mean(x>50))
image(kb, values = p50bayes)

## comparando as predicoes
image(kc.ml)
image(kb)

plot(kc.ml, kb$predictive$pred)


## Simulating data
ex.data <- grf(50, cov.pars=c(10, .25))
as.geodata(ex.data)
##
## Basic usage
##

## a basic and simple call to the function
ex.post <- krige.bayes(ex.data)
ex.post
names(ex.post)

## different input and output options
ex1 <- krige.bayes(ex.data, prior = list(phi.prior = "fixed", phi = 0.3))
ex1 <- krige.bayes(ex.data, model = list(cov.model="spherical"))
ex1 <- krige.bayes(ex.data, output = list(n.posterior = 100))

## now performing prediction
ex.grid <- as.matrix(expand.grid(seq(0,1,l=6), seq(0,1,l=6)))
points(ex.grid)
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, prior =
                        prior.control(phi.discrete=seq(0, 2, l=3),
                                      tausq.rel.discrete=seq(0, 2, l=3)),
                        output=output.control(n.post=100))
names(ex.bayes)

image(ex.bayes)

## some graphical output ...
## ... for the posterior ...
plot(ex.data)
lines(ex.bayes, sum = mean) 
lines(ex.bayes, summ="median", lty=2, post="par")
lines(ex.bayes, summ="mean", lwd=2, lty=2, post="par")

plot(ex.bayes)
## ... and for the predictive
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(mgp = c(2,1,0))
image(ex.bayes, main="predicted values")
image(ex.bayes, val="variance", main="prediction variance")
image(ex.bayes, val= "simulation", number.col=1,
      main="a simulation from the \npredictive distribution")
image(ex.bayes, val= "simulation", number.col=2,
      main="another simulation from \nthe predictive distribution")


##-----------------------------------------------------------------------------##

## Parametros ANO Sigma


## Análise exploratória geoestatística

names(dados)

dim(dados)
Sigma.ANO.geo<-as.geodata(dados,coords.col = 1:2, data.col = 4,borders=T)
Sigma.ANO.geo$borders<-limite

plot(Sigma.ANO.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX



boxcox(Sigma.ANO.geo$data~1) 
boxcox(Sigma.ANO.geo$data~1,lambda=seq(0,1,l=20))
## 
##trans<-boxcox(ANO.sigma~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.sigma.ANO<-with(trans, x[which.max(y)]);lambda.sigma.ANO # Valor máximo de Lambda
## 
##Sigma.ANO.geo$data<-(((ANO.sigma^(lambda.sigma.ANO)) - 1)/lambda.sigma.ANO);Sigma.ANO.geo$data #normalizando - (X^lambda.sigma.ANO)-1/lambda.sigma.ANO
##ANO.sigma.inv<-(((lambda.sigma.ANO*Sigma.ANO.geo$data)+1)^(1/lambda.sigma.ANO))-5;ANO.sigma.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_sigmaANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Sigma.ANO.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(Sigma.ANO.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.sigma.ANO<-variog(Sigma.ANO.geo,max.dist=200,uvec=seq(0, 200, by=10),trend=)
plot(v.sigma.ANO,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)

trend_Sigma.ANO <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.ANO)
ncol(trend_Sigma.ANO)

## Tendência linear
Sigma.ANO <- list()
for(i in 1:ncol(trend_Sigma.ANO)){
Sigma.ANO[[i]] <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend_Sigma.ANO[,i],ini=c(0.04,23.35), nug=0.02,)
Sigma.ANO[[ncol(trend_Sigma.ANO)+1]] <- likfit(Sigma.ANO.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Sigma.ANO, AIC)

names(trend_Sigma.ANO)
summary(Sigma.ANO[[order(sapply(Sigma.ANO, AIC))[1]]])
trend1_Sigma.ANO<-trend_Sigma.ANO[,order(sapply(Sigma.ANO, AIC))[1]];trend_Sigma.ANO <- trend_Sigma.ANO[,-order(sapply(Sigma.ANO, AIC))[1]]

Sigma.ANO1 <- list()
for(i in 1:ncol(trend_Sigma.ANO)){
Sigma.ANO1[[i]] <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend1_Sigma.ANO+trend_Sigma.ANO[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Sigma.ANO)
sapply(Sigma.ANO1, AIC)
summary(Sigma.ANO1[[order(sapply(Sigma.ANO1, AIC))[1]]])

## Tendência quadrática
trend_Sigma.ANO <- dados[,c(5,6,8:25,1,2)]


Sigma.ANO2 <- list()
for(i in 1:ncol(trend_Sigma.ANO)){
Sigma.ANO2[[i]] <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~poly(trend_Sigma.ANO[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Sigma.ANO2, AIC)
summary(Sigma.ANO2[[order(sapply(Sigma.ANO2, AIC))[1]]])
names(trend_Sigma.ANO)
trend3_Sigma.ANO<-trend_Sigma.ANO[,order(sapply(Sigma.ANO2, AIC))[1]];trend_Sigma.ANO <- trend_Sigma.ANO[,-order(sapply(Sigma.ANO2, AIC))[1]]

Sigma.ANO3 <- list()
for(i in 1:ncol(trend_Sigma.ANO)){
Sigma.ANO3[[i]] <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~poly(trend3_Sigma.ANO,2)+poly(trend_Sigma.ANO[,i],2),ini=c(0.04,23.35), nug=0.02,)
}

names(trend_Sigma.ANO)
sapply(Sigma.ANO3, AIC)
summary(Sigma.ANO3[[order(sapply(Sigma.ANO3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Sigma.ANO <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.ANO)

Sigma.ANO <- list()
Sigma.ANO$lf1 <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend_Sigma.ANO[,13]+trend_Sigma.ANO[,18],ini=c(0.04,23.35), nug=0.02,)
Sigma.ANO$lf2 <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend_Sigma.ANO[,13]+poly(trend_Sigma.ANO[,9],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.ANO$lf3 <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend_Sigma.ANO[,18]+poly(trend_Sigma.ANO[,9],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.ANO$lf4 <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~poly(trend_Sigma.ANO[,9],2)+poly(trend_Sigma.ANO[,20],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.ANO$lf5 <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend_Sigma.ANO[,13]+poly(trend_Sigma.ANO[,20],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.ANO$lf6 <- likfit(Sigma.ANO.geo, cov.model="gau", trend =~trend_Sigma.ANO[,18]+poly(trend_Sigma.ANO[,20],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Sigma.ANO, AIC))
plot(Sigma.ANO.geo,low=T,trend=~poly(trend_Sigma.ANO[,9],2)+poly(trend_Sigma.ANO[,20],2))

postscript("Figuras/Vario-sigma_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.sigma.ANO<-variog(Sigma.ANO.geo,max.dist=180,uvec=seq(0, 180, by=10),trend =~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2))
plot(v.sigma.ANO,xlab="Distância (km)",ylab=expression(gamma(u)))#,ylim = c(5e-3,0.01))
lines(lowess(v.sigma.ANO$u,v.sigma.ANO$v))
dev.off()


##-----------------------------------------------------------------------------##
## uMelhor modelo trend=~Prec_media+Param_ANO.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Sigma.ANOlf4<-list()
Sigma.ANOlf4$exp <- likfit(Sigma.ANO.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2))
Sigma.ANOlf4$gau <- likfit(Sigma.ANO.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2))  
Sigma.ANOlf4$sph <- likfit(Sigma.ANO.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2)) 
Sigma.ANOlf4$cir <- likfit(Sigma.ANO.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2))
Sigma.ANOlf4$kappa1.5 <- likfit(Sigma.ANO.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2))  
Sigma.ANOlf4$kappa2.5 <- likfit(Sigma.ANO.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Sigma.ANO[,20],2)+poly(trend_Sigma.ANO[,9],2)) 


##-----------------------------------------------------------------------------


xtable(data.frame("Máxima verossimilhança"=sapply(sigmalf1, logLik)))
data.frame("Modelos"=sapply(Sigma.ANOlf2, logLik),"AIC"=sapply(Sigma.ANOlf2, AIC))
sort(sapply(Sigma.ANOlf4, logLik),decreasing = T)
sort(sapply(Sigma.ANOlf4, AIC))

summary(Sigma.ANOlf4$gau)
Sigma.ANOlf4$gau$parameters.summary
Sigma.ANOlf4$gau$phi
Sigma.ANOlf4$gau$sigmasq
Sigma.ANOlf4$gau$tausq
  
plot(v.sigma.ANO,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Sigma.ANOlf4$exp,col=1)
lines(Sigma.ANOlf4$gau,col=2)
lines(Sigma.ANOlf4$sph,col=3)
lines(Sigma.ANOlf4$cir,col=1)
lines(Sigma.ANOlf4$kappa1.5,col=5)
lines(Sigma.ANOlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Sigma.ANO.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Sigma.ANOlf4$gau
kc.sigma.ANO <- krige.conv(Sigma.ANO.geo, locations = grid, krige=krige.control(obj=Sigma.ANOlf4$gau))
attributes(kc.sigma.ANO)
##kc.sigma.ANO$predict<-(backtransform.moments(lambda=lambda.sigma.ANO,mean=kc.sigma.ANO$predict,variance=kc.sigma.ANO$krige.var)$mean)-5
kc.sigma.ANO$predict

Sigma_ANO_krige <- cbind(grid*1000,kc.sigma.ANO$predict)
write.table(Sigma_ANO_krige,"SigmaANO.ascii",col.names = F,row.names=F,quote=F)


summary(kc.sigma.ANO$predict)

postscript("Figuras/sigma_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.ANO,val=kc.sigma.ANO$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.ANO, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.sigma.ANO<-xvalid(Sigma.ANO.geo,model=Sigma.ANOlf4$gau)
names(xv.sigma.ANO)

plot(xv.sigma.ANO$predicted,xv.sigma.ANO$data)
abline(0,1)
abline(lm(xv.sigma.ANO$data~xv.sigma.ANO$predicted)$coef[1],lm(xv.sigma.ANO$data~xv.sigma.ANO$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.sigma.ANO)

smry(Sigma.ANO.geo$data)


postscript("Figuras/xv_sigma-ANO.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.sigma.ANO$predicted,xv.sigma.ANO$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.sigma.ANO$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Parâmetros Verão

dados <- read.csv("ParamDJF.csv", head=T,dec=",",sep= ";")[,-c(1,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC.csv",h=T,sep=",",dec = ".");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
flori <- read.csv("Limite-Floripa.csv",h=T,sep=",",dec = ".")[,-3];flori$Longitude<-flori$Longitude/1000;flori$Latitude<-flori$Latitude/1000


## Análise exploratória geoestatística

names(dados)

dim(dados)
Mu.DJF.geo<-as.geodata(dados,coords.col = 1:2, data.col = 3,borders=T)
Mu.DJF.geo$borders<-limite



plot(Mu.DJF.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

DJF.mu <- Mu.DJF.geo$data+5

boxcox(DJF.mu~1) 
boxcox(DJF.mu~1,lambda=seq(0,1,l=20))
## 
##trans<-boxcox(DJF.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.DJF<-with(trans, x[which.max(y)]);lambda.mu.DJF # Valor máximo de Lambda
## 
##Mu.DJF.geo$data<-(((DJF.mu^(lambda.mu.DJF)) - 1)/lambda.mu.DJF);Mu.DJF.geo$data #normalizando - (X^lambda.mu.DJF)-1/lambda.mu.DJF
##DJF.mu.inv<-(((lambda.mu.DJF*Mu.DJF.geo$data)+1)^(1/lambda.mu.DJF))-5;DJF.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_muDJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Mu.DJF.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(DJF.mu~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.mu.DJF<-variog(Mu.DJF.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.mu.DJF,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
(trend_Mu.DJF <- dados[,c(5:25,1,2)])
names(trend_Mu.DJF)
ncol(trend_Mu.DJF)

## Tendência linear
Mu.DJF <- list()
for(i in 1:ncol(trend_Mu.DJF)){
Mu.DJF[[i]] <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend_Mu.DJF[,i],ini=c(0.04,23.35), nug=0.02,)
Mu.DJF[[ncol(trend_Mu.DJF)+1]] <- likfit(Mu.DJF.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Mu.DJF, AIC)

names(trend_Mu.DJF)
summary(Mu.DJF[[order(sapply(Mu.DJF, AIC))[1]]])
trend1_Mu.DJF<-trend_Mu.DJF[,order(sapply(Mu.DJF, AIC))[1]];trend_Mu.DJF <- trend_Mu.DJF[,-order(sapply(Mu.DJF, AIC))[1]]

Mu.DJF1 <- list()
for(i in 1:ncol(trend_Mu.DJF)){
Mu.DJF1[[i]] <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend1_Mu.DJF+trend_Mu.DJF[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.DJF)
sapply(Mu.DJF1, AIC)
summary(Mu.DJF1[[order(sapply(Mu.DJF1, AIC))[1]]])

## Tendência quadrática
trend_Mu.DJF <- dados[,c(5:25,1,2)]
names(trend_Mu.DJF)

Mu.DJF2 <- list()
for(i in 1:ncol(trend_Mu.DJF)){
Mu.DJF2[[i]] <- likfit(Mu.DJF.geo, cov.model="gau", trend =~poly(trend_Mu.DJF[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Mu.DJF2, AIC)
summary(Mu.DJF2[[order(sapply(Mu.DJF2, AIC))[1]]])
trend3_Mu.DJF<-trend_Mu.DJF[,order(sapply(Mu.DJF2, AIC))[1]];trend_Mu.DJF <- trend_Mu.DJF[,-order(sapply(Mu.DJF2, AIC))[1]]

Mu.DJF3 <- list()
for(i in 1:ncol(trend_Mu.DJF)){
Mu.DJF3[[i]] <- likfit(Mu.DJF.geo, cov.model="gau", trend =~poly(trend3_Mu.DJF,2)+poly(trend_Mu.DJF[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.DJF)
sapply(Mu.DJF3, AIC)
summary(Mu.DJF3[[order(sapply(Mu.DJF3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Mu.DJF <- dados[,c(5:25,1,2)]
names(trend_Mu.DJF)

Mu.DJF <- list()
Mu.DJF$lf1 <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend_Mu.DJF[,3]+trend_Mu.DJF[,12],ini=c(0.04,23.35), nug=0.02,)
Mu.DJF$lf2 <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend_Mu.DJF[,3]+poly(trend_Mu.DJF[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.DJF$lf3 <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend_Mu.DJF[,12]+poly(trend_Mu.DJF[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.DJF$lf4 <- likfit(Mu.DJF.geo, cov.model="gau", trend =~poly(trend_Mu.DJF[,3],2)+poly(trend_Mu.DJF[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.DJF$lf5 <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend_Mu.DJF[,12]+poly(trend_Mu.DJF[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.DJF$lf6 <- likfit(Mu.DJF.geo, cov.model="gau", trend =~trend_Mu.DJF[,19]+poly(trend_Mu.DJF[,3],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Mu.DJF, AIC))
plot(Mu.DJF.geo,low=T,trend=~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2))

postscript("Figuras/Vario-mu_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.mu.DJF<-variog(Mu.DJF.geo,max.dist=200,uvec=seq(0, 200, by=10),trend =~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2))
plot(v.mu.DJF,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.mu.DJF$u,v.mu.DJF$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_DJF.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Mu.DJFlf4<-list()
Mu.DJFlf4$exp <- likfit(Mu.DJF.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2))
Mu.DJFlf4$gau <- likfit(Mu.DJF.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2))  
Mu.DJFlf4$sph <- likfit(Mu.DJF.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2)) 
Mu.DJFlf4$cir <- likfit(Mu.DJF.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2))
Mu.DJFlf4$kappa1.5 <- likfit(Mu.DJF.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2))  
Mu.DJFlf4$kappa2.5 <- likfit(Mu.DJF.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Mu.DJF[,19],2)+poly(trend_Mu.DJF[,3],2)) 


##-----------------------------------------------------------------------------


sort(sapply(Mu.DJFlf4, AIC))

summary(Mu.DJFlf4$gau)
Mu.DJFlf4$gau$parameters.summary
Mu.DJFlf4$gau$phi
Mu.DJFlf4$gau$sigmasq
Mu.DJFlf4$gau$tausq
  
plot(v.mu.DJF,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Mu.DJFlf4$exp,col=1)
lines(Mu.DJFlf4$gau,col=2)
lines(Mu.DJFlf4$sph,col=3)
lines(Mu.DJFlf4$cir,col=1)
lines(Mu.DJFlf4$kappa1.5,col=5)
lines(Mu.DJFlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Mu.DJF.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(dados)

## Krigagem
Mu.DJFlf4$gau
kc.mu.DJF <- krige.conv(Mu.DJF.geo, locations = grid, krige=krige.control(obj=Mu.DJFlf4$gau))
attributes(kc.mu.DJF)
##kc.mu.DJF$predict<-(backtransform.moments(lambda=lambda.mu.DJF,mean=kc.mu.DJF$predict,variance=kc.mu.DJF$krige.var)$mean)-5
kc.mu.DJF$predict

Mu_DJF_krige <- cbind(grid*1000,kc.mu.DJF$predict)
write.table(Mu_DJF_krige,"MuDJF.ascii",col.names = F,row.names=F,quote=F)


summary(kc.mu.DJF$predict)

postscript("Figuras/mu_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.mu.DJF,val=kc.mu.DJF$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.mu.DJF, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.mu.DJF<-xvalid(Mu.DJF.geo,model=Mu.DJFlf4$gau)
names(xv.mu.DJF)

plot(xv.mu.DJF$predicted,xv.mu.DJF$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.mu.DJF$data~xv.mu.DJF$predicted)$coef[1],lm(xv.mu.DJF$data~xv.mu.DJF$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.mu.DJF)

smry(Mu.DJF.geo$data)


postscript("Figuras/xv_mu-DJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.mu.DJF$predicted,xv.mu.DJF$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.mu.DJF$std.error,breaks = 10,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Parametros DJF Sigma


## Análise exploratória geoestatística

names(dados)

dim(dados)
Sigma.DJF.geo<-as.geodata(dados,coords.col = 1:2, data.col = 4,borders=T)
Sigma.DJF.geo$borders<-limite

plot(Sigma.DJF.geo,low=T)
summary(Sigma.DJF.geo)

## Valores positivos para BOXCOX



boxcox(Sigma.DJF.geo$data~1) 
boxcox(Sigma.DJF.geo$data~1,lambda=seq(-2,-1,l=20))
## 
trans<-boxcox(Sigma.DJF.geo$data~1,lambda=seq(-1.4,-1.2,l=20));trans # dando zoom para ver qual lambda usar
lambda.sigma.DJF<-with(trans, x[which.max(y)]);lambda.sigma.DJF # Valor máximo de Lambda
 
Sigma.DJF.geo$data<-(((Sigma.DJF.geo$data^(lambda.sigma.DJF)) - 1)/lambda.sigma.DJF);Sigma.DJF.geo$data #normalizando - (X^lambda.sigma.DJF)-1/lambda.sigma.DJF
DJF.sigma.inv<-(((lambda.sigma.DJF*Sigma.DJF.geo$data)+1)^(1/lambda.sigma.DJF))-5;DJF.sigma.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_sigmaDJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Sigma.DJF.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(Sigma.DJF.geo$data~1,ylab = "Log-Verossimilhança",lambda=seq(-5,5)) 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.sigma.DJF<-variog(Sigma.DJF.geo,max.dist=200,uvec=seq(0, 200, by=10),trend=)
plot(v.sigma.DJF,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)

trend_Sigma.DJF <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.DJF)
ncol(trend_Sigma.DJF)

## Tendência linear
Sigma.DJF <- list()
for(i in 1:ncol(trend_Sigma.DJF)){
Sigma.DJF[[i]] <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend_Sigma.DJF[,i],ini=c(0.04,23.35), nug=0.02,)
Sigma.DJF[[ncol(trend_Sigma.DJF)+1]] <- likfit(Sigma.DJF.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Sigma.DJF, AIC)

names(trend_Sigma.DJF)
summary(Sigma.DJF[[order(sapply(Sigma.DJF, AIC))[1]]])
trend1_Sigma.DJF<-trend_Sigma.DJF[,order(sapply(Sigma.DJF, AIC))[1]];trend_Sigma.DJF <- trend_Sigma.DJF[,-order(sapply(Sigma.DJF, AIC))[1]]

Sigma.DJF1 <- list()
for(i in 1:ncol(trend_Sigma.DJF)){
Sigma.DJF1[[i]] <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend1_Sigma.DJF+trend_Sigma.DJF[,i],ini=c(0.04,10.35), nug=0.02,)

}

names(trend_Sigma.DJF)
sapply(Sigma.DJF1, AIC)
summary(Sigma.DJF1[[order(sapply(Sigma.DJF1, AIC))[1]]])

## Tendência quadrática
trend_Sigma.DJF <- dados[,c(5,6,8:25,1,2)]


Sigma.DJF2 <- list()
for(i in 1:ncol(trend_Sigma.DJF)){
Sigma.DJF2[[i]] <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~poly(trend_Sigma.DJF[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Sigma.DJF2, AIC)
summary(Sigma.DJF2[[order(sapply(Sigma.DJF2, AIC))[1]]])
names(trend_Sigma.DJF)
trend3_Sigma.DJF<-trend_Sigma.DJF[,order(sapply(Sigma.DJF2, AIC))[1]];trend_Sigma.DJF <- trend_Sigma.DJF[,-order(sapply(Sigma.DJF2, AIC))[1]]

Sigma.DJF3 <- list()
for(i in 1:ncol(trend_Sigma.DJF)){
Sigma.DJF3[[i]] <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~poly(trend3_Sigma.DJF,2)+poly(trend_Sigma.DJF[,i],2),ini=c(0.04,23.35), nug=0.02,)
}

names(trend_Sigma.DJF)
sapply(Sigma.DJF3, AIC)
summary(Sigma.DJF3[[order(sapply(Sigma.DJF3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Sigma.DJF <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.DJF)

Sigma.DJF <- list()
Sigma.DJF$lf1 <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend_Sigma.DJF[,10]+trend_Sigma.DJF[,12],ini=c(0.04,23.35), nug=0.02,)
Sigma.DJF$lf2 <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend_Sigma.DJF[,10]+poly(trend_Sigma.DJF[,6],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.DJF$lf3 <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend_Sigma.DJF[,10]+poly(trend_Sigma.DJF[,18],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.DJF$lf4 <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.DJF$lf5 <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend_Sigma.DJF[,12]+poly(trend_Sigma.DJF[,6],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.DJF$lf6 <- likfit(Sigma.DJF.geo, cov.model="gau", trend =~trend_Sigma.DJF[,12]+poly(trend_Sigma.DJF[,18],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Sigma.DJF, AIC))
plot(Sigma.DJF.geo,low=T,trend=~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2))

postscript("Figuras/Vario-sigma_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.sigma.DJF<-variog(Sigma.DJF.geo,max.dist=150,uvec=seq(0, 150, by=10),trend =~poly(trend_Sigma.DJF[,19],2)+poly(trend_Sigma.DJF[,7],2))
plot(v.sigma.DJF,xlab="Distância (km)",ylab=expression(gamma(u)),ylim = c(0.025,0.08))
lines(lowess(v.sigma.DJF$u,v.sigma.DJF$v))
dev.off()


##-----------------------------------------------------------------------------##
## uMelhor modelo trend=~Prec_media+Param_DJF.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Sigma.DJFlf4<-list()
Sigma.DJFlf4$exp <- likfit(Sigma.DJF.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2))
Sigma.DJFlf4$gau <- likfit(Sigma.DJF.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2))  
Sigma.DJFlf4$sph <- likfit(Sigma.DJF.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2)) 
Sigma.DJFlf4$cir <- likfit(Sigma.DJF.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2))
Sigma.DJFlf4$kappa1.5 <- likfit(Sigma.DJF.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2))  
Sigma.DJFlf4$kappa2.5 <- likfit(Sigma.DJF.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Sigma.DJF[,18],2)+poly(trend_Sigma.DJF[,6],2)) 


##-----------------------------------------------------------------------------

sort(sapply(Sigma.DJFlf4, AIC))

summary(Sigma.DJFlf4$sph)
Sigma.DJFlf4$sph$parameters.summary
Sigma.DJFlf4$sph$phi
Sigma.DJFlf4$sph$sigmasq
Sigma.DJFlf4$sph$tausq
  
plot(v.sigma.DJF,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Sigma.DJFlf4$exp,col=1)
lines(Sigma.DJFlf4$gau,col=2)
lines(Sigma.DJFlf4$sph,col=3)
lines(Sigma.DJFlf4$cir,col=1)
lines(Sigma.DJFlf4$kappa1.5,col=5)
lines(Sigma.DJFlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Sigma.DJF.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Sigma.DJFlf4$sph
kc.sigma.DJF <- krige.conv(Sigma.DJF.geo, locations = grid, krige=krige.control(obj=Sigma.DJFlf4$sph))
attributes(kc.sigma.DJF)
kc.sigma.DJF$predict<-(backtransform.moments(lambda=lambda.sigma.DJF,mean=kc.sigma.DJF$predict,variance=kc.sigma.DJF$krige.var)$mean)
mean(kc.sigma.DJF$predict)

Sigma_DJF_krige <- cbind(grid*1000,kc.sigma.DJF$predict)
write.table(Sigma_DJF_krige,"SigmaDJF.ascii",col.names = F,row.names=F,quote=F)


summary(kc.sigma.DJF$predict)

postscript("Figuras/sigma_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.DJF,val=kc.sigma.DJF$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.DJF, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.sigma.DJF<-xvalid(Sigma.DJF.geo,model=Sigma.DJFlf4$sph)
names(xv.sigma.DJF)

plot(xv.sigma.DJF$predicted,xv.sigma.DJF$data)
abline(0,1)
abline(lm(xv.sigma.DJF$data~xv.sigma.DJF$predicted)$coef[1],lm(xv.sigma.DJF$data~xv.sigma.DJF$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.sigma.DJF)

smry(Sigma.DJF.geo$data)


postscript("Figuras/xv_sigma-DJF.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.sigma.DJF$predicted,xv.sigma.DJF$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.sigma.DJF$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##


## Parâmetros Outono

dados <- read.csv("ParamMAM.csv", head=T,dec=",",sep= ";")[,-c(1,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC.csv",h=T,sep=",",dec = ".");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
flori <- read.csv("Limite-Floripa.csv",h=T,sep=",",dec = ".")[,-3];flori$Longitude<-flori$Longitude/1000;flori$Latitude<-flori$Latitude/1000

## Parametros MAM mu

## Análise exploratória geoestatística

names(dados)

dim(dados)
Mu.MAM.geo<-as.geodata(dados,coords.col = 1:2, data.col = 3,borders=T)
Mu.MAM.geo$borders<-limite



plot(Mu.MAM.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

MAM.mu <- Mu.MAM.geo$data+5

boxcox(MAM.mu~1) 
boxcox(MAM.mu~1,lambda=seq(0,1,l=20))
## 
##trans<-boxcox(MAM.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.MAM<-with(trans, x[which.max(y)]);lambda.mu.MAM # Valor máximo de Lambda
## 
##Mu.MAM.geo$data<-(((MAM.mu^(lambda.mu.MAM)) - 1)/lambda.mu.MAM);Mu.MAM.geo$data #normalizando - (X^lambda.mu.MAM)-1/lambda.mu.MAM
##MAM.mu.inv<-(((lambda.mu.MAM*Mu.MAM.geo$data)+1)^(1/lambda.mu.MAM))-5;MAM.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_muMAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Mu.MAM.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(MAM.mu~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.mu.MAM<-variog(Mu.MAM.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.mu.MAM,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
(trend_Mu.MAM <- dados[,c(5:25,1,2)])
names(trend_Mu.MAM)
ncol(trend_Mu.MAM)

## Tendência linear
Mu.MAM <- list()
for(i in 1:ncol(trend_Mu.MAM)){
Mu.MAM[[i]] <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend_Mu.MAM[,i],ini=c(0.04,23.35), nug=0.02,)
Mu.MAM[[ncol(trend_Mu.MAM)+1]] <- likfit(Mu.MAM.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Mu.MAM, AIC)

names(trend_Mu.MAM)
summary(Mu.MAM[[order(sapply(Mu.MAM, AIC))[1]]])
trend1_Mu.MAM<-trend_Mu.MAM[,order(sapply(Mu.MAM, AIC))[1]];trend_Mu.MAM <- trend_Mu.MAM[,-order(sapply(Mu.MAM, AIC))[1]]

Mu.MAM1 <- list()
for(i in 1:ncol(trend_Mu.MAM)){
Mu.MAM1[[i]] <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend1_Mu.MAM+trend_Mu.MAM[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.MAM)
sapply(Mu.MAM1, AIC)
summary(Mu.MAM1[[order(sapply(Mu.MAM1, AIC))[1]]])

## Tendência quadrática
trend_Mu.MAM <- dados[,c(5:25,1,2)]
names(trend_Mu.MAM)

Mu.MAM2 <- list()
for(i in 1:ncol(trend_Mu.MAM)){
Mu.MAM2[[i]] <- likfit(Mu.MAM.geo, cov.model="gau", trend =~poly(trend_Mu.MAM[,i],2),ini=c(0.04,50.35), nug=0.02,)

}

sapply(Mu.MAM2, AIC)
summary(Mu.MAM2[[order(sapply(Mu.MAM2, AIC))[1]]])
trend3_Mu.MAM<-trend_Mu.MAM[,order(sapply(Mu.MAM2, AIC))[1]];trend_Mu.MAM <- trend_Mu.MAM[,-order(sapply(Mu.MAM2, AIC))[1]]

Mu.MAM3 <- list()
for(i in 1:ncol(trend_Mu.MAM)){
Mu.MAM3[[i]] <- likfit(Mu.MAM.geo, cov.model="gau", trend =~poly(trend3_Mu.MAM,2)+poly(trend_Mu.MAM[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.MAM)
sapply(Mu.MAM3, AIC)
summary(Mu.MAM3[[order(sapply(Mu.MAM3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Mu.MAM <- dados[,c(5:25,1,2)]
names(trend_Mu.MAM)

Mu.MAM <- list()
Mu.MAM$lf1 <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12],ini=c(0.04,23.35), nug=0.02,)
Mu.MAM$lf2 <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend_Mu.MAM[,3]+poly(trend_Mu.MAM[,12],2),ini=c(0.04,23.35), nug=0.02,)
Mu.MAM$lf3 <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend_Mu.MAM[,12]+poly(trend_Mu.MAM[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.MAM$lf4 <- likfit(Mu.MAM.geo, cov.model="gau", trend =~poly(trend_Mu.MAM[,3],2)+poly(trend_Mu.MAM[,12],2),ini=c(0.04,23.35), nug=0.02,)
Mu.MAM$lf5 <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend_Mu.MAM[,12]+poly(trend_Mu.MAM[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.MAM$lf6 <- likfit(Mu.MAM.geo, cov.model="gau", trend =~trend_Mu.MAM[,19]+poly(trend_Mu.MAM[,3],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Mu.MAM, AIC))
plot(Mu.MAM.geo,low=T,trend=~trend_Mu.MAM[,3]+trend_Mu.MAM[,12])

postscript("Figuras/Vario-mu_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.mu.MAM<-variog(Mu.MAM.geo,max.dist=150,uvec=seq(0, 150, by=15),trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12])
plot(v.mu.MAM,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.mu.MAM$u,v.mu.MAM$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_MAM.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Mu.MAMlf1<-list()
Mu.MAMlf1$exp <- likfit(Mu.MAM.geo, ini=c(0.04,23.35), nug=0.02,trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12])
Mu.MAMlf1$gau <- likfit(Mu.MAM.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12])  
Mu.MAMlf1$sph <- likfit(Mu.MAM.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12]) 
Mu.MAMlf1$cir <- likfit(Mu.MAM.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12])
Mu.MAMlf1$kappa1.5 <- likfit(Mu.MAM.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~trend_Mu.MAM[,3]+trend_Mu.MAM[,12])  
Mu.MAMlf1$kappa2.5 <- likfit(Mu.MAM.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~trend_Mu.MAM[,3]+trend_Mu.MAM[,12]) 


##-----------------------------------------------------------------------------


sort(sapply(Mu.MAMlf1, AIC))

summary(Mu.MAMlf1$gau)
Mu.MAMlf1$gau$parameters.summary
Mu.MAMlf1$gau$phi
Mu.MAMlf1$gau$sigmasq
Mu.MAMlf1$gau$tausq
  
plot(v.mu.MAM,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Mu.MAMlf1$exp,col=1)
lines(Mu.MAMlf1$gau,col=2)
lines(Mu.MAMlf1$sph,col=3)
lines(Mu.MAMlf1$cir,col=1)
lines(Mu.MAMlf1$kappa1.5,col=5)
lines(Mu.MAMlf1$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Mu.MAM.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Mu.MAMlf1$gau
kc.mu.MAM <- krige.conv(Mu.MAM.geo, locations = grid, krige=krige.control(obj=Mu.MAMlf1$gau))
attributes(kc.mu.MAM)
##kc.mu.MAM$predict<-(backtransform.moments(lambda=lambda.mu.MAM,mean=kc.mu.MAM$predict,variance=kc.mu.MAM$krige.var)$mean)-5
kc.mu.MAM$predict

Mu_MAM_krige <- cbind(grid*1000,kc.mu.MAM$predict)
write.table(Mu_MAM_krige,"MuMAM.ascii",col.names = F,row.names=F,quote=F)


summary(kc.mu.MAM$predict)

postscript("Figuras/mu_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.mu.MAM,val=kc.mu.MAM$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.mu.MAM, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.mu.MAM<-xvalid(Mu.MAM.geo,model=Mu.MAMlf1$gau)
names(xv.mu.MAM)

plot(xv.mu.MAM$predicted,xv.mu.MAM$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.mu.MAM$data~xv.mu.MAM$predicted)$coef[1],lm(xv.mu.MAM$data~xv.mu.MAM$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.mu.MAM)

smry(Mu.MAM.geo$data)


postscript("Figuras/xv_mu-MAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.mu.MAM$predicted,xv.mu.MAM$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.mu.MAM$std.error,breaks = 10,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Parametros MAM Sigma


## Análise exploratória geoestatística

names(dados)

dim(dados)
Sigma.MAM.geo<-as.geodata(dados,coords.col = 1:2, data.col = 4,borders=T)
Sigma.MAM.geo$borders<-limite

plot(Sigma.MAM.geo,low=T)
summary(Sigma.MAM.geo)

## Valores positivos para BOXCOX



boxcox(Sigma.MAM.geo$data~1) 
boxcox(Sigma.MAM.geo$data~1,lambda=seq(-2,-1,l=20))
## 
trans<-boxcox(Sigma.MAM.geo$data~1,lambda=seq(-1.4,-1.2,l=20));trans # dando zoom para ver qual lambda usar
lambda.sigma.MAM<-with(trans, x[which.max(y)]);lambda.sigma.MAM # Valor máximo de Lambda
 
Sigma.MAM.geo$data<-(((Sigma.MAM.geo$data^(lambda.sigma.MAM)) - 1)/lambda.sigma.MAM);Sigma.MAM.geo$data #normalizando - (X^lambda.sigma.MAM)-1/lambda.sigma.MAM
MAM.sigma.inv<-(((lambda.sigma.MAM*Sigma.MAM.geo$data)+1)^(1/lambda.sigma.MAM))-5;MAM.sigma.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_sigmaMAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Sigma.MAM.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(Sigma.MAM.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.sigma.MAM<-variog(Sigma.MAM.geo,max.dist=200,uvec=seq(0, 200, by=10),trend=)
plot(v.sigma.MAM,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)

trend_Sigma.MAM <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.MAM)
ncol(trend_Sigma.MAM)

## Tendência linear
Sigma.MAM <- list()
for(i in 1:ncol(trend_Sigma.MAM)){
Sigma.MAM[[i]] <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend_Sigma.MAM[,i],ini=c(0.04,23.35), nug=0.02,)
Sigma.MAM[[ncol(trend_Sigma.MAM)+1]] <- likfit(Sigma.MAM.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Sigma.MAM, AIC)

names(trend_Sigma.MAM)
summary(Sigma.MAM[[order(sapply(Sigma.MAM, AIC))[1]]])
trend1_Sigma.MAM<-trend_Sigma.MAM[,order(sapply(Sigma.MAM, AIC))[1]];trend_Sigma.MAM <- trend_Sigma.MAM[,-order(sapply(Sigma.MAM, AIC))[1]]

Sigma.MAM1 <- list()
for(i in 1:ncol(trend_Sigma.MAM)){
Sigma.MAM1[[i]] <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend1_Sigma.MAM+trend_Sigma.MAM[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Sigma.MAM)
sapply(Sigma.MAM1, AIC)
summary(Sigma.MAM1[[order(sapply(Sigma.MAM1, AIC))[1]]])

## Tendência quadrática
trend_Sigma.MAM <- dados[,c(5,6,8:25,1,2)]


Sigma.MAM2 <- list()
for(i in 1:ncol(trend_Sigma.MAM)){
Sigma.MAM2[[i]] <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~poly(trend_Sigma.MAM[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Sigma.MAM2, AIC)
summary(Sigma.MAM2[[order(sapply(Sigma.MAM2, AIC))[1]]])
names(trend_Sigma.MAM)
trend3_Sigma.MAM<-trend_Sigma.MAM[,order(sapply(Sigma.MAM2, AIC))[1]];trend_Sigma.MAM <- trend_Sigma.MAM[,-order(sapply(Sigma.MAM2, AIC))[1]]

Sigma.MAM3 <- list()
for(i in 1:ncol(trend_Sigma.MAM)){
Sigma.MAM3[[i]] <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~poly(trend3_Sigma.MAM,2)+poly(trend_Sigma.MAM[,i],2),ini=c(0.04,23.35), nug=0.02,)
}

names(trend_Sigma.MAM)
sapply(Sigma.MAM3, AIC)
summary(Sigma.MAM3[[order(sapply(Sigma.MAM3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Sigma.MAM <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.MAM)

Sigma.MAM <- list()
Sigma.MAM$lf1 <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend_Sigma.MAM[,13]+trend_Sigma.MAM[,12],ini=c(0.04,23.35), nug=0.02,)
Sigma.MAM$lf2 <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend_Sigma.MAM[,13]+poly(trend_Sigma.MAM[,8],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.MAM$lf3 <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend_Sigma.MAM[,12]+poly(trend_Sigma.MAM[,21],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.MAM$lf4 <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~poly(trend_Sigma.MAM[,8],2)+poly(trend_Sigma.MAM[,21],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.MAM$lf5 <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend_Sigma.MAM[,12]+poly(trend_Sigma.MAM[,21],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.MAM$lf6 <- likfit(Sigma.MAM.geo, cov.model="gau", trend =~trend_Sigma.MAM[,13]+poly(trend_Sigma.MAM[,12],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Sigma.MAM, AIC))
plot(Sigma.MAM.geo,low=T,trend=~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2))

postscript("Figuras/Vario-sigma_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.sigma.MAM<-variog(Sigma.MAM.geo,max.dist=150,uvec=seq(0, 150, by=10),trend =~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2))
plot(v.sigma.MAM,xlab="Distância (km)",ylab=expression(gamma(u)),ylim = c(0.01,0.035))
lines(lowess(v.sigma.MAM$u,v.sigma.MAM$v))
dev.off()


##-----------------------------------------------------------------------------##
## uMelhor modelo trend=~Prec_media+Param_MAM.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Sigma.MAMlf4<-list()
Sigma.MAMlf4$exp <- likfit(Sigma.MAM.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2))
Sigma.MAMlf4$gau <- likfit(Sigma.MAM.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2))  
Sigma.MAMlf4$sph <- likfit(Sigma.MAM.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2)) 
Sigma.MAMlf4$cir <- likfit(Sigma.MAM.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2))
Sigma.MAMlf4$kappa1.5 <- likfit(Sigma.MAM.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2))  
Sigma.MAMlf4$kappa2.5 <- likfit(Sigma.MAM.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Sigma.MAM[,21],2)+poly(trend_Sigma.MAM[,8],2)) 


##-----------------------------------------------------------------------------

sort(sapply(Sigma.MAMlf4, AIC))

summary(Sigma.MAMlf4$sph)
Sigma.MAMlf4$sph$parameters.summary
Sigma.MAMlf4$sph$phi
Sigma.MAMlf4$sph$sigmasq
Sigma.MAMlf4$sph$tausq
  
plot(v.sigma.MAM,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Sigma.MAMlf4$exp,col=1)
lines(Sigma.MAMlf4$gau,col=2)
lines(Sigma.MAMlf4$sph,col=3)
lines(Sigma.MAMlf4$cir,col=1)
lines(Sigma.MAMlf4$kappa1.5,col=5)
lines(Sigma.MAMlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Sigma.MAM.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Sigma.MAMlf4$sph
kc.sigma.MAM <- krige.conv(Sigma.MAM.geo, locations = grid, krige=krige.control(obj=Sigma.MAMlf4$sph))
attributes(kc.sigma.MAM)
##kc.sigma.MAM$predict<-(backtransform.moments(lambda=lambda.sigma.MAM,mean=kc.sigma.MAM$predict,variance=kc.sigma.MAM$krige.var)$mean)
##kc.sigma.MAM$predict

Sigma_MAM_krige <- cbind(grid*1000,kc.sigma.MAM$predict)
write.table(Sigma_MAM_krige,"SigmaMAM.ascii",col.names = F,row.names=F,quote=F)


summary(kc.sigma.MAM$predict)

postscript("Figuras/sigma_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.MAM,val=kc.sigma.MAM$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.MAM, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.sigma.MAM<-xvalid(Sigma.MAM.geo,model=Sigma.MAMlf4$sph)
names(xv.sigma.MAM)

plot(xv.sigma.MAM$predicted,xv.sigma.MAM$data)
abline(0,1)
abline(lm(xv.sigma.MAM$data~xv.sigma.MAM$predicted)$coef[1],lm(xv.sigma.MAM$data~xv.sigma.MAM$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.sigma.MAM)

smry(Sigma.MAM.geo$data)


postscript("Figuras/xv_sigma-MAM.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.sigma.MAM$predicted,xv.sigma.MAM$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.sigma.MAM$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##


## Parâmetros Inverno

dados <- read.csv("ParamJJA.csv", head=T,dec=",",sep= ";")[,-c(1,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC.csv",h=T,sep=",",dec = ".");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
flori <- read.csv("Limite-Floripa.csv",h=T,sep=",",dec = ".")[,-3];flori$Longitude<-flori$Longitude/1000;flori$Latitude<-flori$Latitude/1000

## Parametros JJA mu

## Análise exploratória geoestatística

names(dados)

dim(dados)
Mu.JJA.geo<-as.geodata(dados,coords.col = 1:2, data.col = 3,borders=T)
Mu.JJA.geo$borders<-limite



plot(Mu.JJA.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

JJA.mu <- Mu.JJA.geo$data+5

boxcox(JJA.mu~1) 
boxcox(JJA.mu~1,lambda=seq(0,1,l=20))
## 
##trans<-boxcox(JJA.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.JJA<-with(trans, x[which.max(y)]);lambda.mu.JJA # Valor máximo de Lambda
## 
##Mu.JJA.geo$data<-(((JJA.mu^(lambda.mu.JJA)) - 1)/lambda.mu.JJA);Mu.JJA.geo$data #normalizando - (X^lambda.mu.JJA)-1/lambda.mu.JJA
##JJA.mu.inv<-(((lambda.mu.JJA*Mu.JJA.geo$data)+1)^(1/lambda.mu.JJA))-5;JJA.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_muJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Mu.JJA.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(JJA.mu~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.mu.JJA<-variog(Mu.JJA.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.mu.JJA,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
(trend_Mu.JJA <- dados[,c(5:25,1,2)])
names(trend_Mu.JJA)
ncol(trend_Mu.JJA)

## Tendência linear
Mu.JJA <- list()
for(i in 1:ncol(trend_Mu.JJA)){
Mu.JJA[[i]] <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend_Mu.JJA[,i],ini=c(0.04,23.35), nug=0.02,)
Mu.JJA[[ncol(trend_Mu.JJA)+1]] <- likfit(Mu.JJA.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Mu.JJA, AIC)

names(trend_Mu.JJA)
summary(Mu.JJA[[order(sapply(Mu.JJA, AIC))[1]]])
trend1_Mu.JJA<-trend_Mu.JJA[,order(sapply(Mu.JJA, AIC))[1]];trend_Mu.JJA <- trend_Mu.JJA[,-order(sapply(Mu.JJA, AIC))[1]]

Mu.JJA1 <- list()
for(i in 1:ncol(trend_Mu.JJA)){
Mu.JJA1[[i]] <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend1_Mu.JJA+trend_Mu.JJA[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.JJA)
sapply(Mu.JJA1, AIC)
summary(Mu.JJA1[[order(sapply(Mu.JJA1, AIC))[2]]])

## Tendência quadrática
trend_Mu.JJA <- dados[,c(5:25,1,2)]
names(trend_Mu.JJA)

Mu.JJA2 <- list()
for(i in 1:ncol(trend_Mu.JJA)){
Mu.JJA2[[i]] <- likfit(Mu.JJA.geo, cov.model="gau", trend =~poly(trend_Mu.JJA[,i],2),ini=c(0.04,50.35), nug=0.02,)

}

sapply(Mu.JJA2, AIC)
summary(Mu.JJA2[[order(sapply(Mu.JJA2, AIC))[1]]])
trend3_Mu.JJA<-trend_Mu.JJA[,order(sapply(Mu.JJA2, AIC))[1]];trend_Mu.JJA <- trend_Mu.JJA[,-order(sapply(Mu.JJA2, AIC))[1]]

Mu.JJA3 <- list()
for(i in 1:ncol(trend_Mu.JJA)){
Mu.JJA3[[i]] <- likfit(Mu.JJA.geo, cov.model="gau", trend =~poly(trend3_Mu.JJA,2)+poly(trend_Mu.JJA[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.JJA)
sapply(Mu.JJA3, AIC)
summary(Mu.JJA3[[order(sapply(Mu.JJA3, AIC))[2]]])

## Testar combinações lineares e quadráticas
trend_Mu.JJA <- dados[,c(5:25,1,2)]
names(trend_Mu.JJA)

Mu.JJA <- list()
Mu.JJA$lf1 <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14],ini=c(0.04,23.35), nug=0.02,)
Mu.JJA$lf2 <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend_Mu.JJA[,3]+poly(trend_Mu.JJA[,14],2),ini=c(0.04,23.35), nug=0.02,)
Mu.JJA$lf3 <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend_Mu.JJA[,14]+poly(trend_Mu.JJA[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.JJA$lf4 <- likfit(Mu.JJA.geo, cov.model="gau", trend =~poly(trend_Mu.JJA[,3],2)+poly(trend_Mu.JJA[,14],2),ini=c(0.04,23.35), nug=0.02,)
Mu.JJA$lf5 <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend_Mu.JJA[,14]+poly(trend_Mu.JJA[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.JJA$lf6 <- likfit(Mu.JJA.geo, cov.model="gau", trend =~trend_Mu.JJA[,19]+poly(trend_Mu.JJA[,3],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Mu.JJA, AIC))
plot(Mu.JJA.geo,low=T,trend=~trend_Mu.JJA[,3]+trend_Mu.JJA[,14])

postscript("Figuras/Vario-mu_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.mu.JJA<-variog(Mu.JJA.geo,max.dist=150,uvec=seq(0, 150, by=13),trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14])
plot(v.mu.JJA,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.mu.JJA$u,v.mu.JJA$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_JJA.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Mu.JJAlf1<-list()
Mu.JJAlf1$exp <- likfit(Mu.JJA.geo, ini=c(0.04,23.35), nug=0.02,trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14])
Mu.JJAlf1$gau <- likfit(Mu.JJA.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14])  
Mu.JJAlf1$sph <- likfit(Mu.JJA.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14]) 
Mu.JJAlf1$cir <- likfit(Mu.JJA.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14])
Mu.JJAlf1$kappa1.5 <- likfit(Mu.JJA.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~trend_Mu.JJA[,3]+trend_Mu.JJA[,14])  
Mu.JJAlf1$kappa2.5 <- likfit(Mu.JJA.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~trend_Mu.JJA[,3]+trend_Mu.JJA[,14]) 


##-----------------------------------------------------------------------------


sort(sapply(Mu.JJAlf1, AIC))

summary(Mu.JJAlf1$sph)
Mu.JJAlf1$sph$parameters.summary
Mu.JJAlf1$sph$phi
Mu.JJAlf1$sph$sigmasq
Mu.JJAlf1$sph$tausq
  
plot(v.mu.JJA,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Mu.JJAlf1$exp,col=1)
lines(Mu.JJAlf1$gau,col=2)
lines(Mu.JJAlf1$sph,col=3)
lines(Mu.JJAlf1$cir,col=1)
lines(Mu.JJAlf1$kappa1.5,col=5)
lines(Mu.JJAlf1$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Mu.JJA.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Mu.JJAlf1$sph
kc.mu.JJA <- krige.conv(Mu.JJA.geo, locations = grid, krige=krige.control(obj=Mu.JJAlf1$sph))
attributes(kc.mu.JJA)
##kc.mu.JJA$predict<-(backtransform.moments(lambda=lambda.mu.JJA,mean=kc.mu.JJA$predict,variance=kc.mu.JJA$krige.var)$mean)-5
kc.mu.JJA$predict

Mu_JJA_krige <- cbind(grid*1000,kc.mu.JJA$predict)
write.table(Mu_JJA_krige,"MuJJA.ascii",col.names = F,row.names=F,quote=F)

summary(kc.mu.JJA$predict)

postscript("Figuras/mu_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.mu.JJA,val=kc.mu.JJA$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.mu.JJA, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.mu.JJA<-xvalid(Mu.JJA.geo,model=Mu.JJAlf1$sph)
names(xv.mu.JJA)

plot(xv.mu.JJA$predicted,xv.mu.JJA$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.mu.JJA$data~xv.mu.JJA$predicted)$coef[1],lm(xv.mu.JJA$data~xv.mu.JJA$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.mu.JJA)

smry(Mu.JJA.geo$data)


postscript("Figuras/xv_mu-JJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.mu.JJA$predicted,xv.mu.JJA$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.mu.JJA$std.error,breaks = 10,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Parametros JJA Sigma


## Análise exploratória geoestatística

names(dados)

dim(dados)
Sigma.JJA.geo<-as.geodata(dados,coords.col = 1:2, data.col = 4,borders=T)
Sigma.JJA.geo$borders<-limite

plot(Sigma.JJA.geo,low=T)
summary(Sigma.JJA.geo)

## Valores positivos para BOXCOX



boxcox(Sigma.JJA.geo$data~1) 
boxcox(Sigma.JJA.geo$data~1,lambda=seq(-2,-1,l=20))
## 
trans<-boxcox(Sigma.JJA.geo$data~1,lambda=seq(-1.4,-1.2,l=20));trans # dando zoom para ver qual lambda usar
lambda.sigma.JJA<-with(trans, x[which.max(y)]);lambda.sigma.JJA # Valor máximo de Lambda
 
Sigma.JJA.geo$data<-(((Sigma.JJA.geo$data^(lambda.sigma.JJA)) - 1)/lambda.sigma.JJA);Sigma.JJA.geo$data #normalizando - (X^lambda.sigma.JJA)-1/lambda.sigma.JJA
JJA.sigma.inv<-(((lambda.sigma.JJA*Sigma.JJA.geo$data)+1)^(1/lambda.sigma.JJA))-5;JJA.sigma.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_sigmaJJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Sigma.JJA.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(Sigma.JJA.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.sigma.JJA<-variog(Sigma.JJA.geo,max.dist=200,uvec=seq(0, 200, by=10),trend=)
plot(v.sigma.JJA,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)

trend_Sigma.JJA <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.JJA)
ncol(trend_Sigma.JJA)

## Tendência linear
Sigma.JJA <- list()
for(i in 1:ncol(trend_Sigma.JJA)){
Sigma.JJA[[i]] <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend_Sigma.JJA[,i],ini=c(0.04,23.35), nug=0.02,)
Sigma.JJA[[ncol(trend_Sigma.JJA)+1]] <- likfit(Sigma.JJA.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Sigma.JJA, AIC)

names(trend_Sigma.JJA)
summary(Sigma.JJA[[order(sapply(Sigma.JJA, AIC))[1]]])
trend1_Sigma.JJA<-trend_Sigma.JJA[,order(sapply(Sigma.JJA, AIC))[1]];trend_Sigma.JJA <- trend_Sigma.JJA[,-order(sapply(Sigma.JJA, AIC))[1]]

Sigma.JJA1 <- list()
for(i in 1:ncol(trend_Sigma.JJA)){
Sigma.JJA1[[i]] <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend1_Sigma.JJA+trend_Sigma.JJA[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Sigma.JJA)
sapply(Sigma.JJA1, AIC)
summary(Sigma.JJA1[[order(sapply(Sigma.JJA1, AIC))[1]]])

## Tendência quadrática
trend_Sigma.JJA <- dados[,c(5,6,8:25,1,2)]


Sigma.JJA2 <- list()
for(i in 1:ncol(trend_Sigma.JJA)){
Sigma.JJA2[[i]] <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~poly(trend_Sigma.JJA[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Sigma.JJA2, AIC)
summary(Sigma.JJA2[[order(sapply(Sigma.JJA2, AIC))[1]]])
names(trend_Sigma.JJA)
trend3_Sigma.JJA<-trend_Sigma.JJA[,order(sapply(Sigma.JJA2, AIC))[1]];trend_Sigma.JJA <- trend_Sigma.JJA[,-order(sapply(Sigma.JJA2, AIC))[1]]

Sigma.JJA3 <- list()
for(i in 1:ncol(trend_Sigma.JJA)){
Sigma.JJA3[[i]] <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~poly(trend3_Sigma.JJA,2)+poly(trend_Sigma.JJA[,i],2),ini=c(0.04,23.35), nug=0.02,)
}

names(trend_Sigma.JJA)
sapply(Sigma.JJA3, AIC)
summary(Sigma.JJA3[[order(sapply(Sigma.JJA3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Sigma.JJA <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.JJA)

Sigma.JJA <- list()
Sigma.JJA$lf1 <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend_Sigma.JJA[,13]+trend_Sigma.JJA[,18],ini=c(0.04,23.35), nug=0.02,)
Sigma.JJA$lf2 <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend_Sigma.JJA[,13]+poly(trend_Sigma.JJA[,18],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.JJA$lf3 <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend_Sigma.JJA[,18]+poly(trend_Sigma.JJA[,8],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.JJA$lf4 <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~poly(trend_Sigma.JJA[,8],2)+poly(trend_Sigma.JJA[,18],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.JJA$lf5 <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend_Sigma.JJA[,18]+poly(trend_Sigma.JJA[,8],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.JJA$lf6 <- likfit(Sigma.JJA.geo, cov.model="gau", trend =~trend_Sigma.JJA[,13]+poly(trend_Sigma.JJA[,8],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Sigma.JJA, AIC))
plot(Sigma.JJA.geo,low=T,trend=~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2))

postscript("Figuras/Vario-sigma_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.sigma.JJA<-variog(Sigma.JJA.geo,max.dist=150,uvec=seq(0, 150, by=10),trend =~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2))
plot(v.sigma.JJA,xlab="Distância (km)",ylab=expression(gamma(u)),ylim = c(0.01,0.035))
lines(lowess(v.sigma.JJA$u,v.sigma.JJA$v))
dev.off()


##-----------------------------------------------------------------------------##
## uMelhor modelo trend=~Prec_media+Param_JJA.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Sigma.JJAlf4<-list()
Sigma.JJAlf4$exp <- likfit(Sigma.JJA.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2))
Sigma.JJAlf4$gau <- likfit(Sigma.JJA.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2))  
Sigma.JJAlf4$sph <- likfit(Sigma.JJA.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2)) 
Sigma.JJAlf4$cir <- likfit(Sigma.JJA.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2))
Sigma.JJAlf4$kappa1.5 <- likfit(Sigma.JJA.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2))  
Sigma.JJAlf4$kappa2.5 <- likfit(Sigma.JJA.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Sigma.JJA[,18],2)+poly(trend_Sigma.JJA[,8],2)) 


##-----------------------------------------------------------------------------

sort(sapply(Sigma.JJAlf4, AIC))

summary(Sigma.JJAlf4$sph)
Sigma.JJAlf4$sph$parameters.summary
Sigma.JJAlf4$sph$phi
Sigma.JJAlf4$sph$sigmasq
Sigma.JJAlf4$sph$tausq
  
plot(v.sigma.JJA,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Sigma.JJAlf4$exp,col=1)
lines(Sigma.JJAlf4$gau,col=2)
lines(Sigma.JJAlf4$sph,col=3)
lines(Sigma.JJAlf4$cir,col=1)
lines(Sigma.JJAlf4$kappa1.5,col=5)
lines(Sigma.JJAlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Sigma.JJA.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Sigma.JJAlf4$sph
kc.sigma.JJA <- krige.conv(Sigma.JJA.geo, locations = grid, krige=krige.control(obj=Sigma.JJAlf4$sph))
attributes(kc.sigma.JJA)
##kc.sigma.JJA$predict<-(backtransform.moments(lambda=lambda.sigma.JJA,mean=kc.sigma.JJA$predict,variance=kc.sigma.JJA$krige.var)$mean)
##kc.sigma.JJA$predict

Sigma_JJA_krige <- cbind(grid*1000,kc.sigma.JJA$predict)
write.table(Sigma_JJA_krige,"SigmaJJA.ascii",col.names = F,row.names=F,quote=F)


summary(kc.sigma.JJA$predict)

postscript("Figuras/sigma_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.JJA,val=kc.sigma.JJA$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.JJA, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.sigma.JJA<-xvalid(Sigma.JJA.geo,model=Sigma.JJAlf4$sph)
names(xv.sigma.JJA)

plot(xv.sigma.JJA$predicted,xv.sigma.JJA$data)
abline(0,1)
abline(lm(xv.sigma.JJA$data~xv.sigma.JJA$predicted)$coef[1],lm(xv.sigma.JJA$data~xv.sigma.JJA$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.sigma.JJA)

smry(Sigma.JJA.geo$data)


postscript("Figuras/xv_sigma-JJA.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.sigma.JJA$predicted,xv.sigma.JJA$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos",xlim = c(0.5,1.5))
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.sigma.JJA$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##


## Parâmetros Primavera

dados <- read.csv("ParamSON.csv", head=T,dec=",",sep= ";")[,-c(1,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
limite<- read.csv("Limite_SC.csv",h=T,sep=",",dec = ".");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
flori <- read.csv("Limite-Floripa.csv",h=T,sep=",",dec = ".")[,-3];flori$Longitude<-flori$Longitude/1000;flori$Latitude<-flori$Latitude/1000

## Parametros SON mu

## Análise exploratória geoestatística

names(dados)

dim(dados)
Mu.SON.geo<-as.geodata(dados,coords.col = 1:2, data.col = 3,borders=T)
Mu.SON.geo$borders<-limite


plot(Mu.SON.geo,low=T)
summary(limite)

## Valores positivos para BOXCOX

SON.mu <- Mu.SON.geo$data+5

boxcox(SON.mu~1) 
boxcox(SON.mu~1,lambda=seq(0,1,l=20))
## 
##trans<-boxcox(SON.mu~1,lambda=seq(0.4,0.6,l=20));trans # dando zoom para ver qual lambda usar
##lambda.mu.SON<-with(trans, x[which.max(y)]);lambda.mu.SON # Valor máximo de Lambda
## 
##Mu.SON.geo$data<-(((SON.mu^(lambda.mu.SON)) - 1)/lambda.mu.SON);Mu.SON.geo$data #normalizando - (X^lambda.mu.SON)-1/lambda.mu.SON
##SON.mu.inv<-(((lambda.mu.SON*Mu.SON.geo$data)+1)^(1/lambda.mu.SON))-5;SON.mu.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_muSON.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Mu.SON.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(SON.mu~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.mu.SON<-variog(Mu.SON.geo,max.dist=200,uvec=seq(0, 200, by=18),trend=)
plot(v.mu.SON,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)
(trend_Mu.SON <- dados[,c(5:25,1,2)])
names(trend_Mu.SON)
ncol(trend_Mu.SON)

## Tendência linear
Mu.SON <- list()
for(i in 1:ncol(trend_Mu.SON)){
Mu.SON[[i]] <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend_Mu.SON[,i],ini=c(0.04,23.35), nug=0.02,)
Mu.SON[[ncol(trend_Mu.SON)+1]] <- likfit(Mu.SON.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Mu.SON, AIC)

names(trend_Mu.SON)
summary(Mu.SON[[order(sapply(Mu.SON, AIC))[1]]])
trend1_Mu.SON<-trend_Mu.SON[,order(sapply(Mu.SON, AIC))[1]];trend_Mu.SON <- trend_Mu.SON[,-order(sapply(Mu.SON, AIC))[1]]

Mu.SON1 <- list()
for(i in 1:ncol(trend_Mu.SON)){
Mu.SON1[[i]] <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend1_Mu.SON+trend_Mu.SON[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.SON)
sapply(Mu.SON1, AIC)
summary(Mu.SON1[[order(sapply(Mu.SON1, AIC))[1]]])

## Tendência quadrática
trend_Mu.SON <- dados[,c(5:25,1,2)]
names(trend_Mu.SON)

Mu.SON2 <- list()
for(i in 1:ncol(trend_Mu.SON)){
Mu.SON2[[i]] <- likfit(Mu.SON.geo, cov.model="gau", trend =~poly(trend_Mu.SON[,i],2),ini=c(0.04,50.35), nug=0.02,)

}

sapply(Mu.SON2, AIC)
summary(Mu.SON2[[order(sapply(Mu.SON2, AIC))[1]]])
trend3_Mu.SON<-trend_Mu.SON[,order(sapply(Mu.SON2, AIC))[1]];trend_Mu.SON <- trend_Mu.SON[,-order(sapply(Mu.SON2, AIC))[1]]

Mu.SON3 <- list()
for(i in 1:ncol(trend_Mu.SON)){
Mu.SON3[[i]] <- likfit(Mu.SON.geo, cov.model="gau", trend =~poly(trend3_Mu.SON,2)+poly(trend_Mu.SON[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Mu.SON)
sapply(Mu.SON3, AIC)
summary(Mu.SON3[[order(sapply(Mu.SON3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Mu.SON <- dados[,c(5:25,1,2)]
names(trend_Mu.SON)

Mu.SON <- list()
Mu.SON$lf1 <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend_Mu.SON[,3]+trend_Mu.SON[,2],ini=c(0.04,23.35), nug=0.02,)
Mu.SON$lf2 <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend_Mu.SON[,3]+poly(trend_Mu.SON[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.SON$lf3 <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend_Mu.SON[,2]+poly(trend_Mu.SON[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.SON$lf4 <- likfit(Mu.SON.geo, cov.model="gau", trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2),ini=c(0.04,23.35), nug=0.02,)
Mu.SON$lf5 <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend_Mu.SON[,2]+poly(trend_Mu.SON[,3],2),ini=c(0.04,23.35), nug=0.02,)
Mu.SON$lf6 <- likfit(Mu.SON.geo, cov.model="gau", trend =~trend_Mu.SON[,12]+poly(trend_Mu.SON[,3],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Mu.SON, AIC))
plot(Mu.SON.geo,low=T,trend=~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2))

postscript("Figuras/Vario-mu_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.mu.SON<-variog(Mu.SON.geo,max.dist=150,uvec=seq(0, 150, by=15),trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2))
plot(v.mu.SON,xlab="Distância (km)",ylab=expression(gamma(u)))
lines(lowess(v.mu.SON$u,v.mu.SON$v))
dev.off()


##-----------------------------------------------------------------------------##
## Melhor modelo trend=~Prec_media+Param_SON.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Mu.SONlf4<-list()
Mu.SONlf4$exp <- likfit(Mu.SON.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2))
Mu.SONlf4$gau <- likfit(Mu.SON.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2))  
Mu.SONlf4$sph <- likfit(Mu.SON.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2)) 
Mu.SONlf4$cir <- likfit(Mu.SON.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2))
Mu.SONlf4$kappa1.5 <- likfit(Mu.SON.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2))  
Mu.SONlf4$kappa2.5 <- likfit(Mu.SON.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Mu.SON[,3],2)+poly(trend_Mu.SON[,19],2)) 


##-----------------------------------------------------------------------------


sort(sapply(Mu.SONlf4, AIC))

summary(Mu.SONlf4$sph)
Mu.SONlf4$sph$parameters.summary
Mu.SONlf4$sph$phi
Mu.SONlf4$sph$sigmasq
Mu.SONlf4$sph$tausq
  
plot(v.mu.SON,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Mu.SONlf4$exp,col=1)
lines(Mu.SONlf4$gau,col=2)
lines(Mu.SONlf4$sph,col=3)
lines(Mu.SONlf4$cir,col=1)
lines(Mu.SONlf4$kappa1.5,col=5)
lines(Mu.SONlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Mu.SON.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Mu.SONlf4$sph
kc.mu.SON <- krige.conv(Mu.SON.geo, locations = grid, krige=krige.control(obj=Mu.SONlf4$sph))
attributes(kc.mu.SON)
##kc.mu.SON$predict<-(backtransform.moments(lambda=lambda.mu.SON,mean=kc.mu.SON$predict,variance=kc.mu.SON$krige.var)$mean)-5
kc.mu.SON$predict

Mu_SON_krige <- cbind(grid*1000,kc.mu.SON$predict)
write.table(Mu_SON_krige,"MuSON.ascii",col.names = F,row.names=F,quote=F)


summary(kc.mu.SON$predict)

postscript("Figuras/mu_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.mu.SON,val=kc.mu.SON$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.mu.SON, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.mu.SON<-xvalid(Mu.SON.geo,model=Mu.SONlf4$sph)
names(xv.mu.SON)

plot(xv.mu.SON$predicted,xv.mu.SON$data,xlim=c(-5,-3),ylim=c(-5,-3))
abline(0,1)
abline(lm(xv.mu.SON$data~xv.mu.SON$predicted)$coef[1],lm(xv.mu.SON$data~xv.mu.SON$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.mu.SON)

smry(Mu.SON.geo$data)


postscript("Figuras/xv_mu-SON.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.mu.SON$predicted,xv.mu.SON$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos")
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.mu.SON$std.error,breaks = 10,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##

## Parametros SON Sigma


## Análise exploratória geoestatística

names(dados)

dim(dados)
Sigma.SON.geo<-as.geodata(dados,coords.col = 1:2, data.col = 4,borders=T)
Sigma.SON.geo$borders<-limite

plot(Sigma.SON.geo,low=T)
summary(Sigma.SON.geo)

## Valores positivos para BOXCOX



boxcox(Sigma.SON.geo$data~1) 
boxcox(Sigma.SON.geo$data~1,lambda=seq(-2,-1,l=20))
## 
trans<-boxcox(Sigma.SON.geo$data~1,lambda=seq(-1.4,-1.2,l=20));trans # dando zoom para ver qual lambda usar
lambda.sigma.SON<-with(trans, x[which.max(y)]);lambda.sigma.SON # Valor máximo de Lambda
 
Sigma.SON.geo$data<-(((Sigma.SON.geo$data^(lambda.sigma.SON)) - 1)/lambda.sigma.SON);Sigma.SON.geo$data #normalizando - (X^lambda.sigma.SON)-1/lambda.sigma.SON
SON.sigma.inv<-(((lambda.sigma.SON*Sigma.SON.geo$data)+1)^(1/lambda.sigma.SON))-5;SON.sigma.inv #inverter e diminuir com 5 para encontrar a varivel original


postscript("Figuras/BoxCox_sigmaSON.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")
par(mfrow = c(1,2))
hist(Sigma.SON.geo$data,main="", xlab="",ylab = "Frequência")
boxcox(Sigma.SON.geo$data~1,ylab = "Log-Verossimilhança") 
dev.off()

#-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
## Conclusão: 

## Conclusão: Dado normalizado e sem tendência
##-----------------------------------------------------------------------------##

## Variograma

v.sigma.SON<-variog(Sigma.SON.geo,max.dist=200,uvec=seq(0, 200, by=10),trend=)
plot(v.sigma.SON,xlab="Distância",ylab=expression(gamma(u)))

## Interpolação espacial (krigagem)

trend_Sigma.SON <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.SON)
ncol(trend_Sigma.SON)

## Tendência linear
Sigma.SON <- list()
for(i in 1:ncol(trend_Sigma.SON)){
Sigma.SON[[i]] <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend_Sigma.SON[,i],ini=c(0.04,23.35), nug=0.02,)
Sigma.SON[[ncol(trend_Sigma.SON)+1]] <- likfit(Sigma.SON.geo, cov.model="gau",ini=c(0.04,23.35), nug=0.02,)
}


sapply(Sigma.SON, AIC)

names(trend_Sigma.SON)
summary(Sigma.SON[[order(sapply(Sigma.SON, AIC))[1]]])
trend1_Sigma.SON<-trend_Sigma.SON[,order(sapply(Sigma.SON, AIC))[1]];trend_Sigma.SON <- trend_Sigma.SON[,-order(sapply(Sigma.SON, AIC))[1]]

Sigma.SON1 <- list()
for(i in 1:ncol(trend_Sigma.SON)){
Sigma.SON1[[i]] <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend1_Sigma.SON+trend_Sigma.SON[,i],ini=c(0.04,23.35), nug=0.02,)

}

names(trend_Sigma.SON)
sapply(Sigma.SON1, AIC)
summary(Sigma.SON1[[order(sapply(Sigma.SON1, AIC))[1]]])

## Tendência quadrática
trend_Sigma.SON <- dados[,c(5,6,8:25,1,2)]


Sigma.SON2 <- list()
for(i in 1:ncol(trend_Sigma.SON)){
Sigma.SON2[[i]] <- likfit(Sigma.SON.geo, cov.model="gau", trend =~poly(trend_Sigma.SON[,i],2),ini=c(0.04,23.35), nug=0.02,)

}

sapply(Sigma.SON2, AIC)
summary(Sigma.SON2[[order(sapply(Sigma.SON2, AIC))[1]]])
names(trend_Sigma.SON)
trend3_Sigma.SON<-trend_Sigma.SON[,order(sapply(Sigma.SON2, AIC))[1]];trend_Sigma.SON <- trend_Sigma.SON[,-order(sapply(Sigma.SON2, AIC))[1]]

Sigma.SON3 <- list()
for(i in 1:ncol(trend_Sigma.SON)){
Sigma.SON3[[i]] <- likfit(Sigma.SON.geo, cov.model="gau", trend =~poly(trend3_Sigma.SON,2)+poly(trend_Sigma.SON[,i],2),ini=c(0.04,23.35), nug=0.02,)
}

names(trend_Sigma.SON)
sapply(Sigma.SON3, AIC)
summary(Sigma.SON3[[order(sapply(Sigma.SON3, AIC))[1]]])

## Testar combinações lineares e quadráticas
trend_Sigma.SON <- dados[,c(5,6,8:25,1,2)]
names(trend_Sigma.SON)

Sigma.SON <- list()
Sigma.SON$lf1 <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend_Sigma.SON[,13]+trend_Sigma.SON[,18],ini=c(0.04,23.35), nug=0.02,)
Sigma.SON$lf2 <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend_Sigma.SON[,13]+poly(trend_Sigma.SON[,8],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.SON$lf3 <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend_Sigma.SON[,13]+poly(trend_Sigma.SON[,9],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.SON$lf4 <- likfit(Sigma.SON.geo, cov.model="gau", trend =~poly(trend_Sigma.SON[,8],2)+poly(trend_Sigma.SON[,9],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.SON$lf5 <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend_Sigma.SON[,18]+poly(trend_Sigma.SON[,8],2),ini=c(0.04,23.35), nug=0.02,)
Sigma.SON$lf6 <- likfit(Sigma.SON.geo, cov.model="gau", trend =~trend_Sigma.SON[,18]+poly(trend_Sigma.SON[,9],2),ini=c(0.04,23.35), nug=0.02,)

sort(sapply(Sigma.SON, AIC))
plot(Sigma.SON.geo,low=T,trend=~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2))

postscript("Figuras/Vario-sigma_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")
v.sigma.SON<-variog(Sigma.SON.geo,max.dist=150,uvec=seq(0, 150, by=15),trend =~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2))
plot(v.sigma.SON,xlab="Distância (km)",ylab=expression(gamma(u)))#,ylim = c(0.01,0.035))
lines(lowess(v.sigma.SON$u,v.sigma.SON$v))
dev.off()


##-----------------------------------------------------------------------------##
## uMelhor modelo trend=~Prec_media+Param_SON.sigma
## Conclusão: vale a pena usar o modelo com retirada de tendência "1st" pq p-valor < 0,05 e AIC menor
##-----------------------------------------------------------------------------

Sigma.SONlf4<-list()
Sigma.SONlf4$exp <- likfit(Sigma.SON.geo, ini=c(0.04,23.35), nug=0.02,trend =~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2))
Sigma.SONlf4$gau <- likfit(Sigma.SON.geo, cov.model = "gau", ini=c(0.05,17.84), nug=0.02,trend =~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2))  
Sigma.SONlf4$sph <- likfit(Sigma.SON.geo, cov.model = "sph", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2)) 
Sigma.SONlf4$cir <- likfit(Sigma.SON.geo, cov.model = "cir", ini=c(0.04,47.57), nug=0.02,trend =~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2))
Sigma.SONlf4$kappa1.5 <- likfit(Sigma.SON.geo, cov.model = "mat", kappa= 1.5, ini=c(0.04,11.89), nug=0.02,trend =~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2))  
Sigma.SONlf4$kappa2.5 <- likfit(Sigma.SON.geo, cov.model = "mat", kappa= 2.5, ini=c(0.04,11.89), nug=0.02,trend=~poly(trend_Sigma.SON[,9],2)+poly(trend_Sigma.SON[,8],2)) 


##-----------------------------------------------------------------------------

sort(sapply(Sigma.SONlf4, AIC))

summary(Sigma.SONlf4$gau)
Sigma.SONlf4$gau$parameters.summary
Sigma.SONlf4$gau$phi
Sigma.SONlf4$gau$sigmasq
Sigma.SONlf4$gau$tausq
  
plot(v.sigma.SON,xlab="Distance (km)",ylab=expression(gamma(u)))#,ylim =c(0,4e-6)
lines(Sigma.SONlf4$exp,col=1)
lines(Sigma.SONlf4$gau,col=2)
lines(Sigma.SONlf4$sph,col=3)
lines(Sigma.SONlf4$cir,col=1)
lines(Sigma.SONlf4$kappa1.5,col=5)
lines(Sigma.SONlf4$kappa2.5,col=6)

##-----------------------------------------------------------------------------
## Predição na área

grid <- pred_grid(c(200,800),c(6700,7200), by=1)

points(Sigma.SON.geo)
points(grid, pch=19, cex=0.25, col=2)
gr <- locations.inside(grid, limite)
points(gr, pch=19, cex=0.24,col=4)
summary(gr)

## Krigagem
Sigma.SONlf4$gau
kc.sigma.SON <- krige.conv(Sigma.SON.geo, locations = grid, krige=krige.control(obj=Sigma.SONlf4$gau))
attributes(kc.sigma.SON)
##kc.sigma.SON$predict<-(backtransform.moments(lambda=lambda.sigma.SON,mean=kc.sigma.SON$predict,variance=kc.sigma.SON$krige.var)$mean)
##kc.sigma.SON$predict

Sigma_SON_krige <- cbind(grid*1000,kc.sigma.SON$predict)
write.table(Sigma_SON_krige,"SigmaSON.ascii",col.names = F,row.names=F,quote=F)


summary(kc.sigma.SON$predict)

postscript("Figuras/sigma_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.SON,val=kc.sigma.SON$predict,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.SON, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()



##-----------------------------------------------------------------------------
## Validação cruzada
xv.sigma.SON<-xvalid(Sigma.SON.geo,model=Sigma.SONlf4$gau)
names(xv.sigma.SON)

plot(xv.sigma.SON$predicted,xv.sigma.SON$data)
abline(0,1)
abline(lm(xv.sigma.SON$data~xv.sigma.SON$predicted)$coef[1],lm(xv.sigma.SON$data~xv.sigma.SON$predicted)$coef[2])


par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))

plot(xv.sigma.SON)

smry(Sigma.SON.geo$data)


postscript("Figuras/xv_sigma-SON.eps",onefile = T, horizontal = F, width=25/2.54, height=15/2.54,paper = "special")


par(mfrow = c(1,2))


plot(xv.sigma.SON$predicted,xv.sigma.SON$std.error,ylim = c(-3,3), ylab="Resíduos Padronizados", xlab = "Valores Preditos",xlim = c(0.5,1.5))
abline(h=0)
#mtext("(a)",cex=1.2,adj=0,line=1)

hist(xv.sigma.SON$std.error,breaks = ,freq = F,xlim=c(-3,3), main="", xlab="Resíduos Padronizados",ylab = "Densidade")
#mtext("(b)",cex=1.2,adj=0,line=1)

dev.off()

##-----------------------------------------------------------------------------##


## Mapas de vazões

Qp <- function(mu,sigma,area,p){
    Qp<-(exp(mu+(sigma*qnorm(1-p))))*area
}

Vr <- function(mu, sigma, area, T,Qf){

    pqfi <-1-pnorm((log(Qf/area)-mu)/sigma) # probabilidade inferior na CP
    pqff <- 1-(1/T) # probabilidade superior na CP

##Vr <- ((((pqff-pqfi)*Qf)-simpson(function(p) Qp(mu,sigma,area,p),pqfi,pqff,n=100000))*(60*60*24*365))/10^6
Vr <- ((((pqff-pqfi)*Qf)-integrate(function(p) Qp(mu,sigma,area,p),pqfi,pqff)$value)*(60*60*24*365))/10^6
return(Vr)    
}


Q98_ANO <- Qp(kc.mu.ANO$predict,kc.sigma.ANO$predict,1,0.98)*1000
Q98_DJF <- Qp(kc.mu.DJF$predict,kc.sigma.DJF$predict,1,0.98)*1000
Q98_MAM <- Qp(kc.mu.MAM$predict,kc.sigma.MAM$predict,1,0.98)*1000
Q98_JJA <- Qp(kc.mu.JJA$predict,kc.sigma.JJA$predict,1,0.98)*1000
Q98_SON <- Qp(kc.mu.SON$predict,kc.sigma.SON$predict,1,0.98)*1000

kc.sigma.ANO$Q98_ANO<- Q98_ANO
kc.sigma.DJF$Q98_DJF<- Q98_DJF
kc.sigma.MAM$Q98_MAM<- Q98_MAM
kc.sigma.JJA$Q98_JJA<- Q98_JJA
kc.sigma.SON$Q98_SON<- Q98_SON

str(kc.sigma.ANO)
attr(kc.sigma.ANO,"prediction.locations")

Q98_ANO_krige <- cbind(grid*1000,kc.sigma.ANO$Q98_ANO)
write.table(Q98_ANO_krige,"Q98ANO.ascii",col.names = F,row.names=F,quote=F)

Q98_DJF_krige <- cbind(grid*1000,kc.sigma.DJF$Q98_DJF)
write.table(Q98_DJF_krige,"Q98DJF.ascii",col.names = F,row.names=F,quote=F)

Q98_MAM_krige <- cbind(grid*1000,kc.sigma.MAM$Q98_MAM)
write.table(Q98_MAM_krige,"Q98MAM.ascii",col.names = F,row.names=F,quote=F)

Q98_JJA_krige <- cbind(grid*1000,kc.sigma.JJA$Q98_JJA)
write.table(Q98_JJA_krige,"Q98JJA.ascii",col.names = F,row.names=F,quote=F)

Q98_SON_krige <- cbind(grid*1000,kc.sigma.SON$Q98_SON)
write.table(Q98_SON_krige,"Q98SON.ascii",col.names = F,row.names=F,quote=F)


postscript("Figuras/Q98_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.ANO,val=kc.sigma.ANO$Q98_ANO,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.ANO,val=kc.sigma.ANO$Q98_ANO, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/Q98_DJF.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.DJF,val=kc.sigma.DJF$Q98_DJF,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.DJF,val=kc.sigma.DJF$Q98_DJF, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/Q98_MAM.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.MAM,val=kc.sigma.MAM$Q98_MAM,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.MAM,val=kc.sigma.MAM$Q98_MAM, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/Q98_JJA.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.JJA,val=kc.sigma.JJA$Q98_JJA,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.JJA,val=kc.sigma.JJA$Q98_JJA, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()

postscript("Figuras/Q98_SON.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.SON,val=kc.sigma.SON$Q98_SON,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.SON,val=kc.sigma.SON$Q98_SON, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()


## Vazão média plurianual
p <- seq(0.01,0.99,l=100)
    Qm_ANO<-c()
    for(i in 1:length(kc.mu.ANO$predict)){
        Qm_ANO[i]<-integrate(function(p) Qp(kc.mu.ANO$predict[i],kc.sigma.ANO$predict[i],1,p),0,1)$value
    
    }

Qf_max <- Qm_ANO*0.8

summary(Qf_max)
length(kc.mu.ANO$predict)

summary(kc.mu.ANO$predict)
summary(kc.sigma.ANO$predict)

## Volumes regularizáveis
Vr.Qfmax <- c()
for(i in 1:length(Qf_max)){
    Vr.Qfmax[i] <- Vr(kc.mu.ANO$predict[i],kc.sigma.ANO$predict[i],1,10,Qf_max[i])
}

class(kc.sigma.ANO)
kc.sigma.ANO$Vr <- Vr.Qfmax*10
kc.sigma.ANO$Qm_ANO<- Qm_ANO*1000

Vr.Qfmax_krige <- cbind(grid*1000,kc.sigma.ANO$Vr)
write.table(Vr.Qfmax_krige,"VrQfmax.ascii",col.names = F,row.names=F,quote=F)

Qm.ANO_krige <- cbind(grid*1000,kc.sigma.ANO$Qm_ANO)
write.table(Qm.ANO_krige,"QmANO.ascii",col.names = F,row.names=F,quote=F)



postscript("Figuras/Qm_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.ANO,val=kc.sigma.ANO$Qm_ANO,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.ANO,val=kc.sigma.ANO$Qm_ANO, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)

dev.off()


postscript("Figuras/Vr_ANO.eps",onefile = T,horizontal = F, width=15/2.54, height=15/2.54,paper = "special")

image(kc.sigma.ANO,val=kc.sigma.ANO$Vr,col=brewer.pal(9,"Blues"),x.leg=c(250,450),y.leg=c(6850,6860),xlim = c(200,800),ylim = c(6725,7200), xlab = "Longitude (km)",ylab = "Latitude (km)")
contour(kc.sigma.ANO,val=kc.sigma.ANO$Vr, add=T, nlev=21)
northarrow(loc=c(300,7150),size = 20,cex =0.8)
scalebar(loc=c(250,6800),length=200)
text(400,6750,labels = "Sistema de Coordenadas Projetadas UTM \n Datum: SAD/69 Zona 22S",cex = 0.9)


dev.off()


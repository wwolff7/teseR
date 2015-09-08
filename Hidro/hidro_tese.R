# 1.0 CARREGAR LIVRARIA/PACOTES

library(hydroTSM)
library(hydroGOF)
library(vegan)
library(Kendall)
library(fitdistrplus) 
library(VGAM)
library(actuar)
library(ADGofTest)
library(EcoHydRology)
citation("fitdistrplus")
options(OutDec=",",digits = 7)

## função simpson, usada para calcular integrais numericamente
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
#Fluvi1<-read.csv("CA.csv",sep=";",dec=",", head=F)
#Fluvi.ma<-as.matrix(Fluvi1)
#Fluvi.vec<-as.vector(t(Fluvi.ma))



# ABRIR ARQUIVO/BANCO DE DADOS

dir()
setwd("~/MEGA/Doutorado/Rotinas R/Tese/Hidro/")

Peixe<-read.zoo("Qimp_peixe.csv",sep=",",dec=".", head=T)
data.peixe<-seq(as.Date("1985/01/01"),as.Date("2014/08/31"),"day")
areas.peixe<-c(420,3710,2010,801,3660,2800,1620)
PeixeDJF <- as.data.frame(extractzoo(Peixe,trgt="DJF"))
PeixeMAM <- as.data.frame(extractzoo(Peixe,trgt="MAM"))
PeixeJJA <- as.data.frame(extractzoo(Peixe,trgt="JJA"))
PeixeSON <- as.data.frame(extractzoo(Peixe,trgt="SON"))

head(Peixe)
summary(Peixe)
x_peixe <- c(1:nrow(Peixe))
x_peixeDJF <- c(1:nrow(PeixeDJF)) 
x_peixeMAM <- c(1:nrow(PeixeMAM)) 
x_peixeJJA <- c(1:nrow(PeixeJJA)) 
x_peixeSON <- c(1:nrow(PeixeSON)) 

## Vazão de base
Qb.Peixe <- lapply(as.data.frame(Peixe),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeixeDJF <- lapply(PeixeDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeixeMAM <- lapply(PeixeMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeixeJJA <- lapply(PeixeJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeixeSON <- lapply(PeixeSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(PeixeMAM)[,1],streamflow2=Qb.PeixeMAM$E72870000,timeSeries = data.peixe)

## Splines
Q_peixe_sp <- lapply(Peixe,FUN = function(x) splinefun(x_peixe,x))
Q_peixeDJF_sp <- lapply(PeixeDJF,FUN = function(x) splinefun(x_peixeDJF,x))
Q_peixeMAM_sp <- lapply(PeixeMAM,FUN = function(x) splinefun(x_peixeMAM,x))
Q_peixeJJA_sp <- lapply(PeixeJJA,FUN = function(x) splinefun(x_peixeJJA,x))
Q_peixeSON_sp <- lapply(PeixeSON,FUN = function(x) splinefun(x_peixeSON,x))
str(Q_peixe_sp)

## Volume total
VT_peixe<- c()
for(i in 1:ncol(Peixe)){
    VT_peixe[i] <- simpson(function(x) Q_peixe_sp[[i]](x),min(x_peixe),max(x_peixe),n=100000)
}
VT_peixeDJF<- c()
for(i in 1:ncol(PeixeDJF)){
    VT_peixeDJF[i] <- simpson(function(x) Q_peixeDJF_sp[[i]](x),min(x_peixeDJF),max(x_peixeDJF),n=100000)
}
VT_peixeMAM<- c()
for(i in 1:ncol(PeixeMAM)){
    VT_peixeMAM[i] <- simpson(function(x) Q_peixeMAM_sp[[i]](x),min(x_peixeMAM),max(x_peixeMAM),n=100000)
}
VT_peixeJJA<- c()
for(i in 1:ncol(PeixeJJA)){
    VT_peixeJJA[i] <- simpson(function(x) Q_peixeJJA_sp[[i]](x),min(x_peixeJJA),max(x_peixeJJA),n=100000)
}
VT_peixeSON<- c()
for(i in 1:ncol(PeixeSON)){
    VT_peixeSON[i] <- simpson(function(x) Q_peixeSON_sp[[i]](x),min(x_peixeSON),max(x_peixeSON),n=100000)
}
VT_peixe

## Spline vazão de base
Qb_peixe_sp <- lapply(Qb.Peixe,FUN = function(x) splinefun(x_peixe,x))
Qb_peixeDJF_sp <- lapply(Qb.PeixeDJF,FUN = function(x) splinefun(x_peixeDJF,x))
Qb_peixeMAM_sp <- lapply(Qb.PeixeMAM,FUN = function(x) splinefun(x_peixeMAM,x))
Qb_peixeJJA_sp <- lapply(Qb.PeixeJJA,FUN = function(x) splinefun(x_peixeJJA,x))
Qb_peixeSON_sp <- lapply(Qb.PeixeSON,FUN = function(x) splinefun(x_peixeSON,x))

str(Qb_peixe_sp)


## Volume de base
VB_peixe <- c()
for(i in 1:ncol(Peixe)){
    VB_peixe[i] <- simpson(function(x) Qb_peixe_sp[[i]](x),min(x_peixe),max(x_peixe),n=100000)
}

VB_peixeDJF <- c()
for(i in 1:ncol(PeixeDJF)){
    VB_peixeDJF[i] <- simpson(function(x) Qb_peixeDJF_sp[[i]](x),min(x_peixeDJF),max(x_peixeDJF),n=100000)
}

VB_peixeMAM <- c()
for(i in 1:ncol(PeixeMAM)){
    VB_peixeMAM[i] <- simpson(function(x) Qb_peixeMAM_sp[[i]](x),min(x_peixeMAM),max(x_peixeMAM),n=100000)
}

VB_peixeJJA <- c()
for(i in 1:ncol(PeixeJJA)){
    VB_peixeJJA[i] <- simpson(function(x) Qb_peixeJJA_sp[[i]](x),min(x_peixeJJA),max(x_peixeJJA),n=100000)
}

VB_peixeSON <- c()
for(i in 1:ncol(Peixe)){
    VB_peixeSON[i] <- simpson(function(x) Qb_peixeSON_sp[[i]](x),min(x_peixeSON),max(x_peixeSON),n=100000)
}

VB_peixe
## Indice de escoamento de base

IEB_peixe <- cbind(VB_peixe/VT_peixe,VB_peixeDJF/VT_peixeDJF,VB_peixeMAM/VT_peixeMAM,VB_peixeJJA/VT_peixeJJA,VB_peixeSON/VT_peixeSON)
colnames(IEB_peixe) <- c("ANO","DJF","MAM","JJA","SON");IEB_peixe
    
(Q2_peixe<- apply(Peixe,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_peixe<- apply(Peixe,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_peixe<- apply(Peixe,2,mean))
(Qmin_peixe<- apply(Peixe,2,min))

##-----------------------------------------------------------------------------##

Pelotas<-read.zoo("Qimp_pelotas.csv",sep=",",dec=".", head=T)
data.pelotas<-seq(as.Date("1977/01/01"),as.Date("2014/07/31"),"day")
areas.pelotas<-c(550,1170,2820,533,1820,1120)
PelotasDJF <- as.data.frame(extractzoo(Pelotas,trgt="DJF"))
PelotasMAM <- as.data.frame(extractzoo(Pelotas,trgt="MAM"))
PelotasJJA <- as.data.frame(extractzoo(Pelotas,trgt="JJA"))
PelotasSON <- as.data.frame(extractzoo(Pelotas,trgt="SON"))

head(Pelotas)
summary(Pelotas)
x_pelotas <- c(1:nrow(Pelotas))
x_pelotasDJF <- c(1:nrow(PelotasDJF)) 
x_pelotasMAM <- c(1:nrow(PelotasMAM)) 
x_pelotasJJA <- c(1:nrow(PelotasJJA)) 
x_pelotasSON <- c(1:nrow(PelotasSON)) 

## Vazão de base
Qb.Pelotas <- lapply(as.data.frame(Pelotas),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PelotasDJF <- lapply(PelotasDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PelotasMAM <- lapply(PelotasMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PelotasJJA <- lapply(PelotasJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PelotasSON <- lapply(PelotasSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(PelotasMAM)[,1],streamflow2=Qb.PelotasMAM$E72870000,timeSeries = data.pelotas)

## Splines
Q_pelotas_sp <- lapply(Pelotas,FUN = function(x) splinefun(x_pelotas,x))
Q_pelotasDJF_sp <- lapply(PelotasDJF,FUN = function(x) splinefun(x_pelotasDJF,x))
Q_pelotasMAM_sp <- lapply(PelotasMAM,FUN = function(x) splinefun(x_pelotasMAM,x))
Q_pelotasJJA_sp <- lapply(PelotasJJA,FUN = function(x) splinefun(x_pelotasJJA,x))
Q_pelotasSON_sp <- lapply(PelotasSON,FUN = function(x) splinefun(x_pelotasSON,x))
str(Q_pelotas_sp)

## Volume total
VT_pelotas<- c()
for(i in 1:ncol(Pelotas)){
    VT_pelotas[i] <- simpson(function(x) Q_pelotas_sp[[i]](x),min(x_pelotas),max(x_pelotas),n=100000)
}
VT_pelotasDJF<- c()
for(i in 1:ncol(PelotasDJF)){
    VT_pelotasDJF[i] <- simpson(function(x) Q_pelotasDJF_sp[[i]](x),min(x_pelotasDJF),max(x_pelotasDJF),n=100000)
}
VT_pelotasMAM<- c()
for(i in 1:ncol(PelotasMAM)){
    VT_pelotasMAM[i] <- simpson(function(x) Q_pelotasMAM_sp[[i]](x),min(x_pelotasMAM),max(x_pelotasMAM),n=100000)
}
VT_pelotasJJA<- c()
for(i in 1:ncol(PelotasJJA)){
    VT_pelotasJJA[i] <- simpson(function(x) Q_pelotasJJA_sp[[i]](x),min(x_pelotasJJA),max(x_pelotasJJA),n=100000)
}
VT_pelotasSON<- c()
for(i in 1:ncol(PelotasSON)){
    VT_pelotasSON[i] <- simpson(function(x) Q_pelotasSON_sp[[i]](x),min(x_pelotasSON),max(x_pelotasSON),n=100000)
}
VT_pelotas

## Spline vazão de base
Qb_pelotas_sp <- lapply(Qb.Pelotas,FUN = function(x) splinefun(x_pelotas,x))
Qb_pelotasDJF_sp <- lapply(Qb.PelotasDJF,FUN = function(x) splinefun(x_pelotasDJF,x))
Qb_pelotasMAM_sp <- lapply(Qb.PelotasMAM,FUN = function(x) splinefun(x_pelotasMAM,x))
Qb_pelotasJJA_sp <- lapply(Qb.PelotasJJA,FUN = function(x) splinefun(x_pelotasJJA,x))
Qb_pelotasSON_sp <- lapply(Qb.PelotasSON,FUN = function(x) splinefun(x_pelotasSON,x))

str(Qb_pelotas_sp)


## Volume de base
VB_pelotas <- c()
for(i in 1:ncol(Pelotas)){
    VB_pelotas[i] <- simpson(function(x) Qb_pelotas_sp[[i]](x),min(x_pelotas),max(x_pelotas),n=100000)
}

VB_pelotasDJF <- c()
for(i in 1:ncol(PelotasDJF)){
    VB_pelotasDJF[i] <- simpson(function(x) Qb_pelotasDJF_sp[[i]](x),min(x_pelotasDJF),max(x_pelotasDJF),n=100000)
}

VB_pelotasMAM <- c()
for(i in 1:ncol(PelotasMAM)){
    VB_pelotasMAM[i] <- simpson(function(x) Qb_pelotasMAM_sp[[i]](x),min(x_pelotasMAM),max(x_pelotasMAM),n=100000)
}

VB_pelotasJJA <- c()
for(i in 1:ncol(PelotasJJA)){
    VB_pelotasJJA[i] <- simpson(function(x) Qb_pelotasJJA_sp[[i]](x),min(x_pelotasJJA),max(x_pelotasJJA),n=100000)
}

VB_pelotasSON <- c()
for(i in 1:ncol(Pelotas)){
    VB_pelotasSON[i] <- simpson(function(x) Qb_pelotasSON_sp[[i]](x),min(x_pelotasSON),max(x_pelotasSON),n=100000)
}

VB_pelotas
## Indice de escoamento de base

IEB_pelotas <- cbind(VB_pelotas/VT_pelotas,VB_pelotasDJF/VT_pelotasDJF,VB_pelotasMAM/VT_pelotasMAM,VB_pelotasJJA/VT_pelotasJJA,VB_pelotasSON/VT_pelotasSON)
colnames(IEB_pelotas) <- c("ANO","DJF","MAM","JJA","SON");IEB_pelotas
    
(Q2_pelotas<- apply(Pelotas,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_pelotas<- apply(Pelotas,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_pelotas<- apply(Pelotas,2,mean))
(Qmin_pelotas<- apply(Pelotas,2,min))

##-----------------------------------------------------------------------------##

Canoas<-read.zoo("Qimp_canoas.csv",sep=",",dec=".", head=T)
data.canoas<-seq(as.Date("1986/01/01"),as.Date("2014/07/31"),"day")
areas.canoas<-c(10000,3680,4610,3230,2000,489,1010)
CanoasDJF <- as.data.frame(extractzoo(Canoas,trgt="DJF"))
CanoasMAM <- as.data.frame(extractzoo(Canoas,trgt="MAM"))
CanoasJJA <- as.data.frame(extractzoo(Canoas,trgt="JJA"))
CanoasSON <- as.data.frame(extractzoo(Canoas,trgt="SON"))

head(Canoas)
summary(Canoas)
x_canoas <- c(1:nrow(Canoas))
x_canoasDJF <- c(1:nrow(CanoasDJF)) 
x_canoasMAM <- c(1:nrow(CanoasMAM)) 
x_canoasJJA <- c(1:nrow(CanoasJJA)) 
x_canoasSON <- c(1:nrow(CanoasSON)) 

## Vazão de base
Qb.Canoas <- lapply(as.data.frame(Canoas),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.CanoasDJF <- lapply(CanoasDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.CanoasMAM <- lapply(CanoasMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.CanoasJJA <- lapply(CanoasJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.CanoasSON <- lapply(CanoasSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(CanoasMAM)[,1],streamflow2=Qb.CanoasMAM$E72870000,timeSeries = data.canoas)

## Splines
Q_canoas_sp <- lapply(Canoas,FUN = function(x) splinefun(x_canoas,x))
Q_canoasDJF_sp <- lapply(CanoasDJF,FUN = function(x) splinefun(x_canoasDJF,x))
Q_canoasMAM_sp <- lapply(CanoasMAM,FUN = function(x) splinefun(x_canoasMAM,x))
Q_canoasJJA_sp <- lapply(CanoasJJA,FUN = function(x) splinefun(x_canoasJJA,x))
Q_canoasSON_sp <- lapply(CanoasSON,FUN = function(x) splinefun(x_canoasSON,x))
str(Q_canoas_sp)

## Volume total
VT_canoas<- c()
for(i in 1:ncol(Canoas)){
    VT_canoas[i] <- simpson(function(x) Q_canoas_sp[[i]](x),min(x_canoas),max(x_canoas),n=100000)
}
VT_canoasDJF<- c()
for(i in 1:ncol(CanoasDJF)){
    VT_canoasDJF[i] <- simpson(function(x) Q_canoasDJF_sp[[i]](x),min(x_canoasDJF),max(x_canoasDJF),n=100000)
}
VT_canoasMAM<- c()
for(i in 1:ncol(CanoasMAM)){
    VT_canoasMAM[i] <- simpson(function(x) Q_canoasMAM_sp[[i]](x),min(x_canoasMAM),max(x_canoasMAM),n=100000)
}
VT_canoasJJA<- c()
for(i in 1:ncol(CanoasJJA)){
    VT_canoasJJA[i] <- simpson(function(x) Q_canoasJJA_sp[[i]](x),min(x_canoasJJA),max(x_canoasJJA),n=100000)
}
VT_canoasSON<- c()
for(i in 1:ncol(CanoasSON)){
    VT_canoasSON[i] <- simpson(function(x) Q_canoasSON_sp[[i]](x),min(x_canoasSON),max(x_canoasSON),n=100000)
}
VT_canoas

## Spline vazão de base
Qb_canoas_sp <- lapply(Qb.Canoas,FUN = function(x) splinefun(x_canoas,x))
Qb_canoasDJF_sp <- lapply(Qb.CanoasDJF,FUN = function(x) splinefun(x_canoasDJF,x))
Qb_canoasMAM_sp <- lapply(Qb.CanoasMAM,FUN = function(x) splinefun(x_canoasMAM,x))
Qb_canoasJJA_sp <- lapply(Qb.CanoasJJA,FUN = function(x) splinefun(x_canoasJJA,x))
Qb_canoasSON_sp <- lapply(Qb.CanoasSON,FUN = function(x) splinefun(x_canoasSON,x))

str(Qb_canoas_sp)


## Volume de base
VB_canoas <- c()
for(i in 1:ncol(Canoas)){
    VB_canoas[i] <- simpson(function(x) Qb_canoas_sp[[i]](x),min(x_canoas),max(x_canoas),n=100000)
}

VB_canoasDJF <- c()
for(i in 1:ncol(CanoasDJF)){
    VB_canoasDJF[i] <- simpson(function(x) Qb_canoasDJF_sp[[i]](x),min(x_canoasDJF),max(x_canoasDJF),n=100000)
}

VB_canoasMAM <- c()
for(i in 1:ncol(CanoasMAM)){
    VB_canoasMAM[i] <- simpson(function(x) Qb_canoasMAM_sp[[i]](x),min(x_canoasMAM),max(x_canoasMAM),n=100000)
}

VB_canoasJJA <- c()
for(i in 1:ncol(CanoasJJA)){
    VB_canoasJJA[i] <- simpson(function(x) Qb_canoasJJA_sp[[i]](x),min(x_canoasJJA),max(x_canoasJJA),n=100000)
}

VB_canoasSON <- c()
for(i in 1:ncol(Canoas)){
    VB_canoasSON[i] <- simpson(function(x) Qb_canoasSON_sp[[i]](x),min(x_canoasSON),max(x_canoasSON),n=100000)
}

VB_canoas
## Indice de escoamento de base

IEB_canoas <- cbind(VB_canoas/VT_canoas,VB_canoasDJF/VT_canoasDJF,VB_canoasMAM/VT_canoasMAM,VB_canoasJJA/VT_canoasJJA,VB_canoasSON/VT_canoasSON)
colnames(IEB_canoas) <- c("ANO","DJF","MAM","JJA","SON");IEB_canoas
    
(Q2_canoas<- apply(Canoas,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_canoas<- apply(Canoas,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_canoas<- apply(Canoas,2,mean))
(Qmin_canoas<- apply(Canoas,2,min))

##-----------------------------------------------------------------------------## 

Itajai<-read.zoo("Qimp_itajai.csv",sep=",",dec=".", head=T)
data.itajai<-seq(as.Date("1986/01/01"),as.Date("2014/07/31"),"day")
areas.itajai<-c(536,1430,648,717,11803,827,1240,3330,9850,1650,104,5160,286,434,1570,1600,397,9790)
ItajaiDJF <- as.data.frame(extractzoo(Itajai,trgt="DJF"))
ItajaiMAM <- as.data.frame(extractzoo(Itajai,trgt="MAM"))
ItajaiJJA <- as.data.frame(extractzoo(Itajai,trgt="JJA"))
ItajaiSON <- as.data.frame(extractzoo(Itajai,trgt="SON"))

head(Itajai)
summary(Itajai)
x_itajai <- c(1:nrow(Itajai))
x_itajaiDJF <- c(1:nrow(ItajaiDJF)) 
x_itajaiMAM <- c(1:nrow(ItajaiMAM)) 
x_itajaiJJA <- c(1:nrow(ItajaiJJA)) 
x_itajaiSON <- c(1:nrow(ItajaiSON)) 

## Vazão de base
Qb.Itajai <- lapply(as.data.frame(Itajai),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItajaiDJF <- lapply(ItajaiDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItajaiMAM <- lapply(ItajaiMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItajaiJJA <- lapply(ItajaiJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItajaiSON <- lapply(ItajaiSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(ItajaiMAM)[,1],streamflow2=Qb.ItajaiMAM$E72870000,timeSeries = data.itajai)

## Splines
Q_itajai_sp <- lapply(Itajai,FUN = function(x) splinefun(x_itajai,x))
Q_itajaiDJF_sp <- lapply(ItajaiDJF,FUN = function(x) splinefun(x_itajaiDJF,x))
Q_itajaiMAM_sp <- lapply(ItajaiMAM,FUN = function(x) splinefun(x_itajaiMAM,x))
Q_itajaiJJA_sp <- lapply(ItajaiJJA,FUN = function(x) splinefun(x_itajaiJJA,x))
Q_itajaiSON_sp <- lapply(ItajaiSON,FUN = function(x) splinefun(x_itajaiSON,x))
str(Q_itajai_sp)

## Volume total
VT_itajai<- c()
for(i in 1:ncol(Itajai)){
    VT_itajai[i] <- simpson(function(x) Q_itajai_sp[[i]](x),min(x_itajai),max(x_itajai),n=100000)
}
VT_itajaiDJF<- c()
for(i in 1:ncol(ItajaiDJF)){
    VT_itajaiDJF[i] <- simpson(function(x) Q_itajaiDJF_sp[[i]](x),min(x_itajaiDJF),max(x_itajaiDJF),n=100000)
}
VT_itajaiMAM<- c()
for(i in 1:ncol(ItajaiMAM)){
    VT_itajaiMAM[i] <- simpson(function(x) Q_itajaiMAM_sp[[i]](x),min(x_itajaiMAM),max(x_itajaiMAM),n=100000)
}
VT_itajaiJJA<- c()
for(i in 1:ncol(ItajaiJJA)){
    VT_itajaiJJA[i] <- simpson(function(x) Q_itajaiJJA_sp[[i]](x),min(x_itajaiJJA),max(x_itajaiJJA),n=100000)
}
VT_itajaiSON<- c()
for(i in 1:ncol(ItajaiSON)){
    VT_itajaiSON[i] <- simpson(function(x) Q_itajaiSON_sp[[i]](x),min(x_itajaiSON),max(x_itajaiSON),n=100000)
}
VT_itajai

## Spline vazão de base
Qb_itajai_sp <- lapply(Qb.Itajai,FUN = function(x) splinefun(x_itajai,x))
Qb_itajaiDJF_sp <- lapply(Qb.ItajaiDJF,FUN = function(x) splinefun(x_itajaiDJF,x))
Qb_itajaiMAM_sp <- lapply(Qb.ItajaiMAM,FUN = function(x) splinefun(x_itajaiMAM,x))
Qb_itajaiJJA_sp <- lapply(Qb.ItajaiJJA,FUN = function(x) splinefun(x_itajaiJJA,x))
Qb_itajaiSON_sp <- lapply(Qb.ItajaiSON,FUN = function(x) splinefun(x_itajaiSON,x))

str(Qb_itajai_sp)


## Volume de base
VB_itajai <- c()
for(i in 1:ncol(Itajai)){
    VB_itajai[i] <- simpson(function(x) Qb_itajai_sp[[i]](x),min(x_itajai),max(x_itajai),n=100000)
}

VB_itajaiDJF <- c()
for(i in 1:ncol(ItajaiDJF)){
    VB_itajaiDJF[i] <- simpson(function(x) Qb_itajaiDJF_sp[[i]](x),min(x_itajaiDJF),max(x_itajaiDJF),n=100000)
}

VB_itajaiMAM <- c()
for(i in 1:ncol(ItajaiMAM)){
    VB_itajaiMAM[i] <- simpson(function(x) Qb_itajaiMAM_sp[[i]](x),min(x_itajaiMAM),max(x_itajaiMAM),n=100000)
}

VB_itajaiJJA <- c()
for(i in 1:ncol(ItajaiJJA)){
    VB_itajaiJJA[i] <- simpson(function(x) Qb_itajaiJJA_sp[[i]](x),min(x_itajaiJJA),max(x_itajaiJJA),n=100000)
}

VB_itajaiSON <- c()
for(i in 1:ncol(Itajai)){
    VB_itajaiSON[i] <- simpson(function(x) Qb_itajaiSON_sp[[i]](x),min(x_itajaiSON),max(x_itajaiSON),n=100000)
}

VB_itajai
## Indice de escoamento de base

IEB_itajai <- cbind(VB_itajai/VT_itajai,VB_itajaiDJF/VT_itajaiDJF,VB_itajaiMAM/VT_itajaiMAM,VB_itajaiJJA/VT_itajaiJJA,VB_itajaiSON/VT_itajaiSON)
colnames(IEB_itajai) <- c("ANO","DJF","MAM","JJA","SON");IEB_itajai
    
(Q2_itajai<- apply(Itajai,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_itajai<- apply(Itajai,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_itajai<- apply(Itajai,2,mean))
(Qmin_itajai<- apply(Itajai,2,min))

##-----------------------------------------------------------------------------##

Canoinhas<-read.zoo("Qimp_canoinhas.csv",sep=",",dec=".", head=T)
data.Canoinhas<-seq(as.Date("1986/01/01"),as.Date("2014/07/31"),"day")
areas.canoinhas<-772
colnames(Canoinhas[,1]) <- "E65180000"
CanoinhasDJF <- (extractzoo(Canoinhas,trgt="DJF"))
CanoinhasMAM <- (extractzoo(Canoinhas,trgt="MAM"))
CanoinhasJJA <- (extractzoo(Canoinhas,trgt="JJA"))
CanoinhasSON <- (extractzoo(Canoinhas,trgt="SON"))

length(CanoinhasDJF)
dim(x_canoinhasDJF)
head(Canoinhas)
summary(Canoinhas)
x_canoinhas <- c(1:length(Canoinhas))
x_canoinhasDJF <- c(1:length(CanoinhasDJF)) 
x_canoinhasMAM <- c(1:length(CanoinhasMAM)) 
x_canoinhasJJA <- c(1:length(CanoinhasJJA)) 
x_canoinhasSON <- c(1:length(CanoinhasSON)) 

## Vazão e base
Qb.Canoinhas <- BaseflowSeparation(streamflow = as.data.frame(Canoinhas)[,1])[,1]
Qb.CanoinhasDJF <- BaseflowSeparation(streamflow = as.data.frame(CanoinhasDJF)[,1])[,1]
Qb.CanoinhasMAM <- BaseflowSeparation(streamflow = as.data.frame(CanoinhasMAM)[,1])[,1]
Qb.CanoinhasJJA <- BaseflowSeparation(streamflow = as.data.frame(CanoinhasJJA)[,1])[,1]
Qb.CanoinhasSON <- BaseflowSeparation(streamflow = as.data.frame(CanoinhasSON)[,1])[,1]

## Hidrografa
hydrograph(streamflow = as.data.frame(Canoinhas)[,1],streamflow2=Qb.Canoinhas,timeSeries = data.Canoinhas)

## Spline
Q_canoinhas_sp <- splinefun(x_canoinhas,as.data.frame(Canoinhas)[,1])
Q_canoinhasDJF_sp <- splinefun(x_canoinhasDJF,as.data.frame(CanoinhasDJF)[,1])
Q_canoinhasMAM_sp <- splinefun(x_canoinhasMAM,as.data.frame(CanoinhasMAM)[,1])
Q_canoinhasJJA_sp <- splinefun(x_canoinhasJJA,as.data.frame(CanoinhasJJA)[,1])
Q_canoinhasSON_sp <- splinefun(x_canoinhasSON,as.data.frame(CanoinhasSON)[,1])

str(Q_canoinhas_sp)

## Volume total
VT_canoinhas <- simpson(function(x) Q_canoinhas_sp(x),min(x_canoinhas),max(x_canoinhas),n=100000)
VT_canoinhasDJF <- simpson(function(x) Q_canoinhasDJF_sp(x),min(x_canoinhasDJF),max(x_canoinhasDJF),n=100000)
VT_canoinhasMAM <- simpson(function(x) Q_canoinhasMAM_sp(x),min(x_canoinhasMAM),max(x_canoinhasMAM),n=100000)
VT_canoinhasJJA <- simpson(function(x) Q_canoinhasJJA_sp(x),min(x_canoinhasJJA),max(x_canoinhasJJA),n=100000)
VT_canoinhasSON <- simpson(function(x) Q_canoinhasSON_sp(x),min(x_canoinhasSON),max(x_canoinhasSON),n=100000)

VT_canoinhas

## Spline de base
Qb_canoinhas_sp <- splinefun(x_canoinhas,Qb.Canoinhas)
Qb_canoinhasDJF_sp <- splinefun(x_canoinhasDJF,Qb.CanoinhasDJF)
Qb_canoinhasMAM_sp <- splinefun(x_canoinhasMAM,Qb.CanoinhasMAM)
Qb_canoinhasJJA_sp <- splinefun(x_canoinhasJJA,Qb.CanoinhasJJA)
Qb_canoinhasSON_sp <- splinefun(x_canoinhasSON,Qb.CanoinhasSON)

str(Qb_canoinhas_sp)

## Volume de base
VB_canoinhas <- simpson(function(x) Qb_canoinhas_sp(x),min(x_canoinhas),max(x_canoinhas),n=100000)
VB_canoinhasDJF <- simpson(function(x) Qb_canoinhasDJF_sp(x),min(x_canoinhasDJF),max(x_canoinhasDJF),n=100000)
VB_canoinhasMAM <- simpson(function(x) Qb_canoinhasMAM_sp(x),min(x_canoinhasMAM),max(x_canoinhasMAM),n=100000)
VB_canoinhasJJA <- simpson(function(x) Qb_canoinhasJJA_sp(x),min(x_canoinhasJJA),max(x_canoinhasJJA),n=100000)
VB_canoinhasSON <- simpson(function(x) Qb_canoinhasSON_sp(x),min(x_canoinhasSON),max(x_canoinhasSON),n=100000)


## Indice de escoamento de base
IEB_canoinhas <- cbind(VB_canoinhas/VT_canoinhas,VB_canoinhasDJF/VT_canoinhasDJF,VB_canoinhasMAM/VT_canoinhasMAM,VB_canoinhasJJA/VT_canoinhasJJA,VB_canoinhasSON/VT_canoinhasSON)
colnames(IEB_canoinhas) <- c("ANO","DJF","MAM","JJA","SON");IEB_canoinhas

(Q2_canoinhas<- quantile(as.data.frame(Canoinhas)[,1],probs = (1-0.02)))
(Q98_canoinhas<- quantile(as.data.frame(Canoinhas)[,1],probs = (1-0.98)))
(Qm_canoinhas <- mean(as.data.frame(Canoinhas)[,1]))
(Qmin_canoinhas <- min(as.data.frame(Canoinhas)[,1]))
##-----------------------------------------------------------------------------##

Antas<-read.zoo("Qimp_antas.csv",sep=",",dec=".", head=T)
data.antas<-seq(as.Date("1976/01/01"),as.Date("2014/12/31"),"day")
areas.antas<-c(609,300,61900,5340)
AntasDJF <- as.data.frame(extractzoo(Antas,trgt="DJF"))
AntasMAM <- as.data.frame(extractzoo(Antas,trgt="MAM"))
AntasJJA <- as.data.frame(extractzoo(Antas,trgt="JJA"))
AntasSON <- as.data.frame(extractzoo(Antas,trgt="SON"))

head(Antas)
summary(Antas)
x_antas <- c(1:nrow(Antas))
x_antasDJF <- c(1:nrow(AntasDJF)) 
x_antasMAM <- c(1:nrow(AntasMAM)) 
x_antasJJA <- c(1:nrow(AntasJJA)) 
x_antasSON <- c(1:nrow(AntasSON)) 

## Vazão de base
Qb.Antas <- lapply(as.data.frame(Antas),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AntasDJF <- lapply(AntasDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AntasMAM <- lapply(AntasMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AntasJJA <- lapply(AntasJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AntasSON <- lapply(AntasSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(AntasMAM)[,1],streamflow2=Qb.AntasMAM$E72870000,timeSeries = data.antas)

## Splines
Q_antas_sp <- lapply(Antas,FUN = function(x) splinefun(x_antas,x))
Q_antasDJF_sp <- lapply(AntasDJF,FUN = function(x) splinefun(x_antasDJF,x))
Q_antasMAM_sp <- lapply(AntasMAM,FUN = function(x) splinefun(x_antasMAM,x))
Q_antasJJA_sp <- lapply(AntasJJA,FUN = function(x) splinefun(x_antasJJA,x))
Q_antasSON_sp <- lapply(AntasSON,FUN = function(x) splinefun(x_antasSON,x))
str(Q_antas_sp)

## Volume total
VT_antas<- c()
for(i in 1:ncol(Antas)){
    VT_antas[i] <- simpson(function(x) Q_antas_sp[[i]](x),min(x_antas),max(x_antas),n=100000)
}
VT_antasDJF<- c()
for(i in 1:ncol(AntasDJF)){
    VT_antasDJF[i] <- simpson(function(x) Q_antasDJF_sp[[i]](x),min(x_antasDJF),max(x_antasDJF),n=100000)
}
VT_antasMAM<- c()
for(i in 1:ncol(AntasMAM)){
    VT_antasMAM[i] <- simpson(function(x) Q_antasMAM_sp[[i]](x),min(x_antasMAM),max(x_antasMAM),n=100000)
}
VT_antasJJA<- c()
for(i in 1:ncol(AntasJJA)){
    VT_antasJJA[i] <- simpson(function(x) Q_antasJJA_sp[[i]](x),min(x_antasJJA),max(x_antasJJA),n=100000)
}
VT_antasSON<- c()
for(i in 1:ncol(AntasSON)){
    VT_antasSON[i] <- simpson(function(x) Q_antasSON_sp[[i]](x),min(x_antasSON),max(x_antasSON),n=100000)
}
VT_antas

## Spline vazão de base
Qb_antas_sp <- lapply(Qb.Antas,FUN = function(x) splinefun(x_antas,x))
Qb_antasDJF_sp <- lapply(Qb.AntasDJF,FUN = function(x) splinefun(x_antasDJF,x))
Qb_antasMAM_sp <- lapply(Qb.AntasMAM,FUN = function(x) splinefun(x_antasMAM,x))
Qb_antasJJA_sp <- lapply(Qb.AntasJJA,FUN = function(x) splinefun(x_antasJJA,x))
Qb_antasSON_sp <- lapply(Qb.AntasSON,FUN = function(x) splinefun(x_antasSON,x))

str(Qb_antas_sp)


## Volume de base
VB_antas <- c()
for(i in 1:ncol(Antas)){
    VB_antas[i] <- simpson(function(x) Qb_antas_sp[[i]](x),min(x_antas),max(x_antas),n=100000)
}

VB_antasDJF <- c()
for(i in 1:ncol(AntasDJF)){
    VB_antasDJF[i] <- simpson(function(x) Qb_antasDJF_sp[[i]](x),min(x_antasDJF),max(x_antasDJF),n=100000)
}

VB_antasMAM <- c()
for(i in 1:ncol(AntasMAM)){
    VB_antasMAM[i] <- simpson(function(x) Qb_antasMAM_sp[[i]](x),min(x_antasMAM),max(x_antasMAM),n=100000)
}

VB_antasJJA <- c()
for(i in 1:ncol(AntasJJA)){
    VB_antasJJA[i] <- simpson(function(x) Qb_antasJJA_sp[[i]](x),min(x_antasJJA),max(x_antasJJA),n=100000)
}

VB_antasSON <- c()
for(i in 1:ncol(Antas)){
    VB_antasSON[i] <- simpson(function(x) Qb_antasSON_sp[[i]](x),min(x_antasSON),max(x_antasSON),n=100000)
}

VB_antas
## Indice de escoamento de base

IEB_antas <- cbind(VB_antas/VT_antas,VB_antasDJF/VT_antasDJF,VB_antasMAM/VT_antasMAM,VB_antasJJA/VT_antasJJA,VB_antasSON/VT_antasSON)
colnames(IEB_antas) <- c("ANO","DJF","MAM","JJA","SON");IEB_antas


(Q2_antas<- apply(Antas,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_antas<- apply(Antas,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_antas<- apply(Antas,2,mean))
(Qmin_antas<- apply(Antas,2,min))

##-----------------------------------------------------------------------------##

Ararangua<-read.zoo("Qimp_ararangua.csv",sep=",",dec=".", head=T)
data.ararangua<-seq(as.Date("1985/01/01"),as.Date("2004/12/31"),"day")
areas.ararangua<-c(526,355,119,359)
AraranguaDJF <- as.data.frame(extractzoo(Ararangua,trgt="DJF"))
AraranguaMAM <- as.data.frame(extractzoo(Ararangua,trgt="MAM"))
AraranguaJJA <- as.data.frame(extractzoo(Ararangua,trgt="JJA"))
AraranguaSON <- as.data.frame(extractzoo(Ararangua,trgt="SON"))

head(Ararangua)
summary(Ararangua)
x_ararangua <- c(1:nrow(Ararangua))
x_araranguaDJF <- c(1:nrow(AraranguaDJF)) 
x_araranguaMAM <- c(1:nrow(AraranguaMAM)) 
x_araranguaJJA <- c(1:nrow(AraranguaJJA)) 
x_araranguaSON <- c(1:nrow(AraranguaSON)) 

## Vazão de base
Qb.Ararangua <- lapply(as.data.frame(Ararangua),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AraranguaDJF <- lapply(AraranguaDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AraranguaMAM <- lapply(AraranguaMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AraranguaJJA <- lapply(AraranguaJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.AraranguaSON <- lapply(AraranguaSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(AraranguaMAM)[,1],streamflow2=Qb.AraranguaMAM$E72870000,timeSeries = data.ararangua)

## Splines
Q_ararangua_sp <- lapply(Ararangua,FUN = function(x) splinefun(x_ararangua,x))
Q_araranguaDJF_sp <- lapply(AraranguaDJF,FUN = function(x) splinefun(x_araranguaDJF,x))
Q_araranguaMAM_sp <- lapply(AraranguaMAM,FUN = function(x) splinefun(x_araranguaMAM,x))
Q_araranguaJJA_sp <- lapply(AraranguaJJA,FUN = function(x) splinefun(x_araranguaJJA,x))
Q_araranguaSON_sp <- lapply(AraranguaSON,FUN = function(x) splinefun(x_araranguaSON,x))
str(Q_ararangua_sp)

## Volume total
VT_ararangua<- c()
for(i in 1:ncol(Ararangua)){
    VT_ararangua[i] <- simpson(function(x) Q_ararangua_sp[[i]](x),min(x_ararangua),max(x_ararangua),n=100000)
}
VT_araranguaDJF<- c()
for(i in 1:ncol(AraranguaDJF)){
    VT_araranguaDJF[i] <- simpson(function(x) Q_araranguaDJF_sp[[i]](x),min(x_araranguaDJF),max(x_araranguaDJF),n=100000)
}
VT_araranguaMAM<- c()
for(i in 1:ncol(AraranguaMAM)){
    VT_araranguaMAM[i] <- simpson(function(x) Q_araranguaMAM_sp[[i]](x),min(x_araranguaMAM),max(x_araranguaMAM),n=100000)
}
VT_araranguaJJA<- c()
for(i in 1:ncol(AraranguaJJA)){
    VT_araranguaJJA[i] <- simpson(function(x) Q_araranguaJJA_sp[[i]](x),min(x_araranguaJJA),max(x_araranguaJJA),n=100000)
}
VT_araranguaSON<- c()
for(i in 1:ncol(AraranguaSON)){
    VT_araranguaSON[i] <- simpson(function(x) Q_araranguaSON_sp[[i]](x),min(x_araranguaSON),max(x_araranguaSON),n=100000)
}
VT_ararangua

## Spline vazão de base
Qb_ararangua_sp <- lapply(Qb.Ararangua,FUN = function(x) splinefun(x_ararangua,x))
Qb_araranguaDJF_sp <- lapply(Qb.AraranguaDJF,FUN = function(x) splinefun(x_araranguaDJF,x))
Qb_araranguaMAM_sp <- lapply(Qb.AraranguaMAM,FUN = function(x) splinefun(x_araranguaMAM,x))
Qb_araranguaJJA_sp <- lapply(Qb.AraranguaJJA,FUN = function(x) splinefun(x_araranguaJJA,x))
Qb_araranguaSON_sp <- lapply(Qb.AraranguaSON,FUN = function(x) splinefun(x_araranguaSON,x))

str(Qb_ararangua_sp)


## Volume de base
VB_ararangua <- c()
for(i in 1:ncol(Ararangua)){
    VB_ararangua[i] <- simpson(function(x) Qb_ararangua_sp[[i]](x),min(x_ararangua),max(x_ararangua),n=100000)
}

VB_araranguaDJF <- c()
for(i in 1:ncol(AraranguaDJF)){
    VB_araranguaDJF[i] <- simpson(function(x) Qb_araranguaDJF_sp[[i]](x),min(x_araranguaDJF),max(x_araranguaDJF),n=100000)
}

VB_araranguaMAM <- c()
for(i in 1:ncol(AraranguaMAM)){
    VB_araranguaMAM[i] <- simpson(function(x) Qb_araranguaMAM_sp[[i]](x),min(x_araranguaMAM),max(x_araranguaMAM),n=100000)
}

VB_araranguaJJA <- c()
for(i in 1:ncol(AraranguaJJA)){
    VB_araranguaJJA[i] <- simpson(function(x) Qb_araranguaJJA_sp[[i]](x),min(x_araranguaJJA),max(x_araranguaJJA),n=100000)
}

VB_araranguaSON <- c()
for(i in 1:ncol(Ararangua)){
    VB_araranguaSON[i] <- simpson(function(x) Qb_araranguaSON_sp[[i]](x),min(x_araranguaSON),max(x_araranguaSON),n=100000)
}

VB_ararangua
## Indice de escoamento de base

IEB_ararangua <- cbind(VB_ararangua/VT_ararangua,VB_araranguaDJF/VT_araranguaDJF,VB_araranguaMAM/VT_araranguaMAM,VB_araranguaJJA/VT_araranguaJJA,VB_araranguaSON/VT_araranguaSON)
colnames(IEB_ararangua) <- c("ANO","DJF","MAM","JJA","SON");IEB_ararangua

(Q2_ararangua<- apply(Ararangua,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_ararangua<- apply(Ararangua,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_ararangua<- apply(Ararangua,2,mean))
(Qmin_ararangua<- apply(Ararangua,2,min))

##-----------------------------------------------------------------------------##

Chapeco<-read.zoo("Qimp_chapeco.csv",sep=",",dec=".", head=T)
data.chapeco<-seq(as.Date("1979/01/01"),as.Date("2010/12/31"),"day")
areas.chapeco<-c(1660,418,5550,266,1010,642,740,8240,1840)
ChapecoDJF <- as.data.frame(extractzoo(Chapeco,trgt="DJF"))
ChapecoMAM <- as.data.frame(extractzoo(Chapeco,trgt="MAM"))
ChapecoJJA <- as.data.frame(extractzoo(Chapeco,trgt="JJA"))
ChapecoSON <- as.data.frame(extractzoo(Chapeco,trgt="SON"))

head(Chapeco)
summary(Chapeco)
x_chapeco <- c(1:nrow(Chapeco))
x_chapecoDJF <- c(1:nrow(ChapecoDJF)) 
x_chapecoMAM <- c(1:nrow(ChapecoMAM)) 
x_chapecoJJA <- c(1:nrow(ChapecoJJA)) 
x_chapecoSON <- c(1:nrow(ChapecoSON)) 

## Vazão de base
Qb.Chapeco <- lapply(as.data.frame(Chapeco),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ChapecoDJF <- lapply(ChapecoDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ChapecoMAM <- lapply(ChapecoMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ChapecoJJA <- lapply(ChapecoJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ChapecoSON <- lapply(ChapecoSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(ChapecoMAM)[,1],streamflow2=Qb.ChapecoMAM$E72870000,timeSeries = data.chapeco)

## Splines
Q_chapeco_sp <- lapply(Chapeco,FUN = function(x) splinefun(x_chapeco,x))
Q_chapecoDJF_sp <- lapply(ChapecoDJF,FUN = function(x) splinefun(x_chapecoDJF,x))
Q_chapecoMAM_sp <- lapply(ChapecoMAM,FUN = function(x) splinefun(x_chapecoMAM,x))
Q_chapecoJJA_sp <- lapply(ChapecoJJA,FUN = function(x) splinefun(x_chapecoJJA,x))
Q_chapecoSON_sp <- lapply(ChapecoSON,FUN = function(x) splinefun(x_chapecoSON,x))
str(Q_chapeco_sp)

## Volume total
VT_chapeco<- c()
for(i in 1:ncol(Chapeco)){
    VT_chapeco[i] <- simpson(function(x) Q_chapeco_sp[[i]](x),min(x_chapeco),max(x_chapeco),n=100000)
}
VT_chapecoDJF<- c()
for(i in 1:ncol(ChapecoDJF)){
    VT_chapecoDJF[i] <- simpson(function(x) Q_chapecoDJF_sp[[i]](x),min(x_chapecoDJF),max(x_chapecoDJF),n=100000)
}
VT_chapecoMAM<- c()
for(i in 1:ncol(ChapecoMAM)){
    VT_chapecoMAM[i] <- simpson(function(x) Q_chapecoMAM_sp[[i]](x),min(x_chapecoMAM),max(x_chapecoMAM),n=100000)
}
VT_chapecoJJA<- c()
for(i in 1:ncol(ChapecoJJA)){
    VT_chapecoJJA[i] <- simpson(function(x) Q_chapecoJJA_sp[[i]](x),min(x_chapecoJJA),max(x_chapecoJJA),n=100000)
}
VT_chapecoSON<- c()
for(i in 1:ncol(ChapecoSON)){
    VT_chapecoSON[i] <- simpson(function(x) Q_chapecoSON_sp[[i]](x),min(x_chapecoSON),max(x_chapecoSON),n=100000)
}
VT_chapeco

## Spline vazão de base
Qb_chapeco_sp <- lapply(Qb.Chapeco,FUN = function(x) splinefun(x_chapeco,x))
Qb_chapecoDJF_sp <- lapply(Qb.ChapecoDJF,FUN = function(x) splinefun(x_chapecoDJF,x))
Qb_chapecoMAM_sp <- lapply(Qb.ChapecoMAM,FUN = function(x) splinefun(x_chapecoMAM,x))
Qb_chapecoJJA_sp <- lapply(Qb.ChapecoJJA,FUN = function(x) splinefun(x_chapecoJJA,x))
Qb_chapecoSON_sp <- lapply(Qb.ChapecoSON,FUN = function(x) splinefun(x_chapecoSON,x))

str(Qb_chapeco_sp)


## Volume de base
VB_chapeco <- c()
for(i in 1:ncol(Chapeco)){
    VB_chapeco[i] <- simpson(function(x) Qb_chapeco_sp[[i]](x),min(x_chapeco),max(x_chapeco),n=100000)
}

VB_chapecoDJF <- c()
for(i in 1:ncol(ChapecoDJF)){
    VB_chapecoDJF[i] <- simpson(function(x) Qb_chapecoDJF_sp[[i]](x),min(x_chapecoDJF),max(x_chapecoDJF),n=100000)
}

VB_chapecoMAM <- c()
for(i in 1:ncol(ChapecoMAM)){
    VB_chapecoMAM[i] <- simpson(function(x) Qb_chapecoMAM_sp[[i]](x),min(x_chapecoMAM),max(x_chapecoMAM),n=100000)
}

VB_chapecoJJA <- c()
for(i in 1:ncol(ChapecoJJA)){
    VB_chapecoJJA[i] <- simpson(function(x) Qb_chapecoJJA_sp[[i]](x),min(x_chapecoJJA),max(x_chapecoJJA),n=100000)
}

VB_chapecoSON <- c()
for(i in 1:ncol(Chapeco)){
    VB_chapecoSON[i] <- simpson(function(x) Qb_chapecoSON_sp[[i]](x),min(x_chapecoSON),max(x_chapecoSON),n=100000)
}

VB_chapeco
## Indice de escoamento de base

IEB_chapeco <- cbind(VB_chapeco/VT_chapeco,VB_chapecoDJF/VT_chapecoDJF,VB_chapecoMAM/VT_chapecoMAM,VB_chapecoJJA/VT_chapecoJJA,VB_chapecoSON/VT_chapecoSON)
colnames(IEB_chapeco) <- c("ANO","DJF","MAM","JJA","SON");IEB_chapeco

(Q2_chapeco<- apply(Chapeco,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_chapeco<- apply(Chapeco,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_chapeco<- apply(Chapeco,2,mean))
(Qmin_chapeco<- apply(Chapeco,2,min))

##-----------------------------------------------------------------------------##

Irani<-read.zoo("Qimp_irani.csv",sep=",",dec=".", head=T)
data.irani<-seq(as.Date("1970/01/01"),as.Date("2014/12/31"),"day")
areas.irani<-c(1500,654,933)
IraniDJF <- as.data.frame(extractzoo(Irani,trgt="DJF"))
IraniMAM <- as.data.frame(extractzoo(Irani,trgt="MAM"))
IraniJJA <- as.data.frame(extractzoo(Irani,trgt="JJA"))
IraniSON <- as.data.frame(extractzoo(Irani,trgt="SON"))

head(Irani)
summary(Irani)
x_irani <- c(1:nrow(Irani))
x_iraniDJF <- c(1:nrow(IraniDJF)) 
x_iraniMAM <- c(1:nrow(IraniMAM)) 
x_iraniJJA <- c(1:nrow(IraniJJA)) 
x_iraniSON <- c(1:nrow(IraniSON)) 

## Vazão de base
Qb.Irani <- lapply(as.data.frame(Irani),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.IraniDJF <- lapply(IraniDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.IraniMAM <- lapply(IraniMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.IraniJJA <- lapply(IraniJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.IraniSON <- lapply(IraniSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(IraniMAM)[,1],streamflow2=Qb.IraniMAM$E72870000,timeSeries = data.irani)

## Splines
Q_irani_sp <- lapply(Irani,FUN = function(x) splinefun(x_irani,x))
Q_iraniDJF_sp <- lapply(IraniDJF,FUN = function(x) splinefun(x_iraniDJF,x))
Q_iraniMAM_sp <- lapply(IraniMAM,FUN = function(x) splinefun(x_iraniMAM,x))
Q_iraniJJA_sp <- lapply(IraniJJA,FUN = function(x) splinefun(x_iraniJJA,x))
Q_iraniSON_sp <- lapply(IraniSON,FUN = function(x) splinefun(x_iraniSON,x))
str(Q_irani_sp)

## Volume total
VT_irani<- c()
for(i in 1:ncol(Irani)){
    VT_irani[i] <- simpson(function(x) Q_irani_sp[[i]](x),min(x_irani),max(x_irani),n=100000)
}
VT_iraniDJF<- c()
for(i in 1:ncol(IraniDJF)){
    VT_iraniDJF[i] <- simpson(function(x) Q_iraniDJF_sp[[i]](x),min(x_iraniDJF),max(x_iraniDJF),n=100000)
}
VT_iraniMAM<- c()
for(i in 1:ncol(IraniMAM)){
    VT_iraniMAM[i] <- simpson(function(x) Q_iraniMAM_sp[[i]](x),min(x_iraniMAM),max(x_iraniMAM),n=100000)
}
VT_iraniJJA<- c()
for(i in 1:ncol(IraniJJA)){
    VT_iraniJJA[i] <- simpson(function(x) Q_iraniJJA_sp[[i]](x),min(x_iraniJJA),max(x_iraniJJA),n=100000)
}
VT_iraniSON<- c()
for(i in 1:ncol(IraniSON)){
    VT_iraniSON[i] <- simpson(function(x) Q_iraniSON_sp[[i]](x),min(x_iraniSON),max(x_iraniSON),n=100000)
}
VT_irani

## Spline vazão de base
Qb_irani_sp <- lapply(Qb.Irani,FUN = function(x) splinefun(x_irani,x))
Qb_iraniDJF_sp <- lapply(Qb.IraniDJF,FUN = function(x) splinefun(x_iraniDJF,x))
Qb_iraniMAM_sp <- lapply(Qb.IraniMAM,FUN = function(x) splinefun(x_iraniMAM,x))
Qb_iraniJJA_sp <- lapply(Qb.IraniJJA,FUN = function(x) splinefun(x_iraniJJA,x))
Qb_iraniSON_sp <- lapply(Qb.IraniSON,FUN = function(x) splinefun(x_iraniSON,x))

str(Qb_irani_sp)


## Volume de base
VB_irani <- c()
for(i in 1:ncol(Irani)){
    VB_irani[i] <- simpson(function(x) Qb_irani_sp[[i]](x),min(x_irani),max(x_irani),n=100000)
}

VB_iraniDJF <- c()
for(i in 1:ncol(IraniDJF)){
    VB_iraniDJF[i] <- simpson(function(x) Qb_iraniDJF_sp[[i]](x),min(x_iraniDJF),max(x_iraniDJF),n=100000)
}

VB_iraniMAM <- c()
for(i in 1:ncol(IraniMAM)){
    VB_iraniMAM[i] <- simpson(function(x) Qb_iraniMAM_sp[[i]](x),min(x_iraniMAM),max(x_iraniMAM),n=100000)
}

VB_iraniJJA <- c()
for(i in 1:ncol(IraniJJA)){
    VB_iraniJJA[i] <- simpson(function(x) Qb_iraniJJA_sp[[i]](x),min(x_iraniJJA),max(x_iraniJJA),n=100000)
}

VB_iraniSON <- c()
for(i in 1:ncol(Irani)){
    VB_iraniSON[i] <- simpson(function(x) Qb_iraniSON_sp[[i]](x),min(x_iraniSON),max(x_iraniSON),n=100000)
}

VB_irani
## Indice de escoamento de base

IEB_irani <- cbind(VB_irani/VT_irani,VB_iraniDJF/VT_iraniDJF,VB_iraniMAM/VT_iraniMAM,VB_iraniJJA/VT_iraniJJA,VB_iraniSON/VT_iraniSON)
colnames(IEB_irani) <- c("ANO","DJF","MAM","JJA","SON");IEB_irani


(Q2_irani<- apply(Irani,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_irani<- apply(Irani,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_irani<- apply(Irani,2,mean))
(Qmin_irani<- apply(Irani,2,min))

##-----------------------------------------------------------------------------##

Itapocu<-read.zoo("Qimp_itapocu.csv",sep=",",dec=".", head=T)
data.itapocu<-seq(as.Date("1978/01/01"),as.Date("2001/12/31"),"day")
areas.itapocu<-c(182,794,281,358,392)
ItapocuDJF <- as.data.frame(extractzoo(Itapocu,trgt="DJF"))
ItapocuMAM <- as.data.frame(extractzoo(Itapocu,trgt="MAM"))
ItapocuJJA <- as.data.frame(extractzoo(Itapocu,trgt="JJA"))
ItapocuSON <- as.data.frame(extractzoo(Itapocu,trgt="SON"))

head(Itapocu)
summary(Itapocu)
x_itapocu <- c(1:nrow(Itapocu))
x_itapocuDJF <- c(1:nrow(ItapocuDJF)) 
x_itapocuMAM <- c(1:nrow(ItapocuMAM)) 
x_itapocuJJA <- c(1:nrow(ItapocuJJA)) 
x_itapocuSON <- c(1:nrow(ItapocuSON)) 

## Vazão de base
Qb.Itapocu <- lapply(as.data.frame(Itapocu),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItapocuDJF <- lapply(ItapocuDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItapocuMAM <- lapply(ItapocuMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItapocuJJA <- lapply(ItapocuJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.ItapocuSON <- lapply(ItapocuSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(ItapocuMAM)[,1],streamflow2=Qb.ItapocuMAM$E72870000,timeSeries = data.itapocu)

## Splines
Q_itapocu_sp <- lapply(Itapocu,FUN = function(x) splinefun(x_itapocu,x))
Q_itapocuDJF_sp <- lapply(ItapocuDJF,FUN = function(x) splinefun(x_itapocuDJF,x))
Q_itapocuMAM_sp <- lapply(ItapocuMAM,FUN = function(x) splinefun(x_itapocuMAM,x))
Q_itapocuJJA_sp <- lapply(ItapocuJJA,FUN = function(x) splinefun(x_itapocuJJA,x))
Q_itapocuSON_sp <- lapply(ItapocuSON,FUN = function(x) splinefun(x_itapocuSON,x))
str(Q_itapocu_sp)

## Volume total
VT_itapocu<- c()
for(i in 1:ncol(Itapocu)){
    VT_itapocu[i] <- simpson(function(x) Q_itapocu_sp[[i]](x),min(x_itapocu),max(x_itapocu),n=100000)
}
VT_itapocuDJF<- c()
for(i in 1:ncol(ItapocuDJF)){
    VT_itapocuDJF[i] <- simpson(function(x) Q_itapocuDJF_sp[[i]](x),min(x_itapocuDJF),max(x_itapocuDJF),n=100000)
}
VT_itapocuMAM<- c()
for(i in 1:ncol(ItapocuMAM)){
    VT_itapocuMAM[i] <- simpson(function(x) Q_itapocuMAM_sp[[i]](x),min(x_itapocuMAM),max(x_itapocuMAM),n=100000)
}
VT_itapocuJJA<- c()
for(i in 1:ncol(ItapocuJJA)){
    VT_itapocuJJA[i] <- simpson(function(x) Q_itapocuJJA_sp[[i]](x),min(x_itapocuJJA),max(x_itapocuJJA),n=100000)
}
VT_itapocuSON<- c()
for(i in 1:ncol(ItapocuSON)){
    VT_itapocuSON[i] <- simpson(function(x) Q_itapocuSON_sp[[i]](x),min(x_itapocuSON),max(x_itapocuSON),n=100000)
}
VT_itapocu

## Spline vazão de base
Qb_itapocu_sp <- lapply(Qb.Itapocu,FUN = function(x) splinefun(x_itapocu,x))
Qb_itapocuDJF_sp <- lapply(Qb.ItapocuDJF,FUN = function(x) splinefun(x_itapocuDJF,x))
Qb_itapocuMAM_sp <- lapply(Qb.ItapocuMAM,FUN = function(x) splinefun(x_itapocuMAM,x))
Qb_itapocuJJA_sp <- lapply(Qb.ItapocuJJA,FUN = function(x) splinefun(x_itapocuJJA,x))
Qb_itapocuSON_sp <- lapply(Qb.ItapocuSON,FUN = function(x) splinefun(x_itapocuSON,x))

str(Qb_itapocu_sp)


## Volume de base
VB_itapocu <- c()
for(i in 1:ncol(Itapocu)){
    VB_itapocu[i] <- simpson(function(x) Qb_itapocu_sp[[i]](x),min(x_itapocu),max(x_itapocu),n=100000)
}

VB_itapocuDJF <- c()
for(i in 1:ncol(ItapocuDJF)){
    VB_itapocuDJF[i] <- simpson(function(x) Qb_itapocuDJF_sp[[i]](x),min(x_itapocuDJF),max(x_itapocuDJF),n=100000)
}

VB_itapocuMAM <- c()
for(i in 1:ncol(ItapocuMAM)){
    VB_itapocuMAM[i] <- simpson(function(x) Qb_itapocuMAM_sp[[i]](x),min(x_itapocuMAM),max(x_itapocuMAM),n=100000)
}

VB_itapocuJJA <- c()
for(i in 1:ncol(ItapocuJJA)){
    VB_itapocuJJA[i] <- simpson(function(x) Qb_itapocuJJA_sp[[i]](x),min(x_itapocuJJA),max(x_itapocuJJA),n=100000)
}

VB_itapocuSON <- c()
for(i in 1:ncol(Itapocu)){
    VB_itapocuSON[i] <- simpson(function(x) Qb_itapocuSON_sp[[i]](x),min(x_itapocuSON),max(x_itapocuSON),n=100000)
}

VB_itapocu
## Indice de escoamento de base

IEB_itapocu <- cbind(VB_itapocu/VT_itapocu,VB_itapocuDJF/VT_itapocuDJF,VB_itapocuMAM/VT_itapocuMAM,VB_itapocuJJA/VT_itapocuJJA,VB_itapocuSON/VT_itapocuSON)
colnames(IEB_itapocu) <- c("ANO","DJF","MAM","JJA","SON");IEB_itapocu


(Q2_itapocu<- apply(Itapocu,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_itapocu<- apply(Itapocu,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_itapocu<- apply(Itapocu,2,mean))
(Qmin_itapocu<- apply(Itapocu,2,min))

##-----------------------------------------------------------------------------##

Negro<-read.zoo("Qimp_negro.csv",sep=",",dec=".", head=T)
data.negro<-seq(as.Date("1977/01/01"),as.Date("2010/12/31"),"day")
areas.negro<-c(960,2610,605)
NegroDJF <- as.data.frame(extractzoo(Negro,trgt="DJF"))
NegroMAM <- as.data.frame(extractzoo(Negro,trgt="MAM"))
NegroJJA <- as.data.frame(extractzoo(Negro,trgt="JJA"))
NegroSON <- as.data.frame(extractzoo(Negro,trgt="SON"))

head(Negro)
summary(Negro)
x_negro <- c(1:nrow(Negro))
x_negroDJF <- c(1:nrow(NegroDJF)) 
x_negroMAM <- c(1:nrow(NegroMAM)) 
x_negroJJA <- c(1:nrow(NegroJJA)) 
x_negroSON <- c(1:nrow(NegroSON)) 

## Vazão de base
Qb.Negro <- lapply(as.data.frame(Negro),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.NegroDJF <- lapply(NegroDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.NegroMAM <- lapply(NegroMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.NegroJJA <- lapply(NegroJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.NegroSON <- lapply(NegroSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(NegroMAM)[,1],streamflow2=Qb.NegroMAM$E72870000,timeSeries = data.negro)

## Splines
Q_negro_sp <- lapply(Negro,FUN = function(x) splinefun(x_negro,x))
Q_negroDJF_sp <- lapply(NegroDJF,FUN = function(x) splinefun(x_negroDJF,x))
Q_negroMAM_sp <- lapply(NegroMAM,FUN = function(x) splinefun(x_negroMAM,x))
Q_negroJJA_sp <- lapply(NegroJJA,FUN = function(x) splinefun(x_negroJJA,x))
Q_negroSON_sp <- lapply(NegroSON,FUN = function(x) splinefun(x_negroSON,x))
str(Q_negro_sp)

## Volume total
VT_negro<- c()
for(i in 1:ncol(Negro)){
    VT_negro[i] <- simpson(function(x) Q_negro_sp[[i]](x),min(x_negro),max(x_negro),n=100000)
}
VT_negroDJF<- c()
for(i in 1:ncol(NegroDJF)){
    VT_negroDJF[i] <- simpson(function(x) Q_negroDJF_sp[[i]](x),min(x_negroDJF),max(x_negroDJF),n=100000)
}
VT_negroMAM<- c()
for(i in 1:ncol(NegroMAM)){
    VT_negroMAM[i] <- simpson(function(x) Q_negroMAM_sp[[i]](x),min(x_negroMAM),max(x_negroMAM),n=100000)
}
VT_negroJJA<- c()
for(i in 1:ncol(NegroJJA)){
    VT_negroJJA[i] <- simpson(function(x) Q_negroJJA_sp[[i]](x),min(x_negroJJA),max(x_negroJJA),n=100000)
}
VT_negroSON<- c()
for(i in 1:ncol(NegroSON)){
    VT_negroSON[i] <- simpson(function(x) Q_negroSON_sp[[i]](x),min(x_negroSON),max(x_negroSON),n=100000)
}
VT_negro

## Spline vazão de base
Qb_negro_sp <- lapply(Qb.Negro,FUN = function(x) splinefun(x_negro,x))
Qb_negroDJF_sp <- lapply(Qb.NegroDJF,FUN = function(x) splinefun(x_negroDJF,x))
Qb_negroMAM_sp <- lapply(Qb.NegroMAM,FUN = function(x) splinefun(x_negroMAM,x))
Qb_negroJJA_sp <- lapply(Qb.NegroJJA,FUN = function(x) splinefun(x_negroJJA,x))
Qb_negroSON_sp <- lapply(Qb.NegroSON,FUN = function(x) splinefun(x_negroSON,x))

str(Qb_negro_sp)


## Volume de base
VB_negro <- c()
for(i in 1:ncol(Negro)){
    VB_negro[i] <- simpson(function(x) Qb_negro_sp[[i]](x),min(x_negro),max(x_negro),n=100000)
}

VB_negroDJF <- c()
for(i in 1:ncol(NegroDJF)){
    VB_negroDJF[i] <- simpson(function(x) Qb_negroDJF_sp[[i]](x),min(x_negroDJF),max(x_negroDJF),n=100000)
}

VB_negroMAM <- c()
for(i in 1:ncol(NegroMAM)){
    VB_negroMAM[i] <- simpson(function(x) Qb_negroMAM_sp[[i]](x),min(x_negroMAM),max(x_negroMAM),n=100000)
}

VB_negroJJA <- c()
for(i in 1:ncol(NegroJJA)){
    VB_negroJJA[i] <- simpson(function(x) Qb_negroJJA_sp[[i]](x),min(x_negroJJA),max(x_negroJJA),n=100000)
}

VB_negroSON <- c()
for(i in 1:ncol(Negro)){
    VB_negroSON[i] <- simpson(function(x) Qb_negroSON_sp[[i]](x),min(x_negroSON),max(x_negroSON),n=100000)
}

VB_negro
## Indice de escoamento de base

IEB_negro <- cbind(VB_negro/VT_negro,VB_negroDJF/VT_negroDJF,VB_negroMAM/VT_negroMAM,VB_negroJJA/VT_negroJJA,VB_negroSON/VT_negroSON)
colnames(IEB_negro) <- c("ANO","DJF","MAM","JJA","SON");IEB_negro


(Q2_negro<- apply(Negro,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_negro<- apply(Negro,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_negro<- apply(Negro,2,mean))
(Qmin_negro<- apply(Negro,2,min))

##-----------------------------------------------------------------------------##

Peperi<-read.zoo("Qimp_peperi.csv",sep=",",dec=".", head=T)
data.peperi<-seq(as.Date("1965/01/01"),as.Date("2014/12/31"),"day")
areas.peperi<-c(2020,1540)
PeperiDJF <- as.data.frame(extractzoo(Peperi,trgt="DJF"))
PeperiMAM <- as.data.frame(extractzoo(Peperi,trgt="MAM"))
PeperiJJA <- as.data.frame(extractzoo(Peperi,trgt="JJA"))
PeperiSON <- as.data.frame(extractzoo(Peperi,trgt="SON"))

head(Peperi)
summary(Peperi)
x_peperi <- c(1:nrow(Peperi))
x_peperiDJF <- c(1:nrow(PeperiDJF)) 
x_peperiMAM <- c(1:nrow(PeperiMAM)) 
x_peperiJJA <- c(1:nrow(PeperiJJA)) 
x_peperiSON <- c(1:nrow(PeperiSON)) 

## Vazão de base
Qb.Peperi <- lapply(as.data.frame(Peperi),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeperiDJF <- lapply(PeperiDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeperiMAM <- lapply(PeperiMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeperiJJA <- lapply(PeperiJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.PeperiSON <- lapply(PeperiSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(PeperiMAM)[,1],streamflow2=Qb.PeperiMAM$E72870000,timeSeries = data.peperi)

## Splines
Q_peperi_sp <- lapply(Peperi,FUN = function(x) splinefun(x_peperi,x))
Q_peperiDJF_sp <- lapply(PeperiDJF,FUN = function(x) splinefun(x_peperiDJF,x))
Q_peperiMAM_sp <- lapply(PeperiMAM,FUN = function(x) splinefun(x_peperiMAM,x))
Q_peperiJJA_sp <- lapply(PeperiJJA,FUN = function(x) splinefun(x_peperiJJA,x))
Q_peperiSON_sp <- lapply(PeperiSON,FUN = function(x) splinefun(x_peperiSON,x))
str(Q_peperi_sp)

## Volume total
VT_peperi<- c()
for(i in 1:ncol(Peperi)){
    VT_peperi[i] <- simpson(function(x) Q_peperi_sp[[i]](x),min(x_peperi),max(x_peperi),n=100000)
}
VT_peperiDJF<- c()
for(i in 1:ncol(PeperiDJF)){
    VT_peperiDJF[i] <- simpson(function(x) Q_peperiDJF_sp[[i]](x),min(x_peperiDJF),max(x_peperiDJF),n=100000)
}
VT_peperiMAM<- c()
for(i in 1:ncol(PeperiMAM)){
    VT_peperiMAM[i] <- simpson(function(x) Q_peperiMAM_sp[[i]](x),min(x_peperiMAM),max(x_peperiMAM),n=100000)
}
VT_peperiJJA<- c()
for(i in 1:ncol(PeperiJJA)){
    VT_peperiJJA[i] <- simpson(function(x) Q_peperiJJA_sp[[i]](x),min(x_peperiJJA),max(x_peperiJJA),n=100000)
}
VT_peperiSON<- c()
for(i in 1:ncol(PeperiSON)){
    VT_peperiSON[i] <- simpson(function(x) Q_peperiSON_sp[[i]](x),min(x_peperiSON),max(x_peperiSON),n=100000)
}
VT_peperi

## Spline vazão de base
Qb_peperi_sp <- lapply(Qb.Peperi,FUN = function(x) splinefun(x_peperi,x))
Qb_peperiDJF_sp <- lapply(Qb.PeperiDJF,FUN = function(x) splinefun(x_peperiDJF,x))
Qb_peperiMAM_sp <- lapply(Qb.PeperiMAM,FUN = function(x) splinefun(x_peperiMAM,x))
Qb_peperiJJA_sp <- lapply(Qb.PeperiJJA,FUN = function(x) splinefun(x_peperiJJA,x))
Qb_peperiSON_sp <- lapply(Qb.PeperiSON,FUN = function(x) splinefun(x_peperiSON,x))

str(Qb_peperi_sp)


## Volume de base
VB_peperi <- c()
for(i in 1:ncol(Peperi)){
    VB_peperi[i] <- simpson(function(x) Qb_peperi_sp[[i]](x),min(x_peperi),max(x_peperi),n=100000)
}

VB_peperiDJF <- c()
for(i in 1:ncol(PeperiDJF)){
    VB_peperiDJF[i] <- simpson(function(x) Qb_peperiDJF_sp[[i]](x),min(x_peperiDJF),max(x_peperiDJF),n=100000)
}

VB_peperiMAM <- c()
for(i in 1:ncol(PeperiMAM)){
    VB_peperiMAM[i] <- simpson(function(x) Qb_peperiMAM_sp[[i]](x),min(x_peperiMAM),max(x_peperiMAM),n=100000)
}

VB_peperiJJA <- c()
for(i in 1:ncol(PeperiJJA)){
    VB_peperiJJA[i] <- simpson(function(x) Qb_peperiJJA_sp[[i]](x),min(x_peperiJJA),max(x_peperiJJA),n=100000)
}

VB_peperiSON <- c()
for(i in 1:ncol(Peperi)){
    VB_peperiSON[i] <- simpson(function(x) Qb_peperiSON_sp[[i]](x),min(x_peperiSON),max(x_peperiSON),n=100000)
}

VB_peperi
## Indice de escoamento de base

IEB_peperi <- cbind(VB_peperi/VT_peperi,VB_peperiDJF/VT_peperiDJF,VB_peperiMAM/VT_peperiMAM,VB_peperiJJA/VT_peperiJJA,VB_peperiSON/VT_peperiSON)
colnames(IEB_peperi) <- c("ANO","DJF","MAM","JJA","SON");IEB_peperi


(Q2_peperi<- apply(Peperi,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_peperi<- apply(Peperi,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_peperi<- apply(Peperi,2,mean))
(Qmin_peperi<- apply(Peperi,2,min))

##-----------------------------------------------------------------------------##

Tijucas<-read.zoo("Qimp_tijucas.csv",sep=",",dec=".", head=T)
data.tijucas<-seq(as.Date("1983/01/01"),as.Date("2004/12/31"),"day")
areas.tijucas<-c(1042,598,1964)
TijucasDJF <- as.data.frame(extractzoo(Tijucas,trgt="DJF"))
TijucasMAM <- as.data.frame(extractzoo(Tijucas,trgt="MAM"))
TijucasJJA <- as.data.frame(extractzoo(Tijucas,trgt="JJA"))
TijucasSON <- as.data.frame(extractzoo(Tijucas,trgt="SON"))

head(Tijucas)
summary(Tijucas)
x_tijucas <- c(1:nrow(Tijucas))
x_tijucasDJF <- c(1:nrow(TijucasDJF)) 
x_tijucasMAM <- c(1:nrow(TijucasMAM)) 
x_tijucasJJA <- c(1:nrow(TijucasJJA)) 
x_tijucasSON <- c(1:nrow(TijucasSON)) 

## Vazão de base
Qb.Tijucas <- lapply(as.data.frame(Tijucas),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TijucasDJF <- lapply(TijucasDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TijucasMAM <- lapply(TijucasMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TijucasJJA <- lapply(TijucasJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TijucasSON <- lapply(TijucasSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(TijucasMAM)[,1],streamflow2=Qb.TijucasMAM$E72870000,timeSeries = data.tijucas)

## Splines
Q_tijucas_sp <- lapply(Tijucas,FUN = function(x) splinefun(x_tijucas,x))
Q_tijucasDJF_sp <- lapply(TijucasDJF,FUN = function(x) splinefun(x_tijucasDJF,x))
Q_tijucasMAM_sp <- lapply(TijucasMAM,FUN = function(x) splinefun(x_tijucasMAM,x))
Q_tijucasJJA_sp <- lapply(TijucasJJA,FUN = function(x) splinefun(x_tijucasJJA,x))
Q_tijucasSON_sp <- lapply(TijucasSON,FUN = function(x) splinefun(x_tijucasSON,x))
str(Q_tijucas_sp)

## Volume total
VT_tijucas<- c()
for(i in 1:ncol(Tijucas)){
    VT_tijucas[i] <- simpson(function(x) Q_tijucas_sp[[i]](x),min(x_tijucas),max(x_tijucas),n=100000)
}
VT_tijucasDJF<- c()
for(i in 1:ncol(TijucasDJF)){
    VT_tijucasDJF[i] <- simpson(function(x) Q_tijucasDJF_sp[[i]](x),min(x_tijucasDJF),max(x_tijucasDJF),n=100000)
}
VT_tijucasMAM<- c()
for(i in 1:ncol(TijucasMAM)){
    VT_tijucasMAM[i] <- simpson(function(x) Q_tijucasMAM_sp[[i]](x),min(x_tijucasMAM),max(x_tijucasMAM),n=100000)
}
VT_tijucasJJA<- c()
for(i in 1:ncol(TijucasJJA)){
    VT_tijucasJJA[i] <- simpson(function(x) Q_tijucasJJA_sp[[i]](x),min(x_tijucasJJA),max(x_tijucasJJA),n=100000)
}
VT_tijucasSON<- c()
for(i in 1:ncol(TijucasSON)){
    VT_tijucasSON[i] <- simpson(function(x) Q_tijucasSON_sp[[i]](x),min(x_tijucasSON),max(x_tijucasSON),n=100000)
}
VT_tijucas

## Spline vazão de base
Qb_tijucas_sp <- lapply(Qb.Tijucas,FUN = function(x) splinefun(x_tijucas,x))
Qb_tijucasDJF_sp <- lapply(Qb.TijucasDJF,FUN = function(x) splinefun(x_tijucasDJF,x))
Qb_tijucasMAM_sp <- lapply(Qb.TijucasMAM,FUN = function(x) splinefun(x_tijucasMAM,x))
Qb_tijucasJJA_sp <- lapply(Qb.TijucasJJA,FUN = function(x) splinefun(x_tijucasJJA,x))
Qb_tijucasSON_sp <- lapply(Qb.TijucasSON,FUN = function(x) splinefun(x_tijucasSON,x))

str(Qb_tijucas_sp)


## Volume de base
VB_tijucas <- c()
for(i in 1:ncol(Tijucas)){
    VB_tijucas[i] <- simpson(function(x) Qb_tijucas_sp[[i]](x),min(x_tijucas),max(x_tijucas),n=100000)
}

VB_tijucasDJF <- c()
for(i in 1:ncol(TijucasDJF)){
    VB_tijucasDJF[i] <- simpson(function(x) Qb_tijucasDJF_sp[[i]](x),min(x_tijucasDJF),max(x_tijucasDJF),n=100000)
}

VB_tijucasMAM <- c()
for(i in 1:ncol(TijucasMAM)){
    VB_tijucasMAM[i] <- simpson(function(x) Qb_tijucasMAM_sp[[i]](x),min(x_tijucasMAM),max(x_tijucasMAM),n=100000)
}

VB_tijucasJJA <- c()
for(i in 1:ncol(TijucasJJA)){
    VB_tijucasJJA[i] <- simpson(function(x) Qb_tijucasJJA_sp[[i]](x),min(x_tijucasJJA),max(x_tijucasJJA),n=100000)
}

VB_tijucasSON <- c()
for(i in 1:ncol(Tijucas)){
    VB_tijucasSON[i] <- simpson(function(x) Qb_tijucasSON_sp[[i]](x),min(x_tijucasSON),max(x_tijucasSON),n=100000)
}

VB_tijucas
## Indice de escoamento de base

IEB_tijucas <- cbind(VB_tijucas/VT_tijucas,VB_tijucasDJF/VT_tijucasDJF,VB_tijucasMAM/VT_tijucasMAM,VB_tijucasJJA/VT_tijucasJJA,VB_tijucasSON/VT_tijucasSON)
colnames(IEB_tijucas) <- c("ANO","DJF","MAM","JJA","SON");IEB_tijucas


(Q2_tijucas<- apply(Tijucas,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_tijucas<- apply(Tijucas,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_tijucas<- apply(Tijucas,2,mean))
(Qmin_tijucas<- apply(Tijucas,2,min))

##-----------------------------------------------------------------------------##

Tubarao<-read.zoo("Qimp_tubarão.csv",sep=",",dec=".", head=T)
data.tubarao<-seq(as.Date("1987/01/01"),as.Date("2004/12/31"),"day")
areas.tubarao<-c(1515,380,599,822,2740,379,676,1690,620)
TubaraoDJF <- as.data.frame(extractzoo(Tubarao,trgt="DJF"))
TubaraoMAM <- as.data.frame(extractzoo(Tubarao,trgt="MAM"))
TubaraoJJA <- as.data.frame(extractzoo(Tubarao,trgt="JJA"))
TubaraoSON <- as.data.frame(extractzoo(Tubarao,trgt="SON"))

head(Tubarao)
summary(Tubarao)
x_tubarao <- c(1:nrow(Tubarao))
x_tubaraoDJF <- c(1:nrow(TubaraoDJF)) 
x_tubaraoMAM <- c(1:nrow(TubaraoMAM)) 
x_tubaraoJJA <- c(1:nrow(TubaraoJJA)) 
x_tubaraoSON <- c(1:nrow(TubaraoSON)) 

## Vazão de base
Qb.Tubarao <- lapply(as.data.frame(Tubarao),FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TubaraoDJF <- lapply(TubaraoDJF,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TubaraoMAM <- lapply(TubaraoMAM,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TubaraoJJA <- lapply(TubaraoJJA,FUN = function(x) BaseflowSeparation(x)[,1])
Qb.TubaraoSON <- lapply(TubaraoSON,FUN = function(x) BaseflowSeparation(x)[,1])

## Hidrografa
hydrograph(streamflow = as.data.frame(TubaraoMAM)[,1],streamflow2=Qb.TubaraoMAM$E72870000,timeSeries = data.tubarao)

## Splines
Q_tubarao_sp <- lapply(Tubarao,FUN = function(x) splinefun(x_tubarao,x))
Q_tubaraoDJF_sp <- lapply(TubaraoDJF,FUN = function(x) splinefun(x_tubaraoDJF,x))
Q_tubaraoMAM_sp <- lapply(TubaraoMAM,FUN = function(x) splinefun(x_tubaraoMAM,x))
Q_tubaraoJJA_sp <- lapply(TubaraoJJA,FUN = function(x) splinefun(x_tubaraoJJA,x))
Q_tubaraoSON_sp <- lapply(TubaraoSON,FUN = function(x) splinefun(x_tubaraoSON,x))
str(Q_tubarao_sp)

## Volume total
VT_tubarao<- c()
for(i in 1:ncol(Tubarao)){
    VT_tubarao[i] <- simpson(function(x) Q_tubarao_sp[[i]](x),min(x_tubarao),max(x_tubarao),n=100000)
}
VT_tubaraoDJF<- c()
for(i in 1:ncol(TubaraoDJF)){
    VT_tubaraoDJF[i] <- simpson(function(x) Q_tubaraoDJF_sp[[i]](x),min(x_tubaraoDJF),max(x_tubaraoDJF),n=100000)
}
VT_tubaraoMAM<- c()
for(i in 1:ncol(TubaraoMAM)){
    VT_tubaraoMAM[i] <- simpson(function(x) Q_tubaraoMAM_sp[[i]](x),min(x_tubaraoMAM),max(x_tubaraoMAM),n=100000)
}
VT_tubaraoJJA<- c()
for(i in 1:ncol(TubaraoJJA)){
    VT_tubaraoJJA[i] <- simpson(function(x) Q_tubaraoJJA_sp[[i]](x),min(x_tubaraoJJA),max(x_tubaraoJJA),n=100000)
}
VT_tubaraoSON<- c()
for(i in 1:ncol(TubaraoSON)){
    VT_tubaraoSON[i] <- simpson(function(x) Q_tubaraoSON_sp[[i]](x),min(x_tubaraoSON),max(x_tubaraoSON),n=100000)
}
VT_tubarao

## Spline vazão de base
Qb_tubarao_sp <- lapply(Qb.Tubarao,FUN = function(x) splinefun(x_tubarao,x))
Qb_tubaraoDJF_sp <- lapply(Qb.TubaraoDJF,FUN = function(x) splinefun(x_tubaraoDJF,x))
Qb_tubaraoMAM_sp <- lapply(Qb.TubaraoMAM,FUN = function(x) splinefun(x_tubaraoMAM,x))
Qb_tubaraoJJA_sp <- lapply(Qb.TubaraoJJA,FUN = function(x) splinefun(x_tubaraoJJA,x))
Qb_tubaraoSON_sp <- lapply(Qb.TubaraoSON,FUN = function(x) splinefun(x_tubaraoSON,x))

str(Qb_tubarao_sp)


## Volume de base
VB_tubarao <- c()
for(i in 1:ncol(Tubarao)){
    VB_tubarao[i] <- simpson(function(x) Qb_tubarao_sp[[i]](x),min(x_tubarao),max(x_tubarao),n=100000)
}

VB_tubaraoDJF <- c()
for(i in 1:ncol(TubaraoDJF)){
    VB_tubaraoDJF[i] <- simpson(function(x) Qb_tubaraoDJF_sp[[i]](x),min(x_tubaraoDJF),max(x_tubaraoDJF),n=100000)
}

VB_tubaraoMAM <- c()
for(i in 1:ncol(TubaraoMAM)){
    VB_tubaraoMAM[i] <- simpson(function(x) Qb_tubaraoMAM_sp[[i]](x),min(x_tubaraoMAM),max(x_tubaraoMAM),n=100000)
}

VB_tubaraoJJA <- c()
for(i in 1:ncol(TubaraoJJA)){
    VB_tubaraoJJA[i] <- simpson(function(x) Qb_tubaraoJJA_sp[[i]](x),min(x_tubaraoJJA),max(x_tubaraoJJA),n=100000)
}

VB_tubaraoSON <- c()
for(i in 1:ncol(Tubarao)){
    VB_tubaraoSON[i] <- simpson(function(x) Qb_tubaraoSON_sp[[i]](x),min(x_tubaraoSON),max(x_tubaraoSON),n=100000)
}

VB_tubarao
## Indice de escoamento de base

IEB_tubarao <- cbind(VB_tubarao/VT_tubarao,VB_tubaraoDJF/VT_tubaraoDJF,VB_tubaraoMAM/VT_tubaraoMAM,VB_tubaraoJJA/VT_tubaraoJJA,VB_tubaraoSON/VT_tubaraoSON)
colnames(IEB_tubarao) <- c("ANO","DJF","MAM","JJA","SON");IEB_tubarao


(Q2_tubarao<- apply(Tubarao,2,function(x) quantile(x,probs = (1-0.02))))
(Q98_tubarao<- apply(Tubarao,2,function(x) quantile(x,probs = (1-0.98))))
(Qm_tubarao<- apply(Tubarao,2,mean))
(Qmin_tubarao<- apply(Tubarao,2,min))

##-----------------------------------------------------------------------------##

Mampituba<-read.zoo("Qimp_mampituba.csv",sep=",",dec=".", head=T)
data.Mampituba<-seq(as.Date("1986/01/01"),as.Date("2006/12/31"),"day")
areas.mampituba<-339
colnames(Mampituba[,1]) <- "E84970000"
MampitubaDJF <- (extractzoo(Mampituba,trgt="DJF"))
MampitubaMAM <- (extractzoo(Mampituba,trgt="MAM"))
MampitubaJJA <- (extractzoo(Mampituba,trgt="JJA"))
MampitubaSON <- (extractzoo(Mampituba,trgt="SON"))

length(MampitubaDJF)
dim(x_mampitubaDJF)
head(Mampituba)
summary(Mampituba)
x_mampituba <- c(1:length(Mampituba))
x_mampitubaDJF <- c(1:length(MampitubaDJF)) 
x_mampitubaMAM <- c(1:length(MampitubaMAM)) 
x_mampitubaJJA <- c(1:length(MampitubaJJA)) 
x_mampitubaSON <- c(1:length(MampitubaSON)) 

## Vazão e base
Qb.Mampituba <- BaseflowSeparation(streamflow = as.data.frame(Mampituba)[,1])[,1]
Qb.MampitubaDJF <- BaseflowSeparation(streamflow = as.data.frame(MampitubaDJF)[,1])[,1]
Qb.MampitubaMAM <- BaseflowSeparation(streamflow = as.data.frame(MampitubaMAM)[,1])[,1]
Qb.MampitubaJJA <- BaseflowSeparation(streamflow = as.data.frame(MampitubaJJA)[,1])[,1]
Qb.MampitubaSON <- BaseflowSeparation(streamflow = as.data.frame(MampitubaSON)[,1])[,1]

## Hidrografa
hydrograph(streamflow = as.data.frame(Mampituba)[,1],streamflow2=Qb.Mampituba,timeSeries = data.Mampituba)

## Spline
Q_mampituba_sp <- splinefun(x_mampituba,as.data.frame(Mampituba)[,1])
Q_mampitubaDJF_sp <- splinefun(x_mampitubaDJF,as.data.frame(MampitubaDJF)[,1])
Q_mampitubaMAM_sp <- splinefun(x_mampitubaMAM,as.data.frame(MampitubaMAM)[,1])
Q_mampitubaJJA_sp <- splinefun(x_mampitubaJJA,as.data.frame(MampitubaJJA)[,1])
Q_mampitubaSON_sp <- splinefun(x_mampitubaSON,as.data.frame(MampitubaSON)[,1])

str(Q_mampituba_sp)

## Volume total
VT_mampituba <- simpson(function(x) Q_mampituba_sp(x),min(x_mampituba),max(x_mampituba),n=100000)
VT_mampitubaDJF <- simpson(function(x) Q_mampitubaDJF_sp(x),min(x_mampitubaDJF),max(x_mampitubaDJF),n=100000)
VT_mampitubaMAM <- simpson(function(x) Q_mampitubaMAM_sp(x),min(x_mampitubaMAM),max(x_mampitubaMAM),n=100000)
VT_mampitubaJJA <- simpson(function(x) Q_mampitubaJJA_sp(x),min(x_mampitubaJJA),max(x_mampitubaJJA),n=100000)
VT_mampitubaSON <- simpson(function(x) Q_mampitubaSON_sp(x),min(x_mampitubaSON),max(x_mampitubaSON),n=100000)

VT_mampituba

## Spline de base
Qb_mampituba_sp <- splinefun(x_mampituba,Qb.Mampituba)
Qb_mampitubaDJF_sp <- splinefun(x_mampitubaDJF,Qb.MampitubaDJF)
Qb_mampitubaMAM_sp <- splinefun(x_mampitubaMAM,Qb.MampitubaMAM)
Qb_mampitubaJJA_sp <- splinefun(x_mampitubaJJA,Qb.MampitubaJJA)
Qb_mampitubaSON_sp <- splinefun(x_mampitubaSON,Qb.MampitubaSON)

str(Qb_mampituba_sp)

## Volume de base
VB_mampituba <- simpson(function(x) Qb_mampituba_sp(x),min(x_mampituba),max(x_mampituba),n=100000)
VB_mampitubaDJF <- simpson(function(x) Qb_mampitubaDJF_sp(x),min(x_mampitubaDJF),max(x_mampitubaDJF),n=100000)
VB_mampitubaMAM <- simpson(function(x) Qb_mampitubaMAM_sp(x),min(x_mampitubaMAM),max(x_mampitubaMAM),n=100000)
VB_mampitubaJJA <- simpson(function(x) Qb_mampitubaJJA_sp(x),min(x_mampitubaJJA),max(x_mampitubaJJA),n=100000)
VB_mampitubaSON <- simpson(function(x) Qb_mampitubaSON_sp(x),min(x_mampitubaSON),max(x_mampitubaSON),n=100000)


## Indice de escoamento de base
IEB_mampituba <- cbind(VB_mampituba/VT_mampituba,VB_mampitubaDJF/VT_mampitubaDJF,VB_mampitubaMAM/VT_mampitubaMAM,VB_mampitubaJJA/VT_mampitubaJJA,VB_mampitubaSON/VT_mampitubaSON)
colnames(IEB_mampituba) <- c("ANO","DJF","MAM","JJA","SON");IEB_mampituba


(Q2_mampituba<- quantile(as.data.frame(Mampituba)[,1],probs = (1-0.02)))
(Q98_mampituba<- quantile(as.data.frame(Mampituba)[,1],probs = (1-0.98)))
(Qm_mampituba <- mean(as.data.frame(Mampituba)[,1]))
(Qmin_mampituba <- min(as.data.frame(Mampituba)[,1]))
##-----------------------------------------------------------------------------##


CubataoS<-read.zoo("Qimp_cubataoS.csv",sep=",",dec=".", head=T)
data.CubataoS<-seq(as.Date("1951/01/01"),as.Date("2000/12/31"),"day")
areas.cubataoS<-400
colnames(CubataoS[,1]) <- "E84100000"
CubataoSDJF <- (extractzoo(CubataoS,trgt="DJF"))
CubataoSMAM <- (extractzoo(CubataoS,trgt="MAM"))
CubataoSJJA <- (extractzoo(CubataoS,trgt="JJA"))
CubataoSSON <- (extractzoo(CubataoS,trgt="SON"))

length(CubataoSDJF)
dim(x_cubataoSDJF)
head(CubataoS)
summary(CubataoS)
x_cubataoS <- c(1:length(CubataoS))
x_cubataoSDJF <- c(1:length(CubataoSDJF)) 
x_cubataoSMAM <- c(1:length(CubataoSMAM)) 
x_cubataoSJJA <- c(1:length(CubataoSJJA)) 
x_cubataoSSON <- c(1:length(CubataoSSON)) 

## Vazão e base
Qb.CubataoS <- BaseflowSeparation(streamflow = as.data.frame(CubataoS)[,1])[,1]
Qb.CubataoSDJF <- BaseflowSeparation(streamflow = as.data.frame(CubataoSDJF)[,1])[,1]
Qb.CubataoSMAM <- BaseflowSeparation(streamflow = as.data.frame(CubataoSMAM)[,1])[,1]
Qb.CubataoSJJA <- BaseflowSeparation(streamflow = as.data.frame(CubataoSJJA)[,1])[,1]
Qb.CubataoSSON <- BaseflowSeparation(streamflow = as.data.frame(CubataoSSON)[,1])[,1]

## Hidrografa
hydrograph(streamflow = as.data.frame(CubataoS)[,1],streamflow2=Qb.CubataoS,timeSeries = data.CubataoS)

## Spline
Q_cubataoS_sp <- splinefun(x_cubataoS,as.data.frame(CubataoS)[,1])
Q_cubataoSDJF_sp <- splinefun(x_cubataoSDJF,as.data.frame(CubataoSDJF)[,1])
Q_cubataoSMAM_sp <- splinefun(x_cubataoSMAM,as.data.frame(CubataoSMAM)[,1])
Q_cubataoSJJA_sp <- splinefun(x_cubataoSJJA,as.data.frame(CubataoSJJA)[,1])
Q_cubataoSSON_sp <- splinefun(x_cubataoSSON,as.data.frame(CubataoSSON)[,1])

str(Q_cubataoS_sp)

## Volume total
VT_cubataoS <- simpson(function(x) Q_cubataoS_sp(x),min(x_cubataoS),max(x_cubataoS),n=100000)
VT_cubataoSDJF <- simpson(function(x) Q_cubataoSDJF_sp(x),min(x_cubataoSDJF),max(x_cubataoSDJF),n=100000)
VT_cubataoSMAM <- simpson(function(x) Q_cubataoSMAM_sp(x),min(x_cubataoSMAM),max(x_cubataoSMAM),n=100000)
VT_cubataoSJJA <- simpson(function(x) Q_cubataoSJJA_sp(x),min(x_cubataoSJJA),max(x_cubataoSJJA),n=100000)
VT_cubataoSSON <- simpson(function(x) Q_cubataoSSON_sp(x),min(x_cubataoSSON),max(x_cubataoSSON),n=100000)

VT_cubataoS

## Spline de base
Qb_cubataoS_sp <- splinefun(x_cubataoS,Qb.CubataoS)
Qb_cubataoSDJF_sp <- splinefun(x_cubataoSDJF,Qb.CubataoSDJF)
Qb_cubataoSMAM_sp <- splinefun(x_cubataoSMAM,Qb.CubataoSMAM)
Qb_cubataoSJJA_sp <- splinefun(x_cubataoSJJA,Qb.CubataoSJJA)
Qb_cubataoSSON_sp <- splinefun(x_cubataoSSON,Qb.CubataoSSON)

str(Qb_cubataoS_sp)

## Volume de base
VB_cubataoS <- simpson(function(x) Qb_cubataoS_sp(x),min(x_cubataoS),max(x_cubataoS),n=100000)
VB_cubataoSDJF <- simpson(function(x) Qb_cubataoSDJF_sp(x),min(x_cubataoSDJF),max(x_cubataoSDJF),n=100000)
VB_cubataoSMAM <- simpson(function(x) Qb_cubataoSMAM_sp(x),min(x_cubataoSMAM),max(x_cubataoSMAM),n=100000)
VB_cubataoSJJA <- simpson(function(x) Qb_cubataoSJJA_sp(x),min(x_cubataoSJJA),max(x_cubataoSJJA),n=100000)
VB_cubataoSSON <- simpson(function(x) Qb_cubataoSSON_sp(x),min(x_cubataoSSON),max(x_cubataoSSON),n=100000)


## Indice de escoamento de base
IEB_cubataoS <- cbind(VB_cubataoS/VT_cubataoS,VB_cubataoSDJF/VT_cubataoSDJF,VB_cubataoSMAM/VT_cubataoSMAM,VB_cubataoSJJA/VT_cubataoSJJA,VB_cubataoSSON/VT_cubataoSSON)
colnames(IEB_cubataoS) <- c("ANO","DJF","MAM","JJA","SON");IEB_cubataoS


(Q2_cubataoS<- quantile(as.data.frame(CubataoS)[,1],probs = (1-0.02)))
(Q98_cubataoS<- quantile(as.data.frame(CubataoS)[,1],probs = (1-0.98)))
(Qm_cubataoS <- mean(as.data.frame(CubataoS)[,1]))
(Qmin_cubataoS <- min(as.data.frame(CubataoS)[,1]))

##-----------------------------------------------------------------------------##

CubataoN<-read.zoo("Qimp_cubataoN.csv",sep=",",dec=".", head=T)
data.CubataoN<-seq(as.Date("1986/01/01"),as.Date("2010/12/31"),"day")
areas.cubataoN<-374
colnames(CubataoN[,1]) <- "E82270050"
CubataoNDJF <- (extractzoo(CubataoN,trgt="DJF"))
CubataoNMAM <- (extractzoo(CubataoN,trgt="MAM"))
CubataoNJJA <- (extractzoo(CubataoN,trgt="JJA"))
CubataoNSON <- (extractzoo(CubataoN,trgt="SON"))

length(CubataoNDJF)
dim(x_cubataoNDJF)
head(CubataoN)
summary(CubataoN)
x_cubataoN <- c(1:length(CubataoN))
x_cubataoNDJF <- c(1:length(CubataoNDJF)) 
x_cubataoNMAM <- c(1:length(CubataoNMAM)) 
x_cubataoNJJA <- c(1:length(CubataoNJJA)) 
x_cubataoNSON <- c(1:length(CubataoNSON)) 

## Vazão e base
Qb.CubataoN <- BaseflowSeparation(streamflow = as.data.frame(CubataoN)[,1])[,1]
Qb.CubataoNDJF <- BaseflowSeparation(streamflow = as.data.frame(CubataoNDJF)[,1])[,1]
Qb.CubataoNMAM <- BaseflowSeparation(streamflow = as.data.frame(CubataoNMAM)[,1])[,1]
Qb.CubataoNJJA <- BaseflowSeparation(streamflow = as.data.frame(CubataoNJJA)[,1])[,1]
Qb.CubataoNSON <- BaseflowSeparation(streamflow = as.data.frame(CubataoNSON)[,1])[,1]

## Hidrografa
hydrograph(streamflow = as.data.frame(CubataoN)[,1],streamflow2=Qb.CubataoN,timeSeries = data.CubataoN)

## Spline
Q_cubataoN_sp <- splinefun(x_cubataoN,as.data.frame(CubataoN)[,1])
Q_cubataoNDJF_sp <- splinefun(x_cubataoNDJF,as.data.frame(CubataoNDJF)[,1])
Q_cubataoNMAM_sp <- splinefun(x_cubataoNMAM,as.data.frame(CubataoNMAM)[,1])
Q_cubataoNJJA_sp <- splinefun(x_cubataoNJJA,as.data.frame(CubataoNJJA)[,1])
Q_cubataoNSON_sp <- splinefun(x_cubataoNSON,as.data.frame(CubataoNSON)[,1])

str(Q_cubataoN_sp)

## Volume total
VT_cubataoN <- simpson(function(x) Q_cubataoN_sp(x),min(x_cubataoN),max(x_cubataoN),n=100000)
VT_cubataoNDJF <- simpson(function(x) Q_cubataoNDJF_sp(x),min(x_cubataoNDJF),max(x_cubataoNDJF),n=100000)
VT_cubataoNMAM <- simpson(function(x) Q_cubataoNMAM_sp(x),min(x_cubataoNMAM),max(x_cubataoNMAM),n=100000)
VT_cubataoNJJA <- simpson(function(x) Q_cubataoNJJA_sp(x),min(x_cubataoNJJA),max(x_cubataoNJJA),n=100000)
VT_cubataoNSON <- simpson(function(x) Q_cubataoNSON_sp(x),min(x_cubataoNSON),max(x_cubataoNSON),n=100000)

VT_cubataoN

## Spline de base
Qb_cubataoN_sp <- splinefun(x_cubataoN,Qb.CubataoN)
Qb_cubataoNDJF_sp <- splinefun(x_cubataoNDJF,Qb.CubataoNDJF)
Qb_cubataoNMAM_sp <- splinefun(x_cubataoNMAM,Qb.CubataoNMAM)
Qb_cubataoNJJA_sp <- splinefun(x_cubataoNJJA,Qb.CubataoNJJA)
Qb_cubataoNSON_sp <- splinefun(x_cubataoNSON,Qb.CubataoNSON)

str(Qb_cubataoN_sp)

## Volume de base
VB_cubataoN <- simpson(function(x) Qb_cubataoN_sp(x),min(x_cubataoN),max(x_cubataoN),n=100000)
VB_cubataoNDJF <- simpson(function(x) Qb_cubataoNDJF_sp(x),min(x_cubataoNDJF),max(x_cubataoNDJF),n=100000)
VB_cubataoNMAM <- simpson(function(x) Qb_cubataoNMAM_sp(x),min(x_cubataoNMAM),max(x_cubataoNMAM),n=100000)
VB_cubataoNJJA <- simpson(function(x) Qb_cubataoNJJA_sp(x),min(x_cubataoNJJA),max(x_cubataoNJJA),n=100000)
VB_cubataoNSON <- simpson(function(x) Qb_cubataoNSON_sp(x),min(x_cubataoNSON),max(x_cubataoNSON),n=100000)


## Indice de escoamento de base
IEB_cubataoN <- cbind(VB_cubataoN/VT_cubataoN,VB_cubataoNDJF/VT_cubataoNDJF,VB_cubataoNMAM/VT_cubataoNMAM,VB_cubataoNJJA/VT_cubataoNJJA,VB_cubataoNSON/VT_cubataoNSON)
colnames(IEB_cubataoN) <- c("ANO","DJF","MAM","JJA","SON");IEB_cubataoN


(Q2_cubataoN<- quantile(as.data.frame(CubataoN)[,1],probs = (1-0.02)))
(Q98_cubataoN<- quantile(as.data.frame(CubataoN)[,1],probs = (1-0.98)))
(Qm_cubataoN <- mean(as.data.frame(CubataoN)[,1]))
(Qmin_cubataoN <- min(as.data.frame(CubataoN)[,1]))
##-----------------------------------------------------------------------------##
IEB
IEB <- rbind(IEB_peixe,IEB_pelotas,IEB_canoas,IEB_itajai,IEB_antas,IEB_ararangua,IEB_chapeco,IEB_irani,IEB_itapocu,IEB_negro,IEB_peperi,IEB_tijucas,IEB_tubarao,IEB_canoinhas,IEB_cubataoN,IEB_cubataoS,IEB_mampituba)

nomes <- c(names(Peixe),names(Pelotas),names(Canoas),names(Itajai),names(Antas),names(Ararangua),names(Chapeco),names(Irani),names(Itapocu),names(Negro),names(Peperi),names(Tijucas),names(Tubarao),"E65180000","E82270050","E84100000","E84970000")


Q98 <- c(as.vector(Q98_peixe),as.vector(Q98_pelotas),as.vector(Q98_canoas),as.vector(Q98_itajai),as.vector(Q98_antas),as.vector(Q98_ararangua),as.vector(Q98_chapeco),as.vector(Q98_irani),as.vector(Q98_itapocu),as.vector(Q98_negro),as.vector(Q98_peperi),as.vector(Q98_tijucas),as.vector(Q98_tubarao),as.vector(Q98_canoinhas),as.vector(Q98_cubataoN),as.vector(Q98_cubataoS),as.vector(Q98_mampituba))

Q2 <- c(as.vector(Q2_peixe),as.vector(Q2_pelotas),as.vector(Q2_canoas),as.vector(Q2_itajai),as.vector(Q2_antas),as.vector(Q2_ararangua),as.vector(Q2_chapeco),as.vector(Q2_irani),as.vector(Q2_itapocu),as.vector(Q2_negro),as.vector(Q2_peperi),as.vector(Q2_tijucas),as.vector(Q2_tubarao),as.vector(Q2_canoinhas),as.vector(Q2_cubataoN),as.vector(Q2_cubataoS),as.vector(Q2_mampituba))

Qmed <- c(as.vector(Qm_peixe),as.vector(Qm_pelotas),as.vector(Qm_canoas),as.vector(Qm_itajai),as.vector(Qm_antas),as.vector(Qm_ararangua),as.vector(Qm_chapeco),as.vector(Qm_irani),as.vector(Qm_itapocu),as.vector(Qm_negro),as.vector(Qm_peperi),as.vector(Qm_tijucas),as.vector(Qm_tubarao),as.vector(Qm_canoinhas),as.vector(Qm_cubataoN),as.vector(Qm_cubataoS),as.vector(Qm_mampituba))

Qmin <- c(as.vector(Qmin_peixe),as.vector(Qmin_pelotas),as.vector(Qmin_canoas),as.vector(Qmin_itajai),as.vector(Qmin_antas),as.vector(Qmin_ararangua),as.vector(Qmin_chapeco),as.vector(Qmin_irani),as.vector(Qmin_itapocu),as.vector(Qmin_negro),as.vector(Qmin_peperi),as.vector(Qmin_tijucas),as.vector(Qmin_tubarao),as.vector(Qmin_canoinhas),as.vector(Qmin_cubataoN),as.vector(Qmin_cubataoS),as.vector(Qmin_mampituba))

areas <- c(areas.peixe,areas.pelotas,areas.canoas,areas.itajai,areas.antas,areas.ararangua,areas.chapeco,areas.irani,areas.itapocu,areas.negro,areas.peperi,areas.tijucas,areas.tubarao,areas.canoinhas,areas.cubataoN,areas.cubataoS,areas.mampituba)

IEB
(tab <- data.frame(nomes,IEB,Q98,Q2,Qmed,Qmin,areas))
dim(tab)
write.csv(tab[,c(1:6,11)],"Param/IEB.csv")
(tab1 <- cbind(tab[2:7],Param.ANO))

## Gráfico para demonstrar nenhuma relação das vazões específicas com a área
par(mfrow = c(2,3))
for(i in 1:5){
plot(tab$areas,tab[,i+1])
}

cor(tab1)

plot(tab[,2:6])
cor(tab[,2:7])
hydropairs(tab[,2:7],method = )
summary(Pelotas)
IEB_pelotas


##-----------------------------------------------------------------------------##

matrixplot(dwi(Itajai, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio do Itajaí-Açú")
matrixplot(dwi(Canoas, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio Canoas")
matrixplot(dwi(Canoinhas, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio Canoinhas")
matrixplot(dwi(Peixe, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio do Peixe")
matrixplot(dwi(Pelotas, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio Pelotas")

## Vazão específica média vs área de drenagem


## 3.0 ANÁLISE EXPLORATÓRIA  

plot(Peixe,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Pelotas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Canoas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Canoinhas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Itajai,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)


##-----------------------------------------------------------------------------##
## Curvas de Permanência
## Rodar quanto for escrever

par(mfrow = c(3,2))
fdc(Peixe, lQ.thr="",hQ.thr="", plot=T,thr.shw=F,log="y",leg.pos="bottomleft", main= "", xlab="% Permanência",mgp=c(2.3,1,0), ylab=expression(Q~(m^3.~s^{-1})))
fdc(Canoas, lQ.thr="",hQ.thr="", plot=T,thr.shw=F,log="y",leg.pos="bottomleft", main= "", xlab="% Permanência",mgp=c(2.3,1,0), ylab=expression(Q~(m^3.~s^{-1})))
fdc(Pelotas, lQ.thr="",hQ.thr="", plot=T,thr.shw=F,log="y",leg.pos="bottomleft", main= "", xlab="% Permanência",mgp=c(2.3,1,0), ylab=expression(Q~(m^3.~s^{-1})))
fdc(Itajai, lQ.thr="",hQ.thr="", plot=T,thr.shw=F,log="y",leg.pos="bottomleft", main= "", xlab="% Permanência",mgp=c(2.3,1,0), ylab=expression(Q~(m^3.~s^{-1})))
fdc(Canoinhas, lQ.thr="",hQ.thr="", plot=T,thr.shw=F,log="y",leg.pos="bottomleft", main= "", xlab="% Permanência",mgp=c(2.3,1,0), ylab=expression(Q~(m^3.~s^{-1})))


##-----------------------------------------------------------------------------##
## Ajuste das Qesp para distribuição log-normal
##-----------------------------------------------------------------------------##


Qparam <- function(Estac){
fit <- list()
for(i in 1:ncol(Estac)){
fit$E[i] <- fitdist(Estac[,i], "lnorm")}
x <- do.call(cbind.data.frame, fit)
names(x) <- colnames(Estac)
return(t(x))
}


##-----------------------------------------------------------------------------##
## Anual
Qesp_itajai <- as.data.frame(Itajai)
names(Qesp_itajai)
summary(Qesp_itajai)
dim(Itajai)

Qesp_peixe <- as.data.frame(Peixe)
names(Qesp_peixe)
summary(Qesp_peixe)
dim(Peixe)

Qesp_canoas <- as.data.frame(Canoas)
names(Qesp_canoas)
summary(Qesp_canoas)
dim(Canoas)

Qesp_canoinhas <- as.data.frame(Canoinhas)
names(Qesp_canoinhas) <- "E65180000" 
summary(Qesp_canoinhas)
length(Canoinhas)

Qesp_pelotas <- as.data.frame(Pelotas)
names(Qesp_pelotas)
summary(Qesp_pelotas)
dim(Pelotas)

Qesp_antas <- as.data.frame(Antas)
names(Qesp_antas)
summary(Qesp_antas)
dim(Antas)

Qesp_ararangua <- as.data.frame(Ararangua)
names(Qesp_ararangua)
summary(Qesp_ararangua)
dim(Ararangua)

Qesp_chapeco <- as.data.frame(Chapeco)
names(Qesp_chapeco)
summary(Qesp_chapeco)
dim(Chapeco)

Qesp_irani <- as.data.frame(Irani)
names(Qesp_irani)
summary(Qesp_irani)
dim(Irani)

Qesp_itapocu <- as.data.frame(Itapocu)
names(Qesp_itapocu)
summary(Qesp_itapocu)
dim(Itapocu)

Qesp_negro <- as.data.frame(Negro)
names(Qesp_negro)
summary(Qesp_negro)
dim(Negro)

Qesp_peperi <- as.data.frame(Peperi)
names(Qesp_peperi)
summary(Qesp_peperi)
dim(Peperi)

Qesp_tijucas <- as.data.frame(Tijucas)
names(Qesp_tijucas)
summary(Qesp_tijucas)
dim(Tijucas)

Qesp_tubarao <- as.data.frame(Tubarao)
names(Qesp_tubarao)
summary(Qesp_tubarao)
dim(Tubarao)

Qesp_cubataoN <- as.data.frame(CubataoN)
names(Qesp_cubataoN) <- "E82270050" 
summary(Qesp_cubataoN)
length(CubataoN)

Qesp_cubataoS <- as.data.frame(CubataoS)
names(Qesp_cubataoS) <- "E84100000" 
summary(Qesp_cubataoS)
length(CubataoS)

Qesp_mampituba <- as.data.frame(Mampituba)
names(Qesp_mampituba) <- "E84970000" 
summary(Qesp_mampituba)
length(Mampituba)


(Param_itajai <- Qparam(Qesp_itajai))
(Param_peixe <- Qparam(Qesp_peixe))
(Param_canoas <- Qparam(Qesp_canoas))
(Param_canoinhas <- Qparam(Qesp_canoinhas))
(Param_pelotas <- Qparam(Qesp_pelotas))
(Param_antas <- Qparam(Qesp_antas))
(Param_ararangua <- Qparam(Qesp_ararangua))
(Param_chapeco <- Qparam(Qesp_chapeco))
(Param_cubataoN <- Qparam(Qesp_cubataoN))
(Param_cubataoS <- Qparam(Qesp_cubataoS))
(Param_irani <- Qparam(Qesp_irani))
(Param_itapocu <- Qparam(Qesp_itapocu))
(Param_mampituba <- Qparam(Qesp_mampituba))
(Param_negro <- Qparam(Qesp_negro))
(Param_peperi <- Qparam(Qesp_peperi))
(Param_tijucas <- Qparam(Qesp_tijucas))
(Param_tubarao <- Qparam(Qesp_tubarao))

(Param.ANO <- rbind(Param_peixe,Param_pelotas,Param_canoas,Param_itajai,Param_antas,Param_ararangua,Param_chapeco,Param_irani,Param_itapocu,Param_negro,Param_peperi,Param_tijucas,Param_tubarao,Param_canoinhas,Param_cubataoN,Param_cubataoS,Param_mampituba))
dim(Param.ANO)
##write.csv(Param.ANO,"Param/Param_ANO.csv")

##-----------------------------------------------------------------------------##

## Sazonal Verão 12/1/2
Qesp_itajaiDJF <- as.data.frame(extractzoo(Itajai,trgt="DJF"))
Qesp_peixeDJF <- as.data.frame(extractzoo(Peixe,trgt="DJF"))
Qesp_canoasDJF <- as.data.frame(extractzoo(Canoas,trgt="DJF"))
Qesp_canoinhasDJF <- as.data.frame(extractzoo(Canoinhas,trgt="DJF"))
names(Qesp_canoinhasDJF) <- "E65180000" 
Qesp_pelotasDJF <- as.data.frame(extractzoo(Pelotas,trgt="DJF"))
Qesp_itajaiDJF <- as.data.frame(extractzoo(Itajai,trgt="DJF"))
Qesp_antasDJF <- as.data.frame(extractzoo(Antas,trgt="DJF"))
Qesp_araranguaDJF <- as.data.frame(extractzoo(Ararangua,trgt="DJF"))
Qesp_chapecoDJF <- as.data.frame(extractzoo(Chapeco,trgt="DJF"))
Qesp_iraniDJF <- as.data.frame(extractzoo(Irani,trgt="DJF"))
Qesp_itapocuDJF <- as.data.frame(extractzoo(Itapocu,trgt="DJF"))
Qesp_negroDJF <- as.data.frame(extractzoo(Negro,trgt="DJF"))
Qesp_peperiDJF <- as.data.frame(extractzoo(Peperi,trgt="DJF"))
Qesp_tijucasDJF <- as.data.frame(extractzoo(Tijucas,trgt="DJF"))
Qesp_tubaraoDJF <- as.data.frame(extractzoo(Tubarao,trgt="DJF"))
Qesp_cubataoNDJF <- as.data.frame(extractzoo(CubataoN,trgt="DJF"))
names(Qesp_cubataoNDJF) <- "E82270050" 
Qesp_cubataoSDJF <- as.data.frame(extractzoo(CubataoS,trgt="DJF"))
names(Qesp_cubataoSDJF) <- "E84100000" 
Qesp_mampitubaDJF <- as.data.frame(extractzoo(Mampituba,trgt="DJF"))
names(Qesp_mampitubaDJF) <- "E84970000" 



summary(Qesp_itajai.DJF)
dim(Qesp_itajai.DJF)[1]

(ParamDJF_itajai <- Qparam(Qesp_itajaiDJF))
(ParamDJF_peixe <- Qparam(Qesp_peixeDJF))
(ParamDJF_canoas <- Qparam(Qesp_canoasDJF))
(ParamDJF_canoinhas <- Qparam(Qesp_canoinhasDJF))
(ParamDJF_pelotas <- Qparam(Qesp_pelotasDJF))
(ParamDJF_antas <- Qparam(Qesp_antasDJF))
(ParamDJF_ararangua <- Qparam(Qesp_araranguaDJF))
(ParamDJF_chapeco <- Qparam(Qesp_chapecoDJF))
(ParamDJF_cubataoN <- Qparam(Qesp_cubataoNDJF))
(ParamDJF_cubataoS <- Qparam(Qesp_cubataoSDJF))
(ParamDJF_irani <- Qparam(Qesp_iraniDJF))
(ParamDJF_itapocu <- Qparam(Qesp_itapocuDJF))
(ParamDJF_mampituba <- Qparam(Qesp_mampitubaDJF))
(ParamDJF_negro <- Qparam(Qesp_negroDJF))
(ParamDJF_peperi <- Qparam(Qesp_peperiDJF))
(ParamDJF_tijucas <- Qparam(Qesp_tijucasDJF))
(ParamDJF_tubarao <- Qparam(Qesp_tubaraoDJF))

(Param.DJF <- rbind(ParamDJF_peixe,ParamDJF_pelotas,ParamDJF_canoas,ParamDJF_itajai,ParamDJF_antas,ParamDJF_ararangua,ParamDJF_chapeco,ParamDJF_irani,ParamDJF_itapocu,ParamDJF_negro,ParamDJF_peperi,ParamDJF_tijucas,ParamDJF_tubarao,ParamDJF_canoinhas,ParamDJF_cubataoN,ParamDJF_cubataoS,ParamDJF_mampituba))
dim(Param.DJF)
##write.csv(Param.DJF,"Param/Param_DJF.csv")

##-----------------------------------------------------------------------------##

## Sazonal Outono 3/4/5
Qesp_itajaiMAM <- as.data.frame(extractzoo(Itajai,trgt="MAM"))
Qesp_peixeMAM <- as.data.frame(extractzoo(Peixe,trgt="MAM"))
Qesp_canoasMAM <- as.data.frame(extractzoo(Canoas,trgt="MAM"))
Qesp_canoinhasMAM <- as.data.frame(extractzoo(Canoinhas,trgt="MAM"))
names(Qesp_canoinhasMAM) <- "E65180000" 
Qesp_pelotasMAM <- as.data.frame(extractzoo(Pelotas,trgt="MAM"))
Qesp_itajaiMAM <- as.data.frame(extractzoo(Itajai,trgt="MAM"))
Qesp_antasMAM <- as.data.frame(extractzoo(Antas,trgt="MAM"))
Qesp_araranguaMAM <- as.data.frame(extractzoo(Ararangua,trgt="MAM"))
Qesp_chapecoMAM <- as.data.frame(extractzoo(Chapeco,trgt="MAM"))
Qesp_iraniMAM <- as.data.frame(extractzoo(Irani,trgt="MAM"))
Qesp_itapocuMAM <- as.data.frame(extractzoo(Itapocu,trgt="MAM"))
Qesp_negroMAM <- as.data.frame(extractzoo(Negro,trgt="MAM"))
Qesp_peperiMAM <- as.data.frame(extractzoo(Peperi,trgt="MAM"))
Qesp_tijucasMAM <- as.data.frame(extractzoo(Tijucas,trgt="MAM"))
Qesp_tubaraoMAM <- as.data.frame(extractzoo(Tubarao,trgt="MAM"))
Qesp_cubataoNMAM <- as.data.frame(extractzoo(CubataoN,trgt="MAM"))
names(Qesp_cubataoNMAM) <- "E82270050" 
Qesp_cubataoSMAM <- as.data.frame(extractzoo(CubataoS,trgt="MAM"))
names(Qesp_cubataoSMAM) <- "E84100000" 
Qesp_mampitubaMAM <- as.data.frame(extractzoo(Mampituba,trgt="MAM"))
names(Qesp_mampitubaMAM) <- "E84970000" 



summary(Qesp_itajai.MAM)
dim(Qesp_itajai.MAM)[1]

(ParamMAM_itajai <- Qparam(Qesp_itajaiMAM))
(ParamMAM_peixe <- Qparam(Qesp_peixeMAM))
(ParamMAM_canoas <- Qparam(Qesp_canoasMAM))
(ParamMAM_canoinhas <- Qparam(Qesp_canoinhasMAM))
(ParamMAM_pelotas <- Qparam(Qesp_pelotasMAM))
(ParamMAM_antas <- Qparam(Qesp_antasMAM))
(ParamMAM_ararangua <- Qparam(Qesp_araranguaMAM))
(ParamMAM_chapeco <- Qparam(Qesp_chapecoMAM))
(ParamMAM_cubataoN <- Qparam(Qesp_cubataoNMAM))
(ParamMAM_cubataoS <- Qparam(Qesp_cubataoSMAM))
(ParamMAM_irani <- Qparam(Qesp_iraniMAM))
(ParamMAM_itapocu <- Qparam(Qesp_itapocuMAM))
(ParamMAM_mampituba <- Qparam(Qesp_mampitubaMAM))
(ParamMAM_negro <- Qparam(Qesp_negroMAM))
(ParamMAM_peperi <- Qparam(Qesp_peperiMAM))
(ParamMAM_tijucas <- Qparam(Qesp_tijucasMAM))
(ParamMAM_tubarao <- Qparam(Qesp_tubaraoMAM))

(Param.MAM <- rbind(ParamMAM_peixe,ParamMAM_pelotas,ParamMAM_canoas,ParamMAM_itajai,ParamMAM_antas,ParamMAM_ararangua,ParamMAM_chapeco,ParamMAM_irani,ParamMAM_itapocu,ParamMAM_negro,ParamMAM_peperi,ParamMAM_tijucas,ParamMAM_tubarao,ParamMAM_canoinhas,ParamMAM_cubataoN,ParamMAM_cubataoS,ParamMAM_mampituba))
dim(Param.MAM)
##write.csv(Param.MAM,"Param/Param_MAM.csv")

##-----------------------------------------------------------------------------##

## Sazonal Outono 6/7/8
Qesp_itajaiJJA <- as.data.frame(extractzoo(Itajai,trgt="JJA"))
Qesp_peixeJJA <- as.data.frame(extractzoo(Peixe,trgt="JJA"))
Qesp_canoasJJA <- as.data.frame(extractzoo(Canoas,trgt="JJA"))
Qesp_canoinhasJJA <- as.data.frame(extractzoo(Canoinhas,trgt="JJA"))
names(Qesp_canoinhasJJA) <- "E65180000" 
Qesp_pelotasJJA <- as.data.frame(extractzoo(Pelotas,trgt="JJA"))
Qesp_itajaiJJA <- as.data.frame(extractzoo(Itajai,trgt="JJA"))
Qesp_antasJJA <- as.data.frame(extractzoo(Antas,trgt="JJA"))
Qesp_araranguaJJA <- as.data.frame(extractzoo(Ararangua,trgt="JJA"))
Qesp_chapecoJJA <- as.data.frame(extractzoo(Chapeco,trgt="JJA"))
Qesp_iraniJJA <- as.data.frame(extractzoo(Irani,trgt="JJA"))
Qesp_itapocuJJA <- as.data.frame(extractzoo(Itapocu,trgt="JJA"))
Qesp_negroJJA <- as.data.frame(extractzoo(Negro,trgt="JJA"))
Qesp_peperiJJA <- as.data.frame(extractzoo(Peperi,trgt="JJA"))
Qesp_tijucasJJA <- as.data.frame(extractzoo(Tijucas,trgt="JJA"))
Qesp_tubaraoJJA <- as.data.frame(extractzoo(Tubarao,trgt="JJA"))
Qesp_cubataoNJJA <- as.data.frame(extractzoo(CubataoN,trgt="JJA"))
names(Qesp_cubataoNJJA) <- "E82270050" 
Qesp_cubataoSJJA <- as.data.frame(extractzoo(CubataoS,trgt="JJA"))
names(Qesp_cubataoSJJA) <- "E84100000" 
Qesp_mampitubaJJA <- as.data.frame(extractzoo(Mampituba,trgt="JJA"))
names(Qesp_mampitubaJJA) <- "E84970000" 



summary(Qesp_itajai.JJA)
dim(Qesp_itajai.JJA)[1]

(ParamJJA_itajai <- Qparam(Qesp_itajaiJJA))
(ParamJJA_peixe <- Qparam(Qesp_peixeJJA))
(ParamJJA_canoas <- Qparam(Qesp_canoasJJA))
(ParamJJA_canoinhas <- Qparam(Qesp_canoinhasJJA))
(ParamJJA_pelotas <- Qparam(Qesp_pelotasJJA))
(ParamJJA_antas <- Qparam(Qesp_antasJJA))
(ParamJJA_ararangua <- Qparam(Qesp_araranguaJJA))
(ParamJJA_chapeco <- Qparam(Qesp_chapecoJJA))
(ParamJJA_cubataoN <- Qparam(Qesp_cubataoNJJA))
(ParamJJA_cubataoS <- Qparam(Qesp_cubataoSJJA))
(ParamJJA_irani <- Qparam(Qesp_iraniJJA))
(ParamJJA_itapocu <- Qparam(Qesp_itapocuJJA))
(ParamJJA_mampituba <- Qparam(Qesp_mampitubaJJA))
(ParamJJA_negro <- Qparam(Qesp_negroJJA))
(ParamJJA_peperi <- Qparam(Qesp_peperiJJA))
(ParamJJA_tijucas <- Qparam(Qesp_tijucasJJA))
(ParamJJA_tubarao <- Qparam(Qesp_tubaraoJJA))

(Param.JJA <- rbind(ParamJJA_peixe,ParamJJA_pelotas,ParamJJA_canoas,ParamJJA_itajai,ParamJJA_antas,ParamJJA_ararangua,ParamJJA_chapeco,ParamJJA_irani,ParamJJA_itapocu,ParamJJA_negro,ParamJJA_peperi,ParamJJA_tijucas,ParamJJA_tubarao,ParamJJA_canoinhas,ParamJJA_cubataoN,ParamJJA_cubataoS,ParamJJA_mampituba))
dim(Param.JJA)
##write.csv(Param.JJA,"Param/Param_JJA.csv")
##-----------------------------------------------------------------------------##

## Sazonal Primavera 9/10/11
Qesp_itajaiSON <- as.data.frame(extractzoo(Itajai,trgt="SON"))
Qesp_peixeSON <- as.data.frame(extractzoo(Peixe,trgt="SON"))
Qesp_canoasSON <- as.data.frame(extractzoo(Canoas,trgt="SON"))
Qesp_canoinhasSON <- as.data.frame(extractzoo(Canoinhas,trgt="SON"))
names(Qesp_canoinhasSON) <- "E65180000" 
Qesp_pelotasSON <- as.data.frame(extractzoo(Pelotas,trgt="SON"))
Qesp_itajaiSON <- as.data.frame(extractzoo(Itajai,trgt="SON"))
Qesp_antasSON <- as.data.frame(extractzoo(Antas,trgt="SON"))
Qesp_araranguaSON <- as.data.frame(extractzoo(Ararangua,trgt="SON"))
Qesp_chapecoSON <- as.data.frame(extractzoo(Chapeco,trgt="SON"))
Qesp_iraniSON <- as.data.frame(extractzoo(Irani,trgt="SON"))
Qesp_itapocuSON <- as.data.frame(extractzoo(Itapocu,trgt="SON"))
Qesp_negroSON <- as.data.frame(extractzoo(Negro,trgt="SON"))
Qesp_peperiSON <- as.data.frame(extractzoo(Peperi,trgt="SON"))
Qesp_tijucasSON <- as.data.frame(extractzoo(Tijucas,trgt="SON"))
Qesp_tubaraoSON <- as.data.frame(extractzoo(Tubarao,trgt="SON"))
Qesp_cubataoNSON <- as.data.frame(extractzoo(CubataoN,trgt="SON"))
names(Qesp_cubataoNSON) <- "E82270050" 
Qesp_cubataoSSON <- as.data.frame(extractzoo(CubataoS,trgt="SON"))
names(Qesp_cubataoSSON) <- "E84100000" 
Qesp_mampitubaSON <- as.data.frame(extractzoo(Mampituba,trgt="SON"))
names(Qesp_mampitubaSON) <- "E84970000" 



summary(Qesp_itajai.SON)
dim(Qesp_itajai.SON)[1]

(ParamSON_itajai <- Qparam(Qesp_itajaiSON))
(ParamSON_peixe <- Qparam(Qesp_peixeSON))
(ParamSON_canoas <- Qparam(Qesp_canoasSON))
(ParamSON_canoinhas <- Qparam(Qesp_canoinhasSON))
(ParamSON_pelotas <- Qparam(Qesp_pelotasSON))
(ParamSON_antas <- Qparam(Qesp_antasSON))
(ParamSON_ararangua <- Qparam(Qesp_araranguaSON))
(ParamSON_chapeco <- Qparam(Qesp_chapecoSON))
(ParamSON_cubataoN <- Qparam(Qesp_cubataoNSON))
(ParamSON_cubataoS <- Qparam(Qesp_cubataoSSON))
(ParamSON_irani <- Qparam(Qesp_iraniSON))
(ParamSON_itapocu <- Qparam(Qesp_itapocuSON))
(ParamSON_mampituba <- Qparam(Qesp_mampitubaSON))
(ParamSON_negro <- Qparam(Qesp_negroSON))
(ParamSON_peperi <- Qparam(Qesp_peperiSON))
(ParamSON_tijucas <- Qparam(Qesp_tijucasSON))
(ParamSON_tubarao <- Qparam(Qesp_tubaraoSON))

(Param.SON <- rbind(ParamSON_peixe,ParamSON_pelotas,ParamSON_canoas,ParamSON_itajai,ParamSON_antas,ParamSON_ararangua,ParamSON_chapeco,ParamSON_irani,ParamSON_itapocu,ParamSON_negro,ParamSON_peperi,ParamSON_tijucas,ParamSON_tubarao,ParamSON_canoinhas,ParamSON_cubataoN,ParamSON_cubataoS,ParamSON_mampituba))
dim(Param.SON)
##write.csv(Param.SON,"Param/Param_SON.csv")
##-----------------------------------------------------------------------------##

## Juntar todos os parâmetros para sazonal
(Param <- cbind(Param.ANO,Param.DJF,Param.MAM,Param.JJA,Param.SON))

write.csv(Param,"Param/Param_todos.csv")

##-----------------------------------------------------------------------------##

itajai(dimnames_dimnames$Param$1[[warnings]])
cbind()
Param_itajai$names

fit <- list()
for(i in 1:ncol(Qesp_itajai)){
    fit$E[i] <- fitdist(Qesp_itajai[,i], "lnorm")}
return(fit)

x <- do.call(cbind.data.frame, fit)
names(x) <- colnames(Qesp_itajai)
t(x)-param.itajai


df <- data.frame(matrix(unlist(fit), nrow=2, byrow=F))
t(df)
names(df)
as.matrix(fit)
names(fit)
as.vector(colnames(Qesp_itajai)[1])
as.name(paste("fit.","bacia",sep="")) <- list()

        <- list()


fit.itajaiE1 <- fitdist(Qesp_itajai[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE2 <- fitdist(Qesp_itajai[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE3 <- fitdist(Qesp_itajai[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE4 <- fitdist(Qesp_itajai[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")

cdfcomp(fit.itajaiE1)
fit.itajaiE1

fit.itajaiE1 <- fitdist(Qesp_itajai[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE2 <- fitdist(Qesp_itajai[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE3 <- fitdist(Qesp_itajai[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE4 <- fitdist(Qesp_itajai[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE5 <- fitdist(Qesp_itajai[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE6 <- fitdist(Qesp_itajai[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE7 <- fitdist(Qesp_itajai[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE8 <- fitdist(Qesp_itajai[,8], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE9 <- fitdist(Qesp_itajai[,9], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE10 <- fitdist(Qesp_itajai[,10], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE11 <- fitdist(Qesp_itajai[,11], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE12 <- fitdist(Qesp_itajai[,12], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE13 <- fitdist(Qesp_itajai[,13], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE14 <- fitdist(Qesp_itajai[,14], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE15 <- fitdist(Qesp_itajai[,15], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE16 <- fitdist(Qesp_itajai[,16], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE17 <- fitdist(Qesp_itajai[,17], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiE18 <- fitdist(Qesp_itajai[,18], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


param.itajai <- t(data.frame(fit.itajaiE1$estimate,fit.itajaiE2$estimate,fit.itajaiE3$estimate,fit.itajaiE4$estimate,fit.itajaiE5$estimate,fit.itajaiE6$estimate,fit.itajaiE7$estimate,fit.itajaiE8$estimate,fit.itajaiE9$estimate,fit.itajaiE10$estimate,fit.itajaiE11$estimate,fit.itajaiE12$estimate,fit.itajaiE13$estimate,fit.itajaiE14$estimate,fit.itajaiE15$estimate,fit.itajaiE16$estimate,fit.itajaiE17$estimate,fit.itajaiE18$estimate));rownames(param.itajai) <- colnames(Itajai);param.itajai

##-----------------------------------------------------------------------------##
## Curva de Permanência sazonal
## meses 12/1/2

fdc(extractzoo(Itajai,trgt="DJF"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_itajai.DJF <- as.data.frame(extractzoo(Itajai,trgt="DJF"))
summary(Qesp_itajai.DJF)
dim(Qesp_itajai.DJF)
plot(Qesp_itajai.DJF)

fit.itajaiDJF.E1 <- fitdist(Qesp_itajai.DJF[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E2 <- fitdist(Qesp_itajai.DJF[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E3 <- fitdist(Qesp_itajai.DJF[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E4 <- fitdist(Qesp_itajai.DJF[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E5 <- fitdist(Qesp_itajai.DJF[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E6 <- fitdist(Qesp_itajai.DJF[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E7 <- fitdist(Qesp_itajai.DJF[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E8 <- fitdist(Qesp_itajai.DJF[,8], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E9 <- fitdist(Qesp_itajai.DJF[,9], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E10 <- fitdist(Qesp_itajai.DJF[,10], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E11 <- fitdist(Qesp_itajai.DJF[,11], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E12 <- fitdist(Qesp_itajai.DJF[,12], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E13 <- fitdist(Qesp_itajai.DJF[,13], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E14 <- fitdist(Qesp_itajai.DJF[,14], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E15 <- fitdist(Qesp_itajai.DJF[,15], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E16 <- fitdist(Qesp_itajai.DJF[,16], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E17 <- fitdist(Qesp_itajai.DJF[,17], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiDJF.E18 <- fitdist(Qesp_itajai.DJF[,18], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")




param.itajai.DJF <-t(data.frame(fit.itajaiDJF.E1$estimate,fit.itajaiDJF.E2$estimate,fit.itajaiDJF.E3$estimate,fit.itajaiDJF.E4$estimate,fit.itajaiDJF.E5$estimate,fit.itajaiDJF.E6$estimate,fit.itajaiDJF.E7$estimate,fit.itajaiDJF.E8$estimate,fit.itajaiDJF.E9$estimate,fit.itajaiDJF.E10$estimate,fit.itajaiDJF.E11$estimate,fit.itajaiDJF.E12$estimate,fit.itajaiDJF.E13$estimate,fit.itajaiDJF.E14$estimate,fit.itajaiDJF.E15$estimate,fit.itajaiDJF.E16$estimate,fit.itajaiDJF.E17$estimate,fit.itajaiDJF.E18$estimate));rownames(param.itajai.DJF) <- colnames(Itajai);param.itajai.DJF

## meses 3/4/5

fdc(extractzoo(Itajai,trgt="MAM"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_itajai.MAM <- as.data.frame(extractzoo(Itajai,trgt="MAM"))
summary(Qesp_itajai.MAM)
dim(Qesp_itajai.MAM)
plot(Qesp_itajai.MAM)


fit.itajaiMAM.E1 <- fitdist(Qesp_itajai.MAM[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E2 <- fitdist(Qesp_itajai.MAM[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E3 <- fitdist(Qesp_itajai.MAM[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E4 <- fitdist(Qesp_itajai.MAM[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E5 <- fitdist(Qesp_itajai.MAM[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E6 <- fitdist(Qesp_itajai.MAM[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E7 <- fitdist(Qesp_itajai.MAM[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E8 <- fitdist(Qesp_itajai.MAM[,8], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E9 <- fitdist(Qesp_itajai.MAM[,9], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E10 <- fitdist(Qesp_itajai.MAM[,10], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E11 <- fitdist(Qesp_itajai.MAM[,11], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E12 <- fitdist(Qesp_itajai.MAM[,12], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E13 <- fitdist(Qesp_itajai.MAM[,13], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E14 <- fitdist(Qesp_itajai.MAM[,14], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E15 <- fitdist(Qesp_itajai.MAM[,15], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E16 <- fitdist(Qesp_itajai.MAM[,16], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E17 <- fitdist(Qesp_itajai.MAM[,17], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiMAM.E18 <- fitdist(Qesp_itajai.MAM[,18], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


param.itajai.MAM <-t(data.frame(fit.itajaiMAM.E1$estimate,fit.itajaiMAM.E2$estimate,fit.itajaiMAM.E3$estimate,fit.itajaiMAM.E4$estimate,fit.itajaiMAM.E5$estimate,fit.itajaiMAM.E6$estimate,fit.itajaiMAM.E7$estimate,fit.itajaiMAM.E8$estimate,fit.itajaiMAM.E9$estimate,fit.itajaiMAM.E10$estimate,fit.itajaiMAM.E11$estimate,fit.itajaiMAM.E12$estimate,fit.itajaiMAM.E13$estimate,fit.itajaiMAM.E14$estimate,fit.itajaiMAM.E15$estimate,fit.itajaiMAM.E16$estimate,fit.itajaiMAM.E17$estimate,fit.itajaiMAM.E18$estimate));rownames(param.itajai.MAM) <- colnames(Itajai);param.itajai.MAM

## meses 6/7/8

fdc(extractzoo(Itajai,trgt="JJA"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_itajai.JJA <- as.data.frame(extractzoo(Itajai,trgt="JJA"))
summary(Qesp_itajai.JJA)
dim(Qesp_itajai.JJA)
plot(Qesp_itajai.JJA)


fit.itajaiJJA.E1 <- fitdist(Qesp_itajai.JJA[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E2 <- fitdist(Qesp_itajai.JJA[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E3 <- fitdist(Qesp_itajai.JJA[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E4 <- fitdist(Qesp_itajai.JJA[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E5 <- fitdist(Qesp_itajai.JJA[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E6 <- fitdist(Qesp_itajai.JJA[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E7 <- fitdist(Qesp_itajai.JJA[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E8 <- fitdist(Qesp_itajai.JJA[,8], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E9 <- fitdist(Qesp_itajai.JJA[,9], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E10 <- fitdist(Qesp_itajai.JJA[,10], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E11 <- fitdist(Qesp_itajai.JJA[,11], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E12 <- fitdist(Qesp_itajai.JJA[,12], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E13 <- fitdist(Qesp_itajai.JJA[,13], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E14 <- fitdist(Qesp_itajai.JJA[,14], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E15 <- fitdist(Qesp_itajai.JJA[,15], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E16 <- fitdist(Qesp_itajai.JJA[,16], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E17 <- fitdist(Qesp_itajai.JJA[,17], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiJJA.E18 <- fitdist(Qesp_itajai.JJA[,18], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


param.itajai.JJA <-t(data.frame(fit.itajaiJJA.E1$estimate,fit.itajaiJJA.E2$estimate,fit.itajaiJJA.E3$estimate,fit.itajaiJJA.E4$estimate,fit.itajaiJJA.E5$estimate,fit.itajaiJJA.E6$estimate,fit.itajaiJJA.E7$estimate,fit.itajaiJJA.E8$estimate,fit.itajaiJJA.E9$estimate,fit.itajaiJJA.E10$estimate,fit.itajaiJJA.E11$estimate,fit.itajaiJJA.E12$estimate,fit.itajaiJJA.E13$estimate,fit.itajaiJJA.E14$estimate,fit.itajaiJJA.E15$estimate,fit.itajaiJJA.E16$estimate,fit.itajaiJJA.E17$estimate,fit.itajaiJJA.E18$estimate));rownames(param.itajai.JJA) <- colnames(Itajai);param.itajai.JJA

## meses 9/10/11

fdc(extractzoo(Itajai,trgt="SON"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_itajai.SON <- as.data.frame(extractzoo(Itajai,trgt="SON"))
summary(Qesp_itajai.SON)
dim(Qesp_itajai.SON)
plot(Qesp_itajai.SON)


fit.itajaiSON.E1 <- fitdist(Qesp_itajai.SON[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E2 <- fitdist(Qesp_itajai.SON[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E3 <- fitdist(Qesp_itajai.SON[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E4 <- fitdist(Qesp_itajai.SON[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E5 <- fitdist(Qesp_itajai.SON[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E6 <- fitdist(Qesp_itajai.SON[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E7 <- fitdist(Qesp_itajai.SON[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E8 <- fitdist(Qesp_itajai.SON[,8], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E9 <- fitdist(Qesp_itajai.SON[,9], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E10 <- fitdist(Qesp_itajai.SON[,10], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E11 <- fitdist(Qesp_itajai.SON[,11], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E12 <- fitdist(Qesp_itajai.SON[,12], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E13 <- fitdist(Qesp_itajai.SON[,13], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E14 <- fitdist(Qesp_itajai.SON[,14], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E15 <- fitdist(Qesp_itajai.SON[,15], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E16 <- fitdist(Qesp_itajai.SON[,16], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E17 <- fitdist(Qesp_itajai.SON[,17], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.itajaiSON.E18 <- fitdist(Qesp_itajai.SON[,18], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


param.itajai.SON <-t(data.frame(fit.itajaiSON.E1$estimate,fit.itajaiSON.E2$estimate,fit.itajaiSON.E3$estimate,fit.itajaiSON.E4$estimate,fit.itajaiSON.E5$estimate,fit.itajaiSON.E6$estimate,fit.itajaiSON.E7$estimate,fit.itajaiSON.E8$estimate,fit.itajaiSON.E9$estimate,fit.itajaiSON.E10$estimate,fit.itajaiSON.E11$estimate,fit.itajaiSON.E12$estimate,fit.itajaiSON.E13$estimate,fit.itajaiSON.E14$estimate,fit.itajaiSON.E15$estimate,fit.itajaiSON.E16$estimate,fit.itajaiSON.E17$estimate,fit.itajaiSON.E18$estimate));rownames(param.itajai.SON) <- colnames(Itajai);param.itajai.SON

##-----------------------------------------------------------------------------##
## Peixe

Qesp_peixe <- as.data.frame(Peixe)
names(Qesp_peixe)
summary(Qesp_peixe)

fit.peixeE1 <- fitdist(Qesp_peixe[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeE2 <- fitdist(Qesp_peixe[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeE3 <- fitdist(Qesp_peixe[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeE4 <- fitdist(Qesp_peixe[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeE5 <- fitdist(Qesp_peixe[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeE6 <- fitdist(Qesp_peixe[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeE7 <- fitdist(Qesp_peixe[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


param.peixe <- t(data.frame(fit.peixeE1$estimate,fit.peixeE2$estimate,fit.peixeE3$estimate,fit.peixeE4$estimate,fit.peixeE5$estimate,fit.peixeE6$estimate,fit.peixeE7$estimate));rownames(param.peixe) <- colnames(Peixe);param.peixe

##-----------------------------------------------------------------------------##
## Curva de Permanência sazonal

## meses 12/1/2

fdc(extractzoo(Peixe,trgt="DJF"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_peixe.DJF <- as.data.frame(extractzoo(Peixe,trgt="DJF"))
summary(Qesp_peixe.DJF)
dim(Qesp_peixe.DJF)
plot(Qesp_peixe.DJF)

fit.peixeDJF.E1 <- fitdist(Qesp_peixe.DJF[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeDJF.E2 <- fitdist(Qesp_peixe.DJF[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeDJF.E3 <- fitdist(Qesp_peixe.DJF[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeDJF.E4 <- fitdist(Qesp_peixe.DJF[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeDJF.E5 <- fitdist(Qesp_peixe.DJF[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeDJF.E6 <- fitdist(Qesp_peixe.DJF[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeDJF.E7 <- fitdist(Qesp_peixe.DJF[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")

param.peixe.DJF <-t(data.frame(fit.peixeDJF.E1$estimate,fit.peixeDJF.E2$estimate,fit.peixeDJF.E3$estimate,fit.peixeDJF.E4$estimate,fit.peixeDJF.E5$estimate,fit.peixeDJF.E6$estimate,fit.peixeDJF.E7$estimate));rownames(param.peixe.DJF) <- colnames(Peixe);param.peixe.DJF

## meses 3/4/5

fdc(extractzoo(Peixe,trgt="MAM"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_peixe.MAM <- as.data.frame(extractzoo(Peixe,trgt="MAM"))
summary(Qesp_peixe.MAM)
dim(Qesp_peixe.MAM)
plot(Qesp_peixe.MAM)

fit.peixeMAM.E1 <- fitdist(Qesp_peixe.MAM[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeMAM.E2 <- fitdist(Qesp_peixe.MAM[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeMAM.E3 <- fitdist(Qesp_peixe.MAM[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeMAM.E4 <- fitdist(Qesp_peixe.MAM[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeMAM.E5 <- fitdist(Qesp_peixe.MAM[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeMAM.E6 <- fitdist(Qesp_peixe.MAM[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeMAM.E7 <- fitdist(Qesp_peixe.MAM[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")

param.peixe.MAM <- t(data.frame(fit.peixeMAM.E1$estimate,fit.peixeMAM.E2$estimate,fit.peixeMAM.E3$estimate,fit.peixeMAM.E4$estimate,fit.peixeMAM.E5$estimate,fit.peixeMAM.E6$estimate,fit.peixeMAM.E7$estimate));rownames(param.peixe.MAM) <- colnames(Peixe);param.peixe.MAM

## meses 6/7/8

fdc(extractzoo(Peixe,trgt="JJA"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_peixe.JJA <- as.data.frame(extractzoo(Peixe,trgt="JJA"))
summary(Qesp_peixe.JJA)
dim(Qesp_peixe.JJA)
plot(Qesp_peixe.JJA)

fit.peixeJJA.E1 <- fitdist(Qesp_peixe.JJA[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeJJA.E2 <- fitdist(Qesp_peixe.JJA[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeJJA.E3 <- fitdist(Qesp_peixe.JJA[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeJJA.E4 <- fitdist(Qesp_peixe.JJA[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeJJA.E5 <- fitdist(Qesp_peixe.JJA[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeJJA.E6 <- fitdist(Qesp_peixe.JJA[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeJJA.E7 <- fitdist(Qesp_peixe.JJA[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")

param.peixe.JJA <-t(data.frame(fit.peixeJJA.E1$estimate,fit.peixeJJA.E2$estimate,fit.peixeJJA.E3$estimate,fit.peixeJJA.E4$estimate,fit.peixeJJA.E5$estimate,fit.peixeJJA.E6$estimate,fit.peixeJJA.E7$estimate));rownames(param.peixe.JJA) <- colnames(Peixe);param.peixe.JJA

## meses 9/10/11

fdc(extractzoo(Peixe,trgt="SON"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="y", main= "", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

Qesp_peixe.SON <- as.data.frame(extractzoo(Peixe,trgt="SON"))
summary(Qesp_peixe.SON)
dim(Qesp_peixe.SON)
plot(Qesp_peixe.SON)

fit.peixeSON.E1 <- fitdist(Qesp_peixe.SON[,1], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeSON.E2 <- fitdist(Qesp_peixe.SON[,2], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeSON.E3 <- fitdist(Qesp_peixe.SON[,3], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeSON.E4 <- fitdist(Qesp_peixe.SON[,4], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeSON.E5 <- fitdist(Qesp_peixe.SON[,5], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeSON.E6 <- fitdist(Qesp_peixe.SON[,6], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fit.peixeSON.E7 <- fitdist(Qesp_peixe.SON[,7], "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")

param.peixe.SON <- t(data.frame(fit.peixeSON.E1$estimate,fit.peixeSON.E2$estimate,fit.peixeSON.E3$estimate,fit.peixeSON.E4$estimate,fit.peixeSON.E5$estimate,fit.peixeSON.E6$estimate,fit.peixeSON.E7$estimate));rownames(param.peixe.SON) <- colnames(Peixe);param.peixe.SON

##-----------------------------------------------------------------------------##




Qp.DJF1<-quantile(extractzoo(Peixe[,1],trgt="DJF"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.DJF2<-quantile(extractzoo(Peixe[,2],trgt="DJF"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.DJF3<-quantile(extractzoo(Peixe[,3],trgt="DJF"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.DJF4<-quantile(extractzoo(Peixe[,4],trgt="DJF"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.DJF5<-quantile(extractzoo(Peixe[,5],trgt="DJF"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.DJF6<-quantile(extractzoo(Peixe[,6],trgt="DJF"),prob=seq(0,1,by=0.05),na.rm=T)

Qp.DJF<-cbind(Qp.DJF1,Qp.DJF2,Qp.DJF3,Qp.DJF4,Qp.DJF5,Qp.DJF6)
Qpesp.DJF<-cbind(Qp.DJF1/areas[1],Qp.DJF2/areas[2],Qp.DJF3/areas[3],Qp.DJF4/areas[4],Qp.DJF5/areas[5],Qp.DJF6/areas[6])

# meses 3/4/5

fdc(extractzoo(Q,trgt="MAM"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="", main= "Qp Outono", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))
Qp.MAM1<-quantile(extractzoo(Q[,1],trgt="MAM"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.MAM2<-quantile(extractzoo(Q[,2],trgt="MAM"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.MAM3<-quantile(extractzoo(Q[,3],trgt="MAM"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.MAM4<-quantile(extractzoo(Q[,4],trgt="MAM"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.MAM5<-quantile(extractzoo(Q[,5],trgt="MAM"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.MAM6<-quantile(extractzoo(Q[,6],trgt="MAM"),prob=seq(0,1,by=0.05),na.rm=T)

Qp.MAM<-cbind(Qp.MAM1,Qp.MAM2,Qp.MAM3,Qp.MAM4,Qp.MAM5,Qp.MAM6)
Qpesp.MAM<-cbind(Qp.MAM1/areas[1],Qp.MAM2/areas[2],Qp.MAM3/areas[3],Qp.MAM4/areas[4],Qp.MAM5/areas[5],Qp.MAM6/areas[6])

# meses 6/7/8

fdc(extractzoo(Q,trgt="JJA"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="", main= "Qp Inverno", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))
Qp.JJA1<-quantile(extractzoo(Q[,1],trgt="JJA"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.JJA2<-quantile(extractzoo(Q[,2],trgt="JJA"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.JJA3<-quantile(extractzoo(Q[,3],trgt="JJA"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.JJA4<-quantile(extractzoo(Q[,4],trgt="JJA"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.JJA5<-quantile(extractzoo(Q[,5],trgt="JJA"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.JJA6<-quantile(extractzoo(Q[,6],trgt="JJA"),prob=seq(0,1,by=0.05),na.rm=T)

Qp.JJA<-cbind(Qp.JJA1,Qp.JJA2,Qp.JJA3,Qp.JJA4,Qp.JJA5,Qp.JJA6)
Qpesp.JJA<-cbind(Qp.JJA1/areas[1],Qp.JJA2/areas[2],Qp.JJA3/areas[3],Qp.JJA4/areas[4],Qp.JJA5/areas[5],Qp.JJA6/areas[6])

# meses 9/10/11

fdc(extractzoo(Q,trgt="SON"), lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="", main= "Qp Inverno", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))
Qp.SON1<-quantile(extractzoo(Q[,1],trgt="SON"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.SON2<-quantile(extractzoo(Q[,2],trgt="SON"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.SON3<-quantile(extractzoo(Q[,3],trgt="SON"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.SON4<-quantile(extractzoo(Q[,4],trgt="SON"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.SON5<-quantile(extractzoo(Q[,5],trgt="SON"),prob=seq(0,1,by=0.05),na.rm=T)
Qp.SON6<-quantile(extractzoo(Q[,6],trgt="SON"),prob=seq(0,1,by=0.05),na.rm=T)

Qp.SON<-cbind(Qp.SON1,Qp.SON2,Qp.SON3,Qp.SON4,Qp.SON5,Qp.SON6)
Qpesp.SON<-cbind(Qp.SON1/areas[1],Qp.SON2/areas[2],Qp.SON3/areas[3],Qp.SON4/areas[4],Qp.SON5/areas[5],Qp.SON6/areas[6])


# Ano hidrológico
Qhy<-data.frame(data,Q)
Qhy.esp<-data.frame(data,Qimp)

summary(Qhy)
summary(Qhy.esp)
head(Qhy.esp)

## Criar um vetor que começa no ano hidrológico
(breaks <- seq(as.Date("1960-10-01"), length=54,by="year"))
hydroYear <-cut(Qhy$data, breaks, labels=1960:2012)
ano<-seq(as.Date("1960/01/01"),as.Date("2012/12/31"),"year")
class(Qhy)
head(Qhy)

##write.table(Qhy.2,file = "Qhydr.txt",dec = ",")


## vazão máxima pra cada ano hidrológico
Qmaxhy<-zoo(aggregate(Qhy, list(hydroYear), max)[,-c(1:2)],ano)
plot(Qmaxhy,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
head(Qmaxhy)
summary(Qmaxhy)

Qmaxhy <- as.data.frame(Qmaxhy)
names(Qmaxhy)

## Q7 pra cada ano hidrológico

Q.7<-zoo(aggregate(Qhy, list(hydroYear), function(x) min(rollmean(x,7,align = "left")))[,-c(1:2)],ano)
head(Q.7)
plot(Q.7,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

Q.7 <- as.data.frame(Q.7)

fdc(Q.7, lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="", main="Q7 Anual", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))

##Q.7 <- as.data.frame.ts(Q.7)
)##fdc <- ecdf(Q.7[,1])

Qimp <- as.data.frame(Qimp)
colnames(Qimp) <- c("Q1","Q2","Q3","Q4")
names(Qimp)


fitE4bs <- fitdist(Qimp$Q4, "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fitE3bs <- fitdist(Qimp$Q3, "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fitE2bs <- fitdist(Qimp$Q2, "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")
fitE1bs <- fitdist(Qimp$Q1, "lnorm")#,start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


gofstat(fitE1bs)
plot(fitE1bs)
summary(fitE1bs)
(param <- data.frame(fitE1bs$estimate,fitE2bs$estimate,fitE3bs$estimate,fitE4bs$estimate))

plot(param)

qlnorm(1-0.98,fitE1bs$estimate[1],fitE1bs$estimate[2])*areas[1]

quantile(Qimp$Q1,0.9)
seq(0,0.99,l=10)
bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
fit <- vglm(Qimp$Q1 ~ 1, bisa, data = bdata2, trace = TRUE)
with(bdata2, hist(Qimp$Q1, prob = TRUE, ylim = c(0, 0.5), col = "lightblue"))
coef(fit, matrix = TRUE)
with(bdata2, mean(Qimp$Q1))


fitE3<-list()
    
fitE3$we <- fitdist(Qimp$Q3, "weibull")
fitE3$ga <- fitdist(Qimp$Q3, "gamma")
fitE3$ln <- fitdist(Qimp$Q3, "lnorm")
fitE3$no <- fitdist(Qimp$Q3, "norm")
fitE3$ca <- fitdist(Qimp$Q3, "cauchy")
fitE3$lo <- fitdist(Qimp$Q3, "logis")
fitE3$gu <- fitdist(Qimp$Q3, "gumbel",start = list(loc=10,scale=10))
fitE3$lg <- fitdist(Qimp$Q3, "lgamma",start = c(shapelog = 10, ratelog=10) ,lower =-Inf, upper = Inf,optim.method="BFGS")
fitE3$bu <- fitdist(Qimp$Q3, "burr",start = c(shape1 = 10, shape2= 10, scale=10),lower =0.001, upper = Inf)
fitE3$pe <- fitdist(Qimp$Q3, "trbeta",start = c(shape1=10, shape2=10, shape3=10,
                                 scale = 10),lower = 0.001, upper = Inf) # Pearson 6
fitE3$ge <- fitdist(Qimp$Q3, "gev",start = c(loc=10, scale=10, shape=10),lower =-Inf, upper = Inf,optim.method="BFGS")
fitE3$bs <- fitdist(Qimp$Q3, "bisa",start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


sort(gofstat(fitE2)$aic)
sort(gofstat(fitE1)$aic)
sort(gofstat(fitE3)$aic)
sort(gofstat(fitE4)$aic)

dbisa(
plot(fitE1$no)
plot(fitE1$ln,aps=1)
plot(fitE1$bs,asp=1)
summary(fitE1$ln)
plot(fitE4$bu)
plot(fitE1$bu, demp = TRUE)
plot(fitE1$bu, histo = FALSE, demp = F)


##=============================================================================##
## Ajuste das vazões máximas


## Estação 1
fitmaxE1<-list()

fitmaxE1$we <- fitdist(Qmaxhy$Q1, "weibull")
fitmaxE1$ga <- fitdist(Qmaxhy$Q1, "gamma")
fitmaxE1$ln <- fitdist(Qmaxhy$Q1, "lnorm")
fitmaxE1$no <- fitdist(Qmaxhy$Q1, "norm")
fitmaxE1$ca <- fitdist(Qmaxhy$Q1, "cauchy")
fitmaxE1$lo <- fitdist(Qmaxhy$Q1, "logis")
fitmaxE1$gu <- fitdist(Qmaxhy$Q1, "gumbel",start = list(loc=10,scale=10))
fitmaxE1$lg <- fitdist(Qmaxhy$Q1, "lgamma",start = c(shapelog = 10, ratelog=10), lower =0.0001 , upper = Inf)
fitmaxE1$bu <- fitdist(Qmaxhy$Q1, "burr",start = c(shape1 = 10, shape2= 10, scale=10),lower =0.001, upper = Inf)
fitmaxE1$pe <- fitdist(Qmaxhy$Q1, "trbeta",start = c(shape1=10, shape2=10, shape3=10,
                                               scale = 10),lower = 0.001, upper = Inf) # Pearson 6
fitmaxE1$ge <- fitdist(Qmaxhy$Q1, "gev",start = c(loc=10, scale=10, shape=10),lower =-Inf, upper = Inf,optim.method="BFGS")
fitmaxE1$bs <- fitdist(Qmaxhy$Q1, "bisa",start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")


# Escolher o menor valor
sort(gofstat(fitmaxE1)$aic)
sort(gofstat(fitmaxE1)$bic)

plot(fitmaxE1$bs)
plot(fitmaxE1$lg)
plot(fitmaxE1$ln)

summary(fitmaxE1$ge)
summary(fitmaxE1$lg)

fitln <- fitdist(Qhy$Q2, "lnorm")
summary(fitln)


    
## Estação 2
fitmaxE2<-list()

fitmaxE2$we <- fitdist(Qmaxhy$Q2, "weibull")
fitmaxE2$ga <- fitdist(Qmaxhy$Q2, "gamma")
fitmaxE2$ln <- fitdist(Qmaxhy$Q2, "lnorm")
fitmaxE2$no <- fitdist(Qmaxhy$Q2, "norm")
fitmaxE2$ca <- fitdist(Qmaxhy$Q2, "cauchy")
fitmaxE2$lo <- fitdist(Qmaxhy$Q2, "logis")
fitmaxE2$gu <- fitdist(Qmaxhy$Q2, "gumbel",start = list(location=10,scale=10),lower =0.001, upper = Inf)
fitmaxE2$lg <- fitdist(Qmaxhy$Q2, "lgamma",start = c(location = 10, scale= 10,shape=10), lower =0.0001 , upper = Inf)
fitmaxE2$bu <- fitdist(Qmaxhy$Q2, "burr",start = c(shape1 = 10, shape2= 10, scale=10),lower =-Inf, upper = Inf)
fitmaxE2$pe <- fitdist(Qmaxhy$Q2, "trbeta",start = c(shape1=10, shape2=10, shape3=10,
                                             scale = 10),lower = 0.001, upper = Inf) # Pearson 6
fitmaxE2$ge <- fitdist(Qmaxhy$Q2, "gev",start = c(location=10, scale=10, shape=10),lower =-Inf, upper = Inf,optim.method="BFGS")
fitmaxE2$bs <- fitdist(Qmaxhy$Q2, "bisa",start = c(scale=10, shape=10),lower = 0.0001, upper = Inf)#,optim.method="BFGS")

dgumbel(
# Escolher o menor valor
sort(gofstat(fitmaxE2)$aic)
sort(gofstat(fitmaxE2)$bic)

plot(fitmaxE2$lg)
plot(fitmaxE2$gu)
plot(fitmaxE2$ln)


## Estação 4
fitmaxE4<-list()

fitmaxE4$we <- fitdist(Qmaxhy$Q4, "weibull")
fitmaxE4$ga <- fitdist(Qmaxhy$Q4, "gamma")
fitmaxE4$ln <- fitdist(Qmaxhy$Q4, "lnorm")
fitmaxE4$no <- fitdist(Qmaxhy$Q4, "norm")
fitmaxE4$ca <- fitdist(Qmaxhy$Q4, "cauchy")
fitmaxE4$lo <- fitdist(Qmaxhy$Q4, "logis")
fitmaxE4$gu <- fitdist(Qmaxhy$Q4, "gumbel",start = list(a=10,b=10))



# Escolher o menor valor
sort(gofstat(fitmaxE4)$aic)
sort(gofstat(fitmaxE4)$bic)

plot(fitmaxE4$ln)

##=============================================================================##

##=============================================================================##
## Ajuste das Q7


## Estação 1
fitQ7E1<-list()

fitQ7E1$we <- fitdist(Q.7$Q1, "weibull")
fitQ7E1$ga <- fitdist(Q.7$Q1, "gamma")
fitQ7E1$ln <- fitdist(Q.7$Q1, "lnorm")
fitQ7E1$no <- fitdist(Q.7$Q1, "norm")
fitQ7E1$ca <- fitdist(Q.7$Q1, "cauchy")
fitQ7E1$lo <- fitdist(Q.7$Q1, "logis")
fitQ7E1$gu <- fitdist(Q.7$Q1, "gumbel",start = list(a=10,b=10))


# Escolher o menor valor
sort(gofstat(fitQ7E1)$aic)
sort(gofstat(fitQ7E1)$bic)

plot(fitQ7E1$we)
summary(fitQ7E1$we)

## Estação 2
fitE2<-list()

fitE2$we <- fitdist(Q.7$Q2, "weibull")
fitE2$ga <- fitdist(Q.7$Q2, "gamma")
fitE2$ln <- fitdist(Q.7$Q2, "lnorm")
fitE2$no <- fitdist(Q.7$Q2, "norm")
fitE2$ca <- fitdist(Q.7$Q2, "cauchy")
fitE2$lo <- fitdist(Q.7$Q2, "logis")
fitE2$gu <- fitdist(Q.7$Q2, "gumbel",start = list(a=10,b=10))


# Escolher o menor valor
sort(gofstat(fitE2)$aic)
sort(gofstat(fitE2)$bic)

plot(fitE2$gu)
plot(fitE2$ln)


## Estação 4
fitE4<-list()

fitE4$we <- fitdist(Q.7$Q4, "weibull")
fitE4$ga <- fitdist(Q.7$Q4, "gamma")
fitE4$ln <- fitdist(Q.7$Q4, "lnorm")
fitE4$no <- fitdist(Q.7$Q4, "norm")
fitE4$ca <- fitdist(Q.7$Q4, "cauchy")
fitE4$lo <- fitdist(Q.7$Q4, "logis")
fitE4$gu <- fitdist(Q.7$Q4, "gumbel",start = list(a=10,b=10))


# Escolher o menor valor
sort(gofstat(fitE4)$aic)
sort(gofstat(fitE4)$bic)

plot(fitE4$ln)

##=============================================================================##







cdfcomp(fitE1, legendtext=c("Weibull", "gamma", "lognormal","normal","cauchy","exponencial","logistica","geométrica"))
denscomp(fitE1, legendtext=c("Weibull", "gamma", "lognormal","normal","cauchy","exponencial","logistica","geométrica"))
qqcomp(fitE1, legendtext=c("Weibull", "gamma", "lognormal","normal","cauchy","exponencial","logistica","geométrica"))
ppcomp(fitE1, legendtext=c("Weibull", "gamma", "lognormal","normal","cauchy","exponencial","logistica","geométrica"))

## Ajuste das vazões mínimas Q7

fitW<-list()









Qt<-cbind(Qesp[-c(1:(366)),1],Qesp[,2:4])
colnames(Qt)<-c("Qt1","Qt2","Qt3","Qt4")
names(Qt)

attach(as.data.frame(Qt))
hist(log(Qt4))
dwi(Qt)
f <- ~Qt1+Qt2+Qt3+Qt4
i <- mnimput(f,Qt,eps=1e-3,ts=TRUE, method="spline",sp.control=list(df=c(7,7,7,7,7)))
plot(i)
Qtimp<-zoo(predict(i),data)

length(Qt[,1])-366
par(new=T,mgp=c(2.5,1,0),mfrow=c(2,1))
plot(Qesp.n[1:(366),1])
lines(Qtimp[1:(366),1],col="red")
plot(log(Qesp[1:(366),1]),Qtimp[1:(366),1])
abline(0,1,lty=2)
d(Qtimp[1:(366),1],Qesp[1:(366),1])
cor(Qtimp[1:(366),1],Qesp[1:(366),1])
ggof(Qtimp[1:(366),1]*areas[1],Qesp[1:(366),1]*areas[1])


# 3.0 HIDROLOGIA 
# Q7
Q.7 <- rollapply(data=Q, width=7, FUN=mean, fill=NA, partial= TRUE,align="left")
(Q_7<-daily2annual(Q.7, FUN=min, na.rm=T))
fdc(Q_7, lQ.thr="",hQ.thr="", plot=TRUE,thr.shw=F,log="", main="Q7 Anual", xlab="% Permanência",mgp=c(2.5,1,0), ylab=expression(Q~(m^3.~s^{-1})))
data<-seq(as.Date("1960/1/1"),as.Date("2013/12/31"),"year")
plot.xts(Q_7[,1])
lines(lowess(data,Q_7[,1]), col="blue", lwd=2)

for (i in 1:4){
plot(Q_7[,i])
lines(lowess(data,Q_7[,i]), col="blue", lwd=2)
}

for (i in 1:4){
  plot(Q.7[,i])
  lines(lowess(data,Q.7[,i]), col="blue", lwd=2)
}

for (i in 1:4){
acf(as.data.frame(Q_7[,i]))
pacf(as.data.frame(Q_7[,i]))
}

# Q7 sazonal

Q7.DJF<-dm2seasonal(Q.7, season="DJF",FUN=min, na.rm=TRUE)
Q7.MAM<-dm2seasonal(Q.7, season="MAM",FUN=min, na.rm=TRUE)
Q7.JJA<-dm2seasonal(Q.7, season="JJA",FUN=min, na.rm=TRUE)
Q7.SON<-dm2seasonal(Q.7, season="SON",FUN=min, na.rm=TRUE)

# Vazão Máxima Qmax
Qmax<-daily2annual(Q, FUN=max, na.rm=F)
plot(Qmax)
for (i in 1:4){
  plot(Qmax[,i])
  lines(lowess(data,Qmax[,i]), col="blue", lwd=2)
}
for (i in 1:4){
  acf(as.data.frame(Qmax[,i]))
  pacf(as.data.frame(Qmax[,i]))
}
Qmax

# Qmax sazonal
Qmax.DJF<-dm2seasonal(Q, season="DJF",FUN=max, na.rm=TRUE)
Qmax.MAM<-dm2seasonal(Q, season="MAM",FUN=max, na.rm=TRUE)
Qmax.JJA<-dm2seasonal(Q, season="DJF",FUN=max, na.rm=TRUE)
Qmax.SON<-dm2seasonal(Q, season="SON",FUN=max, na.rm=TRUE)

# Vazão Mínima Qmin
Qmin<-daily2annual(Q, FUN=min, na.rm=F)
plot(Qmin)

# Qmin sazonal
Qmin.DJF<-dm2seasonal(Q, season="DJF",FUN=min, na.rm=TRUE)
Qmin.MAM<-dm2seasonal(Q, season="MAM",FUN=min, na.rm=TRUE)
Qmin.JJA<-dm2seasonal(Q, season="DJF",FUN=min, na.rm=TRUE)
Qmin.SON<-dm2seasonal(Q, season="SON",FUN=min, na.rm=TRUE)

# Vazão Média Qmed
Qmed<-daily2annual(Q, FUN=mean, na.rm=F)
plot(Qmed)

# Qmed sazonal
Qmed.DJF<-dm2seasonal(Q, season="DJF",FUN=mean, na.rm=TRUE)
Qmed.MAM<-dm2seasonal(Q, season="MAM",FUN=mean, na.rm=TRUE)
Qmed.JJA<-dm2seasonal(Q, season="DJF",FUN=mean, na.rm=TRUE)
Qmed.SON<-dm2seasonal(Q, season="SON",FUN=mean, na.rm=TRUE)
plot(Qmed.DJF)

# TESTE DE TENDENCIA

#Q.7
MannKendall(Q_7[,1])
MannKendall(Q_7[,2])
MannKendall(Q_7[,3])
MannKendall(Q_7[,4])

#Qmax
MannKendall(Qmax[,1])
MannKendall(Qmax[,2])
MannKendall(Qmax[,3])
MannKendall(Qmax[,4])

#Qmin
MannKendall(Qmin[,1])
MannKendall(Qmin[,2])
MannKendall(Qmin[,3])
MannKendall(Qmin[,4])

#Qmed
MannKendall(Qmed[,1])
MannKendall(Qmed[,2])
MannKendall(Qmed[,3])
MannKendall(Qmed[,4])

# Tendencias sazonais

# Q7.DJF
MannKendall(Q7.DJF[,1])
MannKendall(Q7.DJF[,2])
MannKendall(Q7.DJF[,3])
MannKendall(Q7.DJF[,4])

# Q7.MAM
MannKendall(Q7.MAM[,1])
MannKendall(Q7.MAM[,2])
MannKendall(Q7.MAM[,3])
MannKendall(Q7.MAM[,4])

# Q7.JJA
MannKendall(Q7.JJA[,1])
MannKendall(Q7.JJA[,2])
MannKendall(Q7.JJA[,3])
MannKendall(Q7.JJA[,4])

# Q7.SON
MannKendall(Q7.SON[,1])
MannKendall(Q7.SON[,2])
MannKendall(Q7.SON[,3])
MannKendall(Q7.SON[,4])

# Qmax.DJF
MannKendall(Qmax.DJF[,1])
MannKendall(Qmax.DJF[,2])
MannKendall(Qmax.DJF[,3])
MannKendall(Qmax.DJF[,4])


# Qmax.MAM
MannKendall(Qmax.MAM[,1])
MannKendall(Qmax.MAM[,2])
MannKendall(Qmax.MAM[,3])
MannKendall(Qmax.MAM[,4])

# Qmax.JJA
MannKendall(Qmax.JJA[,1])
MannKendall(Qmax.JJA[,2])
MannKendall(Qmax.JJA[,3])
MannKendall(Qmax.JJA[,4])

# Qmax.SON
MannKendall(Qmax.SON[,1])
MannKendall(Qmax.SON[,2])
MannKendall(Qmax.SON[,3])
MannKendall(Qmax.SON[,4])


# Qmin.DJF
MannKendall(Qmin.DJF[,2])
MannKendall(Qmin.DJF[,2])
MannKendall(Qmin.DJF[,3])
MannKendall(Qmin.DJF[,4])

# Qmin.MAM
MannKendall(Qmin.MAM[,1])
MannKendall(Qmin.MAM[,2])
MannKendall(Qmin.MAM[,3])
MannKendall(Qmin.MAM[,4])

# Qmin.JJA
MannKendall(Qmin.JJA[,1])
MannKendall(Qmin.JJA[,2])
MannKendall(Qmin.JJA[,3])
MannKendall(Qmin.JJA[,4])

# Qmin.SON
MannKendall(Qmin.SON[,1])
MannKendall(Qmin.SON[,2])
MannKendall(Qmin.SON[,3])
MannKendall(Qmin.SON[,4])

# Qmed.DJF
MannKendall(Qmed.DJF[,2])
MannKendall(Qmed.DJF[,2])
MannKendall(Qmed.DJF[,3])
MannKendall(Qmed.DJF[,4])

# Qmed.MAM
MannKendall(Qmed.MAM[,1])
MannKendall(Qmed.MAM[,2])
MannKendall(Qmed.MAM[,3])
MannKendall(Qmed.MAM[,4])

# Qmed.JJA
MannKendall(Qmed.JJA[,1])
MannKendall(Qmed.JJA[,2])
MannKendall(Qmed.JJA[,3])
MannKendall(Qmed.JJA[,4])

# Qmed.SON
MannKendall(Qmed.SON[,1])
MannKendall(Qmed.SON[,2])
MannKendall(Qmed.SON[,3])
MannKendall(Qmed.SON[,4])



## TABELAS

Table1<-data.frame("DJF"=Q7.10.DJF,"MAM"=Q7.10.MAM,"JJA"=Q7.10.JJA,"SON"=Q7.10.SON,"Anual"=Q7.10.ANO)
Table2<-data.frame("DJF"=Qp.DJF,"MAM"=Qp.MAM,"JJA"=Qp.JJA,"SON"=Qp.SON, "Anual"=Qp)

data(PrecipGL)
plot(PrecipGL)
lines(lowess(time(PrecipGL),PrecipGL),lwd=3, col=2)
acf(PrecipGL)
pacf(PrecipGL)
t<-MannKendall(PrecipGL)
print(t)
summary(t)
library(boot)
data(PrecipGL)
MKtau<-function(z) MannKendall(z)$tau
tsboot(PrecipGL, MKtau, R=500, l=5, sim="fixed")
PrecipGL

library(boot) 
data(manaus)

library(fitdistrplus) 


 (1) fit of a gamma distribution by maximum likelihood estimation
#
data(groundbeef)
(serving <- groundbeef$serving)
fitg <- fitdist(serving, "gamma")
summary(fitg)
plot(fitg)
plot(fitg, demp = TRUE)
plot(fitg, histo = FALSE, demp = TRUE)
cdfcomp(fitg, addlegend=FALSE)
denscomp(fitg, addlegend=FALSE)
ppcomp(fitg, addlegend=FALSE)
qqcomp(fitg, addlegend=FALSE)
# (2) use the moment matching estimation (using a closed formula)
#
fitgmme <- fitdist(serving, "gamma", method="mme")
summary(fitgmme)
# (3) Comparison of various fits
#
fit<-list()


fit$W <- fitdist(serving, "weibull")
fit$g <- fitdist(serving, "gamma")
fit$ln <- fitdist(serving, "lnorm")



summary(fitW)
summary(fitg)
summary(fitln)

cdfcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
denscomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
qqcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
ppcomp(list(fitW, fitg, fitln), legendtext=c("Weibull", "gamma", "lognormal"))
gofstat(list(fitW, fitg, fitln), fitnames=c("Weibull", "gamma", "lognormal"))

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

fitgumbel <- fitdist(serving, "gumbel", start=list(a=mean(serving), b=var(serving)))
summary(fitgumbel)
plot(fitgumbel)
attributes(qgumbel)
attributes(fitgumbel)
qgumbel(c(0.9,0.8),fitgumbel$estimate[1],fitgumbel$estimate[2])

fun <- function(x) x^2 + x - 1                   # create function
curve(fun, xlim=c(-2, 1))                        # plot f unction
( res <- optimize(qwakeby, interval=c(0.001, 0.999)) )
points(res$minimum, res$objective)

v
fitexp <- fitdist(serving, "exp")
plot(fitexp)
summary(fitexp)
lower=c(-Inf, 0.01, 0.01, 0.01), method="L-BFGS-B"

gumbel <- function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
mledist(serving,"exponential",start=list(a=0.0135,b=0),lower=c(-Inf,0.01),optim.method ="L-BFGS-B")
mledist(serving,"gumbel",start=list(a=0,b=10),optim.method ="L-BFGS-B" )



write.csv(serving,file="serving.csv")
dir()

plot(fit.peixe$E2)
cdfcomp(fit.peixe$E4, addlegend=FALSE)
denscomp(fit.peixe$E4, addlegend=FALSE)
ppcomp(fit.peixe$E4, addlegend=FALSE)
qqcomp(fit.peixe$E2, addlegend=FALSE)

fit.peixe

plotdist(fit.peixe$E1,discrete=TRUE)

sd(log(Qesp_peixe[,1]))
mean(log(Qesp_peixe[,1]))

ad.test(Qesp_peixe[,7],plnorm,mean=mean(log(Qesp_peixe[,7])),sd=sd(log(Qesp_peixe[,7])))



gofstat(fit.peixe$E7)$ad

plot(fit.peixe$E7)
Dagum()
ks.test(fit.peixe$E1)
plot(fitE1bs)
summary(fitE1bs)
param <- data.frame(fitE1bs$estimate,fitE2bs$estimate,fitE3bs$estimate,fitE4bs$estimate)

plot(param)
plnorm
qlnorm(1-0.98,fitE1bs$estimate[1],fitE1bs$estimate[2])*areas[1]

Qp <- function(Estac,trgt) {
    Qp <- numeric(0)
    for(i in 1:ncol(Estac)){
    Qp[i]<-quantile(extractzoo(Estac[,i],trgt=trgt),prob=seq(0,1,by=0.05),na.rm=T)
}
    return(Qp)
}

Qp(Peixe,"DJF")

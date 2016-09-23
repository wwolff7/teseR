# 1.0 CARREGAR LIVRARIA/PACOTES

library(hydroTSM)
library(mtsdi)
library(hydroGOF)
library(extrafont)
loadfonts()


##rm(list = ls())
citation("mtsdi")
options(OutDec=",",digits = 7)

# 2.0 ABRIR ARQUIVO/BANCO DE DADOS

getwd()
setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Preen")

Peixe<-read.zoo("Peixe-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.peixe<-seq(as.Date("1985/01/01"),as.Date("2014/08/31"),"day")
areas.peixe<-c(420,3710,2010,801,3660,2800,1620)
head(Peixe)
summary(Peixe)
hydropairs(as.data.frame(Peixe),met="spear")

Pelotas<-read.zoo("Pelotas-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-6]
data.pelotas<-seq(as.Date("1977/01/01"),as.Date("2014/07/31"),"day")
areas.pelotas<-c(550,1170,2820,533,1820,1120)[,-6]
head(Pelotas)
summary(Pelotas)
hydropairs(as.data.frame(Pelotas),met="spear")

Canoas<-read.zoo("Canoas-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.canoas<-seq(as.Date("1986/01/01"),as.Date("2014/07/31"),"day")
areas.canoas<-c(10000,3680,4610,3230,2000,489,1010)
head(Canoas)
summary(Canoas)
hydropairs(as.data.frame(Canoas),met="spear")

Itajai<-read.zoo("Itajaí-Açú-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.itajai<-seq(as.Date("1986/01/01"),as.Date("2014/07/31"),"day")
areas.itajai<-c(536,1430,648,717,11803,827,1240,3330,9850,1650,104,5160,286,434,1570,1600,397,9790)
head(Itajai)
summary(Itajai)
hydropairs(as.data.frame(Itajai),met="spear")

Canoinhas<-read.zoo("Canoinhas-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.Canoinhas<-seq(as.Date("1986/01/01"),as.Date("2014/07/31"),"day")
areas.Canoinhas<-772
head(Canoinhas)
summary(Canoinhas)
hydropairs(as.data.frame(Canoinhas),met="spear")

Antas<-read.zoo("Antas-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-1]
data.Antas<-seq(as.Date("1976/01/01"),as.Date("2010/12/31"),"day")
areas.Antas<-c(272,609,300,61900,5340)[-1]
head(Antas)
summary(Antas)
hydropairs(as.data.frame(Antas),met="spear")

Negro <-read.zoo("Negro-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.Negro<-seq(as.Date("1977/01/01"),as.Date("2010/12/31"),"day")
areas.Negro<-c(960,2610,605)
head(Negro)
summary(Negro)
hydropairs(as.data.frame(Negro),met="spear")

Ararangua <-read.zoo("Araranguá-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-1]
data.Ararangua<-seq(as.Date("1985/01/01"),as.Date("2004/12/31"),"day")
areas.Ararangua<-c(863,526,355,119,359)[-1]
head(Ararangua)
summary(Ararangua)
hydropairs(as.data.frame(Ararangua),met="spear")

Peperi <-read.zoo("Peperi-Guaçu-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.Peperi<-seq(as.Date("1965/01/01"),as.Date("2014/12/31"),"day")
areas.Peperi<-c(2020,1540)
head(Peperi)
summary(Peperi)
hydropairs(as.data.frame(Peperi),met="spear")

Tubarão <-read.zoo("Tubarão-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-1]
data.Tubarão<-seq(as.Date("1987/01/01"),as.Date("2004/12/31"),"day")
areas.Tubarão<-c(770,1515,380,599,822,2740,379,676,1690,620)[-1]
head(Tubarão)
summary(Tubarão)
hydropairs(as.data.frame(Tubarão),met="spear")

CubataoN <-read.zoo("Cubatao-Norte-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.CubataoN<-seq(as.Date("1986/01/01"),as.Date("2009/31/31"),"day")
areas.CubataoN <-c(374)
plot(CubataoN)
summary(CubataoN)
hydropairs(as.data.frame(CubataoN),met="spear")
plot(CubataoN)

Itapocu <-read.zoo("Itapocu-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-5]
data.Itapocu<-seq(as.Date("1978/01/01"),as.Date("2001/12/31"),"day")
areas.Itapocu <-c(182,794,281,358,30,392)[-5]
head(Itapocu)
summary(Itapocu)
hydropairs(as.data.frame(Itapocu),met="spear")
plot(Itapocu)

CubataoS <-read.zoo("Cubatao-Sul-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-1]
CubataoS <- window(CubataoS,end=as.Date("2000/12/31"))
data.CubataoS<-seq(as.Date("1951/01/01"),as.Date("2000/12/31"),"day")
areas.CubataoS <-c(522,400)[-1]
head(CubataoS)

summary(CubataoS)
hydropairs(as.data.frame(CubataoS),met="spear")
plot(CubataoS)

Irani <-read.zoo("Irani-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.Irani<-seq(as.Date("1970/01/01"),as.Date("2014/12/31"),"day")
areas.Irani <-c(1500,654,933)
head(Irani)
summary(Irani)
hydropairs(as.data.frame(Irani),met="spear")

Mampituba <-read.zoo("Mampituba-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.Mampituba<-seq(as.Date("1986/01/01"),as.Date("2006/12/31"),"day")
areas.Mampituba <-c(339)
head(Mampituba)
summary(Mampituba)
hydropairs(as.data.frame(Mampituba),met="spear")
plot(Mampituba)

Tijucas <-read.zoo("Tijucas-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
data.Tijucas<-seq(as.Date("1983/01/01"),as.Date("2004/12/31"),"day")
areas.Tijucas <-c(1042,598,1964)
head(Tijucas)
summary(Tijucas)
hydropairs(as.data.frame(Tijucas),met="spear")

Chapeco <-read.zoo("Chapecó-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)[,-1]
data.Chapeco<-seq(as.Date("1979/01/01"),as.Date("2010/12/31"),"day")
areas.Chapeco <-c(65,1660,418,5550,266,1010,642,740,8240,1840)[-1]
head(Chapeco)
summary(Chapeco)
hydropairs(as.data.frame(Chapeco),met="spear")

Todas <- cbind(Canoas,Pelotas,Peixe,Itajai,Itapocu,Canoinhas,Antas,Negro,Ararangua,
               Peperi,Tubarão,CubataoN,CubataoS,Irani,Mampituba,Tijucas,Chapeco)

dim(Todas)
names(Todas)[c(43,66,67,71)] <- c("E65180000","E84100000","E82270050","E84970000")

pdf("Figuras/Estat_dados.pdf",onefile = T, width=27/2.54, height=35/2.54,paper = "special",colormodel="grey",family
    = "CM Roman")

#par(mar=c(2.3,2.3,.5,.5))
mat <- matrixplot(dwi(Todas, var.type="Days"),ColorRamp="Temperature",aspect=1.5,
           colorkey=list(labels=list(cex=1.5)))

str(mat)

data <- as.character(seq(1929,2009,5))
space <- matrix(rep("",4*length(data)),nrow=length(data),ncol=4)

mat$x.scales$labels <- c(as.vector(t(cbind(data,space))),2014)
mat$x.scales$cex <- c(1.2,1.2)
mat$y.scales$cex <- c(1.2,1.2)

plot(mat)

dev.off()
embed_fonts("Figuras/Estat_dados.pdf",outfile = "Figuras/Estat_dados.pdf")


par(mar=c(2.3,2.3,.5,.5))

postscript("Figuras/Brutos_Itajai.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Itajai, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio do Itajaí-Açú")

dev.off()

postscript("Figuras/Brutos_Canoas.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Canoas, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Canoas")

dev.off()

matrixplot(dwi(Canoinhas, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Canoinhas")

postscript("Figuras/Brutos_Peixe.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Peixe, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio do Peixe")

dev.off()

postscript("Figuras/Brutos_Pelotas.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Pelotas, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Pelotas")

dev.off()

postscript("Figuras/Brutos_Negro.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Negro, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Negro")

dev.off()
postscript("Figuras/Brutos_Antas.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Antas, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio das Antas")

dev.off()
postscript("Figuras/Brutos_Ararangua.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Ararangua, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Araranguá")

dev.off()
postscript("Figuras/Brutos_Peperi.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Peperi, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Peperi-Guaçú")

dev.off()
postscript("Figuras/Brutos_Tubarão.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Tubarão, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Tubarão")

dev.off()
postscript("Figuras/Brutos_Itapocu.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Itapocu, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Itapocu")

dev.off()
postscript("Figuras/Brutos_CubataoS.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(CubataoS, var.type="Days"),ColorRamp="Days",main = "Bacia do Cubatao Sul")

dev.off()
postscript("Figuras/Brutos_Irani.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Irani, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Irani")

dev.off()
postscript("Figuras/Brutos_Tijucas.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Tijucas, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Tijucas")

dev.off()
postscript("Figuras/Brutos_Chapeco.eps",onefile = T,horizontal = T, width=30/2.54, height=10/2.54,paper = "special")

matrixplot(dwi(Chapeco, var.type="Days"),ColorRamp="Days",main = "Bacia do Rio Chapecó")

dev.off()

Pelotas1 <- window(Pelotas, start=as.Date("1985-01-01"))
Canoas1 <- window(Canoas, start=as.Date("1985-01-01"))
Canoinhas1 <- window(Canoinhas, start=as.Date("1985-01-01"))
Peixe1 <- window(Peixe, start=as.Date("1985-01-01"))
Itajai1 <- window(Itajai, start=as.Date("1985-01-01"),end=as.Date("2004-12-31"))
Itajai0 <- window(Itajai, end=as.Date("2005-12-31"))


## 3.0 ANÁLISE EXPLORATÓRIA  E PREENCHIMENTO DE FALHAS

plot(Peixe1,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Pelotas1,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Canoas1,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Canoinhas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Itajai1,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)
plot(Itajai0,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="",cex=3)

## Itajaí-Açú - Preenchimento 1 

Qesp_itajai7<-Itajai[,7]/areas.itajai[7]  
Qesp_itajai8<-Itajai[,8]/areas.itajai[8]
Qesp_itajai9<-Itajai[,9]/areas.itajai[9]
Qesp_itajai10<-Itajai[,10]/areas.itajai[10]  

length(areas.itajai)
dim(Itajai)

Qesp_itajai<-window(cbind.zoo(Qesp_itajai7,Qesp_itajai8,Qesp_itajai9,Qesp_itajai10),end=as.Date("2004-12-31"))

matrixplot(dwi(Qesp_itajai, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio do Itajaí-Açú")
matrixplot(dwi(Itajai, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio do Itajaí-Açú")
hydropairs(as.data.frame(Qesp_itajai))
           
summary(Qesp_itajai)
colnames(Qesp_itajai)<-colnames(Itajai)
names(Qesp_itajai)

plot(Qesp_itajai,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")


## Preenchimento de falhas
f <- ~Qesp_itajai7+Qesp_itajai8+Qesp_itajai9+Qesp_itajai10
i <- mnimput(f,Qesp_itajai,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))
## plot(i)

Qimp_itajai<-zoo(predict(i),seq(as.Date("1929/01/01"),as.Date("2009/12/31"),"day"))
##write.csv(Qimp_peixe,"Q_peixe.csv")

colnames(Qimp_itajai)
plot(Qimp_itajai)
summary(Qimp_itajai)

names(Itajai)

Itajai[,7] <-Qimp_itajai[,1]*areas.itajai[7]
Itajai[,8] <-Qimp_itajai[,2]*areas.itajai[8]
Itajai[,9] <-Qimp_itajai[,3]*areas.itajai[9]  
Itajai[,10] <-Qimp_itajai[,4]*areas.itajai[10]
##-----------------------------------------------------------------------------

## Calibrar 
Qesp_itajai1<-window(cbind.zoo(Qesp_itajai7,Qesp_itajai8,Qesp_itajai9,Qesp_itajai10),start = as.Date("1931-01-01"),end=as.Date("1963-12-31"))

hydropairs(as.data.frame(Qesp_itajai),met="spear")
plot(Qesp_itajai1)
summary(Qesp_itajai1)
dim(Qesp_itajai1)
names(Qesp_itajai1)
areas.itajai[7:10]

fun.imput <- splinefun(Qesp_itajai1[,2],Qesp_itajai1[,1])


##----------------------------------------------------------------------------- 

notNA   <- which(!is.na(Qesp_itajai1), arr.ind=T)
set.seed(333); sel.pos <- notNA[sample(nrow(notNA), 300),]
obs0 <- as.data.frame(Qesp_itajai1[sel.pos[,1],sel.pos[,2]])
head(sel.pos[,1])

Qesp_itajai2 <- Qesp_itajai1
Qesp_itajai2[sel.pos[,1],sel.pos[,2]] <- NA

f <- ~Qesp_itajai2[,1]+Qesp_itajai2[,2]+Qesp_itajai2[,3]+Qesp_itajai2[,4]

i0 <- mnimput(f,Qesp_itajai2,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))

(sim0 <- data.frame(predict(i0))[sel.pos[,1],sel.pos[,2]])

ggof(sim0,obs0)
hydroPSO()
fun.cal <- splinefun(sim0,obs0)
sim0cal <- fun.cal01(sim0)

ggof(sim0cal,obs0)

d(sim0,obs0)*cor(sim0,obs0)

##-----------------------------------------------------------------------------
Qesp_itajai1<-window(cbind.zoo(Qesp_itajai7,Qesp_itajai8,Qesp_itajai9,Qesp_itajai10),start = as.Date("1931-01-01"),end=as.Date("1963-12-31"))

obs0.1 <- Qesp_itajai1[c(2000:(2000+(3*365))),1]
Qesp_itajai1[c(2000:(2000+(3*365))),1] = NA
summary(Qesp_itajai1)
dim(Qesp_itajai1)
i0.1 <- mnimput(f,Qesp_itajai1,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))

(sim0.1 <- data.frame(predict(i0.1))[c(2000:(2000+(3*365))),1])

length(sim0.1)

ggof(sim0.1,obs0.1)

d(sim0.1,obs0.1)*cor(sim0.1,obs0.1)

fun.cal01 <- hydropso()approxfun(sim0.1,obs0.1)
sim0.1cal <- fun.cal01(sim0.1)

ggof(sim0.1cal,obs0.1)



##-----------------------------------------------------------------------------
Qesp_itajai1<-window(cbind.zoo(Qesp_itajai7,Qesp_itajai8,Qesp_itajai9,Qesp_itajai10),start = as.Date("1931-01-01"),end=as.Date("1963-12-31"))

obs2 <- Qesp_itajai1[c(300:(300+1*365)),1]
Qesp_itajai1[c(300:(300+1*365)),1] = NA

i2 <- mnimput(f,Qesp_itajai1,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))

(sim2 <- data.frame(predict(i2))[c(300:(300+1*365)),1])

ggof(sim2,obs2)
d(sim2,obs2)*cor(sim2,obs2)

##-----------------------------------------------------------------------------
Qesp_itajai1<-window(cbind.zoo(Qesp_itajai7,Qesp_itajai8,Qesp_itajai9,Qesp_itajai10),start = as.Date("1931-01-01"),end=as.Date("1963-12-31"))

obs3 <- Qesp_itajai1[c(300:(300+2*365)),1]
Qesp_itajai1[c(300:(300+2*365)),1] = NA

i3 <- mnimput(f,Qesp_itajai1,eps=1e-3,ts=F,log=T, method="arima")#,sp.control=list(df=c(7,7,7,7)))

(sim3 <- data.frame(predict(i3))[c(300:(300+2*365)),1])

ggof(sim3,obs3)
d(sim3,obs3)*cor(sim3,obs3)

##-----------------------------------------------------------------------------

obs4 <- Qesp_itajai1[c(300:(300+3*365)),1]
obs5 <- Qesp_itajai1[c(300:(300+4*365)),1]
obs6 <- Qesp_itajai1[c(300:(300+5*365)),1]


summary(i)
attributes(i)
plot(predict(i))
i$filled.dataset

plot(i)



## Itajaí-Açú - Preenchimento 2

Itajai <- window(Itajai,start=as.Date("1985-01-01"),end=as.Date("2004-12-31"))
## rm(Itajai)


Qesp_itajai1<-Itajai[,1]/areas.itajai[1]
Qesp_itajai2<-Itajai[,2]/areas.itajai[2]
Qesp_itajai3<-Itajai[,3]/areas.itajai[3]  
Qesp_itajai4<-Itajai[,4]/areas.itajai[4]
Qesp_itajai5<-Itajai[,5]/areas.itajai[5]  
Qesp_itajai6<-Itajai[,6]/areas.itajai[6]
Qesp_itajai7<-Itajai[,7]/areas.itajai[7]  
Qesp_itajai8<-Itajai[,8]/areas.itajai[8]
Qesp_itajai9<-Itajai[,9]/areas.itajai[9]
Qesp_itajai10<-Itajai[,10]/areas.itajai[10]  
Qesp_itajai11<-Itajai[,11]/areas.itajai[11]
Qesp_itajai12<-Itajai[,12]/areas.itajai[12]  
Qesp_itajai13<-Itajai[,13]/areas.itajai[13]
Qesp_itajai14<-Itajai[,14]/areas.itajai[14]  
Qesp_itajai15<-Itajai[,15]/areas.itajai[15]  
Qesp_itajai16<-Itajai[,16]/areas.itajai[16]
Qesp_itajai17<-Itajai[,17]/areas.itajai[17]  
Qesp_itajai18<-Itajai[,18]/areas.itajai[18]  


length(areas.itajai)
dim(Itajai)

Qesp_itajai<-cbind.zoo(Qesp_itajai1,Qesp_itajai2,Qesp_itajai3,Qesp_itajai4,Qesp_itajai5,Qesp_itajai6,Qesp_itajai7,Qesp_itajai8,Qesp_itajai9,Qesp_itajai10,Qesp_itajai11,Qesp_itajai12,Qesp_itajai13,Qesp_itajai14,Qesp_itajai15,Qesp_itajai16,Qesp_itajai17,Qesp_itajai18)
colnames(Qesp_itajai)<-colnames(Itajai)


matrixplot(dwi(Qesp_itajai, var.type="Days"),ColorRamp="Days",main = "Macrobacia do Rio do Itajaí-Açú")
hydropairs(as.data.frame(Qesp_itajai))

summary(Qesp_itajai)
colnames(Qesp_itajai)<-colnames(Itajai)
names(Qesp_itajai)

plot(Qesp_itajai,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
example(image3D)

## Preenchimento de falhas
f <- ~Qesp_itajai1+Qesp_itajai2+Qesp_itajai3+Qesp_itajai4+Qesp_itajai5+Qesp_itajai6+Qesp_itajai7+Qesp_itajai8+Qesp_itajai9+Qesp_itajai10+Qesp_itajai11+Qesp_itajai12+Qesp_itajai13+Qesp_itajai14+Qesp_itajai15+Qesp_itajai16+Qesp_itajai17+Qesp_itajai18
i <- mnimput(f,Qesp_itajai,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_itajai<-zoo(p,seq(as.Date("1985/01/01"),as.Date("2004/12/31"),"day"));colnames(Qimp_itajai)<-colnames(Itajai)
write.csv(Qimp_itajai,"Estac.Preen/Qimp_itajai.csv")
dir()

plot(Qimp_itajai[,7])
summary(Qimp_itajai)
colnames(Itajai)

names(Itajai)

Itajai[,1] <-Qimp_itajai[,1]*areas.itajai[1]
Itajai[,2] <-Qimp_itajai[,2]*areas.itajai[2]
Itajai[,3] <-Qimp_itajai[,3]*areas.itajai[3]  
Itajai[,4] <-Qimp_itajai[,4]*areas.itajai[4]
Itajai[,5] <-Qimp_itajai[,5]*areas.itajai[5]
Itajai[,6] <-Qimp_itajai[,6]*areas.itajai[6]
Itajai[,7] <-Qimp_itajai[,7]*areas.itajai[7]  
Itajai[,8] <-Qimp_itajai[,8]*areas.itajai[8]
Itajai[,9] <-Qimp_itajai[,9]*areas.itajai[9]
Itajai[,10] <-Qimp_itajai[,10]*areas.itajai[10]
Itajai[,11] <-Qimp_itajai[,11]*areas.itajai[11]  
Itajai[,12] <-Qimp_itajai[,12]*areas.itajai[12]
Itajai[,13] <-Qimp_itajai[,13]*areas.itajai[13]
Itajai[,14] <-Qimp_itajai[,14]*areas.itajai[14]
Itajai[,15] <-Qimp_itajai[,15]*areas.itajai[15]  
Itajai[,16] <-Qimp_itajai[,16]*areas.itajai[16]
Itajai[,17] <-Qimp_itajai[,17]*areas.itajai[17]  
Itajai[,18] <-Qimp_itajai[,18]*areas.itajai[18]

plot(Itajai)

##-----------------------------------------------------------------------------##

## Peixe 

areas.peixe
summary(Peixe)
plot(Peixe,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Peixe, var.type="Days"),ColorRamp="Days",main="Peixe")
Peixe <- window(Peixe,start=as.Date("1985-01-01"))
## rm(Itajai)
rm(Peixe)


## Vazão específica

Qesp_peixe1<-Peixe[,1]/areas.peixe[1]
Qesp_peixe2<-Peixe[,2]/areas.peixe[2]
Qesp_peixe3<-Peixe[,3]/areas.peixe[3]  
Qesp_peixe4<-Peixe[,4]/areas.peixe[4]
Qesp_peixe5<-Peixe[,5]/areas.peixe[5]  
Qesp_peixe6<-Peixe[,6]/areas.peixe[6]
Qesp_peixe7<-Peixe[,7]/areas.peixe[7]  


Qesp_peixe<-cbind.zoo(Qesp_peixe1,Qesp_peixe2,Qesp_peixe3,Qesp_peixe4,Qesp_peixe5,Qesp_peixe6,Qesp_peixe7)
summary(Qesp_peixe)
colnames(Qesp_peixe)<-colnames(Peixe)
names(Qesp_peixe)

plot(Qesp_peixe,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_peixe1+Qesp_peixe2+Qesp_peixe3+Qesp_peixe4+Qesp_peixe5+Qesp_peixe6+Qesp_peixe7
i <- mnimput(f,Qesp_peixe,eps=1e-3,ts=T,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_peixe<-zoo(p,seq(as.Date("1985/01/01"),as.Date("2014/08/31"),"day"));colnames(Qimp_peixe)<-colnames(Peixe)
write.csv(Qimp_peixe,"Estac.Preen/Qimp_peixe.csv")

colnames(Qimp_peixe)
plot(Qimp_peixe[,2])
dwi(Qimp_peixe)
summary(Qimp_peixe)

Peixe[,1]<-Qimp_peixe[,1]*areas.peixe[1]
Peixe[,2]<-Qimp_peixe[,2]*areas.peixe[2]
Peixe[,3]<-Qimp_peixe[,3]*areas.peixe[3]  
Peixe[,4]<-Qimp_peixe[,4]*areas.peixe[4]
Peixe[,5]<-Qimp_peixe[,5]*areas.peixe[5]
Peixe[,6]<-Qimp_peixe[,6]*areas.peixe[6]  
Peixe[,7]<-Qimp_peixe[,7]*areas.peixe[7]

plot(Peixe)

matrixplot(dwi(Qimp_peixe, var.type="Days"),ColorRamp="Days")

##-----------------------------------------------------------------------------##


## Canoas 


summary(Canoas)
plot(Canoas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Canoas, var.type="Days"),ColorRamp="Days")
Canoas <- window(Canoas,start=as.Date("1985-01-01"))
## rm(Itajai)
rm(Canoas)


## Vazão específica

Qesp_canoas1<-Canoas[,1]/areas.canoas[1]
Qesp_canoas2<-Canoas[,2]/areas.canoas[2]
Qesp_canoas3<-Canoas[,3]/areas.canoas[3]  
Qesp_canoas4<-Canoas[,4]/areas.canoas[4]
Qesp_canoas5<-Canoas[,5]/areas.canoas[5]  
Qesp_canoas6<-Canoas[,6]/areas.canoas[6]
Qesp_canoas7<-Canoas[,7]/areas.canoas[7]  


Qesp_canoas<-cbind.zoo(Qesp_canoas1,Qesp_canoas2,Qesp_canoas3,Qesp_canoas4,Qesp_canoas5,Qesp_canoas6,Qesp_canoas7)
summary(Qesp_canoas)
colnames(Qesp_canoas)<-colnames(Canoas)
names(Qesp_canoas)

plot(Qesp_canoas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_canoas1+Qesp_canoas2+Qesp_canoas3+Qesp_canoas4+Qesp_canoas5+Qesp_canoas6+Qesp_canoas7
i <- mnimput(f,Qesp_canoas,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_canoas<-zoo(p,seq(as.Date("1985/01/01"),as.Date("2014/08/31"),"day"));colnames(Qimp_canoas)<-colnames(Canoas)
write.csv(Qimp_canoas,"Estac.Preen/Qimp_canoas.csv")


colnames(Qimp_canoas)
plot(Qimp_canoas[,2])
dwi(Qimp_canoas)
summary(Qimp_canoas)

Canoas[,1]<-Qimp_canoas[,1]*areas.canoas[1]
Canoas[,2]<-Qimp_canoas[,2]*areas.canoas[2]
Canoas[,3]<-Qimp_canoas[,3]*areas.canoas[3]  
Canoas[,4]<-Qimp_canoas[,4]*areas.canoas[4]
Canoas[,5]<-Qimp_canoas[,5]*areas.canoas[5]
Canoas[,6]<-Qimp_canoas[,6]*areas.canoas[6]  
Canoas[,7]<-Qimp_canoas[,7]*areas.canoas[7]

plot(Canoas)

##-----------------------------------------------------------------------------##


## Pelotas 

summary(Pelotas)
plot(Pelotas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Pelotas, var.type="Days"),ColorRamp="Days")
Pelotas <- window(Pelotas,start=as.Date("1977-01-01"))


## Vazão específica

Qesp_pelotas1<-Pelotas[,1]/areas.pelotas[1]
Qesp_pelotas2<-Pelotas[,2]/areas.pelotas[2]
Qesp_pelotas3<-Pelotas[,3]/areas.pelotas[3]  
Qesp_pelotas4<-Pelotas[,4]/areas.pelotas[4]
Qesp_pelotas5<-Pelotas[,5]/areas.pelotas[5]  
Qesp_pelotas6<-Pelotas[,6]/areas.pelotas[6]



Qesp_pelotas<-cbind.zoo(Qesp_pelotas1,Qesp_pelotas2,Qesp_pelotas3,Qesp_pelotas4,Qesp_pelotas5,Qesp_pelotas6)
summary(Qesp_pelotas)
colnames(Qesp_pelotas)<-colnames(Pelotas)
names(Qesp_pelotas)

plot(Qesp_pelotas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_pelotas1+Qesp_pelotas2+Qesp_pelotas3+Qesp_pelotas4+Qesp_pelotas5+Qesp_pelotas6
i <- mnimput(f,Qesp_pelotas,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_pelotas<-zoo(p,seq(as.Date("1977/01/01"),as.Date("2014/07/31"),"day"));colnames(Qimp_pelotas)<-colnames(Pelotas)
write.csv(Qimp_pelotas,"Estac.Preen/Qimp_pelotas.csv")

colnames(Qimp_pelotas)
plot(Qimp_pelotas[,2])
dwi(Qimp_pelotas)
summary(Qimp_pelotas)

Pelotas[,1]<-Qimp_pelotas[,1]*areas.pelotas[1]
Pelotas[,2]<-Qimp_pelotas[,2]*areas.pelotas[2]
Pelotas[,3]<-Qimp_pelotas[,3]*areas.pelotas[3]  
Pelotas[,4]<-Qimp_pelotas[,4]*areas.pelotas[4]
Pelotas[,5]<-Qimp_pelotas[,5]*areas.pelotas[5]
Pelotas[,6]<-Qimp_pelotas[,6]*areas.pelotas[6]  

plot(Pelotas)

##-----------------------------------------------------------------------------##

## Canoinhas 
Canoinhas<-read.zoo("Canoinhas-R.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
areas.Canoinhas<-772
summary(Canoinhas)
data.Canoinhas<-seq(as.Date("1981/01/01"),as.Date("2010/12/31"),"day")
plot(Canoinhas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

Canoinhas <- as.data.frame(Canoinhas)/areas.Canoinhas
Canoinhas[,1][is.na(Canoinhas[,1])] <- mean(Canoinhas[,1],na.rm = T)

Qimp_canoinhas<-zoo(Canoinhas,data.Canoinhas);colnames(Qimp_canoinhas) <- "E65180000" 
write.csv(Qimp_canoinhas,"Estac.Preen/Qimp_canoinhas.csv")

head(Canoinhas)
summary(Qimp_canoinhas)
plot(Qimp_canoinhas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

##-----------------------------------------------------------------------------##

## Negro 

summary(Negro)
plot(Negro,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Negro, var.type="Days"),ColorRamp="Days")
Negro <- window(Negro,start=as.Date("1977-01-01"),end=as.Date("2010-12-31"))


## Vazão específica

Qesp_negro1<-Negro[,1]/areas.Negro[1]
Qesp_negro2<-Negro[,2]/areas.Negro[2]
Qesp_negro3<-Negro[,3]/areas.Negro[3]  

Qesp_negro<-cbind.zoo(Qesp_negro1,Qesp_negro2,Qesp_negro3)
summary(Qesp_negro)
colnames(Qesp_negro)<-colnames(Negro)
names(Qesp_negro)

plot(Qesp_negro,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_negro1+Qesp_negro2+Qesp_negro3
i <- mnimput(f,Qesp_negro,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_negro<-zoo(p,seq(as.Date("1977/01/01"),as.Date("2010/12/31"),"day"));colnames(Qimp_negro)<-colnames(Negro)
write.csv(Qimp_negro,"Estac.Preen/Qimp_negro.csv")

colnames(Qimp_negro)
plot(Qimp_negro[,2])
dwi(Qimp_negro)
summary(Qimp_negro)

Negro[,1]<-Qimp_negro[,1]*areas.Negro[1]
Negro[,2]<-Qimp_negro[,2]*areas.Negro[2]
Negro[,3]<-Qimp_negro[,3]*areas.Negro[3]  

plot(Negro)

##-----------------------------------------------------------------------------##

## Antas

summary(Antas)
plot(Antas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Antas, var.type="Days"),ColorRamp="Days")
Antas<- window(Antas,start=as.Date("1976-01-01"),end=as.Date("2014-12-31"))


## Vazão específica

Qesp_antas1<-Antas[,1]/areas.Antas[1]
Qesp_antas2<-Antas[,2]/areas.Antas[2]
Qesp_antas3<-Antas[,3]/areas.Antas[3]  
Qesp_antas4<-Antas[,4]/areas.Antas[4]

Qesp_antas<-cbind.zoo(Qesp_antas1,Qesp_antas2,Qesp_antas3,Qesp_antas4)
summary(Qesp_antas)
colnames(Qesp_antas)<-colnames(Antas)
names(Qesp_antas)

plot(Qesp_antas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_antas1+Qesp_antas2+Qesp_antas3+Qesp_antas4
i <- mnimput(f,Qesp_antas,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_antas<-zoo(p,seq(as.Date("1976/01/01"),as.Date("2014/12/31"),"day"));colnames(Qimp_antas)<-colnames(Antas)
write.csv(Qimp_antas,"Estac.Preen/Qimp_antas.csv")

colnames(Qimp_antas)
plot(Qimp_antas[,2])
dwi(Qimp_antas)
summary(Qimp_antas)

Antas[,1]<-Qimp_antas[,1]*areas.Antas[1]
Antas[,2]<-Qimp_antas[,2]*areas.Antas[2]
Antas[,3]<-Qimp_antas[,3]*areas.Antas[3]  
Antas[,4]<-Qimp_antas[,4]*areas.Antas[4]  

plot(Antas)

##-----------------------------------------------------------------------------##

## Araranguá

summary(Ararangua)
plot(Ararangua,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Ararangua, var.type="Days"),ColorRamp="Days")
Ararangua<- window(Ararangua,start=as.Date("1985-01-01"),end=as.Date("2004-12-31"))


## Vazão específica

Qesp_ararangua1<-Ararangua[,1]/areas.Ararangua[1]
Qesp_ararangua2<-Ararangua[,2]/areas.Ararangua[2]
Qesp_ararangua3<-Ararangua[,3]/areas.Ararangua[3]  
Qesp_ararangua4<-Ararangua[,4]/areas.Ararangua[4]

Qesp_ararangua<-cbind.zoo(Qesp_ararangua1,Qesp_ararangua2,Qesp_ararangua3,Qesp_ararangua4)
summary(Qesp_ararangua)
colnames(Qesp_ararangua)<-colnames(Ararangua)
names(Qesp_ararangua)

plot(Qesp_ararangua,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_ararangua1+Qesp_ararangua2+Qesp_ararangua3+Qesp_ararangua4
i <- mnimput(f,Qesp_ararangua,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_ararangua<-zoo(p,seq(as.Date("1985/01/01"),as.Date("2004/12/31"),"day"));colnames(Qimp_ararangua)<-colnames(Ararangua)
write.csv(Qimp_ararangua,"Estac.Preen/Qimp_ararangua.csv")

colnames(Qimp_ararangua)
plot(Qimp_ararangua[,2])
dwi(Qimp_ararangua)
summary(Qimp_ararangua)

Ararangua[,1]<-Qimp_ararangua[,1]*areas.Ararangua[1]
Ararangua[,2]<-Qimp_ararangua[,2]*areas.Ararangua[2]
Ararangua[,3]<-Qimp_ararangua[,3]*areas.Ararangua[3]  
Ararangua[,4]<-Qimp_ararangua[,4]*areas.Ararangua[4]  

plot(Ararangua)

##-----------------------------------------------------------------------------##

## Peperi-guaçú

summary(Peperi)
plot(Peperi,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Peperi, var.type="Days"),ColorRamp="Days")
Peperi<- window(Peperi,start=as.Date("1965-01-01"),end=as.Date("2014-12-31"))


## Vazão específica

Qesp_peperi1<-Peperi[,1]/areas.Peperi[1]
Qesp_peperi2<-Peperi[,2]/areas.Peperi[2]

Qesp_peperi<-cbind.zoo(Qesp_peperi1,Qesp_peperi2)
summary(Qesp_peperi)
colnames(Qesp_peperi)<-colnames(Peperi)
names(Qesp_peperi)

plot(Qesp_peperi,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_peperi1+Qesp_peperi2
i <- mnimput(f,Qesp_peperi,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_peperi<-zoo(p,seq(as.Date("1965/01/01"),as.Date("2014/12/31"),"day"));colnames(Qimp_peperi)<-colnames(Peperi)
write.csv(Qimp_peperi,"Estac.Preen/Qimp_peperi.csv")

colnames(Qimp_peperi)
plot(Qimp_peperi[,2])
dwi(Qimp_peperi)
summary(Qimp_peperi)

Peperi[,1]<-Qimp_peperi[,1]*areas.Peperi[1]
Peperi[,2]<-Qimp_peperi[,2]*areas.Peperi[2]

plot(Peperi)

##-----------------------------------------------------------------------------##

## Tubarão

summary(Tubarão)
plot(Tubarão,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Tubarão, var.type="Days"),ColorRamp="Days")
Tubarão<- window(Tubarão,start=as.Date("1987-01-01"),end=as.Date("2004-12-31"))


## Vazão específica

Qesp_tubarão1<-Tubarão[,1]/areas.Tubarão[1]
Qesp_tubarão2<-Tubarão[,2]/areas.Tubarão[2]
Qesp_tubarão3<-Tubarão[,3]/areas.Tubarão[3]  
Qesp_tubarão4<-Tubarão[,4]/areas.Tubarão[4]
Qesp_tubarão5<-Tubarão[,5]/areas.Tubarão[5]
Qesp_tubarão6<-Tubarão[,6]/areas.Tubarão[6]
Qesp_tubarão7<-Tubarão[,7]/areas.Tubarão[7]  
Qesp_tubarão8<-Tubarão[,8]/areas.Tubarão[8]
Qesp_tubarão9<-Tubarão[,9]/areas.Tubarão[9]

Qesp_tubarão<-cbind.zoo(Qesp_tubarão1,Qesp_tubarão2,Qesp_tubarão3,Qesp_tubarão4,Qesp_tubarão5,Qesp_tubarão6,Qesp_tubarão7,Qesp_tubarão8,Qesp_tubarão9)


summary(Qesp_tubarão)
colnames(Qesp_tubarão)<-colnames(Tubarão)
names(Qesp_tubarão)

plot(Qesp_tubarão,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_tubarão1+Qesp_tubarão2+Qesp_tubarão3+Qesp_tubarão4+Qesp_tubarão5+Qesp_tubarão6+Qesp_tubarão7+Qesp_tubarão8+Qesp_tubarão9
i <- mnimput(f,Qesp_tubarão,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_tubarão<-zoo(p,seq(as.Date("1987/01/01"),as.Date("2004/12/31"),"day"));colnames(Qimp_tubarão)<-colnames(Tubarão)
write.csv(Qimp_tubarão,"Estac.Preen/Qimp_tubarão.csv")

colnames(Qimp_tubarão)
plot(Qimp_tubarão[,2])
dwi(Qimp_tubarão)
summary(Qimp_tubarão)

Tubarão[,1]<-Qimp_tubarão[,1]*areas.Tubarão[1]
Tubarão[,2]<-Qimp_tubarão[,2]*areas.Tubarão[2]
Tubarão[,3]<-Qimp_tubarão[,3]*areas.Tubarão[3]  
Tubarão[,4]<-Qimp_tubarão[,4]*areas.Tubarão[4]  
Tubarão[,5]<-Qimp_tubarão[,5]*areas.Tubarão[5]
Tubarão[,6]<-Qimp_tubarão[,6]*areas.Tubarão[6]
Tubarão[,7]<-Qimp_tubarão[,7]*areas.Tubarão[7]  
Tubarão[,8]<-Qimp_tubarão[,8]*areas.Tubarão[8]  
Tubarão[,9]<-Qimp_tubarão[,9]*areas.Tubarão[9]  

plot(Tubarão)

##-----------------------------------------------------------------------------##

## Itapocu

summary(Itapocu)
plot(Itapocu,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Itapocu, var.type="Days"),ColorRamp="Days")
Itapocu<- window(Itapocu,start=as.Date("1978-01-01"),end=as.Date("2001-12-31"))


## Vazão específica

Qesp_itapocu1<-Itapocu[,1]/areas.Itapocu[1]
Qesp_itapocu2<-Itapocu[,2]/areas.Itapocu[2]
Qesp_itapocu3<-Itapocu[,3]/areas.Itapocu[3]  
Qesp_itapocu4<-Itapocu[,4]/areas.Itapocu[4]
Qesp_itapocu5<-Itapocu[,5]/areas.Itapocu[5]


Qesp_itapocu<-cbind.zoo(Qesp_itapocu1,Qesp_itapocu2,Qesp_itapocu3,Qesp_itapocu4,Qesp_itapocu5)


summary(Qesp_itapocu)
colnames(Qesp_itapocu)<-colnames(Itapocu)
names(Qesp_itapocu)

plot(Qesp_itapocu,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_itapocu1+Qesp_itapocu2+Qesp_itapocu3+Qesp_itapocu4+Qesp_itapocu5
i <- mnimput(f,Qesp_itapocu,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_itapocu<-zoo(p,seq(as.Date("1978/01/01"),as.Date("2001/12/31"),"day"));colnames(Qimp_itapocu)<-colnames(Itapocu)
write.csv(Qimp_itapocu,"Estac.Preen/Qimp_itapocu.csv")

colnames(Qimp_itapocu)
plot(Qimp_itapocu[,2])
dwi(Qimp_itapocu)
summary(Qimp_itapocu)

Itapocu[,1]<-Qimp_itapocu[,1]*areas.Itapocu[1]
Itapocu[,2]<-Qimp_itapocu[,2]*areas.Itapocu[2]
Itapocu[,3]<-Qimp_itapocu[,3]*areas.Itapocu[3]  
Itapocu[,4]<-Qimp_itapocu[,4]*areas.Itapocu[4]  
Itapocu[,5]<-Qimp_itapocu[,5]*areas.Itapocu[5]


plot(Itapocu)

##-----------------------------------------------------------------------------##

## Irani

summary(Irani)
plot(Irani,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Irani, var.type="Days"),ColorRamp="Days")
Irani<- window(Irani,start=as.Date("1970-01-01"),end=as.Date("2014-12-31"))


## Vazão específica

Qesp_irani1<-Irani[,1]/areas.Irani[1]
Qesp_irani2<-Irani[,2]/areas.Irani[2]
Qesp_irani3<-Irani[,3]/areas.Irani[3]  

Qesp_irani<-cbind.zoo(Qesp_irani1,Qesp_irani2,Qesp_irani3)

summary(Qesp_irani)
colnames(Qesp_irani)<-colnames(Irani)
names(Qesp_irani)

plot(Qesp_irani,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_irani1+Qesp_irani2+Qesp_irani3
i <- mnimput(f,Qesp_irani,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_irani<-zoo(p,seq(as.Date("1970/01/01"),as.Date("2014/12/31"),"day"));colnames(Qimp_irani)<-colnames(Irani)
write.csv(Qimp_irani,"Estac.Preen/Qimp_irani.csv")

colnames(Qimp_irani)
plot(Qimp_irani[,2])
dwi(Qimp_irani)
summary(Qimp_irani)

Irani[,1]<-Qimp_irani[,1]*areas.Irani[1]
Irani[,2]<-Qimp_irani[,2]*areas.Irani[2]
Irani[,3]<-Qimp_irani[,3]*areas.Irani[3]  

plot(Irani)
##-----------------------------------------------------------------------------##

## Tijucas

summary(Tijucas)
plot(Tijucas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Tijucas, var.type="Days"),ColorRamp="Days")
Tijucas<- window(Tijucas,start=as.Date("1983-01-01"),end=as.Date("2004-12-31"))


## Vazão específica

Qesp_tijucas1<-Tijucas[,1]/areas.Tijucas[1]
Qesp_tijucas2<-Tijucas[,2]/areas.Tijucas[2]
Qesp_tijucas3<-Tijucas[,3]/areas.Tijucas[3]  

Qesp_tijucas<-cbind.zoo(Qesp_tijucas1,Qesp_tijucas2,Qesp_tijucas3)

summary(Qesp_tijucas)
colnames(Qesp_tijucas)<-colnames(Tijucas)
names(Qesp_tijucas)

plot(Qesp_tijucas,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_tijucas1+Qesp_tijucas2+Qesp_tijucas3
i <- mnimput(f,Qesp_tijucas,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_tijucas<-zoo(p,seq(as.Date("1983/01/01"),as.Date("2004/12/31"),"day"));colnames(Qimp_tijucas)<-colnames(Tijucas)
write.csv(Qimp_tijucas,"Estac.Preen/Qimp_tijucas.csv")

colnames(Qimp_tijucas)
plot(Qimp_tijucas[,2])
dwi(Qimp_tijucas)
summary(Qimp_tijucas)

Tijucas[,1]<-Qimp_tijucas[,1]*areas.Tijucas[1]
Tijucas[,2]<-Qimp_tijucas[,2]*areas.Tijucas[2]
Tijucas[,3]<-Qimp_tijucas[,3]*areas.Tijucas[3]  

plot(Tijucas)
##-----------------------------------------------------------------------------##

## Chapeco

summary(Chapeco)
plot(Chapeco[,1],ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
matrixplot(dwi(Chapeco, var.type="Days"),ColorRamp="Days")
Chapeco<- window(Chapeco,start=as.Date("1979-01-01"),end=as.Date("2010-12-31"))

## Vazão específica

Qesp_chapeco1<-Chapeco[,1]/areas.Chapeco[1]
Qesp_chapeco2<-Chapeco[,2]/areas.Chapeco[2]
Qesp_chapeco3<-Chapeco[,3]/areas.Chapeco[3]  
Qesp_chapeco4<-Chapeco[,4]/areas.Chapeco[4]
Qesp_chapeco5<-Chapeco[,5]/areas.Chapeco[5]
Qesp_chapeco6<-Chapeco[,6]/areas.Chapeco[6]
Qesp_chapeco7<-Chapeco[,7]/areas.Chapeco[7]
Qesp_chapeco8<-Chapeco[,8]/areas.Chapeco[8]  
Qesp_chapeco9<-Chapeco[,9]/areas.Chapeco[9]

Qesp_chapeco<-cbind.zoo(Qesp_chapeco1,Qesp_chapeco2,Qesp_chapeco3,Qesp_chapeco4,Qesp_chapeco5,Qesp_chapeco6,Qesp_chapeco7,Qesp_chapeco8,Qesp_chapeco9)


summary(Qesp_chapeco)
colnames(Qesp_chapeco)<-colnames(Chapeco)
names(Qesp_chapeco)

plot(Qesp_chapeco,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_chapeco1+Qesp_chapeco2+Qesp_chapeco3+Qesp_chapeco4+Qesp_chapeco5+Qesp_chapeco6+Qesp_chapeco7+Qesp_chapeco8+Qesp_chapeco9
i <- mnimput(f,Qesp_chapeco,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_chapeco<-zoo(p,seq(as.Date("1979/01/01"),as.Date("2010/12/31"),"day"));colnames(Qimp_chapeco)<-colnames(Chapeco)
write.csv(Qimp_chapeco,"Estac.Preen/Qimp_chapeco.csv")

colnames(Qimp_chapeco)
plot(Qimp_chapeco[,1])
dwi(Qimp_chapeco)
summary(Qimp_chapeco)

Chapeco[,1]<-Qimp_chapeco[,1]*areas.Chapeco[1]
Chapeco[,2]<-Qimp_chapeco[,2]*areas.Chapeco[2]
Chapeco[,3]<-Qimp_chapeco[,3]*areas.Chapeco[3]  
Chapeco[,4]<-Qimp_chapeco[,4]*areas.Chapeco[4]  
Chapeco[,5]<-Qimp_chapeco[,5]*areas.Chapeco[5]
Chapeco[,6]<-Qimp_chapeco[,6]*areas.Chapeco[6]
Chapeco[,7]<-Qimp_chapeco[,7]*areas.Chapeco[7]
Chapeco[,8]<-Qimp_chapeco[,8]*areas.Chapeco[8]  
Chapeco[,9]<-Qimp_chapeco[,9]*areas.Chapeco[9]  

plot(Chapeco)

##-----------------------------------------------------------------------------##

## Cubatão Norte

summary(CubataoN)
plot(CubataoN,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
Negro1<- window(Negro,start=as.Date("1986-01-01"),end=as.Date("2010-02-28"))
matrixplot(dwi(Negro1, var.type="Days"),ColorRamp="Days")
CubataoN1 <- cbind(CubataoN,Negro1);colnames(CubataoN1) <- c("E82270050",colnames(Negro1))
hydropairs(as.data.frame(CubataoN1),met="spear")

## Vazão específica

Qesp_cubataoN11<-CubataoN1[,1]/areas.CubataoN[1]
Qesp_cubataoN12<-CubataoN1[,2]/areas.Negro[1]
Qesp_cubataoN13<-CubataoN1[,3]/areas.Negro[2]  
Qesp_cubataoN14<-CubataoN1[,4]/areas.Negro[3]

Qesp_cubataoN<-cbind.zoo(Qesp_cubataoN11,Qesp_cubataoN12,Qesp_cubataoN13,Qesp_cubataoN14)

summary(Qesp_cubataoN)
colnames(Qesp_cubataoN)<-colnames(CubataoN1)
names(Qesp_cubataoN)

plot(Qesp_cubataoN,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_cubataoN11+Qesp_cubataoN12+Qesp_cubataoN13+Qesp_cubataoN14
i <- mnimput(f,Qesp_cubataoN,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_cubataoN<-zoo(p,seq(as.Date("1986/01/01"),as.Date("2010/02/28"),"day"));colnames(Qimp_cubataoN)<-colnames(CubataoN1)
write.csv(Qimp_cubataoN[,1],"Estac.Preen/Qimp_cubataoN.csv")

colnames(Qimp_cubataoN)
plot(Qimp_cubataoN[,1])
dwi(Qimp_cubataoN)
summary(Qimp_cubataoN)

CubataoN1[,1]<-Qimp_cubataoN[,1]*areas.CubataoN[1]
CubataoN1[,2]<-Qimp_cubataoN[,2]*areas.Negro[2]
CubataoN1[,3]<-Qimp_cubataoN[,3]*areas.Negro[3]  
CubataoN1[,4]<-Qimp_cubataoN[,4]*areas.Negro[4]  

plot(CubataoN1)

##-----------------------------------------------------------------------------##

## Cubatão Sul

summary(CubataoS)
plot(CubataoS,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
Tijucas1<- window(Tijucas,start=as.Date("1951-01-01"),end=as.Date("2000-12-31"))[,-3]
matrixplot(dwi(Tijucas1, var.type="Days"),ColorRamp="Days")
CubataoS1 <- cbind(CubataoS,Tijucas1);colnames(CubataoS1) <- c("E84100000",colnames(Tijucas1))
hydropairs(as.data.frame(CubataoS1),met="spear")

## Vazão específica

Qesp_cubataoS11<-CubataoS1[,1]/areas.CubataoS
Qesp_cubataoS12<-CubataoS1[,2]/areas.Tijucas[1]
Qesp_cubataoS13<-CubataoS1[,3]/areas.Tijucas[2]  

Qesp_cubataoS<-cbind.zoo(Qesp_cubataoS11,Qesp_cubataoS12,Qesp_cubataoS13)

summary(Qesp_cubataoS)
colnames(Qesp_cubataoS)<-colnames(CubataoS1)
names(Qesp_cubataoS)

plot(Qesp_cubataoS,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_cubataoS11+Qesp_cubataoS12+Qesp_cubataoS13
i <- mnimput(f,Qesp_cubataoS,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_cubataoS<-zoo(p,seq(as.Date("1951/01/01"),as.Date("2000/12/31"),"day"));colnames(Qimp_cubataoS)<-colnames(CubataoS1)
write.csv(Qimp_cubataoS[,1],"Estac.Preen/Qimp_cubataoS.csv")

colnames(Qimp_cubataoS)
plot(Qimp_cubataoS[,1])
dwi(Qimp_cubataoS)
summary(Qimp_cubataoS)

CubataoS1[,1]<-Qimp_cubataoS[,1]*areas.CubataoS
CubataoS1[,2]<-Qimp_cubataoS[,2]*areas.Tijucas[1]
CubataoS1[,3]<-Qimp_cubataoS[,3]*areas.Tijucas[2]  


plot(CubataoS1)

##-----------------------------------------------------------------------------##

## Mampituba

summary(Mampituba)
summary(Ararangua)
plot(Mampituba,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")
Ararangua1<- window(Ararangua,start=as.Date("1986-01-01"),end=as.Date("2006-12-31"))
matrixplot(dwi(Ararangua1, var.type="Days"),ColorRamp="Days")
Mampituba1 <- cbind(Mampituba,Ararangua1);colnames(Mampituba1) <- c("E84970000",colnames(Ararangua1))
hydropairs(as.data.frame(Mampituba1),met="spear")

## Vazão específica

Qesp_mampituba1<-Mampituba1[,1]/areas.Mampituba
Qesp_mampituba2<-Mampituba1[,2]/areas.Ararangua[1]
Qesp_mampituba3<-Mampituba1[,3]/areas.Ararangua[2]  
Qesp_mampituba4<-Mampituba1[,4]/areas.Ararangua[3]
Qesp_mampituba5<-Mampituba1[,5]/areas.Ararangua[4]  

Qesp_mampituba<-cbind.zoo(Qesp_mampituba1,Qesp_mampituba2,Qesp_mampituba3,Qesp_mampituba4,Qesp_mampituba5)

summary(Qesp_mampituba)
colnames(Qesp_mampituba)<-colnames(Mampituba1)
names(Qesp_mampituba)

plot(Qesp_mampituba,ylab=expression(Q~(m^3.~s^{-1})),xlab="Tempo",main="")

## Preenchimento de falhas
f <- ~Qesp_mampituba1+Qesp_mampituba2+Qesp_mampituba3+Qesp_mampituba4+Qesp_mampituba5
i <- mnimput(f,Qesp_mampituba,eps=1e-3,ts=F,log=T, method="spline",sp.control=list(df=c(7,7,7,7,7)))
## plot(i)
p <- ifelse(predict(i)<0,predict(i)*-1,predict(i))

Qimp_mampituba<-zoo(p,seq(as.Date("1986/01/01"),as.Date("2006/12/31"),"day"));colnames(Qimp_mampituba)<-colnames(Mampituba1)
write.csv(Qimp_mampituba[,1],"Estac.Preen/Qimp_mampituba.csv")

colnames(Qimp_mampituba)
plot(Qimp_mampituba[,1])
dwi(Qimp_mampituba)
summary(Qimp_mampituba)

Mampituba1[,1]<-Qimp_mampituba[,1]*areas.Mampituba
Mampituba1[,2]<-Qimp_mampituba[,2]*areas.Ararangua[1]
Mampituba1[,3]<-Qimp_mampituba[,3]*areas.Ararangua[2]  


plot(Mampituba1)

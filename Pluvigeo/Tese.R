library(mtsdi)
library(hydroTSM)
options(OutDec=",")

#Abrir arquivo

getwd()
setwd("/home/wagner/GDrive/Doutorado/Estações/Fluviométricos/imputação de dados")
dir()
vazao<-read.zoo("EstaçõesR.csv",format="%d/%m/%Y",sep=";",dec=",", head=T)
class(vazao)
vazao
help(plotmath)
matrixplot(dwi(vazao), var.type="Days")

# Casas decimais/normaliza??o
w<-log(vazao)
w
# Histograma
hist(log(vazao[,2]))

# Teste de normalidade 1? coluna [l,c]
ks.test(w[,1],"pnorm",mean(w[,1]),sd(w[,1]))
shapiro.test(w[,1])
qqnorm(w[,2])
abline(0,1)
qqline(w[,1])
qqnorm(w[,1])

summary(w)
summary(vazao)

f <- ~E72810000+E74295000
i <- mnimput(f,w,eps=1e-3,ts=TRUE, method="spline",sp.control=list(df=c(7,7,7,7,7)))
predict(i)
imput<-as.zoo(exp(predict(i)))
imput
class(imput)

m<-cbind.zoo(vazao[,1],imput[,2])
class(m)

hydroplot(m[,1], FUN=max, var.type="flow",ylab= "Q", var.unit = "m3/s")
hydroplot(m[,1], pfreq="o",ptype="ts", FUN=mean, stype="default", ylab="Q (m³/s)", adj=1, mgp=c(2.5,1,0))


# log neperiano da 3? coluna
n<-log(w[,3])
n
a<-w[-seq(1,10),1]
a
plot(i)
w[,2]<-log(w[,2])

log(w)
shapiro.test(w[,3])
hist(w[,3])

library(hydroTSM)
data(OcaEnOnaQts)
hydroplot(OcaEnOnaQts, pfreq="seasonal", FUN=mean, stype="default",mgp=c(2,1,0))

Q.7 <- rollapply(data=serie[,5], width=7, FUN=mean, fill=NA, partial= TRUE,align="left")
hydroplot(Q.7, ptype="ts+boxplot", pfreq="o", var.unit="mm")          
daily2annual(Q.7, FUN=min, na.rm=TRUE) 
fdc(imput[,2])




x <- as.numeric(OcaEnOnaQts)
x
y <- x + rnorm(length(x), mean=10)
y
xx <- data.frame(x=x, y=y)
xx
fdc(xx, thr.shw=T,log="", main= "Curva de Perman?ncia", xlab="% Perman?ncia",mgp=c(2,1,0), ylab=expression(Q~(m3.~s^{-1})))

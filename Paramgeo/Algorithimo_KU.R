## Leitura dos dados

dados <- read.table("http://www.stat.ucla.edu/~nchristo/statistics_c173_c273/kriging_11.txt", header=TRUE)

## dados <- read.csv("Param_geo.csv", head=T,dec=",",sep =
## ";")[-22,];dados$Centroid.Y<-dados$Centroid.Y/1000;dados$Centroid.X<-dados$Centroid.X/1000

## limite<- read.csv("Limite_Canoas.csv",h=T,sep=",");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
## names(dados)
## dim(dados)
## head(dados)

## Definição dos parâmetros tese
## mu <- mean(dados$Param_ANO.sigma)
## sigma2 <- 0.00136
## phi <- 9.8905
## tau2 <- 0.0011

## Definição dos parâmetros teste

sigma2 <- 10
phi <- 3.33
tau2 <- 0.4


## Coordenada do ponto a ser predito


Xo <- c(65,137)

## Calcular a matrix de Distância entre os dados 

d <- dist(dados[,1:2],diag = T,upper=T)
d

## calcular o vetor de distância entre os dados e o ponto a ser predito

do <- dist(rbind(Xo,dados[,1:2]))[1:length(dados$z)]
do           

## Calcular a matrix R de correlação entre os dados

R <- as.matrix((sigma2*(exp(-d/phi)))-tau2);colnames(R) <- NULL
(trend.matrix <- as.matrix(cbind(rep(1,l=nrow(R)),dados[,1:2])));colnames(trend.matrix) <- NULL
diag(R) <- sigma2
R
trend.matrix
R3 <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko
R1 <- cbind(R,trend.matrix);R1 
R2 <- cbind(t(R[,-c(1:length(dados$z))]),matrix(rep(0,l=ncol(trend.matrix)^2),nrow = 3,ncol=3));R2
R3 <- rbind(R1,R2);R3


## Calcular a regressão com as coordenadas
lm <- lm(z~x+y,dados)
lm$residuals
plot(dados$x,dados$z)
lm1 <-lm$coef[1]+ dados$x[1]*lm$coef[2]+dados$y[1]*lm$coef[3]
dados$z[1]-lm1
lm$residuals

## calcular o vetor de correlação entre os dados e o ponto a ser predito

r <- rbind(as.matrix((sigma2*(exp(-do/phi)))-tau2),1);r


## calcular a Matriz inversa de R

R.inv <- solve(R3)
R.inv
(Pesos <- R.inv%*%r)
sum(Pesos[c(1:(length(Pesos)-ncol(trend.matrix)))]) ## Esse deve ser igual a 1 

Pesos1 <- Pesos[c(1:(length(Pesos)-ncol(trend.matrix)))]

## Predição

Yo <- sum(Pesos1*dados$z)
Yo

Yo <- sum(Pesos1*lm$residuals)
lm$residuals
Yo
head(dados)
## No geoR para conferir
require(geoR)
b <-as.geodata(dados,coords.col = 1:2, data.col = 3)
names(dados)



univ_kr <- ksline(b, cov.model="exp", cov.pars=c(10,3.33), nugget=0.4,locations=Xo, m0="kt", trend=1)


univ_kr$predict

kc <- krige.conv(b, locations=Xo,krige=krige.control(type.krige="ok",
cov.model="exp", cov.pars=c(10,3.33), nugget=0.2,trend.l="1st",trend.d="1st"))

(kc)
kc$predict
bmodel <- likfit(b, ini=c(0.004,25.71), nug=0.005,trend=~poly(coords[,1],2)+coords[,2]+Param_SON.mu)


567.6581-0.0011

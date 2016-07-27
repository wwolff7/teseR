## Leitura dos dados
options(OutDec=",")
rm(list = ls())
require(geoR)
require(MASS)
require(matrixcalc)
dados <- read.csv("Param.csv", head=T,dec=",",sep = "");dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000

limite<- read.csv("Limite_Canoas.csv",h=T,sep=",");limite$Longitude<-limite$Longitude/1000;limite$Latitude<-limite$Latitude/1000
names(dados)
dim(dados)
head(dados)

## Coordenada do ponto a ser predito
Xo <- c(670563.231,6831157.59)/1000

## Area
area <- 219.7403073




## Calcular a matrix de Distância entre os dados 
d <- dist(dados[,2:3],diag = T,upper=T)
d
dados
## calcular o vetor de distância entre os dados e o ponto a ser predito
do <- dist(rbind(Xo,dados[,2:3]))[1:nrow(dados)]
do           

## Definição dos parâmetros 
## ANO Mu

tau2.mu.ANO <- 0.020651983740
sigma2.mu.ANO <- 0.037100700650
phi.mu.ANO <-61.550407630000

## Calcular a matrix R de semivariância entre os dados

R <- as.matrix(tau2.mu.ANO+(sigma2.mu.ANO*(1-exp(-d^2/phi.mu.ANO^2))));colnames(R) <- NULL;R
## diag(R) <- tau2.mu.ANO
R
## ou covariância
## R <- as.matrix((sigma2.mu.ANO*(exp(-d^2/phi.mu.ANO^2)))-tau2.mu.ANO);colnames(R) <- NULL;R
## diag(R) <- sigma2.mu.ANO
## R

## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de semivariância entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.mu.ANO+(sigma2.mu.ANO*(1-exp(-do^2/phi.mu.ANO^2)))),1);r0
## ou covariância
## r0 <- rbind(as.matrix((sigma2.mu.ANO*(exp(-do^2/phi.mu.ANO^2)))-tau2.mu.ANO),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
(R.ko.inv <- ginv(R.ko))

## verificando o posto ou rank da matriz Agenomic
qr(R.ko) $ rank
sum(diag(R.ko)) # posto incompleto

## obs: para ser posto completo, é necessario o rank da matriz ser igual ao traço da
## matriz (a Matriz deve ser quadrada, simetrica e idempotente)

## t(r0)%*%R.ko%*%r0>0
## (auto_val <- eigen(R.ko)[[1]])>0 #Se tdos forem > 0 é PD
## sum(auto_val)
## sum(diag(R.ko)) ## traço
##  
## prod(diag(R.ko))
## det(R.ko)

## (R.ko.inv1 <- solve(R.ko))

## R.ko.inv-R.ko.inv1
##  
## diag(R.ko.inv)
## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)
## (W.ko <- solve(R.ko,r0))

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição
(mu.ANO.ko <- sum(dados$mu_ANO*W.ko1))
mu.ANO.ko
## Conferir na GeoR, não fazer em PHP
dados <- read.csv("ParamDJF.csv", head=T,dec=",",sep= ";")[,-c(1,27)];dados$Latitude<-dados$Latitude/1000;dados$Longitude<-dados$Longitude/1000
b <-as.geodata(dados,coords.col = 1:2, data.col = 4)
names(dados)
dir()
summary(b)

tau2.sig.DJF
sigma2.sig.DJF
phi.sig.DJF

trend_Mu.ANO <- dados[,c(5,6,8:25,1,2)]
names(trend_Mu.ANO)

likkc <- likfit(b,cov.model = "gau",ini=c(sigma2.sig.DJF-0.01,phi.sig.DJF),fix.nugget = T, trend = ~poly(trend_Mu.ANO[,8],2)+poly(trend_Mu.ANO[,9],2),nug=tau2.sig.DJF,lambda=1)
likkc$parameters.summary
likkc$trend
likkc$sigmasq

-2*logLik(likkc)[1]+2*5==AIC(likkc)

kc <- krige.conv(b, locations=Xo,krige=krige.control(obj=likkc))
kc$predict
kc$krige.var

sig.DJF.ko

(kc$predict <- backtransform.moments(lambda = lambda.mu.ANO, mean = kc$predict, variance = kc$krige.var)$mean)-5

(kc$predict <- (((lambda.mu.ANO*kc$predict)+1)^(1/lambda.mu.ANO)*(1+((kc$krige.var*(1-lambda.mu.ANO))/(2*((lambda.mu.ANO*kc$predict)+1)^2))))-5)


##-----------------------------------------------------------------------------##

## Definição dos parâmetros 
## ANO sigma

tau2.sig.ANO   <- 0.009944800760
sigma2.sig.ANO <-0.025371118320
phi.sig.ANO    <-102.5468827


## Calcular a matrix R de semivariância entre os dados
R <- as.matrix(tau2.sig.ANO+(sigma2.sig.ANO*(1-exp(-d^2/phi.sig.ANO^2))));colnames(R) <- NULL;R
## diag(R) <- tau2.sig.ANO
## KO das variáveis externas 

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.sig.ANO+(sigma2.sig.ANO*(1-exp(-do^2/phi.sig.ANO^2)))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição
(sigma.ANO.ko <- sum(dados$sigma_ANO*W.ko1))
sigma.ANO.ko

(Q98_ANO <- Qp(mu.ANO.ko,sigma.ANO.ko,area,0.98))
area
## FIM

## Definição dos parâmetros 
## DJF Mu

tau2.mu.DJF <- 0.024194871550
sigma2.mu.DJF <- 0.185078996800
phi.mu.DJF <- 92.469561470000

## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.mu.DJF+(sigma2.mu.DJF*(1-exp(-d^2/phi.mu.DJF^2))));colnames(R) <- NULL;R
## diag(R) <- tau2.mu.DJF

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.mu.DJF+(sigma2.mu.DJF*(1-exp(-do^2/phi.mu.DJF^2)))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(mu.DJF.ko <- sum(dados$mu_DJF*W.ko1))

##-----------------------------------------------------------------------------##


## Definição dos parâmetros 
## DJF sigma

tau2.sig.DJF   <- 0.01308691549
sig2.sig.DJF <-0.01831340988
phi.sig.DJF    <-61.98482893


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.sig.DJF+(sig2.sig.DJF*(1-exp(-d^2/phi.sig.DJF^2))));colnames(R) <- NULL;R

##diag(R) <- tau2.sig.DJF
## KO das variáveis externas 

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.sig.DJF+(sigma2.sig.DJF*(1-exp(-do^2/phi.sig.DJF^2)))),1);r0


## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição
(sig.DJF.ko <- sum(dados$sigma_DJF*W.ko1))


##-----------------------------------------------------------------------------##

## Definição dos parâmetros 
## MAM Mu

tau2.mu.MAM <- 0.023997466370
sigma2.mu.MAM <- 0.055831276960
phi.mu.MAM <-59.443812610000


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.mu.MAM+(sigma2.mu.MAM*(1-exp(-d^2/phi.mu.MAM^2))));colnames(R) <- NULL;R
## diag(R) <- tau2.mu.MAM

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.mu.MAM+(sigma2.mu.MAM*(1-exp(-do^2/phi.mu.MAM^2)))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(mu.MAM.ko <- sum(dados$mu_MAM*W.ko1))
mu.MAM.ko
##-----------------------------------------------------------------------------##


## Definição dos parâmetros 
## MAM sigma

tau2.sig.MAM   <- 0.008345854746
sigma2.sig.MAM <- 0.016656155770
phi.sig.MAM    <-109.024978700000


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.sig.MAM+(sigma2.sig.MAM*((1.5*(d/phi.sig.MAM)-(0.5*(d/phi.sig.MAM)^3)))));colnames(R) <- NULL;R

##diag(R) <- tau2.sig.MAM

## KO das variáveis externas 

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.sig.MAM+(sigma2.sig.MAM*((1.5*(do/phi.sig.MAM)-(0.5*(do/phi.sig.MAM)^3))))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(sigma.MAM.ko <- sum(dados$sigma_MAM*W.ko1))

##-----------------------------------------------------------------------------##

## Definição dos parâmetros 
## JJA Mu

tau2.mu.JJA <- 0.014478953390
sigma2.mu.JJA <- 0.064626601480
phi.mu.JJA <- 161.702330500000


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.mu.JJA+(sigma2.mu.JJA*((1.5*(d/phi.mu.JJA)-(0.5*(d/phi.mu.JJA)^3)))));colnames(R) <- NULL;R
## diag(R) <- tau2.mu.JJA
## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.mu.JJA+(sigma2.mu.JJA*((1.5*(do/phi.mu.JJA)-(0.5*(do/phi.mu.JJA)^3))))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(mu.JJA.ko <- sum(dados$mu_JJA*W.ko1))

##-----------------------------------------------------------------------------##


## Definição dos parâmetros 
## JJA sigma


tau2.sig.JJA   <- 0.008241509970
sigma2.sig.JJA <- 0.021191031240
phi.sig.JJA    <- 245.965267800000


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.sig.JJA+(sigma2.sig.JJA*((1.5*(d/phi.sig.JJA)-(0.5*(d/phi.sig.JJA)^3)))));colnames(R) <- NULL;R
## diag(R) <- tau2.mu.JJA
## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.sig.JJA+(sigma2.sig.JJA*((1.5*(do/phi.sig.JJA)-(0.5*(do/phi.sig.JJA)^3))))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(sigma.JJA.ko <- sum(dados$sigma_JJA*W.ko1))
sigma.JJA.ko
##-----------------------------------------------------------------------------##
## Definição dos parâmetros 
## SON Mu

tau2.mu.SON <- 0.016795539330
sigma2.mu.SON <- 0.040937742080
phi.mu.SON <-142.663965100000


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.mu.SON+(sigma2.mu.SON*((1.5*(d/phi.mu.SON)-(0.5*(d/phi.mu.SON)^3)))));colnames(R) <- NULL;R
## diag(R) <- tau2.mu.SON

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.mu.SON+(sigma2.mu.SON*(1.5*(do/phi.mu.SON)-(0.5*(do/phi.mu.SON)^3)))),1);r0
cov.spatial()
## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(mu.SON.ko <- sum(dados$mu_SON*W.ko1))
mu.SON.ko
##-----------------------------------------------------------------------------##

## Definição dos parâmetros 
## SON sigma

tau2.sig.SON   <- 0.009145455499
sigma2.sig.SON <-0.020674607710
phi.sig.SON    <-77.377897380000


## Calcular a matrix R de semivariância entre os dados
dados
R <- as.matrix(tau2.sig.SON+(sigma2.sig.SON*(1-exp(-d^2/phi.sig.SON^2))));colnames(R) <- NULL;R
## diag(R) <- tau2.sig.SON

## KO das variáveis externas 

## R1 <- cbind(R,rep(1,l=nrow(R))
## Acrescenter uma linha e coluna de 1 "Langrangeano"
R.ko <- rbind(cbind(R,rep(1,l=nrow(R))),c(rep(1,l=nrow(R)),0));R.ko

## calcular o vetor de correlação entre os dados e o ponto a ser predito
r0 <- rbind(as.matrix(tau2.sig.SON+(sigma2.sig.SON*(1-exp(-do^2/phi.sig.SON^2)))),1);r0

## calcular a Matriz inversa de R.ko
is.positive.definite(R.ko) #Se for falso fazer a inversa generalizada
R.ko.inv <- ginv(R.ko)

## Determinar os Pesos
(W.ko <- R.ko.inv%*%r0)

sum(W.ko[-length(r0)]) ## Esse deve ser igual a 1 

## O ultimo valor não é peso, é multiplicador de Lagrange
(W.ko1 <- W.ko[-length(r0)])

## Predição

(sigma.SON.ko <- sum(dados$sigma_SON*W.ko1))


##-----------------------------------------------------------------------------##
## Regionalização da curva de permanência

## Funções necessárias

## Curva de permanência (CP)
Qp <- function(mu,sigma,area,p){
    Qp<-(exp(mu+(sigma*qnorm(1-p))))*area
}

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

## função para determinar o volume do reservatório, somente usada para ANO
Vr <- function(mu, sigma, area, T,Qf){

    pqfi <-1-pnorm((log(Qf/area)-mu)/sigma) # probabilidade inferior na CP
    pqff <- 1-(1/T) # probabilidade superior na CP

Vr <- ((((pqff-pqfi)*Qf)-simpson(function(p) Qp(mu,sigma,area,p),pqfi,pqff,n=100000))*(60*60*24*365))/10^6
##Vr <- ((((pqff-pqfi)*Qf)-integrate(function(p) Qp(mu,sigma,area,p),pqfi,pqff,n=100000))*(60*60*24*365))/10^6
return(Vr)    
}

##-----------------------------------------------------------------------------##

## Dados de entrada
## vetor de probabilidades (p)
p <- seq(0.01,0.99,l=100)

1-1/100
p
qnorm(1-p)
## Área da bacia de drenagem, fornecido pelo usuário 
area <- 123

##-----------------------------------------------------------------------------##
## Resultados

## Vazão média plurianual
## (Qm_ANO_teste<- integrate(f=function(p) Qp(mu.ANO.ko,sigma.ANO.ko,area,p),low=0.00000001,uppe=1,subdivisions = 100000))
(Qm_ANO<- simpson(fun=function(p) Qp(mu.ANO.ko,sigma.ANO.ko,area,p),0.00000001,1,n=100000))
(Qm_DJF<- simpson(fun=function(p) Qp(mu.DJF.ko,sigma.DJF.ko,area,p),0.00000001,1,n=100000))
(Qm_MAM<- simpson(fun=function(p) Qp(mu.MAM.ko,sigma.MAM.ko,area,p),0.00000001,1,n=100000))
(Qm_JJA<- simpson(fun=function(p) Qp(mu.JJA.ko,sigma.JJA.ko,area,p),0.00000001,1,n=100000))
(Qm_SON<- simpson(fun=function(p) Qp(mu.SON.ko,sigma.SON.ko,area,p),0.00000001,1,n=100000))

## Vazão com 90% de permanência no tempo
(Q90_ANO <- Qp(mu.ANO.ko,sigma.ANO.ko,area,0.9))
(Q90_DJF <- Qp(mu.DJF.ko,sigma.DJF.ko,area,0.9))
(Q90_MAM <- Qp(mu.MAM.ko,sigma.MAM.ko,area,0.9))
(Q90_JJA <- Qp(mu.JJA.ko,sigma.JJA.ko,area,0.9))
(Q90_SON <- Qp(mu.SON.ko,sigma.SON.ko,area,0.9))

## Vazão com 95% de permanência no tempo
(Q95_ANO <- Qp(mu.ANO.ko,sigma.ANO.ko,area,0.95))
(Q95_DJF <- Qp(mu.DJF.ko,sigma.DJF.ko,area,0.95))
(Q95_MAM <- Qp(mu.MAM.ko,sigma.MAM.ko,area,0.95))
(Q95_JJA <- Qp(mu.JJA.ko,sigma.JJA.ko,area,0.95))
(Q95_SON <- Qp(mu.SON.ko,sigma.SON.ko,area,0.95))

## Vazão com 98% de permanência no tempo
(Q98_ANO <- Qp(mu.ANO.ko,sigma.ANO.ko,area,0.98))
(Q98_DJF <- Qp(mu.DJF.ko,sigma.DJF.ko,area,0.98))
(Q98_MAM <- Qp(mu.MAM.ko,sigma.MAM.ko,area,0.98))
(Q98_JJA <- Qp(mu.JJA.ko,sigma.JJA.ko,area,0.98))
(Q98_SON <- Qp(mu.SON.ko,sigma.SON.ko,area,0.98))

## Máxima vazão possível de ser outorgada
(Qout_max_ANO <- Q98_ANO*0.5)
(Qout_max_DJF <- Q98_DJF*0.5)
(Qout_max_MAM <- Q98_MAM*0.5)
(Qout_max_JJA <- Q98_JJA*0.5)
(Qout_max_SON <- Q98_SON*0.5)

## Gráfico da curva de permanência
## Valor máximo no eixo y, o [2] significa que eu peguei o elemento 2 do vetor, pois o
## primeiro tende ao infinito
ymax <- (exp(mu.ANO.ko+(mean(sigma.ANO.ko)*qnorm(1-p)))*area)[2]*1.4

## Curvas de permanência

plot(p*100,exp(mu.ANO.ko+(sigma.ANO.ko*qnorm(1-p)))*area,type="l",xlab =
         "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0),ylim = c(0,ymax))
lines(p*100,exp(mu.DJF.ko+(sigma.DJF.ko*qnorm(1-p)))*area,lty=1,col=2)
lines(p*100,exp(mu.MAM.ko+(sigma.MAM.ko*qnorm(1-p)))*area,lty=1,col=3)
lines(p*100,exp(mu.JJA.ko+(sigma.JJA.ko*qnorm(1-p)))*area,lty=1,col=4)
lines(p*100,exp(mu.SON.ko+(sigma.SON.ko*qnorm(1-p)))*area,lty=1,col=5)
legend("topright",c("Anual","Verão","Outono","Inverno","Primavera"),bty = "n",lty=c(1,1,1,1,1),col=c(1,2,3,4,13),lwd = c(1.5,1.5,1.5,1.5,1.5))

exp(mu.ANO.ko+(sigma.ANO.ko*qnorm(1-p)))*area

par(mfrow = c(2,2))#,mar = c(3, 3, 3, 3))
plot(p*100,exp(mu.DJF.ko+(sigma.DJF.ko*qnorm(1-p)))*area,type="l",xlab =
         "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0),main =
             "(a) - Verão")
plot(p*100,exp(mu.MAM.ko+(sigma.MAM.ko*qnorm(1-p)))*area,type="l",xlab =
         "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0),main =
             "(b) - Outono")
plot(p*100,exp(mu.JJA.ko+(sigma.JJA.ko*qnorm(1-p)))*area,type="l",xlab =
         "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0),main =
             "(c) - Inverno")
plot(p*100,exp(mu.SON.ko+(sigma.SON.ko*qnorm(1-p)))*area,type="l",xlab =
         "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0),main =
             "(d) - Primavera")

##-----------------------------------------------------------------------------##
## Regularização de volume SOMENTE PARA ANO

## Vazão demandada
(Qdem <- 1.679)

## Vazão outorgada pela SDS
(Qout <- 0.075)

## vazão a jusante
(Qj <- Q98_ANO - Qout)

## vazão firme, a vazão pontecial de ser regularizada 10% de perdas
(Qf <- (Qdem+Qj)*1.10)

## Tempo de retorno (anos) relacionado a probabilidade de deficit hídrico
T <- 30

## Volume do reservatório(m^3x10^6) para o tempo de retorno em questão
Vr(mu.ANO.ko,sigma.ANO.ko,area,T,Qf)

##-----------------------------------------------------------------------------##
## Gráfico volume do reservatório, obs. não precisa fazer por agora
plot(p*100,exp(mu.ANO.ko+(mean(sigma.ANO.ko)*qnorm(1-p)))*area,type="l",xlab =
         "Permanência (%)",ylab = expression(Q[p]~~(m^3%.%s^{-1})),mgp=c(2.5,1,0))
lines(p*100,rep(Qf,times=length(p)),col=2,lty=2,lwd = 2) # Linha vazão firme

## probabilidade da Qf na curva de permanência 
pqfi <-1-pnorm((log(Qf/area)-mu.ANO.ko)/sigma.ANO.ko)

## probabilidade do T na curva de permanência 
pqff <- 1-(1/T)

## vetor de sequência de pqfi a pqff em %
pqfseq <- seq(pqfi*100,pqff*100,l=1000) 

## trecho da curva de permanência entre pqfi e pqff
Qpseq <- (exp(mu.ANO.ko+(sigma.ANO.ko*qnorm(1-(pqfseq/100)))))*area

## "pintar" polígono do Vr
coords.x <- c(pqfi*100,pqfseq,pqff*100) 
coords.y <- c(Qf,Qpseq,Qf)
polygon(coords.x,coords.y,col=1)

## textos no gráfico
text(0,Qf*1.50,expression(Q[f]),cex=1.5)
text(mean(pqfseq),Qf*1.50,bquote(V[r]~para~T==.(T)~anos),cex = 1.5)

##-----------------------------------------------------------------------------##
## Simulações reservatório (sim)

## Criar vetor com tempos de retorno
(Tsim <- c(2,5,10,20,30,40,50,60,70,80,90,100))

## Volume do reservatório para cada tempo de retorno
Vr.Tsim <- c()
for(i in 1:length(Tsim)){
    Vr.Tsim[i] <- Vr(mu.ANO.ko,sigma.ANO.ko,area,Tsim[i],Qf)
}
Vr.Tsim

## Criar vetor de vazões firmes, intervalo de 10% da Qm até 80% da Qm
(Qfsim <- seq(Qm_ANO*0.1,Qm_ANO*0.80,l=11))

## Volume do reservatório para cada vazão firme, com tempo de retorno (T) fixo, estabelecido
## pelo usuário
Vr.Qfsim <- c()
for(i in 1:length(Qfsim)){
    Vr.Qfsim[i] <- Vr(mu.ANO.ko,sigma.ANO.ko,area,T,Qfsim[i])
}
Vr.Qfsim

## Gráficos T vs Vr e Vr vs Qf
par(mfrow = c(1,2))
plot(Tsim,Vr.Tsim,type = "l",ylab = expression(V[r]~~(m^{3}~10^6)),xlab="T (anos)",mgp=c(2.5,1,0))
plot(Qfsim,Vr.Qfsim,type = "l",ylab =expression(V[r]~~(m^3~10^6)),xlab=expression(Q[f]~~(m^3~s^{-1})) ,mgp=c(2.5,1,0))

##-----------------------------------------------------------------------------##
## Com a curva Vr vs Qf é possível que o usuário insira o volume útil armazenado de água
## do seu reservatório (Vol.util=Vol.total-Vol.morto) e assim é calculado a Qf para este volume

## Reservatório, volume útil (m^3), dado fornecido pelo usuário
vol <- 38000

## Criar função de interpolação Vr vs Qf
fun.spline <- splinefun(Vr.Qfsim, Qfsim)
Vr.Qfsim
Qfsim
## Vazão firme possível com reservatório existente
(Qf.possivel <- fun.spline(vol/10^6))

## Graficos da simulação
plot(Vr.Qfsim, Qfsim,type = "l",ylab = expression(Q[f]~~(m^3~s^{-1})),xlab=expression(V[r]~~(m^3~10^6)),mgp=c(2.5,1,0))
points(vol/10^6,Qf.possivel,pch=19)
##lines(Vr.Qfsim,fun.spline(Vr.Qfsim),col=3,) #Interpolação sobreposta a curva Vr vs Qf

## Relação da vazão regularizada com a máxima possível de ser 90 % da Qm
Qf.possivel/(Qm_ANO)

##-----------------------------------------------------------------------------##
## Final do reservatório, lembrando que os calculos para reservatório é somente para os
## parâmetros ANO




## Essa parte é só para conferir a interpolação dos parâmetros
## No geoR para conferir
require(geoR)
(dados <- cbind(dados,res))
b <-as.geodata(dados,coords.col = 3:4, data.col = 7,covar.col = c(5,8))
names(dados)
(b$data)
summary(b)

tau2.mu.ANO
sigma2.mu.ANO
phi.mu.ANO

likkc <- likfit(b,cov.model = "gau",ini=c(sigma2.mu.ANO,phi.mu.ANO),nug=tau2.mu.ANO,trend =~Prec_media+Param_ANO.sigma)
likkc$parameters.summary
likkc$tausq
kc <- krige.conv(b, locations=Xo,krige=krige.control(obj=likkc))

kc$predict


c <-as.geodata(dados,coords.col = 3:4, data.col = 9 ,covar.col = c(3,10))
names(dados)
(c$data)
likkc <- likfit(c,cov.model = "sph",ini=c(0.04,47.57), nug=0.02,trend=~coords[,1]+Param_DJF.sigma)
likkc
kc <- krige.conv(c, locations=Xo,krige=krige.control(obj=likkc))
kc$predict

mu.inv<-(((lambda*kc$predict)+1)^(1/lambda))-5 #inverter e diminuir com 5 para encontrar a varivel original
mu.inv



summary(likkc)
likkc
t <- trend.spatial(trend=~Prec_media+Param_ANO.sigma,b)
krige=krige.control(obj=Mu.ANOlf14$gau)


svd(
     # The function is currently defined as
     function(X, tol = sqrt(.Machine$double.eps))
     {
     ## Generalized Inverse of a Matrix
       dnx <- dimnames(X)
       if(is.null(dnx)) dnx <- vector("list", 2)
       s <- svd(X)
       nz <- s$d > tol * s$d[1]
       structure(
         if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
         dimnames = dnx[2:1])
     }
     ## End(Not run)

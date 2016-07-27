# ABRIR ARQUIVO/BANCO DE DADOS
options(OutDec = ",",digits = 7)
dir()
setwd("~/MEGA/Doutorado/Rotinas R/Tese/Multivar/")

Vars<-read.csv("Var-BHs.csv",sep = "",dec=",", head=T)[,7:28]
names(Vars)
head(Vars)
summary(Vars)
plot(Vars)

distancia <- dist(Vars, method = "euclidean", diag=F)

agrupamento <- hclust(distancia, method="ward.D")
plot(agrupamento, hang=0, main="", axes=T, xlab="", labels=F, frame.plot=T)

hcd <- as.dendrogram(agrupamento)
par(mfrow=c(1,3))
plot(cut(hcd, h=50000)$lower[[1]],
    main="Corte em h=10") # Primeiro grupo
plot(cut(hcd, h=10)$lower[[2]],
     main="Corte em h=10") # Segundo grupo
plot(cut(hcd, h=10)$lower[[3]],
     main="Corte em h=10") # Terceiro grupo

grupos <- cutree(agrupamento, h=10)
x <- cbind(parametros, grupos)
grupo1 <- subset(x, grupos == 1)
grupo2 <- subset(x, grupos == 2)
grupo3 <- subset(x, grupos == 3)
rbind(colMeans(grupo1), colMeans(grupo2), colMeans(grupo3))
table(Vars[,1])

(Antas <- subset(Vars,Vars[,1]=="Antas"))

Macrobac <- split(Vars,Vars[,1])

corre <- cor(Vars)
str(corre)
corre[,4]

colnames(Vars) <- c("1","2","Mu","Sigma",as.character(seq(5,25)))

colnames(Vars)


resu <- prcomp(Vars,scale. = TRUE)
summ <- summary(resu)
attributes(summ)
summ$sdev


biplot(resu,col = c("gray", "black"),xlim = c(-0.4,0.3),ylim = c(-0.4,0.3))

summary(resu)
str(Vars)
predict(resu)[,3]
resu$x[,3]

Vars1 <- Vars[,c(1,6:15,17,19,21,22)]
resu1 <- prcomp(Vars1,scale. = TRUE)
resu1
biplot(resu1,col = c("gray", "black"),xlim = c(-0.4,0.3),ylim = c(-0.4,0.3))
summary(resu1)

Vars2 <- Vars1[,-c(5,8,9)]
resu2 <- prcomp(Vars2,scale. = TRUE)
resu2
biplot(resu2,col = c("gray", "black"))#,xlim = c(-0.4,0.3),ylim = c(-0.4,0.3))
summary(resu2)

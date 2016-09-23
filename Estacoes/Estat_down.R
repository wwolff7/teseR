## para instalar pacotes 
## carregando pacotes necessários

library(hydroTSM)
library(httr)
library(XML)
library(extrafont)
loadfonts()

getwd()
setwd("/home/wagner/MEGA/Doutorado/Rotinas R/Tese/Estacoes")

## Download automático dos dados
## link da ANA-HIDROWEB
baseurl <-c("http://hidroweb.ana.gov.br/Estacao.asp?Codigo=", "&CriaArq=true&TipoArq=1")

## vetor com código das estações
estacoes <- c("71550000","71498000","71383000","71350001","71300000","71250000",
              "71200000","70500000","70300000","70200000","70100000","86100000",
              "72870000","72849000","72810000","72715000","72680000","72630000",
              "72430000","83675000","83345000","83029900","83660000","83800002",
              "83892998","83900000","83440000","83690000","83250000","83880000",
              "83300200","83892990","83105000","83050000","83677000","83069900",
              "83520000","82320000","82350000","82370000","82549000","82770000",
              "65180000","74320000","74295000","74100000","74270000","65094500",
              "65095000","65135000","84820000","84853000","84800000","84949000",
              "74370000","74470000","84559800","84520000","84249998","84300000",
              "84580000","84551000","84520010","84560000","84598002","84100000",
              "82270050","73350000","73300000","73330000","84970000","84071000",
              "84095000","84095500","65925000","73900000","73770000","73765000",
              "73820000","73780000","73690001","73960000","73600000")             


## looping para download dos arquivos das estações flu
for (est in estacoes){
  r <- POST(url = paste0(baseurl[1], est, baseurl[2]), body = list(cboTipoReg = "9"), encode = "form")

  if (r$status_code == 200) {
      cont <- content(r, as = "text",encoding = "ISO-8859-1")
    # readLines(textConnection(cont))
    arquivo <- unlist(regmatches(cont, gregexpr("ARQ.+/VAZOES.ZIP", cont)))
    }
   
  if(length(arquivo) > 0){ 
    (arq.url <- paste0("http://hidroweb.ana.gov.br/", arquivo))
    download.file(url = arq.url 
                  ,destfile = gsub("VAZOES.ZIP", paste0(est, ".zip"), basename(arq.url))
                  ,mode = "wb")
    cat("Arquivo", est, "salvo com sucesso.\n")
  } else {
    cat("*** Arquivo", est, "sem dados de Vazão.***\n")
  }
}

## Organização dos dados 

files_csv <- grep("csv",dir(),value = T)
Nomes <- gsub(".csv","",files_csv)


Estat <- list()
for(i in 1:length(files_csv)){
    Estat[[i]] <- suppressWarnings(read.zoo(files_csv[i],format="%d/%m/%Y",sep=";",dec=",",head=T))

    ## Verificar e índice de datas repetidas
    ix <- !duplicated(time(Estat[[i]]), fromLast = TRUE)
    anyDuplicated(time(Estat[[i]]))

    ## Completar as datas faltantes e remover as repetidas
    Estat[[i]] <- merge(Estat[[i]][ix],zoo(,seq(start(Estat[[i]]),end(Estat[[i]]),by="month")), all=TRUE)

    ## Características da série
    serie <- seq.Date(start(Estat[[i]]),end(Estat[[i]])+30,by="day")
    year <- as.numeric(unique(format(serie,"%Y")))
    month <- sort(unique(format(serie,"%m")))
    day <- sort(unique(format(serie,"%d")))

    ## Criar datas falsas, todos os meses com 31 dias
    dataf <- sort(paste(as.character(sort(c(rep(year[1],(12-as.numeric(format(serie[1],"%m"))+1)*31), rep(year[2:length(year)],31*12)))),
                   c(sort(rep(sprintf("%02d", format(serie[1],"%m"):12),31)),rep(month,31*(year[length(year)]-year[1]))),
                   c(rep(day,(12-as.numeric(format(serie[1],"%m"))+1)),rep(day,12*(year[length(year)]-year[1]))),
                   sep="-"))

    anyDuplicated(dataf)

    ## Criar zoo objeto para ajustar as datas verdadeiras e excluir as falsas
    Estat[[i]] <- suppressWarnings(zoo(as.vector(t(Estat[[i]])),as.Date(dataf))[-c(length(serie)+1:length(dataf)),])

}

names(Estat) <- Nomes

Todas <- do.call(cbind.zoo, Estat)

pdf("Figuras/Estat_dados.pdf",onefile = T, width=27/2.54, height=35/2.54,paper = "special",colormodel="grey",family
    = "CM Roman")

#par(mar=c(2.3,2.3,.5,.5))
mat <- matrixplot(dwi(Todas, var.type="Days"),ColorRamp="Temperature",aspect=1.5,
           colorkey=list(labels=list(cex=1.5)))

str(mat)

data <- as.character(seq(1928,2016,4))
space <- matrix(rep("",3*length(data)),nrow=length(data),ncol=3)

mat$x.scales$labels <- c(as.vector(t(cbind(data,space))),2016)
mat$x.scales$cex <- c(1.2,1.2)
mat$y.scales$cex <- c(1.2,1.2)

plot(mat)

dev.off()
embed_fonts("Figuras/Estat_dados.pdf",outfile = "Figuras/Estat_dados.pdf")



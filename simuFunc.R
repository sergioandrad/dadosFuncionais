## Simulacao do paper de Sorensen
set.seed(2)
library(fda)
library(tidyverse)
library(pracma)
library(refund)

# Exemplo de base de b-splines cubicos
create.bspline.basis(rangeval= c(0,1), norder = 4, breaks = c(0, .25, .5, .75, 1)) %>% plot()

# Gerando as bases em um grid discreto e criando coeficientes da expansao em bases
arg.t     <- seq(0,10, length.out= 256)
base      <- fda::bsplineS(x = arg.t, breaks = seq(0,10, length.out = 11), norder=4) %>% data.frame() %>% t()
C         <- matrix(rnorm(n = 150*13), nrow = 150, ncol = 13) %*% matrix(runif(n = 13*13), ncol = 13, nrow = 13)

# Gerando amostra da funcao trajetoria a partir da representacoa em bases e gerando ruido
X         <- C %*% base

# Gerando o ruido gaussiano
ruido     <- matrix(rnorm(n = 256*150, mean = 0, sd = .5), nrow = 150, ncol = 256) 

# Gerando amostra da funcao coeficiente
betatFunc <- function(x){-pnorm(x, 2, 0.3) + 3*pnorm(x, 5, 0.4) + pnorm(x, 7.5, 0.5)} # Funcao do paper
betat     <- apply(as.matrix(seq(0,10, length.out= 256), nrow = 256, ncol =1), MARGIN = c(1,2), FUN = betatFunc)

# Calculando o logito pelo produto de x(t) e \beta(t) interno usando a regra do trapezio
Xbeta <- matrix(nrow = 150, ncol = 256)
for(i in 1:150){
  Xbeta[i,1:256] <- X[i,1:256]*betat
}
intXbeta <- apply(Xbeta, 1, pracma::trapz)

# Funcao que dummeriza  as probabilidades atraves do logito com quebra em p=0.5 
binary <- function(x){if(x >= 0){1}else{0}}
# Binarizando a variavel 
Y      <- apply(data.frame(intXbeta), c(1,2), binary)


# Quebrando a amostra em treino e teste e estimando a functional logistic regression
# lf() gera o termo funcional da regressao
binary2  <- function(x){if(x >= 0.5){1}else{0}}

accuracy <- data.frame()
iter     <- 0
metodo   <- 'REML' # Insira aqui o metodo desejado e compativel com o GAM
for(i in 10^seq(from = -1, to = 3, by = 1)){
  iter             <- iter + 1
  for (j in c(1:100)) {
  set.seed(j)
  ruido            <- matrix(rnorm(n = 256*150, mean = 0, sd = .5), nrow = 150, ncol = 256) 
  W                <- X + i*ruido
  Ytrain           <- Y[1:100,]
  Ytest            <- Y[101:150,]
  Wtrain           <- W[1:100,]
  Wtest            <- W[101:150,]
  reg              <- pfr(Ytrain ~ lf(Wtrain), family = binomial(), drop.intercept = FALSE, method = metodo)
  Ytestp           <- predict(reg, list('Wtrain'= Wtest), type='response') %>% as.data.frame()
  Ytestpbin        <- apply(Ytestp, c(1,2), binary2)
  accuracy[iter,j] <- sum(abs(Ytest-Ytestpbin))/50
}}
t(accuracy) %>% boxplot()
apply(accuracy, MARGIN = 1, FUN = mean)

# Verificando as dimensoes das matrizes
C     %>% dim() # 150 repeticoes dos coeficientes das 13 curvas
base  %>% dim() # 13 curvas (K_x) discretizadas em 256 pontos igualmente espaçados
ruido %>% dim() # ruido nos 256 pontos de cada uma das 150 curvas
X     %>% dim() # 256 pontos de cada uma das 150 curvas
betat %>% dim() # Funcao coeficiente avaliada em cada ponto
Xbeta %>% dim()

#Plotando o gráfico das curvas
tX        <- t(X) %>% as.data.frame()
tX[,151]  <- seq(from = 1, to = 256, by = 1) 
curvesX   <- tX %>% rename(t = 151) %>% group_by(t) %>% pivot_longer(cols = -t) %>% rename(key = name)
Ydf       <- Y %>% as.data.frame()
Ydf[,2]   <- paste0('V', seq(from =1, to= 150, by = 1))
Ydf       <- Ydf %>% rename(Y = 1, key = 2) 
fontSize  <- 15

curvesX                      %>% #filter(key %in% c(paste0('V', seq(from =101, to= 150, by = 1)))) %>%
  left_join(y=Ydf, by='key') %>%
  group_by(Y)                %>%
  ggplot() + geom_line(aes(x=t, y = value, color = Y, shape = key),size=1) + 
  theme( axis.text        = element_text(size = fontSize),
         axis.text.x      = element_text(size = fontSize),
         axis.title       = element_text(size = fontSize),
         legend.position  = "none",
         panel.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         plot.title       = element_text(size=1.5*fontSize),
         strip.text       = element_text(size = fontSize)) + ggtitle('') + xlab('t') + ylab('X(t)')

# Utilizando ondaletas
library(wavelets)
library(glmnet)

accuracyWav <- data.frame()
iter     <- 0
for(i in 10^seq(from = -1, to = 3, by = 1)){
  iter             <- iter + 1
  for (j in c(1:100)) {
    set.seed(j)
    ruido            <- matrix(rnorm(n = 256*150, mean = 0, sd = 0.5), nrow = 150, ncol = 256) 
    W                <- X + i*ruido
    Ytrain           <- Y[1:100,]
    Ytest            <- Y[101:150,]
    Wtrain           <- W[1:100,]
    Wtest            <- W[101:150,]
    Xwav <- data.frame()
    for(i2 in c(1:100)){
      waveTrW        <- dwt(Wtrain[i2,], filter="la8", boundary="periodic", fast=TRUE)
      Ww             <- reduce(slot(waveTrW, 'W'), rbind) # Coefs das wavelets
      Wv             <- reduce(slot(waveTrW, 'V'), rbind) # Coefs das funcoes escala
      Xwav[i2,1:480] <- rbind(Ww,Wv) %>% t()
    }
    Xwavmat  <- Xwav %>% as.matrix()
    waveReg  <- glmnet(x = Xwavmat, y=Ytrain, family = 'binomial', nfolds = 1, nlamda = 1, lambda = 0.0001)
    Xwavtest <- data.frame()
    for(i3 in c(1:50)){
      waveTrWtest        <- dwt(Wtest[i3,], filter="la10", boundary="periodic", fast=TRUE)
      Wwtest             <- reduce(slot(waveTrWtest, 'W'), rbind) # Coefs das wavelets
      Wvtest             <- reduce(slot(waveTrWtest, 'V'), rbind) # Coefs das funcoes escala
      Xwavtest[i3,1:480] <- rbind(Wwtest, Wvtest) %>% t()
    }
  Xwavtestmat <- as.matrix(Xwavtest)
  Ytest            <- predict(waveReg, Xwavtestmat, type = 'response') %>% as.data.frame()
  Ytestbin         <- apply(Ytest, c(1,2), binary)
  accuracyWav[iter,j] <- sum(abs(Ytestbin-Ytest))/50
  }
}
accuracyWav
apply(accuracyWav, MARGIN = 1, FUN = mean)

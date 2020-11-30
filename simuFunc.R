## Simulacao do paper de Sorensen
set.seed(1)
library(fda)
library(tidyverse)
library(pracma)
library(refund)

# Avaliando as bases em um grid discreto e criando coeficientes da expansao em bases
arg.t     <- seq(0,10, length.out= 256)
base      <- fda::bsplineS(x = arg.t, breaks = seq(0,10, length.out = 11), norder=4) %>% data.frame() %>% t()
C         <- matrix(rnorm(n = 150*13), nrow = 150, ncol = 13) %*% matrix(runif(n = 13*13), ncol = 13, nrow = 13)

C     %>% dim() # 150 repeticoes dos coeficientes das 13 curvas
base  %>% dim() # 13 curvas (K_x) discretizadas em 256 pontos igualmente espaçados

# Gerando amostra da funcao trajetoria a partir da representacoa em bases e gerando ruido
X         <- C %*% base
ruido     <- matrix(rnorm(n = 256*150, mean = 0, sd = 0.1), nrow = 150, ncol = 256) 

# Gerando amostra da funcao coeficiente
betatFunc <- function(x){-pnorm(x, 2, 0.3) + 3*pnorm(x, 5, 0.4) + pnorm(x, 7.5, 0.5)} # Funcao do paper
betat     <- apply(as.matrix(seq(0,10, length.out= 256), nrow = 256, ncol =1), MARGIN = c(1,2), FUN = betatFunc)

# Verificando as dimensoes das matrizes
ruido %>% dim() # ruido nos 256 pontos de cada uma das 150 curvas
X     %>% dim() # 256 pontos de cada uma das 150 curvas
betat %>% dim() # Funcao coeficiente avaliada em cada ponto

# Calculando o logito por regra do trapezio
Xbeta <- matrix(nrow = 150, ncol = 256)
for(i in 1:150){
  Xbeta[i,1:256] <- X[i,1:256]*betat
}

Xbeta %>% dim()
intXbeta <- apply(Xbeta, 1, pracma::trapz)

# Funcao que dummeriza  as probabilidades atraves do logito com quebra em p=0.5 
binary <- function(x){if(x >= 0){1}else{0}}
# Binarizando a variavel 
Y   <- apply(data.frame(intXbeta), c(1,2), binary)

#Plotando o gráfico das curvas
tX        <- t(X) %>% as.data.frame()
tX[,151]  <- seq(from = 1, to = 256, by = 1) 
curvesX   <- tX %>% rename(t = 151) %>% group_by(t) %>% pivot_longer(cols = -t) %>% rename(key = name)
Ydf       <- Y %>% as.data.frame()
Ydf[,2]   <- paste0('V', seq(from =1, to= 150, by = 1))
Ydf       <- Ydf %>% rename(Y = 1, key = 2)
fontSize  <- 15
curvesX                      %>%
  left_join(y=Ydf, by='key') %>%
  group_by(key)              %>%
  ggplot() + geom_line(aes(x=t, y = value, color = key),size=1) + 
    theme( axis.text        = element_text(size = fontSize),
           axis.text.x      = element_text(size = fontSize),
           axis.title       = element_text(size = fontSize, face = "bold"),
           legend.position  = "none",
           panel.background = element_blank(),
           panel.grid.major = element_blank(),
           plot.title       = element_text(size=1.5*fontSize),
           panel.grid.minor = element_blank(),
           strip.text       = element_text(size = fontSize)) + ggtitle('Curvas geradas por combinações lineares de B-Splines cúbicos') + xlab('t') + ylab('X(t)')

# Estimando a penalized functional logistic regression
# lf() gera o termo funcional da regressao
W <- X
Ytrain <- Y[1:100,]
Ytest  <- Y[101:150,]
Wtrain <- W[1:100,]
Wtest  <- W[101:150,]

reg <- pfr(Ytrain ~ lf(Wtrain), family = binomial(), drop.intercept = TRUE, method = 'ML')
reg               %>% summary()
reg$coefficients
reg$residuals     %>% boxplot()
reg$fitted.values %>% plot()

# Comparando out-of-sample
binary2  <- function(x){if(x >= 0.5){1}else{0}}
Ytest    <- predict(reg, list(Wtrain = Wtest)) %>% as.data.frame()
Ytestbin <- apply(Ytest, c(1,2), binary)
sum(abs(Ytestbin-Y[101:150,]))/50

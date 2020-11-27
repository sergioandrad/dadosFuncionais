## Simulacao do paper de Sorensen
set.seed(1)
library(fda)
library(tidyverse)
library(pracma)

# Criando bases
bsplinesbasis <- fda::create.bspline.basis(rangeval=c(0,10), nbasis=13, norder=4, names="bspl")

# Avaliando as bases em um grid discreto e criando coeficientes da expansao em bases
X         <- fda::bsplineS(x = seq(0,10, length.out= 256), breaks = seq(0,10, length.out = 11), norder=4) %>% data.frame()
C         <- matrix(rnorm(n = 150*13), nrow = 150, ncol = 13) %*% matrix(runif(n = 13*13), ncol = 13, nrow = 13)

t(C)  %>% dim()
X     %>% dim()

# Gerando amostra da funcao trajetoria a partir da representacoa em bases e gerando ruido
XC        <- (as.matrix(X, nrow = 256, ncol = 13) %*% t(C))
ruido     <- as.matrix(rnorm(n = 256*150), nrow = 256, ncol = 150) 

# Gerando amostra da funcao coeficiente
betatFunc <- function(x){-pnorm(x, 2, 0.3) + 3*pnorm(x, 5, 0.4) + pnorm(x, 7.5, 0.5)} # Funcao do paper
betat     <- apply(as.matrix(seq(0,10, length.out= 256), nrow = 256, ncol =1), MARGIN = c(1,2), FUN = betatFunc)

# Verificando as dimensoes das matrizes
ruido %>% dim()
XC    %>% dim()
betat %>% dim()

# Calculando o logito por regra do trapezio
columnwiseTrap <- function(x){trapzmat(X = x, Y = betat)}
logitY         <- apply(XC, 2, columnwiseTrap)

# Funcao que dummeriza  as probabilidades atraves do logito com quebra em p=0.5 
binary <- function(x){if(x >= 0){1}else{0}}

# Atribuindo as covariaveis com ruido e a variavel resposta
W  <- XC + ruido
Y  <- apply(data.frame(logitY), c(1,2), binary)

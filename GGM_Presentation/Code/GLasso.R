# Primal Set-Up 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(pacman)

# Set-Up
p_load(glasso, dplyr)
frobenius = function(M) sum(M^2)


# 
x<-matrix(rnorm(20*50),ncol=20)
s<- var(x)
a <-glasso(s, rho=0)
aa<-glasso(s,rho=.02, w.init=a$w, wi.init=a$wi)

all.equal(a$wi, t(a$wi))
.i = a$w %*% a$wi
diag(.i)
.i[col(.i) != row(.i)] %>% summary()


set.seed(1)
test = glasso(s, rho=0.4)
test2 = glasso(s, rho=0.4)





# 
load("preprocessed_Allen_GIJOE.rdata")
.list = load("preprocessed_Allen_GIJOE.rdata")




x = seq(1, 5, length.out=100)
y = exp(x)
y2 = 4*exp(x)+90
plot(x, y, "l"
     , xlab = "Time (seconds)"
     , ylab = "Love (quettaHearts) / (teraHearts)")
lines(x, y2, "l", col="blue")

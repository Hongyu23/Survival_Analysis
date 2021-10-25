#Question 2 plot
sigmoid <-function(x) 0.9*x-log(exp(0.5*x)+exp(x)+exp(0.4*x)+exp(2*x))-log(exp(0.4*x)+exp(2*x))

curve(sigmoid,-15,10)

#Question 2 Newton
x <- 0

a1 <- exp(0.5*x)+exp(x)+exp(0.4*x)+exp(2*x)
a2 <- exp(0.4*x)+exp(2*x)
a3 <- 0.5*exp(0.5*x)+exp(x)+0.4*exp(0.4*x)+2*exp(2*x)
a4 <- 0.4*exp(0.4*x)+2*exp(2*x)
a5 <- 0.25*exp(0.5*x)+exp(x)+0.16*exp(0.4*x)+4*exp(2*x)
a6 <- 0.16*exp(0.4*x)+4*exp(2*x)

u <- 0.9-(a3/a1)-(a4/a2)
H <- (-a5/a1)+((a3)^2)/((a1)^2)-(a6/a2)+((a4)^2)/((a2)^2)

for (i in 1:10){
  x = x + (u/H)
}


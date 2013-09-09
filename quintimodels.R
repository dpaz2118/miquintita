library("sfsmisc")
library("deSolve")
options(x11updates=0.25)


#tabla de parametros
tmodel <- read.table("params_model.dat",col.names=c("w0","wm","am","dm"),
		     row.names=c("inv1","inv2","sugra","exp2","as","cnr"))
assign("tmodel",tmodel, envir = .GlobalEnv)

#estilo de linea para cada modelo
lsty <- matrix(c(1,2,3,4,5,6),ncol=1)
rownames(lsty)=c("inv1","inv2","sugra","exp2","as","cnr")
lsty <- as.data.frame(lsty)
assign("lsty",lsty, envir = .GlobalEnv)

#parametro de densidad de materia
assign("omega.mat",0.3, envir = .GlobalEnv)


plot.wall <- function(new.dev) {
  if(missing(new.dev)){new.dev=FALSE}
  plot.w("cnr",new.dev)
  plot.w("inv1")
  plot.w("inv2")
  plot.w("sugra")
  plot.w("exp2")
  plot.w("as")

}

plot.xall <- function(new.dev) {
  if(missing(new.dev)){new.dev=FALSE}
  plot.xfun("cnr",new.dev)
  plot.xfun("inv1")
  plot.xfun("inv2")
  plot.xfun("sugra")
  plot.xfun("exp2")
  plot.xfun("as")

}

plot.w <- function(model.name,new.dev) {

 a=seq(-2,0,0.01)
 a=10.0^a

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    plot(a, w.model(a,model.name), type = "l",
	 xlab="a",ylab="w(a)",
	 xaxt = "n",yaxt="n", log = "x",lty=lsty[model.name,])
    eaxis(1,at=c(1E-2,1E-1,1))
    eaxis(2,small.mult=4)
    eaxis(3,labels=FALSE,at=c(1E-2,1E-1,1))
    eaxis(4,labels=FALSE)
 } else {
    lines(a,w.model(a,model.name),lty=lsty[model.name,])
 }
	
}

plot.xfun <- function(model.name,new.dev) {

 a=seq(-2,0,0.01)
 a=10.0^a

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    plot(a, X.fun(a,model.name), type = "l",
	 xlab="a",ylab="X(a)",
	 xaxt = "n",yaxt="n", log = "x",lty=lsty[model.name,])
    eaxis(1,at=c(1E-2,1E-1,1))
    eaxis(2,small.mult=4)
    eaxis(3,labels=FALSE,at=c(1E-2,1E-1,1))
    eaxis(4,labels=FALSE)
 } else {
    lines(a,X.fun(a,model.name),lty=lsty[model.name,])
 }
	
}

X.fun <- function(a,model.name) {
	res=a*0.0
	n <-length(a)
	for(i in 1:n){
		ff=integrate(w.model,lower=a[i],upper=1,model.name)
		res[i]=ff$value
	}

	res=omega.mat/(1-omega.mat)*exp(-3*res)
	return(res)
}

w.model <-function(a,model.name) {
	
   w0<-tmodel[model.name,]$w0
   wm<-tmodel[model.name,]$wm
   am<-tmodel[model.name,]$am
   dm<-tmodel[model.name,]$dm

   c4 = exp(am/dm)
   c3 = exp(1/dm)
   c2 = (wm-w0)*(1+c4)/(1-c3)
   c1 = w0

   x=exp(-a/dm)
   w.phi= c1+c2*(1-c3*x)/(1+c4*x)

   return(w.phi)
}

evol.grow.G <- function(model.name) {
	parms <- c(name = model.name)
        times <- seq(0.01, 1, length = 101)
	xstart <- c(x1 = 0.001, x2 = 0.001) #condiciones iniciales para el factor de 
	                            #crecimiento y sus derivadas
        out <- ode(xstart,times,dgrow.G,parms)

}


dgrow.G <- function(a,yfunc,parms)
{
         with(as.list(c(parms, yfunc)), {
	      #se queja de que x2 es no numerico o que es character hay algo roto
         dx1 <- 0.0#x2
         dx2 <- 0.0#-(7/2-3/2*w.model(a,name)/(1+X.fun(a,name)))*x2/a
	 dx2 <- 0.0# dx2 -3/2*(1-w.model(a,name))/(1+X.fun(a,name))*x1/a/a
         res <- c(dx1, dx2)
         return(list(res))
       })

}

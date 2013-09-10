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

plot.dall <- function(new.dev) {
  if(missing(new.dev)){new.dev=FALSE}
  plot.dgrow("cnr",new.dev)
  plot.dgrow("inv1")
  plot.dgrow("inv2")
  plot.dgrow("sugra")
  plot.dgrow("exp2")
  plot.dgrow("as")

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

plot.dgrow <- function(model.name,new.dev) {

 a <- seq(0.001, 1, length = 301)

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    d <- evol.grow.G(model.name)
    d <- d/d[301]
    plot(a, d, type = "l",
	 xlab="a",ylab="D(a)",
	 xaxt = "n",yaxt="n", lty=lsty[model.name,])
    eaxis(1)#,at=c(1E-2,1E-1,1))
    eaxis(2)#,small.mult=4)
    eaxis(3)#,labels=FALSE,at=c(1E-2,1E-1,1))
    eaxis(4)#,labels=FALSE)
 } else {
    d <- evol.grow.G(model.name)
    d <- d/d[301]
    lines(a,d,lty=lsty[model.name,])
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

#integrado en el growfactor para cada modelo
evol.grow.G <- function(model.name) {
        atimes <- seq(0.001, 1, length = 301)
	xstart <- c(D =0.1, Dprime =1.0) #condiciones iniciales para el factor de 
	                            #crecimiento y sus derivadas
        out <- ode(xstart,atimes,dgrow.G,model.name)
	sld <- as.data.frame(out)
	a=sld$time
	dfact=sld$D*a
	return(dfact)
}

dgrow.G <- function(a,yfunc,model.name)
{
         with(as.list(c(yfunc)), {
         dD <- Dprime
         dDprime <- -(7/2-3/2*w.model(a,model.name)/(1+X.fun(a,model.name)))*Dprime/a
	 dDprime <- dDprime -3/2*(1-w.model(a,model.name))/(1+X.fun(a,model.name))*D/a/a
         res <- c(dD, dDprime)
         return(list(res))
       })
}
  

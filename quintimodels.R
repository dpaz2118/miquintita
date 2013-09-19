library("sfsmisc")
library("deSolve")
options(x11updates=0.25)


#tabla de parametros
tmodel <- read.table("params_model.dat",col.names=c("w0","wm","am","dm"),
		     row.names=c("inv1","inv2","sugra","exp2","as","cnr"))
assign("tmodel",tmodel, envir = .GlobalEnv)

#estilo de linea para cada modelo
lsty <- matrix(c(1,2,3,4,5,6,7),ncol=1)
rownames(lsty)=c("lcdm","inv1","inv2","sugra","exp2","as","cnr")
lsty <- as.data.frame(lsty)
assign("lsty",lsty, envir = .GlobalEnv)

#parametro de densidad de materia
assign("omega.mat",0.26, envir = .GlobalEnv)
assign("H0",76, envir = .GlobalEnv)

#Los H/Hl no dan igual a Ellise por un factor 0.9565217
#parece mejorar cuando se aumenta el omega.mat, por ejemplo 0.32

#Graficos generales #######################
plot.oall <- function(new.dev) 
{
  if(missing(new.dev)){new.dev=FALSE}
  plot.omega("cnr",new.dev)
  plot.omega("inv1")
  plot.omega("inv2")
  plot.omega("sugra")
  plot.omega("exp2")
  plot.omega("as")
  plot.omega("lcdm")

}

plot.wall <- function(new.dev) 
{
  if(missing(new.dev)){new.dev=FALSE}
  plot.w("cnr",new.dev)
  plot.w("inv1")
  plot.w("inv2")
  plot.w("sugra")
  plot.w("exp2")
  plot.w("as")

}

plot.Hall <- function(new.dev) 
{
  if(missing(new.dev)){new.dev=FALSE}
  plot.H("cnr",new.dev)
  plot.H("inv1")
  plot.H("inv2")
  plot.H("sugra")
  plot.H("exp2")
  plot.H("as")

}

plot.xall <- function(new.dev) 
{
  if(missing(new.dev)){new.dev=FALSE}
  plot.xfun("cnr",new.dev)
  plot.xfun("inv1")
  plot.xfun("inv2")
  plot.xfun("sugra")
  plot.xfun("exp2")
  plot.xfun("as")

}

plot.dall <- function(new.dev) 
{
  plot.dgrow("cnr",new.dev)
  plot.dgrow("inv1")
  plot.dgrow("inv2")
  plot.dgrow("sugra")
  plot.dgrow("exp2")
  plot.dgrow("as")

}

plot.dpall <- function(new.dev) 
{
  plot.dprimegrow("cnr",new.dev)
  plot.dprimegrow("inv1")
  plot.dprimegrow("inv2")
  plot.dprimegrow("sugra")
  plot.dprimegrow("exp2")
  plot.dprimegrow("as")

}
#Graficos individuales ######################
plot.w <- function(model.name,new.dev) 
{
 a=seq(-2,0,0.01)
 a=10.0^a

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    plot(a, w.model(log(a),model.name), type = "l",
	 xlab="a",ylab="w(a)",
	 xaxt = "n",yaxt="n", log = "x",lty=lsty[model.name,])
    eaxis(1,at=c(1E-2,1E-1,1))
    eaxis(2,small.mult=4)
    eaxis(3,labels=FALSE,at=c(1E-2,1E-1,1))
    eaxis(4,labels=FALSE)
 } else {
    lines(a,w.model(log(a),model.name),lty=lsty[model.name,])
 }
	
}

plot.xfun <- function(model.name,new.dev) 
{

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

plot.H <- function(model.name,new.dev) 
{

 a=seq(-2,0,0.01)
 a=10.0^a
 z=1/a-1
 hz=hubble.par(a,model.name)
 hl=hubble.par(a,"lcdm")
 #hl=sqrt(omega.mat/a/a/a+(1-omega.mat))*H0 #hubble.par(a,"lcdm")

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    plot(z,hz/hl , type = "l",
	 xlab="z",ylab="H/Hlcdm",
	 xaxt = "n",yaxt="n",lty=lsty[model.name,],ylim=c(0.9,1.3),xlim=c(0,10))
    eaxis(1)
    eaxis(2)
    eaxis(3,labels=FALSE)
    eaxis(4,labels=FALSE)
 } else {
    lines(z,hz/hl,lty=lsty[model.name,])
 }
	
}

plot.omega <- function(model.name,new.dev) 
{
 a=seq(-3,0,0.01)
 a=10.0^a
 om=omega.de(a,model.name)

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    plot(a,om , type = "l",log="xy",
	 xlab="a",ylab="omega.DE",
	 xaxt = "n",yaxt="n",lty=lsty[model.name,],ylim=c(0.001,1.0),xlim=c(0.002,1.0))
    eaxis(1,at=c(1E-2,1E-1,1))
    eaxis(2,at=c(1E-3,1E-2,1E-1,1))
    eaxis(3,labels=FALSE)
    eaxis(4,labels=FALSE)
 } else {
    lines(a,om,lty=lsty[model.name,])
 }
	
}

plot.dgrow <- function(model.name,new.dev)
{

 a <- seq(0.001, 1, length = 301)

 nuevo=TRUE
 if(missing(new.dev)){nuevo=FALSE; new.dev=dev.new}
 if(attr(dev.cur(),"names")== "null device" || nuevo) {
    new.dev()
    par(tcl=1)
    d <- evol.grow.G(a,model.name)
    d <- d/d[301]
    plot(a, d, type = "l",
	 xlab="a",ylab="D(a)",
	 xaxt = "n",yaxt="n", lty=lsty[model.name,])
    eaxis(1)
    eaxis(2)
    eaxis(3,labels=FALSE)
    eaxis(4,labels=FALSE)
 } else {
    d <- evol.grow.G(a,model.name)
    d <- d/d[301]
    lines(a,d,lty=lsty[model.name,])
 }
	
}

plot.dprimegrow <- function(model.name,new.dev)
{

 z <- seq(0., 1.0, length = 101)
 a <- 1/(1+z)
 a <- rev(a)
 z <- rev(z)

 nuevo=TRUE
 dl <- evol.growprime.G(a,"lcdm")
 if(missing(new.dev)){nuevo=FALSE; new.dev=dev.new}
 if(attr(dev.cur(),"names")== "null device" || nuevo) {
    new.dev()
    par(tcl=1)
    d <- evol.growprime.G(a,model.name)
    d <- d/dl
    plot(z, d, type = "l",
	 xlab="z",ylab="dDda/D Normalizado a LCDM",
	 xaxt = "n",yaxt="n", lty=lsty[model.name,],ylim=c(0.7,1.1))
    eaxis(1)
    eaxis(2)
    eaxis(3,labels=FALSE)
    eaxis(4,labels=FALSE)
 } else {
    d <- evol.growprime.G(a,model.name)
    d <- d/dl
    lines(z,d,lty=lsty[model.name,])
 }
	
}

#Funciones utilizadas #######################

omega.de <-function(a,model.name)
{
	fz <- omega.mat/(1-omega.mat)/X.fun(a,model.name)/a/a/a
	Ez2 <- (hubble.par(a,model.name))^2/(H0^2)
	res <- (1-omega.mat)*fz/Ez2
	return(res)
}

hubble.par <-function(a,model.name)
{
	res <- omega.mat/a^3*(1.0+1.0/X.fun(a,model.name))
	res <- H0*sqrt(res)
	return(res)
}

X.fun <- function(a,model.name)
{
	# a viene en lineal
	res=a*0.0
	n <-length(a)
	for(i in 1:n){
		ff=integrate(w.model,lower=log(a[i]),upper=0,model.name)
		res[i]=ff$value
	}

	res=(omega.mat/(1-omega.mat))*exp(-3*res)
	return(res)
}

w.model <-function(a1,model.name)
{
   # a1 viene en logaritmo, por la integral	
   a=exp(a1)
   if(model.name == "lcdm") {
      w.phi=-1.0 + a*0.0
   }else{
      w0<-tmodel[model.name,]$w0
      wm<-tmodel[model.name,]$wm
      am<-tmodel[model.name,]$am
      dm<-tmodel[model.name,]$dm

      q1=(1+exp(am/dm))/(1+exp(-(a-am)/dm))
      q2=(1-exp(-(a-1)/dm))/(1-exp(1/dm))
      
      w.phi=w0+(wm-w0)*q1*q2

      w.phi= w.phi
   }

   return(w.phi)
}

evol.grow.G <- function(atimes,model.name) 
{
	xstart <- c(D =0.1, Dprime =1.0) #condiciones iniciales para el factor de 
	                            #crecimiento y sus derivadas
	at=c(0.001,atimes)
        out <- ode(xstart,at,dgrow.G,model.name,method=rkMethod("rk45f"))
	sld <- as.data.frame(out)
	a=sld$time
	dfact=sld$D*a
	nn=length(dfact)
	dfact=dfact[2:nn]

	return(dfact)
}

evol.growprime.G <- function(atimes,model.name) 
{
	xstart <- c(D =0.1, Dprime =1.0) #condiciones iniciales para el factor de 
	                            #crecimiento y sus derivadas
	at=c(0.001,atimes)
        out <- ode(xstart,at,dgrow.G,model.name,method=rkMethod("rk45f"))
	sld <- as.data.frame(out)
	a <- sld$time
	dfact <- sld$D + a*sld$Dprime
	dfact <- dfact/(sld$D*a)
	nn=length(dfact)
	dfact=dfact[2:nn]
	return(dfact)
}

dgrow.G <- function(a,yfunc,model.name)
{
         with(as.list(c(yfunc)), {
         dD <- Dprime
	 #X.fun come lineal
	 #w.model come log
	 alog=log(a)
         dDprime <- -(7/2-3/2*w.model(alog,model.name)/(1+X.fun(a,model.name)))*Dprime/a
	 dDprime <- dDprime -3/2*(1-w.model(alog,model.name))/(1+X.fun(a,model.name))*D/a/a
         res <- c(dD, dDprime)
         return(list(res))
       })
}
  

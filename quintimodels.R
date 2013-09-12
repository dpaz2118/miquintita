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

 if(missing(new.dev)){new.dev=FALSE}
 if(attr(dev.cur(),"names")== "null device" || new.dev) {
    dev.new()
    par(tcl=1)
    plot(z,hz/hl , type = "l",
	 xlab="a",ylab="X(a)",
	 xaxt = "n",yaxt="n",lty=lsty[model.name,])
    eaxis(1)
    eaxis(2)
    eaxis(3,labels=FALSE)
    eaxis(4,labels=FALSE)
 } else {
    lines(z,hz/hl,lty=lsty[model.name,])
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
    d <- evol.grow.G(model.name)
    d <- d/d[301]
    plot(a, d, type = "l",
	 xlab="a",ylab="D(a)",
	 xaxt = "n",yaxt="n", lty=lsty[model.name,])
    eaxis(1)
    eaxis(2)
    eaxis(3,labels=FALSE)
    eaxis(4,labels=FALSE)
 } else {
    d <- evol.grow.G(model.name)
    d <- d/d[301]
    lines(a,d,lty=lsty[model.name,])
 }
	
}

hubble.par <-function(a,model.name)
{
	res <- omega.mat/a^3 + omega.mat/a^3/X.fun(a,model.name)
	res <- H0*res
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
      w.phi=-1.0
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

#integrado en el growfactor para cada modelo
evol.grow.G <- function(model.name) 
{
        atimes <- seq(0.001, 1, length = 301)
	xstart <- c(D =0.1, Dprime =1.0) #condiciones iniciales para el factor de 
	                            #crecimiento y sus derivadas
        out <- ode(xstart,atimes,dgrow.G,model.name,method=rkMethod("rk45f"))
	sld <- as.data.frame(out)
	a=sld$time
	dfact=sld$D*a
	return(dfact)
}

dgrow.G <- function(a,yfunc,model.name)
{
         with(as.list(c(yfunc)), {
         dD <- Dprime
	 ### guarda!!! muchos cambios en w.model chechear
	 #xfun come lineal
	 #w come log
         dDprime <- -(7/2-3/2*w.model(a,model.name)/(1+X.fun(a,model.name)))*Dprime/a
	 dDprime <- dDprime -3/2*(1-w.model(a,model.name))/(1+X.fun(a,model.name))*D/a/a
         res <- c(dD, dDprime)
         return(list(res))
       })
}
  

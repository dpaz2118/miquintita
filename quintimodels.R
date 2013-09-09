
options(x11updates=0.25)
library("sfsmisc")

#tabla de parametros
tmodel <- read.table("params_model.dat",col.names=c("w0","wm","am","dm"),
		     row.names=c("inv1","inv2","sugra","exp2","as","cnr"))
assign("tmodel",tmodel, envir = .GlobalEnv)

#estilo de linea para cada modelo
lsty <- matrix(c(1,2,3,4,5,6),ncol=1)
rownames(lsty)=c("inv1","inv2","sugra","exp2","as","cnr")
lsty <- as.data.frame(lsty)
assign("lsty",lsty, envir = .GlobalEnv)

w.model <-function(model.name,a) {
	
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

plot.w <- function(model.name) {

 a=seq(-2,0,0.01)
 a=10.0^a

 if(attr(dev.cur(),"names")== "null device") {
    dev.new()
    par(tcl=1)
    plot(a, w.model(model.name,a), type = "l", xaxt = "n",yaxt="n", log = "x",lty=lsty[model.name,])
    eaxis(1,at=c(1E-2,1E-1,1))
    eaxis(2,small.mult=4)
 } else {
    lines(a,w.model(model.name,a),lty=lsty[model.name,])
 }
	
}

require(libamrfile)


cg <- colorRampPalette(c("#003366","#54A9FF","#FFFF00","#FF0000"),space="rgb")(128)


uplot <- function(fab){
  #test plot : color plot of log(|u|,10)
  print(names(fab))
  modu <-sqrt(fab$v[,,1]^2 + fab$v[,,2]^2)
  modu <- ifelse(modu > 5000,5000,modu)
  image(fab$x,fab$y,log(modu,10),add=TRUE,zlim=log(c(1,5000),10),col=cg,useRaster=TRUE)
  rect(min(fab$x),min(fab$y),max(fab$x),max(fab$y),border="black")
}

settoabs <- function(fab)
  {
    #example modification of fab data.
    fab$v <- abs(fab$v)
    for (ic in fab$comp)
      {
        amr.write.fab(fab$amrID, fab$level, fab$ID, fab$v[,,ic], ic, ng=fab$nghost)
      }
  }




pthwID <- amr.load("plot.amundsen.2d.hdf5")


names <- amr.query.compnames(pthwID)
ncomp <- length(names)
comp <- 0:(ncomp-1)
uxc <- comp[names=="xVel"]
uyc <- comp[names=="yVel"]

#multilevel plot
par(xaxs="i",yaxs="i",las=1,mar=c(1,1,1,1),mfrow=c(2,1))
plot(c(0,768e+3),c(0,768e+3),type='n',axes=FALSE)
box()

amr.apply(pthwID, uplot, comp=c(uxc,uyc),maxlev=2)

#modification
amr.apply(pthwID, settoabs, comp=c(uxc,uyc),maxlev=2)

#show the result
plot(c(0,768e+3),c(0,768e+3),type='n',axes=FALSE)
box()
amr.apply(pthwID, uplot, comp=c(uxc,uyc),maxlev=2)

amr.free(pthwID)



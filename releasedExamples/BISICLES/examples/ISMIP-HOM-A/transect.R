require(libamrfile)

transect <- function(fx,fy,amrID,comp){
  n <- length(fx)
  fz <- numeric(n)
  
  maxlev <- amr.query.nlevel(amrID) -1
  maxlev <- 0
  for (lev in 0:maxlev){
    maxfab <- amr.query.nfab(amrID,lev) - 1
    for (fab in 0:maxfab){
      t <- amr.read.fab(amrID,lev,fab,comp,ng=1)
     # print(names(t))
      dx <- t$x[2] - t$x[1]
      dy <- dx
      xr <- range(t$x)#+ c(-dx,dx)
      yr <- range(t$y)#+ c(-dx,dx)
      print (xr)
      print (yr)
      for (i in 1:n){
        if (fx[i] > xr[1] & fx[i] < xr[2] & fy[i] >yr[1] & fy[i] < yr[2]){
                                        #closest cell centers
          ic <- which.min( abs(t$x-fx[i]))
          jc <- which.min( abs(t$y-fy[i]))
          ip <- if(fx[i] > t$x[ic]){1}else{-1}
          jp <- if(fy[i] > t$y[jc]){1}else{-1}
          dcx <- abs(t$x[ic] - fx[i])
          dcy <- abs(t$y[jc] - fy[i])
          wcc <- 0 + 1*(dx - dcx)*(dy - dcy) / (dx*dy)
          wcp <- 0 + 1*(dx - dcx)*(dcy)/ (dx*dy)
          wpc <- 0 + 1*(dcx)*(dy - dcy)/ (dx*dy)
          wpp <- 0 + 1*(dcx)*(dcy)/ (dx*dy)
          fz[i] =  wcc * t$v[ic,jc] +
            wpp * t$v[ic+ip,jc+jp] +
              wcp * t$v[ic,jc+jp] +
                wpc * t$v[ic+ip,jc]
          
        } else {
          #fz[i] =  t$v[which.min( abs(t$x-fx[i])),which.min( abs(t$y-fy[i]))]
        }
      }
    }
  }
  fz
}

pdf("ISMIP-HOM-slippyA%d.pdf",width=7,height=5,onefile=TRUE)
par(mfrow=c(2,3))
par(mar=c(5,5,5,1))

#domain sizes
L <- c(5,10,20,40,80,160)
lower <- c(0.0,0.0,0,0,0,0)
upper <- c(300,300,300,300,300,300)

par(xaxs="i")

for (i in 1:length(L)){

  a <- paste("plot.ISMIP-HOMA.l1l2.0lev.0064.",L[i],"km.000000.2d.hdf5",sep="")
  ID <- amr.load(a)
  nx <- 64
  dx = L[i] / nx
  fx <- seq(dx/2,L[i]-dx/2,by=dx) 
  fy <- rep(3*L[i]/4,nx) #our domains are shifted L/2 in y , for now
  
                                        #lines(fx,fy)
  xfl <- sqrt((fx-max(fx))^2 + (fy-max(fy))^2)
  
                                        #base velocity
  ub <- transect(fx*1e3,fy*1e3,ID,1)
  vb <- transect(fx*1e3,fy*1e3,ID,2)
  umodb = sqrt(ub^2 + vb^2)
  
                                        #surface velocity
  #us <- transect(fx*1e3,fy*1e3,ID,14)
  #vs <- transect(fx*1e3,fy*1e3,ID,15)
  #umods = sqrt(us^2 + vs^2)


  #ua <- transect(fx*1e3,fy*1e3,ID,8)
  #va <- transect(fx*1e3,fy*1e3,ID,9)
  #umoda = sqrt(ua^2 + va^2)
  
  plot(xfl,umodb,type='l',xlab="x (km)",ylim=c(lower[i],upper[i]),
       main=paste("ISMIP-HOM-A, ",L[i], " km"),
       ylab=expression( abs(u) * phantom(0) * (m * phantom(0)*a^{-1})))
  #lines(xfl,umods,col="red")
  #lines(xfl,umoda,col="blue")
  amr.free(ID)
  amr.free.all()

  #legend(x="topleft",leg=c("base","surface","vert. avg."),col=c("black","red","blue"),lwd=1)
  
}

dev.off()

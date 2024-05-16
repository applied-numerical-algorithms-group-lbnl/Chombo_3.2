#require(hdf5)
require(hdf5, lib.loc="/local/R/lib")
cg <- colorRampPalette(c("blue","orange"),space="Lab")(64)

extractLevel <- function(a, l, nGhost, 
                         v = function(i){paste("v",i,sep="")}) {

  #rename list entries for convenience
  names(a[[l]]) <- c("atts","offs","d","boxes","procs")

  nc <- attr(a[[l]]$atts,"comps")
  os <- a[[l]]$offs
  dx <- attr(a[[l]],"dx")
  boxes <- list()

  
  
  for (b in 1:length(a[[l]]$boxes$hi_i)){
    
    ni <- 2*nGhost + 1 + a[[l]]$boxes$hi_i[b] - a[[l]]$boxes$lo_i[b]
    nj <- 2*nGhost + 1 + a[[l]]$boxes$hi_j[b] - a[[l]]$boxes$lo_j[b]
    
    d <- list()
  
   # d$x <- (0.5 + seq(a[[l]]$boxes$lo_i[b],  a[[l]]$boxes$hi_i[b])) * dx
   # d$y <- (0.5 + seq(a[[l]]$boxes$lo_j[b],  a[[l]]$boxes$hi_j[b])) * dx

    d$x <- (a[[l]]$boxes$lo_i[b] + 1:ni - nGhost - 0.5)*dx  / 1000
    d$y <- (a[[l]]$boxes$lo_j[b] + 1:nj - nGhost - 0.5)*dx / 1000
    
    for (i in 1:nc){
    
    s <- os[b] + nj*ni*(i-1) + 1:(nj*ni)
    d[[v(i)]] <- matrix(a[[l]]$d[s],nrow=nj,byrow=TRUE)
  }
    
    boxes[[b]] <- d
    
  }
  list(boxes=boxes, t = time, nc = nc)
}


readChombo <- function(f, nGhost){

  a <- hdf5load(f,load=FALSE)
  
  ln <- sort(names(a)[-length(names(a))])
  levels <- list()
  
  for (l in ln){
    levels[[l]] <- extractLevel(a,l, nGhost)
  }

  levels
  
}
xslice <- function( levels ,ys, v, maxlevel=length(levels)){
  
  xx <- numeric()
  zz <- numeric()

  for (level in levels[1:maxlevel]){
    for (box in level$boxes){
      y <- box$y
      sq <- abs(y - ys) < abs(y[2] - y[1])
      if (any(sq)){
        x <- box$x
        z <- box[[v]][sq,]
        if (is.matrix(z)){
          z <- apply(z,2,mean)
        }
        zz <- c(zz[xx < min(x) | xx > max(x)],z)
        xx <- c(xx[xx < min(x) | xx > max(x)],x)
      }
     
    }
  }

  o <- order(xx)
  list(x=xx[o],z=zz[o])
  
  
}

xsect2D <- function(f , y, vv, col=rep("black",length(vv)), scale = rep(1,length(vv))){
  
     a <- readChombo(f,nGhost=1)
     for (i in 1:length(vv)){
       s <- xslice(a, y, vv[i])
       lines(s$x,scale[i]*s$z,col=col[i])
     }
  
 }

mlplot <- function(levels, v, maxlevel = length(levels),
                   lv=seq(0,1,l=3), r = NULL, f = function(x){x}){

  box0 <- levels[[1]]$boxes[[1]]
  hx <- abs(box0$x[1]-box0$x[2])/2
  hy <- abs(box0$y[1]-box0$y[2])/2
  
  xr <- NULL
  yr <- NULL
  for (box in levels[[1]]$boxes) {
    xr <- range(xr,box$x)
    yr <- range(yr,box$x)
    
  }

  xr <- xr + c(-hx,hx)
  yr <- yr + c(-hy,hy)

  if (is.null(r)){
    r <- range(f(box0[[v]]))
  }

 

  plot(xr, yr , type = 'n',
       xlab = "x (km)", ylab = "y (km)",axes=TRUE)

  #axis(1,at=pretty(xr,n=2))
  #axis(2,at=pretty(yr,n=3))
  #axis(1,at=seq(xr[1],xr[2],l=5))
  #axis(2,at=seq(yr[1],yr[2],l=5))
  for (level in levels[1:maxlevel]){
    for (box in level$boxes){

      x <- box$x
      y <- box$y

      h <- abs(x[2]-x[1])/2
     

      polygon(c(min(x)-h,min(x)-h,max(x)+h,max(x)+h),
              c(min(y)-h,max(y)+h,max(y)+h,min(y)-h),
              border="black",lty=1,lwd=2)
      z <- t(f(box[[v]]))
     
      image(x,y,z, col=cg, zlim=r, add=TRUE)
      polygon( c(min(x)-h,min(x)-h,min(x)+h,min(x)+h),
               c(min(y)-h,min(y)+h,min(y)+h,min(y)-h),
              border="#111111",lty=1)
      box()

     
      
     
    }
  }

}



yslice <- function(levels, xs, v, maxlevel=length(levels)){

  yy <- numeric()
  zz <- numeric()
  
   for (level in levels[1:maxlevel]){
    for (box in level$boxes){
      x <- box$x
      sqx <-  abs(x - xs) < abs(x[2] - x[1])

      if (any(sqx)){
        y <- box$y
        z <- box[[v]][,sqx]
        if (is.matrix(z)){
          z <- apply(z,1,mean)
        }

        zz <- c(zz[yy < min(y) | yy > max(y)],z) 
        yy <- c(yy[yy < min(y) | yy > max(y)],y) 
      }
     
      
    }
  }

  o <- order(yy)
  list(y=yy[o],z=zz[o])
  
}


compare <- function(f,g,v){
#normd of the difference in v between files f and g using
#ChomboCompare
  cat (paste("compare.exactRoot =", f,"\n"),
       paste("compare.computedRoot =", g,"\n"),
       paste("compare.error_var = ",v,"\n"),
       "compare.doPlots = 0\n",
       "compare.HOaverage = 1\n",
       file="inputs.compare")
  cmd="../../../../Chombo/lib/util/ChomboCompare/compare2d.Linux.64.g++.gfortran.OPTHIGH.ex"
  a <- system(paste(cmd,"inputs.compare"),intern=TRUE)
  p <- as.numeric(strsplit(a[5]," |, ")[[1]])
  list(l1=p[2],l2=p[3],max=p[4])
  
}

cseq <- function(v,f,norm="l2"){
  #l2 norms of f[2]-f[1],f[3]-f[2] etc
  if (length(f) > 1){
    a <- compare(f[2],f[1],v)
    c(a[[norm]],cseq(v,f[-1],norm))
  } else {
    NULL
  }
}



eplot <- function(h,e,ylab="y",pch=1, xlab=expression(1/h * " (" * km^{-1} * ")"),add=FALSE){
  n <- length(e)
 
  if (!add){
    plot(range(1.1/h, 0.9/h),
         #range(1.5*e,0.5*e) ,
         c(min(e)/2,2*max(e)),
         log='xy', type='n',
         xlab=xlab,
         ylab=ylab)
    
    lines(1/h,h * e[1]/h[1],col="red",lwd=2)
    lines(1/h,h^2 *e[1]/h[1]^2 ,col="blue",lwd=2)
    lines(1/h,h^2 *e[n]/h[n]^2 ,col="blue",lwd=2)
    
    
    legend(x="bottomleft",leg=expression(C*h,C*h^2),
           lwd=1,col=c("red","blue"))
  }

  points(1/h, e, pch=pch, type = 'p', cex = 2)
  
}

cplot <- function(h,x,ylab="y",pch=1,add=FALSE){
  n <- length(x)
  e <- abs( x[2:n] - x[1:(n-1)]) + 1e-30
  eplot(h[-1],e,ylab="y",pch=1,add=FALSE)
  
}




plots <- function(lbl, model, v, amr = FALSE){



#a256 <- rf(lbl(model,"256"))
#mlplot(a256,"v2",r=c(0,500))

v <- c("v5","v7","v2","v8","v11")
s <- c(.1,.1,1,1,1)
col <- c("black","brown","blue","red","purple")
plot(c(0,160),c(-150,600),type='n',xlab="x (km)",
     ylab = "y (m)")

legend(x="topleft",leg=expression(s/10,b/10,bar(u),u[s],u[b]),
       lwd=1,col=col)

legend(x="bottomleft",leg=model,lwd=1,lty=1:length(model))
par(lty=1,lwd=2)
xsect2D(lbl(model,"128"),40,v,col,s)

if (length(model) > 1){
  for (i in 2:length(model)){
    v <- c("v2","v8","v11")
    col<- c("blue","red","purple")
    s <- c(1,1,1)
    par(lty=i,lwd=3)
    xsect2D(lbl(model[i],"128"),40,v,col,s)
    par(lty = 1,lwd=1)
  }
}



h <- 160 / 2^(5:10)
r <- c("016","032","064","128","256","512","1024")

#r <- r[1:5]
#h <- h[1:6]

l2vxs <- cseq("xbVel",lbl(model[1],r))
eplot(h,l2vxs,ylab=expression(abs(abs(u[b]^h - u[b]^{2*h}))[2]))


l2vxs <- cseq("xsVel",lbl(model[1],r))
eplot(h,l2vxs,ylab=expression(abs(abs(u[s]^h - u[s]^{2*h}))[2]))




}

stucklbl <- function(model,n){
  paste("plot.pattynA.",model,".",n,"00000.2d.hdf5",sep="")
}

stickylbl <- function(model,n){
  paste("plot.stickyA.",model,".",n,"00000.2d.hdf5",sep="")
}

slippylbl <- function(model,n){
  paste("plot.slippyA.",model,".",n,"00000.2d.hdf5",sep="")
}

slippynlbl <- function(model,n,z){
  paste("plot.",z,".slippyA.",model,".",n,"00000.2d.hdf5",sep="")
}

rf <- function(f){

  readChombo(f,nGhost=1)

}




par(mfrow=c(1,1))
w <- 600

printpdf <- FALSE

if (printpdf) {
 png(file="~/Desktop/ismipA%d.png",width=400,height=400,type="cairo",
     bg="white")
 # pdf(file="pattynA%d.pdf",width=6,height=9)
  #par(mfrow=c(3,2))
  
  par(mar=c(5,5,1,1),xaxs="i",yaxs="i")
} else {
  par(mfrow=c(4,3))
}
par(mar=c(5,5,1,1),xaxs="i",yaxs="i",las=1)





plots(stucklbl,c("L1L2","GlensLaw")) 
plots(stickylbl,c("L1L2","GlensLaw"))
plots(slippylbl,c("L1L2","GlensLaw"))

l2n <- cseq("xbVel",slippynlbl("L1L2","064",c("n6","n12","n24","n48","n96","n192")))
h <- 1/c(12,24,48,96,192)
eplot(h,l2n, xlab=expression(1/h * " (" * km^{-1} * ")"),
      ylab=expression(abs(abs(u[b]^h - u[b]^{2*h}))[2]))





if (printpdf) {
  dev.off()
}

v <- c("v1","v5","v6","v7","v2","v8","v11")
col <- c("orange","black","black","brown","blue","red","purple")



#pdf("~/Desktop/pa.pdf")
# par(mar=c(5,5,1,1),xaxs="i",yaxs="i",las=1, bg = "transparent")
#plot(c(0,160),c(-1000,2000),type='n')
#xsect2D(stucklbl("128"),xs,v,col)

#plot(c(0,160),c(-1000,2000),type='n')
#xsect2D(stickylbl("128"),xs,v,col)

#plot(c(0,160),c(0,500),type='n')
#par(lty=1)
#xsect2D(slippylbl("016"),xs,v,col)

#par(lty=2)
#xsect2D(slippylbl("064"),xs,v,col)

#xsect2D(slippylbl("128"),xs,v,col)
#par(lty=2)
#xsect2D(slippylbl("256"),xs,v,col)
#par(lty=3)
#xsect2D(slippylbl("512"),xs,v,col)
#dev.off()
#par(lty=5)
#xsect2D(slippylbl("2048"),xs,v,col)

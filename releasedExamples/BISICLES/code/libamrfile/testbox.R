require(libamrfile)
amrID <- amr.load("plot.amundsen.2d.hdf5")



#read a box of thickness data from level 0
#first, determine the domain corners
dom = amr.query.domaincorners(amrID,lev=0)
b0 <- amr.read.box(amrID,lev=0,lo=dom$lo,hi=dom$hi,comp="thickness")

#read a box of thickness data from level 0. some of this data
#will be copied from level 0, other data will be from finer levels
#choose
b1 <- amr.read.box(amrID,lev=1,lo=c(50,50),hi=c(120,120),comp="thickness"
                   ,interp=0)



#a = 3
#png("libamrfile_amr_read_box.png",width=128*a,height=192*a,pointsize=10,res=72)
par(mar=c(5,5,1,1))
#plot box data 
par(mfrow=c(1,1))
thkzl <- c(1,4000.0) #suitable range for thickness data
thkcol <- topo.colors(128)

#low res data
image(b0$x,b0$y,b0$v,zlim=thkzl,col=thkcol,xlab="x (m)", ylab="y (m)")
contour(b0$x,b0$y,b0$v,add=TRUE,lev=c(0,500,1000,1500,2000))

#paste the higher res data on top
image(b1$x,b1$y,b1$v,add=TRUE,zlim=thkzl,col=thkcol)
contour(b1$x,b1$y,b1$v,add=TRUE,lev=c(0,500,1000,1500,2000))
#draw a border round the high res box
dx = b1$x[2] - b1$x[1]
rect(min(b1$x)-dx/2, min(b1$y)-dx/2,max(b1$x)+dx/2, max(b1$y)+ dx/2,border="pink")

#find the thickness gradient
#d <- amr.readgrad.box(amrID,lev=1,lo=c(50,50),hi=c(120,120),comp=THKID)

#image(d$xf,d$y,d$dvdx,col=heat.colors(25),zlim=c(-.1,.1))
#image(d$x,d$yf,d$dvdy,col=heat.colors(25),zlim=c(-.1,.1))

time=amr.query.time(amrID)
text(min(b0$x),min(b0$y),paste("time = ", time))

#free up memory storing the the amr data
amr.free(amrID)
amr.free.all()




#dev.off()

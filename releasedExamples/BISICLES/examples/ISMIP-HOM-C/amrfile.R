dyn.load("libamrfile2d.Linux.64.g++.gfortran.DEBUG.OPT.so")

amr.load <- function(f)
  {

    #test the file is present first...
    if (is.na(file.info(f)$size)){
      stop(c("file ",f," not found"))
      rc <- -1
    } else {
      r <- .C("amr_read_file_R",
              status=integer(1),
              amrID=integer(1),
              file=as.character(f))
      
      if (r$status == 0)
        {
          rc <- r$amrID
        }
      else
        {
          stop(c("file ",f," not opened"))
          rc <- -1
        }
    }
    rc
  }

amr.free <- function(amrID)
  {
    r <- .C("amr_free",
            status=integer(1),
            amrID=as.integer(amrID))
  }

amr.free.all <- function(amrID)
  {
    r <- .C("amr_free_all")
  }
amr.query.nlevel <- function(amrID)
  {
     r <- .C("amr_query_n_level",
            status=integer(1),
            nlevel=integer(1),
            amrID=as.integer(amrID))

     if (r$status == 0){
       r$nlevel
     } else {
       -1
     }
  }

amr.query.nfab <- function(amrID, level)
  {
    
    r <- .C("amr_query_n_fab",
            status=integer(1),
            nfab=integer(1),
            amrID=as.integer(amrID),
            level=as.integer(level))

     if (r$status == 0){
       r$nfab
     } else {
       -1
     }
    
}

amr.read.fab <- function(amrID, level, fab, comp, ng=0)
  {
    r <- .C("amr_query_fab_dimensions_2d",
           status=integer(1),nx=integer(1),ny=integer(1),
           ncomp=integer(1),amrID=as.integer(amrID),
           level=as.integer(level),fab=as.integer(fab))

    if (r$status == 0) {
  
      s <- .C("amr_read_fab_data_2d",
              status=integer(1),
              v=matrix(0,r$nx+2*ng,r$ny+2*ng),
              x=numeric(r$nx+2*ng),
              y=numeric(r$ny+2*ng),
              amrID=as.integer(amrID),
              level=as.integer(level),
              fab=as.integer(fab),
              comp=as.integer(comp),
              nghost=as.integer(ng))
      if (s$status == 0){
        s
      } else {
        -2
      }
    } else {
      -1
    }
  }

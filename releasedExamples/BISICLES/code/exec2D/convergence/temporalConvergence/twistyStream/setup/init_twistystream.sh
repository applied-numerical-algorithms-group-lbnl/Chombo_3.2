#!/bin/sh
#RUNPREFIX="mpirun -np 16"
RUNPREFIX="mpirun -np 8"
#EXECFILE=../../../driver2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex
EXECFILE=../../../driver2d.Linux.64.mpiCC.mpif90.OPT.MPI.PETSC.ex
COMPAREEXEC=compare2d
INFILE_TEMPLATE=inputs.iceStream.template
INFILE_BASE=inputs.iceStream
COMPARE_TEMPLATE=inputs.compare.template
SCRIPTDIR=../../../scripts
TEMPFILE=temp.out

CRE=L1L2
NLAYER=16
NSMOOTH=32
NPLUS=16
TAGFACTOR=0.25

#function to set dt
getdt()
{
    case $RES in
        0032) BASEDT="0.5";;
        0064) BASEDT="0.25";;
        0128) BASEDT="0.125";;
        0256) BASEDT="0.0625";;
        0512) BASEDT="0.03125";;
        1024) BASEDT="0.015625";;
        2048) BASEDT="0.0078125";;
    esac
}

getdthalf()
{
    case $BASEDT in
        0.5) HALFDT="0.25";;
        0.25) HALFDT="0.125";;
        0.125) HALFDT="0.0625";;
        0.0625) HALFDT="0.03125";;
        0.03125) HALFDT="0.015625";;
        0.015625) HALFDT="0.0078125";;
        0.0078125) HALFDT="0.00390625";;
    esac
}

getdtquarter()
{
    case $BASEDT in
        0.5) QUARTERDT="0.125";;
        0.25) QUARTERDT="0.0625";;
        0.125) QUARTERDT="0.03125";;
        0.0625) QUARTERDT="0.015625";;
        0.03125) QUARTERDT="0.0078125";;
        0.015625) QUARTERDT="0.00390625";;
        0.0078125) QUARTERDT="0.001953125";;
    esac
}

getdteighth()
{
    case $BASEDT in
        0.5) EIGTHDT="0.0625";;
        0.25) EIGTHDT="0.03125";;
        0.125) EIGTHDT="0.015625";;
        0.0625) EIGTHDT="0.0078125";;
        0.03125) EIGTHDT="0.00390625";;
        0.015625) EIGTHDT="0.001953125";;
        0.0078125) EIGTHDT="0.0009765625";;
    esac
}

#function to set tagging values
gettagval()
{
   case $RES in
     0032) TAGVAL="4.0";;
     0064) TAGVAL="1.0";; 
     0128) TAGVAL="0.25";; 
     0256) TAGVAL="0.0625";; 
     0512) TAGVAL="0.015625";; 
     1024) TAGVAL="0.00390625";; 
     2048) TAGVAL="0.0009765625";; 
     *) echo "Unanticipated RES val for Tagging values";;
   esac
}

#function to set tagging values
gettagval3lev()
{
   case $RES in
     0032) TAGVAL="1.0";; 
     0064) TAGVAL="0.25";; 
     0128) TAGVAL="0.0625";; 
     0256) TAGVAL="0.015625";; 
     0512) TAGVAL="0.00390625";; 
     *) echo "Unanticipated RES val for Tagging values";;
   esac
}

TAGSGROW=1
BLOCKFACTOR=4
NREF=2
MAXLEVEL=0
CRSERES=0032
FINESTRES=2048
COMPAREFILE=../doCompare.single
RICHCOMPAREFILE=../doRichardsonCompare.single
RUNFILE=../doRuns.single
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $RICHCOMPAREFILE ]; then
  echo "deleting  $RICHCOMPAREFILE"
  rm  $RICHCOMPAREFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
EXACTDIR=
echo "generating single-level inputs"
for RES in 0032 0064 0128 0256 0512 1024 2048  
do
    RESRUNFILE=../doRuns-$RES.single
    if [ -f $RESRUNFILE ]; then
        echo "deleting $RESRUNFILE"
        rm $RESRUNFILE
    fi
    
    YRES=$RES
    gettagval
    getdt

    echo $RES "-- dt = "  $BASEDT;
    infile=$INFILE_BASE.$CRE.$RES.dt$BASEDT.l$MAXLEVEL
    of=../$infile
    sed  -e s/@DT/$BASEDT/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    
    richcomparein=inputs.compare.Richardson.$CRSERES.l$MAXLEVEL
    rcof=../$richcomparein
#    sed -e s/@CRE/$CRE/ -e s/@FINERES/$RES/ -e s/@CRSERES/$CRSERES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@EXACTDIR/$EXACTDIR/  $COMPARE_TEMPLATE > $rcof
    comparein=inputs.compare.$CRSERES.l$MAXLEVEL
    cof=../$comparein
#    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINESTRES/ -e s/@CRSERES/$CRSERES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@EXACTDIR/$EXACTDIR/  $COMPARE_TEMPLATE > $cof
 
    outfile="run.$CRE.$RES.dt$BASEDT.l$MAXLEVEL"
    innerConvergename="../solverConverge/stream.$CRE.$RES.l$MAXLEVEL"
    outerConvergename="../solverConverge/stream.$CRE.$RES.l$MAXLEVEL.outer"
    runcommand="$RUNPREFIX $EXECFILE $infile > $outfile"
    echo "echo \"doing $RES run with dt = $BASEDT\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
#    if [ $RES!="0032" ]; then
#      richcomparecommand="$COMPAREEXEC $richcomparein"
#      echo $richcomparecommand >> $RICHCOMPAREFILE
#      comparecommand="$COMPAREEXEC $comparein"
#      echo $comparecommand >> $COMPAREFILE
#    fi

    #set up temporal convergence
    # base dt
    echo "echo \"doing $RES run with dt = $BASEDT\" " >> $RESRUNFILE
    echo $runcommand >> $RESRUNFILE    

    # dt/2
    getdthalf
    echo $RES "-- dt = "  $HALFDT;
    infile=$INFILE_BASE.$CRE.$RES.dt$HALFDT.l$MAXLEVEL
    of=../$infile
    sed  -e s/@DT/$HALFDT/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of    

    outfile="run.$CRE.$RES.dt$HALFDT.l$MAXLEVEL"    
    runcommand="$RUNPREFIX $EXECFILE $infile > $outfile"
    echo "echo \"doing $RES run with dt = $HALFDT\" " >> $RESRUNFILE
    echo $runcommand >> $RESRUNFILE

    
    # dt/4
    getdtquarter
    echo $RES "-- dt = "  $QUARTERDT;
    infile=$INFILE_BASE.$CRE.$RES.dt$QUARTERDT.l$MAXLEVEL
    of=../$infile
    sed  -e s/@DT/$QUARTERDT/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of    

    outfile="run.$CRE.$RES.dt$QUARTERDT.l$MAXLEVEL"    
    runcommand="$RUNPREFIX $EXECFILE $infile > $outfile"
    echo "echo \"doing $RES run with dt = $QUARTERDT\" " >> $RESRUNFILE
    echo $runcommand >> $RESRUNFILE
    

    #dt/8
    getdteighth
    echo $RES "-- dt = "  $EIGTHDT;
    infile=$INFILE_BASE.$CRE.$RES.dt$EIGTHDT.l$MAXLEVEL
    of=../$infile
    sed  -e s/@DT/$EIGTHDT/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of    

    outfile="run.$CRE.$RES.dt$EIGTHDT.l$MAXLEVEL"    
    runcommand="$RUNPREFIX $EXECFILE $infile > $outfile"
    echo "echo \"doing $RES run with dt = $EIGTHDT\" " >> $RESRUNFILE
    echo $runcommand >> $RESRUNFILE
        
#    chmod +x $COMPAREFILE
#    chmod +x $RICHCOMPAREFILE
    chmod +x $RUNFILE
    chmod +x $RESRUNFILE

    CRSERES=$RES
done 

exit 0



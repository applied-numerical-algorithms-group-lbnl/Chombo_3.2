#!/bin/sh
#RUNPREFIX="mpirun -np 8"
RUNPREFIX=""
#EXECFILE=../../driver2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex
EXECFILE=../../driver2d.Linux.64.g++.gfortran.OPT.ex 
COMPAREEXEC=compare2d
INFILE_TEMPLATE=inputs.hump.template
INFILE_BASE=inputs.hump
COMPARE_TEMPLATE=inputs.compare.template
CRE=L1L2
NLAYER=16
NSMOOTH=8
NPLUS=16
TAGFACTOR=0.25

#function to set tagging values
gettagval()
{
   case $RES in
     0032) TAGVAL="10";;
     0064) TAGVAL="2.5";; 
     0128) TAGVAL="0.625";; 
     0256) TAGVAL="0.15625";; 
     0512) TAGVAL="0.0390625";; 
     1024) TAGVAL="0.009765625";; 
     2048) TAGVAL="0.00244140625";; 
     4096) TAGVAL="0.0006103515";; 
     8192) TAGVAL="0.0001525878";; 
     *) echo "Unanticipated RES val for Tagging values";;
   esac
}

#function to set tagging values
gettagval3lev()
{
   case $RES in
     0032) TAGVAL="2.5";; 
     0064) TAGVAL="0.625";; 
     0128) TAGVAL="0.15625";; 
     0256) TAGVAL="0.0390625";; 
     0512) TAGVAL="0.009765625";; 
     1024) TAGVAL="0.00244140625";; 
     2048) TAGVAL="0.0006103515";; 
     4096) TAGVAL="0.0001525878";; 
     *) echo "Unanticipated RES val for Tagging values";;
   esac
}

NREF=2
MAXLEVEL=0
CRSERES="0032"
FINESTRES=2048
COMPAREFILE=doCompare.single
RICHCOMPAREFILE=doRichardsonCompare.single
RUNFILE=doRuns.single
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
echo "generating single-level inputs"
for RES in 0032 0064 0128 0256 0512 1024 2048 4096 8192
do
    YRES=$RES
    gettagval 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@NREF1/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    rcof=inputs.compare.Richardson.$CRSERES.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$RES/ -e s/@CRSERES/$CRSERES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $rcof
    cof=inputs.compare.$CRSERES.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINESTRES/ -e s/@CRSERES/$CRSERES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    runcommand="$RUNPREFIX $EXECFILE $of > run.$CRE.$RES.l$MAXLEVEL"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "cp pout.0 pout.$RES" >> $RUNFILE
    if [ $RES!="0032" ]; then
      richcomparecommand="$COMPAREEXEC $rcof"
      echo $richcomparecommand >> $RICHCOMPAREFILE
      comparecommand="$COMPAREEXEC $cof"
      echo $comparecommand >> $COMPAREFILE
    fi
chmod +x $COMPAREFILE
chmod +x $RICHCOMPAREFILE
chmod +x $RUNFILE


CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 

COMPAREFILE=doCompare.2Ref
RUNFILE=doRuns.2Ref
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
FINERES=2048
NREF=2
NSMOOTH=32
MAXLEVEL=1
echo "generating nRef = 2 inputs"
for RES in 0032 0064 0128 0256 0512 1024 2048
do
    YRES=$RES
    gettagval 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@NREF1/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    cof=inputs.compare.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINERES/ -e s/@CRSERES/$RES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    runcommand="$RUNPREFIX $EXECFILE $of > run.$CRE.$RES.r$NREF.l$MAXLEVEL"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "cp pout.0 pout.2Ref.$RES" >> $RUNFILE

    comparecommand="$COMPAREEXEC $cof"
    echo $comparecommand >> $COMPAREFILE
CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 
chmod +x $COMPAREFILE
chmod +x $RUNFILE


COMPAREFILE=doCompare.3Lev
RUNFILE=doRuns.3Lev
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
FINERES=2048
NREF=2
NSMOOTH=32
MAXLEVEL=2
echo "generating three-level inputs" 
for RES in 0032 0064 0128 0256 0512 1024
do
    YRES=$RES
    gettagval3lev 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@NREF1/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    cof=inputs.compare.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINERES/ -e s/@CRSERES/$RES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    runcommand="$RUNPREFIX $EXECFILE $of > run.$CRE.$RES.l$MAXLEVEL"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "cp pout.0 pout.3Lev.$RES" >> $RUNFILE    

    comparecommand="$COMPAREEXEC $cof"
    echo $comparecommand >> $COMPAREFILE

CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 
chmod +x $COMPAREFILE
chmod +x $RUNFILE

COMPAREFILE=doCompare.4Ref
RUNFILE=doRuns.4Ref
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
if [ -f $COMPAREFILE ]; then
  echo "deleting  $COMPAREFILE"
  rm  $COMPAREFILE
fi
FINERES=2048
NREF=4
NSMOOTH=32
MAXLEVEL=1
echo "generating nRef = 4 inputs" 
for RES in 0032 0064 0128 0256 0512 1024
do
    YRES=$RES
    gettagval 
    echo $RES "-- tagval = "  $TAGVAL;
    of=$INFILE_BASE.$CRE.$RES.$NREF.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@NREF/$NREF/ -e s/@NREF1/$NREF/ -e s/@TAGVAL/$TAGVAL/ -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $of
    cof=inputs.compare.$RES.r$NREF.l$MAXLEVEL
    sed -e s/@CRE/$CRE/ -e s/@FINERES/$FINERES/ -e s/@CRSERES/$RES/ -e s/@NREF/$NREF/ -e s/@MAXLEVEL/$MAXLEVEL/  $COMPARE_TEMPLATE > $cof
 
    runcommand="$RUNPREFIX $EXECFILE $of > run.$CRE.$RES.r$NREF.l$MAXLEVEL"
    echo "echo \"doing $RES run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
    echo "cp pout.0 pout.4Ref.$RES" >> $RUNFILE

    comparecommand="$COMPAREEXEC $cof"
    echo $comparecommand >> $COMPAREFILE

CRSERES=$RES
NSMOOTH=$((NSMOOTH + NPLUS))
done 
chmod +x $COMPAREFILE
chmod +x $RUNFILE

exit 0



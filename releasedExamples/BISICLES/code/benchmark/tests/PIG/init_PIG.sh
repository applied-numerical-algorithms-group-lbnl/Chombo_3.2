#!/bin/sh
EXECFILE=../../solverBenchmark2d.Linux.64.g++.gfortran.OPT.ex
SCRIPTDIR=../../../exec2D/scripts
DOCDIR=./doc
INFILE_TEMPLATE=inputs.pig.l1l2.template
INFILE_BASE=inputs.pig
OUTFILE_BASE=run.pig
TEMPFILE=temp.out
RESULTDIR=./doc
VTOPSAFETY=0.9

NLAYER=11
NSMOOTH=4

#functions to convert from names to numbers
getprolongtype()
{
   case $PROLONGTYPE in
     "pc") PROLONG="0";;
     "linear") PROLONG="1";; 
     *) echo "Unanticipated PROLONGTYPE val $PROLONGTYPE for linear solver type";;
   esac
}

getsolvertype()
{
   case $LINSOLVERNAME in
     BICG) LINEARSOLVER="1";;
     MG) LINEARSOLVER="0";; 
     *) echo "Unanticipated LINSOLVERNAME val for linear solver type";;
   esac
}



NREF=2
MAXLEVEL=1
#just do jfnk for now
NONLINSOLVER=1
RUNFILE=doRuns.2lev
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE"
  rm $RUNFILE
fi
echo "generating two-level inputs"

for LINSOLVERTOL in 3 4 5 6 7 8 9 10  
do

for LINSOLVERNAME in BICG MG
do

for PROLONGTYPE in pc linear
do

for NUMMG in 1 2 3 4
do

echo "solver = $LINSOLVERNAME, prolong = $PROLONGTYPE, numMG = $NUMMG, linear solver tol = 1.0e-$LINSOLVERTOL"
    getsolvertype
    getprolongtype

#   echo "PROLONG = $PROLONG, LINEARSOLVER = $LINEARSOLVER"
    infile=$INFILE_BASE.l$MAXLEVEL.$LINSOLVERNAME.$PROLONGTYPE.mg$NUMMG.tol$LINSOLVERTOL;
    outfile=$OUTFILE_BASE.l$MAXLEVEL.$LINSOLVERNAME.$PROLONGTYPE.mg$NUMMG.tol$LINSOLVERTOL;
    innerConvergename=$RESULTDIR/solverConverge.l$MAXLEVEL.$LINSOLVERNAME.$PROLONGTYPE.mg$NUMMG.tol$LINSOLVERTOL;
    outerConvergename=$RESULTDIR/solverConverge.l$MAXLEVEL.$LINSOLVERNAME.$PROLONGTYPE.mg$NUMMG.tol$LINSOLVERTOL.outer; 
    sed -e s/@LINEARSOLVER/$LINEARSOLVER/ -e s/@LINSOLVERTOL/$LINSOLVERTOL/ -e s/@NSMOOTH/$NSMOOTH/ -e s/@NONLINSOLVER/$NONLINSOLVER/ -e s/@NLAYER/$NLAYER/ -e s/@NUMMG/$NUMMG/ -e s/@PROLONG/$PROLONG/ -e s/@NREF/$NREF/  -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@VTOPSAFETY/$VTOPSAFETY/  $INFILE_TEMPLATE > $infile
#    echo $infile;
#    echo $outfile;
    runcommand="$EXECFILE $infile > $outfile"
    echo "echo \"doing $infile run\" " >> $RUNFILE
    echo $runcommand >> $RUNFILE
#run post-processing scripts
    echo "$SCRIPTDIR/innerJFNK.awk < $outfile > $TEMPFILE " >> $RUNFILE
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename" >> $RUNFILE
for NUMSCALE in 2 3 4
do
    innerConvergeScaleName=$innerConvergename.scale$NUMSCALE
    echo "$SCRIPTDIR/a.out $TEMPFILE  $" >> $RUNFILE  $innerConvergeScaleName
done
for NUMSCALE in 2 3 4
do
    innerConvergeScaleName=$innerConvergename.scale$NUMSCALE
    echo "$SCRIPTDIR/a.out $TEMPFILE $innerConvergeScaleName $NUMSCALE" >> $RUNFILE 
done
    echo "$SCRIPTDIR/jfnk.awk < $outfile > $outerConvergename " >>$RUNFILE
#save time.table file if there is one
    echo "mv time.table $RESULTDIR/time.l$MAXLEVEL.$LINSOLVERNAME.$PROLONGTYPE.mg$NUMMG.tol$LINSOLVERTOL.table" >> $RUNFILE
done 

done 

done 

done

chmod +x $RUNFILE

exit 0



#!/bin/sh
EXECFILE1=../../../../driver2d.Linux.64.mpiCC.mpif90.OPT.MPI.PETSC.ex 
EXECFILE2=../../../../driver2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex 
COMPAREEXEC=compare2d
INFILE1_TEMPLATE=inputs.petscBottom.template
INFILE1_BASE=inputs.petsc
INFILE2_TEMPLATE=inputs.BiCGBottom.template
INFILE2_BASE=inputs.BiCGBottom
SCRIPTDIR=../../../../scripts
TEMPFILE=temp.out
MAXDEPTH=-1
NPROC=8

CRSERES="0032"
FINESTRES=2048
RUNFILE1=doRuns.petscBottom
RUNFILE2=doRuns.BiCGBottom
RUNFILE3=doRuns.petscBottom.MGdepth
RUNFILE4=doRuns.BiCGBottom.MGdepth
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE1, $RUNFILE2, $RUNFILE3, $RUNFILE4"
  rm $RUNFILE1
  rm $RUNFILE2
  rm $RUNFILE3
  rm $RUNFILE4
fi
echo "generating single-level inputs"
for RES in 0032 0064 0128 0256 0512 1024 2048  
do
    YRES=$RES
    echo $RES 
    of1=$INFILE1_BASE.$RES
    of2=$INFILE2_BASE.$RES
    sed  -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXDEPTH/$MAXDEPTH/ $INFILE1_TEMPLATE > $of1
    sed  -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXDEPTH/$MAXDEPTH/ $INFILE2_TEMPLATE > $of2

    outfile1="run.petscBottom.$RES"
    innerConvergename1="solverConverge/resid.petscBottom.$RES"
    outerConvergename1="solverConverge/resid.petscBottom.$RES.outer"
    poutname1="pout.petscBottom.$RES"
    runcommand1="mpirun -np $NPROC $EXECFILE1 $of1 > $outfile1"
    echo "echo \"doing $RES run\" " >> $RUNFILE1
    echo $runcommand1 >> $RUNFILE1
    echo "mv pout.0 $poutname1" >> $RUNFILE1
    echo "$SCRIPTDIR/innerPetsc.awk < $outfile1 > $TEMPFILE " >> $RUNFILE1
    echo "$SCRIPTDIR/parseMG $TEMPFILE  $innerConvergename1" >> $RUNFILE1
    echo "$SCRIPTDIR/petsc.awk < $outfile1 > $TEMPFILE" >> $RUNFILE1
    echo "$SCRIPTDIR/parseJFNK $TEMPFILE  $outerConvergename1" >> $RUNFILE1

    outfile2="run.BiCGBottom.$RES"
    innerConvergename2="solverConverge/resid.BiCGBottom.$RES"
    outerConvergename2="solverConverge/resid.BiCGBottom.$RES.outer"
    runcommand2="mpirun -np $NPROC $EXECFILE2 $of2 > $outfile2"
    echo "echo \"doing $RES run\" " >> $RUNFILE2
    poutname2="pout.BiCGBottom.$RES"
    echo $runcommand2 >> $RUNFILE2
    echo "mv pout.0 $poutname2 " >> $RUNFILE2
    echo "$SCRIPTDIR/innerJFNK.awk < $poutname2 > $TEMPFILE " >> $RUNFILE2
    echo "$SCRIPTDIR/parseMG $TEMPFILE  $innerConvergename2" >> $RUNFILE2
    echo "$SCRIPTDIR/jfnk.awk < $poutname2 > $TEMPFILE" >> $RUNFILE2
    echo "$SCRIPTDIR/parseJFNK $TEMPFILE  $outerConvergename2" >> $RUNFILE2

chmod +x $RUNFILE1
chmod +x $RUNFILE2


CRSERES=$RES
done 

echo "generating MG depth test inputs"
for MAXDEPTH in 0 1 2 3 4 
do
    YRES=$RES
    echo "$RES, depth = $MAXDEPTH"
    of1=$INFILE1_BASE.$RES.depth$MAXDEPTH
    of2=$INFILE2_BASE.$RES.depth$MAXDEPTH
    sed  -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXDEPTH/$MAXDEPTH/  $INFILE1_TEMPLATE > $of1
    sed  -e s/@RES/$RES/ -e s/@YRES/$RES/ -e s/@MAXDEPTH/$MAXDEPTH/  $INFILE2_TEMPLATE > $of2

    outfile1="run.petscBottom.$RES.depth$MAXDEPTH"
    innerConvergename1="solverConverge/resid.petscBottom.$RES.depth$MAXDEPTH"
    outerConvergename1="solverConverge/resid.petscBottom.$RES.depth$MAXDEPTH.outer"
    poutname1="pout.petscBottom.$RES.$MAXDEPTH"
    runcommand1="mpirun -np $NPROC $EXECFILE1 $of1 > $outfile1"
    echo "echo \"doing $RES, $MAXDEPTH depth run\" " >> $RUNFILE3
    echo $runcommand1 >> $RUNFILE3
    echo "mv pout.0 $poutname1" >> $RUNFILE3
    echo "$SCRIPTDIR/innerPetsc.awk < $outfile1 > $TEMPFILE " >> $RUNFILE3
    echo "$SCRIPTDIR/parseMG $TEMPFILE  $innerConvergename1" >> $RUNFILE3
    echo "$SCRIPTDIR/petsc.awk < $outfile1 > $TEMPFILE" >> $RUNFILE3
    echo "$SCRIPTDIR/parseJFNK $TEMPFILE  $outerConvergename1" >> $RUNFILE3

    outfile2="run.BiCGBottom.$RES.depth$MAXDEPTH"
    innerConvergename2="solverConverge/resid.BiCGBottom.$RES.depth$MAXDEPTH"
    outerConvergename2="solverConverge/resid.BiCGBottom.$RES.depth$MAXDEPTH.outer"
    runcommand2="mpirun -np $NPROC $EXECFILE2 $of2 > $outfile2"
    echo "echo \"doing $RES run\" " >> $RUNFILE4
    poutname2="pout.BiCGBottom.$RES.$MAXDEPTH"
    echo $runcommand2 >> $RUNFILE4
    echo "mv pout.0 $poutname2 " >> $RUNFILE4
    echo "$SCRIPTDIR/innerJFNK.awk < $poutname2 > $TEMPFILE " >> $RUNFILE4
    echo "$SCRIPTDIR/parseMG $TEMPFILE  $innerConvergename2" >> $RUNFILE4
    echo "$SCRIPTDIR/jfnk.awk < $poutname2 > $TEMPFILE" >> $RUNFILE4
    echo "$SCRIPTDIR/parseJFNK $TEMPFILE  $outerConvergename2" >> $RUNFILE4

chmod +x $RUNFILE3
chmod +x $RUNFILE4


CRSERES=$RES
done 

exit 0



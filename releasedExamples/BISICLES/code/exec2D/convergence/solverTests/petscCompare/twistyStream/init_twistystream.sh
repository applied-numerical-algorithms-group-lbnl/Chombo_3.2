#!/bin/sh
EXECFILE1=../../../../driver2d.Linux.64.mpiCC.mpif90.OPT.MPI.PETSC.ex 
EXECFILE2=../../../../driver2d.Linux.64.mpiCC.gfortran.OPT.MPI.ex 
COMPAREEXEC=compare2d
INFILE1_TEMPLATE=inputs.petsc.template
INFILE1_BASE=inputs.petsc
INFILE2_TEMPLATE=inputs.petscCompare.template
INFILE2_BASE=inputs.petscCompare
SCRIPTDIR=../../../../scripts
TEMPFILE=temp.out
NPROC=1

CRSERES="0032"
FINESTRES=2048
RUNFILE1=doRuns.petsc
RUNFILE2=doRuns.MG-JFNK
if [ -f $RUNFILE ]; then
  echo "deleting $RUNFILE1, $RUNFILE2"
  rm $RUNFILE1
  rm $RUNFILE2
fi
echo "generating single-level inputs"
for RES in 0032 0064 0128 0256 0512 1024 2048  
do
    YRES=$RES
    echo $RES 
    of1=$INFILE1_BASE.$RES
    of2=$INFILE2_BASE.$RES
    sed  -e s/@RES/$RES/ -e s/@YRES/$RES/  $INFILE1_TEMPLATE > $of1
    sed  -e s/@RES/$RES/ -e s/@YRES/$RES/  $INFILE2_TEMPLATE > $of2

    outfile1="run.petsc-gamg.$RES"
    innerConvergename1="solverConverge/resid.petsc-gamg.$RES"
    outerConvergename1="solverConverge/resid.petsc-gamg.$RES.outer"
    poutname1="pout.petsc-gamg.$RES"
    runcommand1="mpirun -np $NPROC $EXECFILE1 $of1 > $outfile1"
    echo "echo \"doing $RES run\" " >> $RUNFILE1
    echo $runcommand1 >> $RUNFILE1
    echo "mv pout.0 $poutname1" >> $RUNFILE1
    echo "$SCRIPTDIR/innerPetsc.awk < $outfile1 > $TEMPFILE " >> $RUNFILE1
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename1" >> $RUNFILE1
    echo "$SCRIPTDIR/petsc.awk < $outfile1 > $outerConvergename1" >> $RUNFILE1

    outfile2="run.MG-JFNK.$RES"
    innerConvergename2="solverConverge/resid.MG-JFNK.$RES"
    outerConvergename2="solverConverge/resid.MG-JFNK.$RES.outer"
    runcommand2="mpirun -np $NPROC $EXECFILE2 $of2 > $outfile2"
    echo "echo \"doing $RES run\" " >> $RUNFILE2
    poutname2="pout.MG-JFNK.$RES"
    echo $runcommand2 >> $RUNFILE2
    echo "mv pout.0 $poutname2 " >> $RUNFILE2
    echo "$SCRIPTDIR/innerJFNK.awk < $poutname2 > $TEMPFILE " >> $RUNFILE2
    echo "$SCRIPTDIR/a.out $TEMPFILE  $innerConvergename2" >> $RUNFILE2
    echo "$SCRIPTDIR/jfnk.awk < $poutname2 > $outerConvergename2" >> $RUNFILE2

chmod +x $RUNFILE1
chmod +x $RUNFILE2


CRSERES=$RES
done 

exit 0



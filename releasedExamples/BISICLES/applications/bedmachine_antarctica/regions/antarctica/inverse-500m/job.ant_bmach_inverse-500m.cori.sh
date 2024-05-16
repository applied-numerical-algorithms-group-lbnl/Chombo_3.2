#!/bin/bash 
#
# job script for nersc cori
#
#SBATCH -A m1041
#SBATCH -J 500m_hypre_bdmach_inverse
##SBATCH --qos=debug
#SBATCH --qos=regular
#SBATCH --time=10:00:00
#SBATCH --nodes=12
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
#SBATCH --mail-type=ALL
#SBATCH --mail-user=DFMartin@lbl.gov


module load python


export PETSC_DIR=/global/common/software/m1041/petsc_install/petsc_haswell_gnu/
export PETSC_ARCH=""
export LD_LIBRARY_PATH=${PYTHON_DIR}/lib:${LD_LIBRARY_PATH}


NCORE=$((32 * $SLURM_JOB_NUM_NODES))
echo "n_node: $SLURM_JOB_NUM_NODES , n_core: $NCORE"
DRIVER=/global/common/software/m1041/BISICLES/haswell/BISICLES/code/exec2D/driver2d.Linux.64.CC.ftn.OPT.MPI.PETSC.ex
NAME=ant_bmach.inverse-500m
RUNDIR=$SLURM_SUBMIT_DIR/$SLURM_JOB_ID
mkdir -p $RUNDIR
INFILEBASE=inputs.$NAME
INFILE=$INFILEBASE"."$SLURM_JOB_ID
cp $SLURM_SUBMIT_DIR/$INFILEBASE $RUNDIR/$INFILE
cp $SLURM_SUBMIT_DIR/.petscrc $RUNDIR/
cd $RUNDIR

#work out what the latest checkpoint file is (if it exists)
if test -n "$(find ../ -maxdepth 1 -name 'chk.ant_bmach_inverse.??????.2d.hdf5' -print -quit)"
    then
    LCHK=`ls -th ../chk.ant_bmach_inverse.??????.2d.hdf5 | head -n 1`
    echo "" >> $INFILE #ensure line break
    echo "amr.restart_file=$LCHK" >> $INFILE
    echo "amr.restart_set_time=false" >> $INFILE
    echo "" >> $INFILE #ensure line break
fi

export CH_TIMER=1
export CH_OUTPUT_INTERVAL=999
export PYTHONPATH=$PWD:$PYTHONPATH

echo "srun $DRIVER $INFILE"
srun $DRIVER $INFILE

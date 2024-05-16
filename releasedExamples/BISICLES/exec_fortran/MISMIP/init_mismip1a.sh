CRE=L1L2
DIFFUSION=implicit
RES=00128
YRES=00032
PI=100
CI=100
NSMOOTH=8
NPLUS=8
for MAXLEVEL in 0 1 2 3 4 5 6;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of

NSMOOTH=$((NSMOOTH + NPLUS))
done 


RES=00256
YRES=00064
PI=200
CI=200
NSMOOTH=16

for MAXLEVEL in 0 1 2 3 4 5;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of
NSMOOTH=$((NSMOOTH + NPLUS))
done 



RES=00512
YRES=00128
PI=400
CI=400
NSMOOTH=24
for MAXLEVEL in 0 1 2 3 4;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of
done 


RES=01024
YRES=00256
PI=800
CI=800
NSMOOTH=32

for MAXLEVEL in 0 1 2 3;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed  -e s/@NSMOOTH/$NSMOOTH/ --e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of
NSMOOTH=$((NSMOOTH + NPLUS))
done 


RES=02048
YRES=00512
PI=1600
CI=1600
NSMOOTH=48

for MAXLEVEL in 0 1 2;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of
done 

RES=04096
YRES=01024
PI=3200
CI=3200
NSMOOTH=64

for MAXLEVEL in 0 1;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of

NSMOOTH=$((NSMOOTH + NPLUS))
done 


RES=08192
YRES=02048
PI=6400
CI=6400
NSMOOTH=72

for MAXLEVEL in 0;
do
    of=inputs.init_mismip1a.$CRE.$DIFFUSION.$RES.r2.l$MAXLEVEL
    sed -e s/@NSMOOTH/$NSMOOTH/ -e s/@CRE/$CRE/ -e s/@DIFFUSION/$DIFFUSION/ -e s/@RES/$RES/ -e s/@YRES/$YRES/ -e s/@MAXLEVEL/$MAXLEVEL/ -e s/@PI/$PI/ -e s/@CI/$CI/ inputs.init_mismip1a.template > $of

NSMOOTH=$((NSMOOTH + NPLUS))
done 





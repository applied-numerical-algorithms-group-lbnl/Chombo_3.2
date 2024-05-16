#!/bin/sh

INFILE_TEMPLATE=inputs.iceStream.template
INFILE_BASE=inputs.iceStream

CRE=L1L2
NLAYER=16
TAGSGROW=1
BLOCKFACTOR=8


echo "generating single-level inputs"
MAXLEVEL=0
TAGVAL=0.0
for XRES in 0032 0064 0128 0256 0512 1024 2048 ; do
    YRES=$XRES
    for SIAONLY in true false; do
	for SIALIMIT in 1e-2 1e-1 1e-0; do
	    NAME=$CRE.siaonly$SIAONLY.sialimit$SIALIMIT.$XRES.$MAXLEVEL"lev"
	    infile=$INFILE_BASE.$NAME
	    
	    sed -e s/@NAME/$NAME/ -e s/@SIALIMIT/$SIALIMIT/ -e s/@SIAONLY/$SIAONLY/ -e s/@CRE/$CRE/ -e s/@NLAYER/$NLAYER/ -e s/@XRES/$XRES/ -e s/@YRES/$YRES/ -e s/@TAGVAL/$TAGVAL/ -e s/@TAGSGROW/$TAGSGROW/ -e s/@BLOCKFACTOR/$BLOCKFACTOR/  -e s/@MAXLEVEL/$MAXLEVEL/  $INFILE_TEMPLATE > $infile
	done
    done
done
exit 0



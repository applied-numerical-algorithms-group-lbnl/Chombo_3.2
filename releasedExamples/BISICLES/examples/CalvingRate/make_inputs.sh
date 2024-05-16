getcre()
{
  case $SMOD in
    l1l2) CRE=L1L2;;
    ssa)  CRE=GlensLaw;;
    l1l2_full) CRE=L1L2;;
    tanh) CRE=GlensLaw;;
    *) echo "unknown  stress model"
  esac

  case $SMOD in
      tanh) NL_SOLVER=5 ;;
      *) NL_SOLVER=1
  esac


}


getcalv()
{

    case $CALVTYPENAME in
	u) CALVTYPE=RateProportionalToSpeedCalvingModel ;;
	*) CALVTYPE=VariableRateCalvingModel
    esac
    }


getbase()
{
   case $BASENW in
	8) BASENAME="crse1000m";;
	16) BASENAME="crse500m";;
        32) BASENAME="crse250m";;
   esac 

   BASENL=$(( $BASENW * 16 ))

   BCW=1
   if [ $SLIP != no ]
   then
       BCW=0
   fi
   
   if [ $DIR != x ]
   then
       BASENX=$BASENW
       BASENY=$BASENL
       BC_LO_Y=1
       BC_LO_X=$BCW
       BC_HI_Y=0
       BC_HI_X=$BCW
       FRONT_HI_Y=1
       FRONT_HI_X=0
       LX=8.0e+3
       LY=128.0e+3
   else
       BASENX=$BASENL
       BASENY=$BASENW
       BC_LO_X=1
       BC_LO_Y=$BCW
       BC_HI_Y=$BCW
       BC_HI_X=0
       FRONT_HI_Y=0
       FRONT_HI_X=1
       LY=8.0e+3
       LX=128.0e+3
   fi

   echo "BASENX="$BASENX",BASENY="$BASENY
}
getsgname()
{
    case $SGN in
	*) SGNAME="sg"$SGN;;
    esac
}

getbflawname()
{

    BFLAW=$BFLAWA
    PLLMODEL=None
    PLLCOEF=-1.0;
    case $BFLAWA in
	tsaiLaw) BFLAW=pressureLimitedLaw;BFLAWNAME=.prlim.; PLLMODEL=Tsai; PLLCOEF=0.5;;
	leguyLaw) BFLAW=pressureLimitedLaw;BFLAWNAME=.leguy.; PLLMODEL=Leguy; PLLCOEF=8.0e12;;
	*) BFLAWNAME=""
    esac
}


mkdir -p scripts

BASEDIR=$PWD
ACOEF=2.0e-17
CFUNC=constfriction
solver=Chombo

BFLAWA=tsaiLaw
TEMPLATE=inputs.cr.template
BASENW=8

for SMOD in ssa tanh
do
for SGN in 0 4
do
    for lev in 0 1 2 3 4 5 6 7
    do
	for DIR in x y
	do
	    for CALVTYPENAME in "" "u"
	    do
		for SLIP in no free
		do
		    tagcap=$(( lev - 1 )) 
		    getbase
		    getcre
		    getsgname
		    getbflawname
		    getcalv
		    NAME=cr_$CALVTYPENAME$DIR.$SLIP-slip.$SMOD$BFLAWNAME$BASENAME.$lev"lev."$SGNAME
		    echo $NAME"..."
		    INFILE=inputs.$NAME
		    sed -e s/@SLIP/$SLIP/ -e s/@CALVTYPE/$CALVTYPE/ -e s/@NL_SOLVER/$NL_SOLVER/ -e s/@SGN/$SGN/ -e s/@LX/$LX/ -e s/@LY/$LY/ -e s/@DIR/$DIR/ -e s/@BC_LO_X/$BC_LO_X/ -e s/@BC_LO_Y/$BC_LO_Y/ -e s/@BC_HI_X/$BC_HI_X/ -e s/@BC_HI_Y/$BC_HI_Y/ -e s/@FRONT_HI_X/$FRONT_HI_X/ -e s/@FRONT_HI_Y/$FRONT_HI_Y/ -e s/@PLLMODEL/$PLLMODEL/ -e s/@PLLCOEF/$PLLCOEF/ -e s/@BFLAW/$BFLAW/ -e s/@BASENX/$BASENX/ -e s/@BASENY/$BASENY/ -e s/#$SMOD// -e s/#$solver// -e s/@SOLVER/$solver/ -e s/@NAME/$NAME/ -e s/@ACOEF/$ACOEF/ -e s/@CFUNC/$CFUNC/ -e s/@ACAB/$ACAB/ -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$SMOD/ -e s/@CRE/$CRE/ $TEMPLATE > scripts/$INFILE
	    
		    echo "... done"
		done
	    done
	done
    done	
done
done

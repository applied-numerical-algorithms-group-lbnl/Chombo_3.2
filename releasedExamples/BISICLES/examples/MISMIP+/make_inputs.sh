getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    l1l2_full) cre=L1L2;;
    *) echo "unknown stress model"
  esac
}

getnodes()
{
  case $lev in
    5)  NODES=1;;
    *) NODES=1
  esac
}

getsolveinterval()
{
    case $lev in
	7) SOLVEINTERVAL=8;;
	6) SOLVEINTERVAL=8;;
	5) SOLVEINTERVAL=8;;
	4) SOLVEINTERVAL=8;;
	3) SOLVEINTERVAL=4;;
	2) SOLVEINTERVAL=2;;
	1) SOLVEINTERVAL=1;;
	*) SOLVEINTERVAL=1
    esac  

    case $MELT in 
	spin) SOLVEINTERVAL=$SOLVEINTERVAL;;
	*) SOLVEINTERVAL=1
    esac

}
getplotinterval()
{
    case $MELT in
	spin) PLOTINTERVAL=1000;;
	*) PLOTINTERVAL=10
    esac  
}

getmaxtime()
{
    case $MELT in
	spin) MAXTIME=50000;;
	*) MAXTIME=1500
    esac  
}

gettaglapname()
{
    case $TAGLAP in
	0) TAGLAPNAME="";;
	1) TAGLAPNAME=".taglap";;
    esac
}

getsgname()
{
    case $SGN in
	0) SGNAME="";;
	*) SGNAME=".sg"$SGN;;
    esac
}

gettaggrowname()
{
    TAGGROWNAME=".grow"$TAGGROW
    case $TAGGROW in
	4) TAGGROWNAME="";;
    esac
}

getbase()
{
   case $BASENY in
	20) BASENAME="";;
	40) BASENAME=".base2k";;
        80) BASENAME=".base1k";;
   esac 

   BASENX=$BASENY"0"
}

getbflawname()
{

    BFLAW=$BFLAWA
    PLLMODEL=None
    PLLCOEF=-1.0;
    case $BFLAWA in
	tsaiLaw) BFLAW=pressureLimitedLaw;BFLAWNAME=.prlim; PLLMODEL=Tsai; PLLCOEF=0.5;;
	leguyLaw) BFLAW=pressureLimitedLaw;BFLAWNAME=.leguy; PLLMODEL=Leguy; PLLCOEF=8.0e12;;
	*) BFLAWNAME=""
    esac
}


mkdir -p scripts
TAGLAP=0
BASEDIR=$PWD
for SGN in 0 4
do
for BFLAWA in powerLaw tsaiLaw leguyLaw
do
    for BASENY in 20 #40 # 80
    do
	for MELT in spin melt0 melt4 melt5 melt42 melt52
	do
	    WIDTH=24
	    GEOM=isomip$WIDTH
	    TOPG=topography$WIDTH
	    THCK=thickness$WIDTH
	    TEMPLATE=inputs.isomip.template
	    ACAB=0.3
	    for ACOEF in 2e-17 2.2e-17 
	    do
		CFUNC=constfriction 	
		for solver in Chombo 
		do
		    for smod in ssa l1l2 l1l2_full;
		    do
			for lev in 0 1 2 3 4; 
			do
			    for TAGGROW in 4 #8 16
			    do
				tagcap=$(( lev - 1 )) 
				getbase
				getcre
				getnodes
				getsolveinterval
				getplotinterval
				gettaglapname
				getsgname
				gettaggrowname
				getbflawname
				getmaxtime
				NAME=$GEOM.$MELT.$smod$BFLAWNAME$BASENAME.l$lev""$TAGGROWNAME.$solver.A$ACOEF.$CFUNC""$SGNAME.a$ACAB
				SPINNAME=$GEOM.spin.$smod$BFLAWNAME$BASENAME.l$lev""$TAGGROWNAME.$solver.A$ACOEF.$CFUNC""$SGNAME.a$ACAB
				INFILE=inputs.$NAME
				sed -e s/@SGN/$SGN/ -e s/@TAGGROW/$TAGGROW/ -e s/@PLLMODEL/$PLLMODEL/ -e s/@PLLCOEF/$PLLCOEF/ -e s/@BFLAW/$BFLAW/ -e s/@BASENX/$BASENX/ -e s/@BASENY/$BASENY/ -e s/@MAXTIME/$MAXTIME/ -e s/@THCK/$THCK/ -e s/@TOPG/$TOPG/ -e s/@SPINNAME/$SPINNAME/ -e s/@TAGLAP/$TAGLAP/ -e s/@PLOTINTERVAL/$PLOTINTERVAL/ -e s/@MELT/$MELT/ -e s/@SOLVEINTERVAL/$SOLVEINTERVAL/  -e s/#$smod// -e s/#$solver// -e s/@SOLVER/$solver/ -e s/@GEOM/$GEOM/ -e s/@NAME/$NAME/ -e s/@ACOEF/$ACOEF/ -e s/@CFUNC/$CFUNC"_"$WIDTH/ -e s/@ACAB/$ACAB/ -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ $TEMPLATE > scripts/$INFILE
				sed -e s/@SPINNAME/$SPINNAME/ -e s/@SOLVEINTERVAL/$SOLVEINTERVAL/ -e s:@BASEDIR:$BASEDIR: -e s/@NAME/$NAME/   -e s/@NODES/$NODES/ -e s/@INFILE/$INFILE/ -e s/@LEV/$lev/  job.newblue.template.sh > scripts/job.newblue.$NAME.sh
			    done
			done
		    done
		done
	    done
	done
    done
done
done


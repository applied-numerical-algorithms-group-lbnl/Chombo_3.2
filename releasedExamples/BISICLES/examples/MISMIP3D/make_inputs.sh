getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}

getinres()
{
    case l$lev in
	l0) inres=128;indx=6250;;
	l1) inres=128;indx=6250;;
	l2) inres=256;indx=3125;;
	l3) inres=512;indx=1562.5;;
	l4) inres=1024;indx=781.25;;
	l5) inres=2048;indx=390.625;;
	*) echo "unknown refinement level";;
    esac
}



for smod in l1l2 ssa;
  do
  for lev in 0 1 2 3 4 5; 
    do
    tagcap=$(( lev - 1 )) 
    getcre
    getinres
    sed -e s/@INRES/$inres/ -e s/@INDX/$indx/ -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.mismip3D.stnd.template > inputs.mismip3D.stnd.$smod.l$lev

    sed -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.mismip3D.p075.template > inputs.mismip3D.p075.$smod.l$lev
  done
done

for smod in l1l2 ssa;
  do
  for lev in 6 7 8; 
    do
    tagcap=$(( lev - 1 )) 
    getcre
    sed -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.mismip3D.p075.template > inputs.mismip3D.p075.$smod.l$lev
  done
done

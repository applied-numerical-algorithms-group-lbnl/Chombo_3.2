getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}

res=0064 
smod=l1l2   
for lev in 0 1 2 3 4 5; do
    tagcap=$(( lev - 1 )) 
    getcre
    NAME=shearMargin.$smod.$res."$lev"lev
    sed -e s/@RESX/$res/ -e s/@RESY/$res/  -e s/@NAME/$NAME/ -e s/@TAGCAP/$tagcap/ -e s/@CRE/$cre/ inputs.shearMargin.template > inputs.$NAME
done




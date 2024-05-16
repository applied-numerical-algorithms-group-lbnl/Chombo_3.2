getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa) cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}

for smod in l1l2;
  do
  for lev in 0 1 2 3 4 5; 
    do
    tagcap=$(( lev - 1 )) 
    getcre
    
    sed  -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.GIS.5km.template > inputs.GIS.5km.$smod.l$lev
#    sed  -e s/@MAXLEVEL/$lev/ -e s/@SMOD/$smod/ pigv5.1km.template.config > pigv5.1km.$smod.l$lev.config
    
  done
done


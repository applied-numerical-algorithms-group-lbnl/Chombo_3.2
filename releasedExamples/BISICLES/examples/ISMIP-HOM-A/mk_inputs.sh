getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}



nlayer=10
for smod in l1l2 ssa;do
    for lev in 0; do
	tagcap=$(( lev - 1 )) 
	getcre
	for km in 5 10 20 30 40 80 160; do
	    domx="$km".0e+3
	    for ncell in 0032 0064 0128 0256; do
		name=$smod."$lev"lev.$ncell."$km"km
		subs="-e s/@NLAYER/$nlayer/ -e s/@TAGCAP/$tagcap/ -e s/@CRE/$cre/ -e s/@NCELLX/$ncell/ -e s/@NCELLY/$ncell/ -e s/@DOMX/$domx/ -e s/@DOMY/$domx/ -e s/@NAME/$name/"
		sed $subs  inputs.ISMIP-HOM-A.template > inputs.ISMIP-HOM-A.$name
		
	    done
	done
    done
done

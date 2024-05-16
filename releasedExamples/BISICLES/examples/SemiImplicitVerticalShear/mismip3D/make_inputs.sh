getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}


getbasey()
{
    case $basex in
	128) basey=16;;
	256) basey=32;; 
        *) echo "unknown basex"
    esac
}

getaddvel()
{
    case $method in 
	none) addvel=false;;
	implicit) addvel=true;;
	explicit) addvel=true;;
	*) echo "unknown method"
    esac
}

for method in implicit none explicit; do
    getaddvel
    for smod in l1l2; do
	for basex in 128 256; do
	    getbasey
	    for lev in 0 1 2 3 4 5; 
	    do
		tagcap=$(( lev - 1 )) 
		getcre
		name=mismip3D.spinup.$method.$smod.$basex.$lev"lev"
		mkdir -p $name
		sed -e s/@METHOD/$method/ -e s/@ADDVEL/$addvel/ -e s/@BASEX/$basex/ -e s/@BASEY/$basey/ -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.mismip3D.spinup.template > $name/inputs.$name
	    done
	done
    done
done


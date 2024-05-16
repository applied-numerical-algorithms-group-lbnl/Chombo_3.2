
for WD in wd0 wd3; do
for TAGCAP in 2; do
    LEV=$((TAGCAP + 1))
    for ENS in 00 06 12 24; do
	for JSPEED in 300 600; do
	    NAME="ase_fwd_ens"$ENS"_uj"$JSPEED"_"$WD"_lev"$LEV
	    echo $NAME
	    sed -e s/@WD/$WD/ -e s/@JSPEED/$JSPEED/ -e s/@NAME/$NAME/ -e s/@ENS/$ENS/ -e s/@TAGCAP/$TAGCAP/ inputs.ase_fwd_ens_template > inputs.$NAME
	done
    done
done
done

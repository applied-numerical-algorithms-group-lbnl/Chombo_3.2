
for LEV in 0 1 2 3 4 5; 
do
    tagcap=$(( LEV - 1 ))  
    sed -e s/@TAGCAP/$tagcap/ -e s/@LEV/$LEV/ inputs.twc.template > "inputs.twc."$LEV"lev"
#    sed -e s/@TAGCAP/$tagcap/ -e s/@LEV/$LEV/ inputs.twc_smooth.template > "inputs.twc_smooth."$LEV"lev"
done


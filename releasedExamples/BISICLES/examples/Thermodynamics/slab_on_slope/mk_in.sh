for NLAYER in 4 8 16 32 64 128 256 512
do
    sed -e s/@NLAYER/$NLAYER/ inputs.slab_on_slope.template > inputs.slab_on_slope_$NLAYER
done

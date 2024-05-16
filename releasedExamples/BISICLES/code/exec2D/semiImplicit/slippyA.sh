res=064
tf=inputs.slippyA.template
ex=../driver2d.Linux.64.g++.gfortran.OPT.ex
#SSA models: none model is unstable at dt = 4, but fine below that

#for diffusion in implicit none;
#do
#    for cfl in 800 400 200 100 50 25;
#    do#
#	of=inputs.slippyA.GlensLaw.$cfl.$diffusion.$res
#	sed -e s/@NSTEP/$(( 128 * 200 / cfl ))/ -e s/@CR/GlensLaw/ -e s/@CFL/$cfl/ -e s/@DIFFUSION/$diffusion/ -e s/@RES/$res/ 
#-e s/@RES/$res/ $tf > $of
#	lf=log.slippyA.GlensLaw.$cfl.$diffusion.$res
#	$ex $of > $lf
#   done
#done




#L1L2 models: the 'none' model is unstable down to dt = 1/64  so only do implicit 
#and compare with explicit with t = 1/64 

for diffusion in implicit;
do
    for cfl in 800 400 200 100 50 25;
    do
	of=inputs.slippyA.L1L2.$cfl.$diffusion.$res
	sed -e s/@NSTEP/$(( 128 * 200 / cfl ))/ -e s/@CR/L1L2/ -e s/@CFL/$cfl/ -e s/@DIFFUSION/$diffusion/ -e s/@RES/$res/ -e s/@RES/$res/ $tf > $of
	lf=log.slippyA.L1L2.$cfl.$diffusion.$res
	$ex $of > $lf
    done
done

for diffusion in explicit;
do
    for cfl in 500
    do
	of=inputs.slippyA.L1L2.$cfl.$diffusion.$res
	sed  -e 's/amr.plot_interval = 1/amr.plot_interval = 128/' -e s/@NSTEP/$(( 8192 ))/ -e s/@CR/L1L2/ -e s/@CFL/$cfl/ -e s/@DIFFUSION/$diffusion/ -e s/@RES/$res/ -e s/@RES/$res/ $tf > $of
	lf=log.slippyA.L1L2.$cfl.$diffusion.$res
	$ex $of > $lf
    done
done
